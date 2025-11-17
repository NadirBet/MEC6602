// FASMultigrid.cpp - WITH DIAGNOSTIC CHECKS
#include "FASMultigrid.hpp"
#include <iomanip> 
#include <iostream>
#include <stdexcept>
#include <algorithm>

extern void advanceExplicitRK5_LTS_MG(
    Mesh& m,
    const Config& cfg,
    BoundaryConditions& bc,
    FluxCalculator& flux,
    std::vector<Conservative>& U,
    std::vector<Conservative>& R,
    const std::vector<double>& dt_local,
    const std::vector<Conservative>* forcing
);
// ============================================================================
// 2x coarsening 
// ============================================================================

namespace {

Mesh coarsen2x(const Mesh& F)
{
    int ni_f = F.ni;
    int nj_f = F.nj;

    if (ni_f % 2 != 0 || nj_f % 2 != 0) {
        throw std::runtime_error("coarsen2x: fine ni,nj must be divisible by 2");
    }

    int ni_c = ni_f / 2;
    int nj_c = nj_f / 2;

    Mesh C;
    C.allocate(ni_c, nj_c);

    int imin_f, imax_f, jmin_f, jmax_f;
    F.getInteriorBounds(imin_f, imax_f, jmin_f, jmax_f);

    int imin_c, imax_c, jmin_c, jmax_c;
    C.getInteriorBounds(imin_c, imax_c, jmin_c, jmax_c);

    // physical node mapping: (iC_phys,jC_phys) → fine (2*iC_phys, 2*jC_phys)
    for (int jC_phys = 0; jC_phys <= nj_c; ++jC_phys) {
        int jC = jmin_c + jC_phys;
        int jF = jmin_f + 2 * jC_phys;

        for (int iC_phys = 0; iC_phys <= ni_c; ++iC_phys) {
            int iC = imin_c + iC_phys;
            int iF = imin_f + 2 * iC_phys;

            int nC = C.nodeIndex(iC, jC);
            int nF = F.nodeIndex(iF, jF);

            C.xNodes[nC] = F.xNodes[nF];
            C.yNodes[nC] = F.yNodes[nF];
        }
    }

    C.fillGhostNodes();
    C.computeCellCenters();
    C.computeCellAreas();
    C.computeFaceGeometry();

    return C;
}

} // anonymous namespace

// ============================================================================
// Constructor with per-level CFL + sweep control
// ============================================================================

FASMultigridSolver::FASMultigridSolver(
    const Config&  baseCfg,
    bool           use_simple_ff,
    bool           freestream_all,
    const std::string& finestMeshFile,
    int            numLevels,
    const LevelControl& defaultCtrl)
    : useSimpleFarfield(use_simple_ff)
    , freestreamEverywhere(freestream_all)
    , nLevels(numLevels)
    , finestMeshFile_(finestMeshFile)
{
    if (nLevels < 1) {
        throw std::runtime_error("FASMultigridSolver: numLevels must be >= 1");
    }

    levels.reserve(nLevels);

    for (int ell = 0; ell < nLevels; ++ell) {
        levels.emplace_back(
            baseCfg,
            useSimpleFarfield,
            freestreamEverywhere,
            ell,
            defaultCtrl
        );

        Level& L = levels.back();

        // ----------------- CFL per level (keep as before) -----------------
        if (ell == 0) {
            // Finest grid
            L.cfg.CFL = 2;
        } else if (ell == nLevels - 1) {
            // Coarsest grid: be more aggressive
            L.cfg.CFL = 4;
        } else {
            // Intermediate grids
            L.cfg.CFL = 3;
        }

        // ----------------- SWEEP CONTROL---------

    if (ell == 0) {
        // Finest grid: moderate smoothing
        L.ctrl.minSweeps = 7;
        L.ctrl.maxSweeps = 10;
    } else if (ell == nLevels - 1) {
        // Coarsest grid: aggressive smoothing (no recursion, so smooth more)
        L.ctrl.minSweeps = 8;
        L.ctrl.maxSweeps = 8;
    } else {
        // Intermediate grids
        L.ctrl.minSweeps = 8;
        L.ctrl.maxSweeps = 8;
    }
        L.ctrl.stallWindow     = 0;    // disable stall detection
        L.ctrl.stallTol        = 0.0;  // unused when stallWindow=0
        L.ctrl.targetReduction = 0.0;  // reduction <= 0.0 is never true
    }
}
// ============================================================================
// initialize
// ============================================================================

void FASMultigridSolver::initialize()
{
    // 1) Read finest mesh once
    Mesh fine;
    //std::cout << "[MG] Reading finest mesh '" << finestMeshFile_ << "'\n";
    fine.readPlot3D(finestMeshFile_);
    /**std::cout << "[MG] Finest mesh: ni x nj = "
              << fine.ni << " x " << fine.nj << "\n";*/

    // 2) Build hierarchy of meshes by 2x coarsening
    levels[0].mesh = fine;

    for (int ell = 1; ell < nLevels; ++ell) {
        levels[ell].mesh = coarsen2x(levels[ell-1].mesh);
    /**    std::cout << "[MG] Level " << ell
                  << " coarsened: ni x nj = "
                  << levels[ell].mesh.ni << " x " << levels[ell].mesh.nj << "\n";*/
    }

    // 3) Initialize state / residuals on each level
    for (int ell = 0; ell < nLevels; ++ell) {
        Level& L = levels[ell];

        L.cfg.initialize();

        const int nCells = L.mesh.niTotal * L.mesh.njTotal;
        Conservative Uinf = L.cfg.getFreestream().toConservative(L.cfg.gamma);

        L.U.assign(nCells, Uinf);
        L.R.assign(nCells, Conservative());
        L.dt.assign(nCells, 0.0);
        L.forcing.assign(nCells, Conservative());
        L.scratch.assign(nCells, Conservative());

        // Finest level has zero forcing; coarse levels will be filled by
        // restrictSolution during the first V-cycle.

        // Initial ghosts
        L.bc.initializeGhosts(L.U, L.mesh);

        // Initial spectral radii
        std::vector<Primitive> W(nCells);
        for (int c = 0; c < nCells; ++c)
            W[c] = Primitive(L.U[c], L.cfg.gamma);
        L.mesh.computeSpectralRadius(W, L.cfg.gamma);

        // Initial residual (includes Q_F if index>0, but Q_F=0 here)
        computeResidual(L);
        L.bestResL2 = L.currentResL2;
        L.resHistory.clear();
        L.resHistory.push_back(L.currentResL2);

    /**     std::cout << "[MG] Level " << ell
                  << " init L2 residual = " << std::scientific
                  << L.currentResL2 << "\n";*/
    }

    /**std::cout << "[MG] Initialized " << nLevels << " level(s).\n";*/
}

// ============================================================================
// runVCycles
// ============================================================================


void FASMultigridSolver::runVCycles(int nCycles, double resTol)
{
    // Open residual file for finest level
    std::ofstream residual_file("residuals_mg.txt");
    if (residual_file.is_open()) {
        residual_file << std::scientific << std::setprecision(12);
        residual_file << "# V-Cycle    L2[rho]         L2[rhou]        L2[rhov]        L2[rhoE]        L2_total\n";
        residual_file << "# " << std::string(100, '-') << "\n";
    } else {
        std::cerr << "[MG] WARNING: Could not open residuals_mg.txt for writing\n";
    }

    // Print header
    std::cout << "[MG] Starting V-cycles (target residual: " << resTol << ")\n";
    std::cout << "[MG] " << std::string(80, '-') << "\n";

    Level& L0 = levels[0];  // Finest level
    bool converged = false;
    int final_cycle = 0;
    
    for (int cyc = 0; cyc < nCycles; ++cyc) {
        std::cout << "[MG] === V-cycle " << cyc << " ===\n";
        
        // Perform V-cycle
        vCycle(0);

        // Recompute finest residual (no forcing on finest)
        computeResidual(L0);
        double R_total = L0.currentResL2;

        // Compute component residuals for file output
        int imin, imax, jmin, jmax;
        L0.mesh.getInteriorBounds(imin, imax, jmin, jmax);
        
        double acc2[4] = {0, 0, 0, 0};
        double volSum = 0.0;
        
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int c = L0.mesh.cellIndex(i, j);
                double V = L0.mesh.cellArea[c];
                const Conservative& R = L0.R[c];
                
                acc2[0] += V * R.rho  * R.rho;
                acc2[1] += V * R.rhou * R.rhou;
                acc2[2] += V * R.rhov * R.rhov;
                acc2[3] += V * R.rhoE * R.rhoE;
                volSum += V;
            }
        }
        
        double L2[4];
        for (int k = 0; k < 4; ++k) {
            L2[k] = std::sqrt(acc2[k] / volSum);
        }

        // Print to console (finest level only)
        std::cout << "[MG] V-cycle " << cyc 
                  << ": L2_total = " << std::scientific << R_total << "\n";

        // Write to residual file
        if (residual_file.is_open()) {
            residual_file << std::setw(8) << cyc << "  "
                          << std::setw(15) << L2[0] << "  "
                          << std::setw(15) << L2[1] << "  "
                          << std::setw(15) << L2[2] << "  "
                          << std::setw(15) << L2[3] << "  "
                          << std::setw(15) << R_total << "\n";
            residual_file.flush();
        }

        // Check convergence
        if (R_total < resTol) {
            converged = true;
            final_cycle = cyc;
            std::cout << "[MG] CONVERGED: Reached tolerance " << resTol
                      << " at V-cycle " << cyc << "\n";
            std::cout << "[MG] " << std::string(80, '-') << "\n";
            break;
        }
        
        final_cycle = cyc;
    }

    if (residual_file.is_open()) {
        residual_file.close();
    }

    // Status message
    if (converged) {
        std::cout << "[MG] Converged at V-cycle " << final_cycle 
                  << " with L2_total = " << L0.currentResL2 << "\n";
    } else {
        std::cout << "[MG] Reached max V-cycles (" << nCycles 
                  << ") with L2_total = " << L0.currentResL2 << "\n";
    }
    
    std::cout << "[MG] Final finest-level L2 residual = " << L0.currentResL2 << "\n";
}

// ============================================================================
// resetFinestFromCurrent
// ============================================================================

void FASMultigridSolver::resetFinestFromCurrent()
{
    Level& L = levels[0];

    // 1) Ghosts
    L.bc.initializeGhosts(L.U, L.mesh);

    // 2) Spectral radii
    {
        std::vector<Primitive> W(L.U.size());
        for (std::size_t c = 0; c < L.U.size(); ++c) {
            W[c] = Primitive(L.U[c], L.cfg.gamma);
        }
        L.mesh.computeSpectralRadius(W, L.cfg.gamma);
    }

    // 3) Residual (no forcing on finest)
    L.flux.computeResidual(L.U, L.R, L.mesh, L.cfg);

    // 4) L2 + history
    L.currentResL2 = computeL2Residual(L);
    L.bestResL2    = L.currentResL2;
    L.resHistory.clear();
    L.resHistory.push_back(L.currentResL2);
}

// ============================================================================
// finestSolution accessors
// ============================================================================

std::vector<Conservative>& FASMultigridSolver::finestSolution()
{
    return finestLevel().U;
}

const std::vector<Conservative>& FASMultigridSolver::finestSolution() const
{
    return finestLevel().U;
}

// ============================================================================
// vCycle
// ============================================================================

void FASMultigridSolver::vCycle(int ell)
{
    Level& L = levels[ell];

    // 1) Pre-smoothing on level ell
    smoothWithControl(L);

    // 2) Compute fresh forced residual on this level
    computeResidual(L);

    // 3) If we are at the coarsest level, stop here
    if (ell == nLevels - 1) {
        return;
    }

    // 4) Build coarse level from this level (FAS restriction)
    restrictSolution(ell, ell + 1);

    // 5) V-cycle on coarser level
    vCycle(ell + 1);

    // 6) Prolongate correction back to this level
    prolongateCorrection(ell + 1, ell);


}

// ============================================================================
// smoothWithControl
// ============================================================================

void FASMultigridSolver::smoothWithControl(Level& L)
{
    // Baseline forced residual at the start of THIS smoothing call
    computeResidual(L);  // sets L.currentResL2 and L.R

    L.bestResL2 = L.currentResL2;
    L.resHistory.clear();
    L.resHistory.push_back(L.currentResL2);

    const double R0 = std::max(L.currentResL2, 1e-300);  // avoid division by zero

    int sweeps = 0;
    while (true) {
        ++sweeps;

        // One smoothing sweep (updates U and recomputes residual)
        smootherSweep(L);   // at exit, computeResidual(L) has been called

        L.bestResL2 = std::min(L.bestResL2, L.currentResL2);

        L.resHistory.push_back(L.currentResL2);
        const std::size_t maxHist = static_cast<std::size_t>(L.ctrl.stallWindow + 1);
        while (L.resHistory.size() > maxHist && maxHist > 0) {
            L.resHistory.pop_front();
        }

        /** // ---- PRINT FOR EVERY LEVEL ----
        if (sweeps == 1) {
            std::cout << "[MG] level " << L.index
                      << " R0 = " << std::scientific << R0 << "\n";
        }
        std::cout << "[MG] level " << L.index
                  << " sweep " << sweeps
                  << "  L2 = " << std::scientific << L.currentResL2
                  << " (rel=" << L.currentResL2 / R0 << ")\n";*/
        // -------------------------------

        const bool reachedMin = (sweeps >= L.ctrl.minSweeps);
        const bool reachedMax = (sweeps >= L.ctrl.maxSweeps);

        bool stalled = false;
        if (L.ctrl.stallWindow > 0) {
            stalled = reachedMin && isStalled(L);
        }

        const double reduction  = L.currentResL2 / R0;  // <1 => good
        const bool   enoughDrop = (reduction <= L.ctrl.targetReduction);

        if (reachedMin && (enoughDrop || stalled || L.ctrl.stallWindow == 0)) {
            // stallWindow==0 → we just do exactly minSweeps == maxSweeps
            if (sweeps >= L.ctrl.minSweeps) break;
        }
        if (reachedMax) {
            break;
        }
    }


    if (L.index == 0) {
        std::cout << "  [MG] Finest level smoothing: sweeps=" << sweeps
                << ", best L2=" << std::scientific << L.bestResL2
                << ", final L2=" << L.currentResL2 << "\n";
    }
}


// ============================================================================
// smootherSweep  — one smoothing sweep
//   - ALL levels: RK5-LTS using the same routine as single-grid.
//   - FAS forcing still handled in computeResidual() (for coarse levels).
// ============================================================================

void FASMultigridSolver::smootherSweep(Level& L)
{
    // ============================================================================
    // CHECK 3: Verify forcing is being applied correctly during time integration
    // ============================================================================
    if (L.index > 0) {
        // Before time stepping, compute pure PDE residual
        L.bc.apply(L.U, L.mesh);
        {
            std::vector<Primitive> W(L.U.size());
            for (std::size_t c = 0; c < L.U.size(); ++c) {
                W[c] = Primitive(L.U[c], L.cfg.gamma);
            }
            L.mesh.computeSpectralRadius(W, L.cfg.gamma);
        }
        
        std::vector<Conservative> R_pure(L.U.size());
        L.flux.computeResidual(L.U, R_pure, L.mesh, L.cfg);
        
        // Now compute what should be the forced residual manually
        int imin, imax, jmin, jmax;
        L.mesh.getInteriorBounds(imin, imax, jmin, jmax);
        
        double max_diff_rho = 0.0;
        double max_diff_rhou = 0.0;
        double max_diff_rhov = 0.0;
        double max_diff_rhoE = 0.0;
        
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int c = L.mesh.cellIndex(i, j);
                
                // What the forced residual SHOULD be
                double expected_rho  = R_pure[c].rho  + L.forcing[c].rho;
                double expected_rhou = R_pure[c].rhou + L.forcing[c].rhou;
                double expected_rhov = R_pure[c].rhov + L.forcing[c].rhov;
                double expected_rhoE = R_pure[c].rhoE + L.forcing[c].rhoE;
                
                // What it actually is (from previous computeResidual call)
                double actual_rho  = L.R[c].rho;
                double actual_rhou = L.R[c].rhou;
                double actual_rhov = L.R[c].rhov;
                double actual_rhoE = L.R[c].rhoE;
                
                max_diff_rho  = std::max(max_diff_rho,  std::abs(expected_rho  - actual_rho));
                max_diff_rhou = std::max(max_diff_rhou, std::abs(expected_rhou - actual_rhou));
                max_diff_rhov = std::max(max_diff_rhov, std::abs(expected_rhov - actual_rhov));
                max_diff_rhoE = std::max(max_diff_rhoE, std::abs(expected_rhoE - actual_rhoE));
            }
        }
        
        /**std::cout << "[FORCING CHECK] Level " << L.index 
                  << " forcing application: max|R_computed - (R_pure + Q_F)| = ["
                  << std::scientific 
                  << max_diff_rho << ", "
                  << max_diff_rhou << ", "
                  << max_diff_rhov << ", "
                  << max_diff_rhoE << "]\n";
        
        if (max_diff_rho > 1e-12 || max_diff_rhou > 1e-12 || 
            max_diff_rhov > 1e-12 || max_diff_rhoE > 1e-12) {
            std::cout << "[FORCING CHECK] WARNING: Forcing not being applied correctly!\n";
        }*/
    }
    // ============================================================================
    
    // Apply boundary conditions
    L.bc.apply(L.U, L.mesh);

    // Compute spectral radii based on current solution
    {
        std::vector<Primitive> W(L.U.size());
        for (std::size_t c = 0; c < L.U.size(); ++c) {
            W[c] = Primitive(L.U[c], L.cfg.gamma);
        }
        L.mesh.computeSpectralRadius(W, L.cfg.gamma);
    }

    // Build local time steps for this level
    buildLocalDt(L);

    // Call the multigrid-aware RK5 routine
    // Pass forcing pointer only for coarse levels (index > 0)
    const std::vector<Conservative>* forcing_ptr = (L.index > 0) ? &L.forcing : nullptr;
    
    advanceExplicitRK5_LTS_MG(
        L.mesh, 
        L.cfg, 
        L.bc, 
        L.flux, 
        L.U, 
        L.R, 
        L.dt,
        forcing_ptr
    );

    // The residual R now contains the forced residual (for coarse levels)
    // or pure PDE residual (for finest level), ready for convergence checking
    L.currentResL2 = computeL2Residual(L);
}

// ============================================================================
// buildLocalDt  (CFL * Area / (lambdaI + lambdaJ))
// ============================================================================

void FASMultigridSolver::buildLocalDt(Level& L)
{
    Mesh&   m = L.mesh;
    Config& c = L.cfg;

    int imin, imax, jmin, jmax;
    m.getInteriorBounds(imin, imax, jmin, jmax);

    const double EPS = 1e-14;

    L.dt.assign(m.niTotal * m.njTotal, 0.0);

    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            const int    cell  = m.cellIndex(i, j);
            const double denom = m.lambdaI[cell] + m.lambdaJ[cell];
            const double omega = m.cellArea[cell];
            L.dt[cell] = c.CFL * omega / std::max(denom, EPS);
        }
    }
}

// ============================================================================
// computeResidual  (BC + spectral radii + flux + FAS forcing on coarse levels)
// ============================================================================

void FASMultigridSolver::computeResidual(Level& L)
{
    // Apply BCs to U
    L.bc.apply(L.U, L.mesh);

    // Update spectral radii
    {
        std::vector<Primitive> W(L.U.size());
        for (std::size_t c = 0; c < L.U.size(); ++c) {
            W[c] = Primitive(L.U[c], L.cfg.gamma);
        }
        L.mesh.computeSpectralRadius(W, L.cfg.gamma);
    }

    // Compute PDE residual
    L.flux.computeResidual(L.U, L.R, L.mesh, L.cfg);

    // Add FAS forcing on coarse levels: R_F = R + Q_F
    if (L.index > 0) {
        int imin, imax, jmin, jmax;
        L.mesh.getInteriorBounds(imin, imax, jmin, jmax);

        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int c = L.mesh.cellIndex(i, j);
                L.R[c].rho  += L.forcing[c].rho;
                L.R[c].rhou += L.forcing[c].rhou;
                L.R[c].rhov += L.forcing[c].rhov;
                L.R[c].rhoE += L.forcing[c].rhoE;
            }
        }
    }

    // L2 norm of (possibly forced) residual
    L.currentResL2 = computeL2Residual(L);
}

// ============================================================================
// computeL2Residual   (volume-weighted L2 over physical cells)
// ============================================================================

double FASMultigridSolver::computeL2Residual(const Level& L) const
{
    const Mesh& m = L.mesh;

    int imin, imax, jmin, jmax;
    m.getInteriorBounds(imin, imax, jmin, jmax);

    double acc2_rho  = 0.0;
    double acc2_rhou = 0.0;
    double acc2_rhov = 0.0;
    double acc2_rhoE = 0.0;
    double volSum    = 0.0;

    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            const int c = m.cellIndex(i, j);
            const double V = m.cellArea[c];
            const Conservative& R = L.R[c];

            acc2_rho  += V * R.rho  * R.rho;
            acc2_rhou += V * R.rhou * R.rhou;
            acc2_rhov += V * R.rhov * R.rhov;
            acc2_rhoE += V * R.rhoE * R.rhoE;
            volSum    += V;
        }
    }

    if (volSum <= 0.0) {
        return 0.0;
    }

    const double invVol = 1.0 / volSum;
    const double L2_rho  = std::sqrt(acc2_rho  * invVol);
    const double L2_rhou = std::sqrt(acc2_rhou * invVol);
    const double L2_rhov = std::sqrt(acc2_rhov * invVol);
    const double L2_rhoE = std::sqrt(acc2_rhoE * invVol);

    return std::sqrt(L2_rho*L2_rho + L2_rhou*L2_rhou +
                     L2_rhov*L2_rhov + L2_rhoE*L2_rhoE);
}

// ============================================================================
// isStalled
// ============================================================================

bool FASMultigridSolver::isStalled(const Level& L) const
{
    // With stallWindow=0 on all levels (for single V-cycle), this is never used.
    if (L.ctrl.stallWindow <= 0) {
        return false;
    }

    const int W = L.ctrl.stallWindow;
    if (L.resHistory.size() < static_cast<std::size_t>(W + 1)) {
        return false; // not enough history yet
    }

    const double r_old = L.resHistory.front();
    const double r_new = L.resHistory.back();

    if (r_new <= 0.0 || !std::isfinite(r_new) || !std::isfinite(r_old)) {
        return false;
    }

    const double ratio = r_old / r_new; // >1 = good decrease, ~1 = flat
    if (ratio <= 1.0) return false;     // increased or unchanged

    const double dropFrac = 1.0 - 1.0/ratio; // e.g. ratio=1.05 → ~4.76% drop
    return (dropFrac < L.ctrl.stallTol);
}

// ============================================================================
// restrictSolution  (fine → coarse, FAS construction of coarse level)
// ============================================================================

void FASMultigridSolver::restrictSolution(int fineIdx, int coarseIdx)
{
    Level& F = levels[fineIdx];
    Level& C = levels[coarseIdx];

    int iFmin, iFmax, jFmin, jFmax;
    int iCmin, iCmax, jCmin, jCmax;
    F.mesh.getInteriorBounds(iFmin, iFmax, jFmin, jFmax);
    C.mesh.getInteriorBounds(iCmin, iCmax, jCmin, jCmax);

    const int nCellsC = C.mesh.niTotal * C.mesh.njTotal;

    // Restricted forced residual on coarse grid
    std::vector<Conservative> R_restr(nCellsC);
    for (int c = 0; c < nCellsC; ++c) {
        R_restr[c] = Conservative();
    }

    C.forcing.assign(nCellsC, Conservative());
    C.R.assign(nCellsC, Conservative());

    // 1) Interpolate solution and restrict residual (as before)
    for (int jC = jCmin; jC < jCmax; ++jC) {
        const int jC_rel = jC - jCmin;
        const int jF     = jFmin + 2 * jC_rel;
        const int jF1    = jF + 1;

        for (int iC = iCmin; iC < iCmax; ++iC) {
            const int iC_rel = iC - iCmin;
            const int iF     = iFmin + 2 * iC_rel;
            const int iF1    = iF + 1;

            const int cF00 = F.mesh.cellIndex(iF,  jF );
            const int cF10 = F.mesh.cellIndex(iF1, jF );
            const int cF01 = F.mesh.cellIndex(iF,  jF1);
            const int cF11 = F.mesh.cellIndex(iF1, jF1);

            const double A00 = F.mesh.cellArea[cF00];
            const double A10 = F.mesh.cellArea[cF10];
            const double A01 = F.mesh.cellArea[cF01];
            const double A11 = F.mesh.cellArea[cF11];
            const double Atot = A00 + A10 + A01 + A11;

            Conservative Uc{};
            if (Atot > 0.0) {
                const Conservative& U00 = F.U[cF00];
                const Conservative& U10 = F.U[cF10];
                const Conservative& U01 = F.U[cF01];
                const Conservative& U11 = F.U[cF11];

                Uc.rho  = (U00.rho  * A00 + U10.rho  * A10 + U01.rho  * A01 + U11.rho  * A11) / Atot;
                Uc.rhou = (U00.rhou * A00 + U10.rhou * A10 + U01.rhou * A01 + U11.rhou * A11) / Atot;
                Uc.rhov = (U00.rhov * A00 + U10.rhov * A10 + U01.rhov * A01 + U11.rhov * A11) / Atot;
                Uc.rhoE = (U00.rhoE * A00 + U10.rhoE * A10 + U01.rhoE * A01 + U11.rhoE * A11) / Atot;
            }

            const int cC = C.mesh.cellIndex(iC, jC);
            C.U[cC] = Uc;

            // Restrict the forced residual from fine grid (F.R already includes forcing on fine)
            const Conservative& RF00 = F.R[cF00];
            const Conservative& RF10 = F.R[cF10];
            const Conservative& RF01 = F.R[cF01];
            const Conservative& RF11 = F.R[cF11];




            if (Atot > 0.0) {
                R_restr[cC].rho  = (RF00.rho  * A00 + RF10.rho  * A10 + 
                                    RF01.rho  * A01 + RF11.rho  * A11) / Atot;
                R_restr[cC].rhou = (RF00.rhou * A00 + RF10.rhou * A10 + 
                                    RF01.rhou * A01 + RF11.rhou * A11) / Atot;
                R_restr[cC].rhov = (RF00.rhov * A00 + RF10.rhov * A10 + 
                                    RF01.rhov * A01 + RF11.rhov * A11) / Atot;
                R_restr[cC].rhoE = (RF00.rhoE * A00 + RF10.rhoE * A10 + 
                                    RF01.rhoE * A01 + RF11.rhoE * A11) / Atot;
            } else {
                R_restr[cC] = Conservative();
            }
        }
    }

    // 1b) Apply BCs and compute spectral radii on coarse grid
    C.bc.initializeGhosts(C.U, C.mesh);
    {
        std::vector<Primitive> Wc(C.U.size());
        for (std::size_t c = 0; c < C.U.size(); ++c) {
            Wc[c] = Primitive(C.U[c], C.cfg.gamma);
        }
        C.mesh.computeSpectralRadius(Wc, C.cfg.gamma);
    }

    // 2) Compute pure PDE residual on coarse grid at restricted solution
    std::vector<Conservative> R2h0(nCellsC);
    C.flux.computeResidual(C.U, R2h0, C.mesh, C.cfg);

    // 3) Build FAS forcing: Q_F = R_restr - R2h0
    for (int jC = jCmin; jC < jCmax; ++jC) {
        for (int iC = iCmin; iC < iCmax; ++iC) {
            int cC = C.mesh.cellIndex(iC, jC);
            C.forcing[cC].rho  = R2h0[cC].rho  - R_restr[cC].rho;
            C.forcing[cC].rhou = R2h0[cC].rhou - R_restr[cC].rhou;
            C.forcing[cC].rhov = R2h0[cC].rhov - R_restr[cC].rhov;
            C.forcing[cC].rhoE = R2h0[cC].rhoE - R_restr[cC].rhoE;
        }
    }

    // ============================================================================
    // VERIFICATION: Check FAS consistency condition
    // At this point, R2h0 + Q_F should equal R_restr
    // ============================================================================
    double max_error_rho = 0.0;
    double max_error_rhou = 0.0;
    double max_error_rhov = 0.0;
    double max_error_rhoE = 0.0;
    
    for (int jC = jCmin; jC < jCmax; ++jC) {
        for (int iC = iCmin; iC < iCmax; ++iC) {
            int cC = C.mesh.cellIndex(iC, jC);
            
            // Compute R2h0 + Q_F
            double reconstructed_rho  = R2h0[cC].rho  + C.forcing[cC].rho;
            double reconstructed_rhou = R2h0[cC].rhou + C.forcing[cC].rhou;
            double reconstructed_rhov = R2h0[cC].rhov + C.forcing[cC].rhov;
            double reconstructed_rhoE = R2h0[cC].rhoE + C.forcing[cC].rhoE;
            
            // Compare with R_restr
            double err_rho  = std::abs(reconstructed_rho  - R_restr[cC].rho);
            double err_rhou = std::abs(reconstructed_rhou - R_restr[cC].rhou);
            double err_rhov = std::abs(reconstructed_rhov - R_restr[cC].rhov);
            double err_rhoE = std::abs(reconstructed_rhoE - R_restr[cC].rhoE);
            
            max_error_rho  = std::max(max_error_rho,  err_rho);
            max_error_rhou = std::max(max_error_rhou, err_rhou);
            max_error_rhov = std::max(max_error_rhov, err_rhov);
            max_error_rhoE = std::max(max_error_rhoE, err_rhoE);
        }
    }
    
    /**std::cout << "[FAS CHECK] Level " << coarseIdx 
              << " FAS consistency: max|R2h0 + Q_F - R_restr| = ["
              << std::scientific 
              << max_error_rho << ", "
              << max_error_rhou << ", "
              << max_error_rhov << ", "
              << max_error_rhoE << "]\n";*/
    
    // These errors should be at machine precision (around 1e-15 or smaller)
    // If they're large, there's a sign error or algebraic mistake
    const double tol = 1e-12;
    if (max_error_rho > tol || max_error_rhou > tol || 
        max_error_rhov > tol || max_error_rhoE > tol) {
        //std::cout << "[FAS CHECK] WARNING: FAS consistency violated! Check sign conventions.\n";
    }
    // ============================================================================

    // 4) Initialize coarse residual with forced residual
    for (int c = 0; c < nCellsC; ++c) {
        C.R[c].rho  = R2h0[c].rho  + C.forcing[c].rho;
        C.R[c].rhou = R2h0[c].rhou + C.forcing[c].rhou;
        C.R[c].rhov = R2h0[c].rhov + C.forcing[c].rhov;
        C.R[c].rhoE = R2h0[c].rhoE + C.forcing[c].rhoE;
    }

    C.currentResL2 = computeL2Residual(C);
    C.bestResL2    = C.currentResL2;
    C.resHistory.clear();
    C.resHistory.push_back(C.currentResL2);
}

// ============================================================================
// prolongateCorrection  (coarse → fine)
// ============================================================================

void FASMultigridSolver::prolongateCorrection(int coarseIdx, int fineIdx)
{
    Level& C = levels[coarseIdx];
    Level& F = levels[fineIdx];

    int iFmin, iFmax, jFmin, jFmax;
    int iCmin, iCmax, jCmin, jCmax;
    F.mesh.getInteriorBounds(iFmin, iFmax, jFmin, jFmax);
    C.mesh.getInteriorBounds(iCmin, iCmax, jCmin, jCmax);

    // 1) Build restricted fine solution on the coarse grid (R U_h)
    std::vector<Conservative> U_restr(C.mesh.niTotal * C.mesh.njTotal);

    for (int jC = jCmin; jC < jCmax; ++jC) {
        int jC_rel = jC - jCmin;
        int jF     = jFmin + 2 * jC_rel;
        int jF1    = jF + 1;

        for (int iC = iCmin; iC < iCmax; ++iC) {
            int iC_rel = iC - iCmin;
            int iF     = iFmin + 2 * iC_rel;
            int iF1    = iF + 1;

            int cF00 = F.mesh.cellIndex(iF,  jF );
            int cF10 = F.mesh.cellIndex(iF1, jF );
            int cF01 = F.mesh.cellIndex(iF,  jF1);
            int cF11 = F.mesh.cellIndex(iF1, jF1);

            double A00 = F.mesh.cellArea[cF00];
            double A10 = F.mesh.cellArea[cF10];
            double A01 = F.mesh.cellArea[cF01];
            double A11 = F.mesh.cellArea[cF11];

            double Atot = A00 + A10 + A01 + A11;

            Conservative Uc{};
            if (Atot > 0.0) {
                const Conservative& U00 = F.U[cF00];
                const Conservative& U10 = F.U[cF10];
                const Conservative& U01 = F.U[cF01];
                const Conservative& U11 = F.U[cF11];

                Uc.rho  = (U00.rho  * A00 + U10.rho  * A10 + U01.rho  * A01 + U11.rho  * A11) / Atot;
                Uc.rhou = (U00.rhou * A00 + U10.rhou * A10 + U01.rhou * A01 + U11.rhou * A11) / Atot;
                Uc.rhov = (U00.rhov * A00 + U10.rhov * A10 + U01.rhov * A01 + U11.rhov * A11) / Atot;
                Uc.rhoE = (U00.rhoE * A00 + U10.rhoE * A10 + U01.rhoE * A01 + U11.rhoE * A11) / Atot;
            }

            int cC = C.mesh.cellIndex(iC, jC);
            U_restr[cC] = Uc;
        }
    }

    // 2) Coarse correction e_H = U_H (smoothed) - R(U_h)
    std::vector<Conservative> eC(C.mesh.niTotal * C.mesh.njTotal);
    for (int jC = jCmin; jC < jCmax; ++jC) {
        for (int iC = iCmin; iC < iCmax; ++iC) {
            int cC = C.mesh.cellIndex(iC, jC);
            eC[cC].rho  = C.U[cC].rho  - U_restr[cC].rho;
            eC[cC].rhou = C.U[cC].rhou - U_restr[cC].rhou;
            eC[cC].rhov = C.U[cC].rhov - U_restr[cC].rhov;
            eC[cC].rhoE = C.U[cC].rhoE - U_restr[cC].rhoE;
        }
    }

    // 3) Prolongate this correction e_H to fine (P e_H) using 9-3-3-1 pattern
    std::fill(F.scratch.begin(), F.scratch.end(), Conservative{0.0,0.0,0.0,0.0});

    for (int jC = jCmin; jC < jCmax; ++jC) {
        int jC_rel = jC - jCmin;
        int jF     = jFmin + 2 * jC_rel;
        int jF1    = jF + 1;

        for (int iC = iCmin; iC < iCmax; ++iC) {
            int iC_rel = iC - iCmin;
            int iF     = iFmin + 2 * iC_rel;
            int iF1    = iF + 1;

            int cC = C.mesh.cellIndex(iC, jC);
            const Conservative& ec = eC[cC];

            const double w00 = 9.0 / 16.0;
            const double w10 = 3.0 / 16.0;
            const double w01 = 3.0 / 16.0;
            const double w11 = 1.0 / 16.0;

            int cF00 = F.mesh.cellIndex(iF,  jF );
            int cF10 = F.mesh.cellIndex(iF1, jF );
            int cF01 = F.mesh.cellIndex(iF,  jF1);
            int cF11 = F.mesh.cellIndex(iF1, jF1);

            auto addWeighted = [&](int cF, double w) {
                F.scratch[cF].rho  += w * ec.rho;
                F.scratch[cF].rhou += w * ec.rhou;
                F.scratch[cF].rhov += w * ec.rhov;
                F.scratch[cF].rhoE += w * ec.rhoE;
            };

            addWeighted(cF00, w00);
            addWeighted(cF10, w10);
            addWeighted(cF01, w01);
            addWeighted(cF11, w11);
        }
    }

    // ============================================================================
    // CHECK 1: Verify correction reduces residual (or at least doesn't increase it)
    // ============================================================================
    
    // Compute residual before correction
    F.bc.apply(F.U, F.mesh);
    {
        std::vector<Primitive> W(F.U.size());
        for (std::size_t c = 0; c < F.U.size(); ++c) {
            W[c] = Primitive(F.U[c], F.cfg.gamma);
        }
        F.mesh.computeSpectralRadius(W, F.cfg.gamma);
    }
    std::vector<Conservative> R_before(F.mesh.niTotal * F.mesh.njTotal);
    F.flux.computeResidual(F.U, R_before, F.mesh, F.cfg);
    
    // Add forcing if this is a coarse level
    if (fineIdx > 0) {
        int imin, imax, jmin, jmax;
        F.mesh.getInteriorBounds(imin, imax, jmin, jmax);
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int c = F.mesh.cellIndex(i, j);
                R_before[c].rho  += F.forcing[c].rho;
                R_before[c].rhou += F.forcing[c].rhou;
                R_before[c].rhov += F.forcing[c].rhov;
                R_before[c].rhoE += F.forcing[c].rhoE;
            }
        }
    }
    
    double L2_before = 0.0;
    {
        int imin, imax, jmin, jmax;
        F.mesh.getInteriorBounds(imin, imax, jmin, jmax);
        double acc2 = 0.0, volSum = 0.0;
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int c = F.mesh.cellIndex(i, j);
                double V = F.mesh.cellArea[c];
                const Conservative& R = R_before[c];
                acc2 += V * (R.rho*R.rho + R.rhou*R.rhou + R.rhov*R.rhov + R.rhoE*R.rhoE);
                volSum += V;
            }
        }
        L2_before = std::sqrt(acc2 / volSum);
    }
    
    // ============================================================================

    // 4) Apply correction on fine grid
    int iFmin2, iFmax2, jFmin2, jFmax2;
    F.mesh.getInteriorBounds(iFmin2, iFmax2, jFmin2, jFmax2);
    

    /** 
    // ============================================================================
    // CHECK 2: Verify prolongation only touches interior cells
    // ============================================================================
    bool touches_boundary = false;
    for (int jF = 0; jF < F.mesh.njTotal; ++jF) {
        for (int iF = 0; iF < F.mesh.niTotal; ++iF) {
            // Check if this is a boundary cell
            bool is_interior = (iF >= iFmin2 && iF < iFmax2 && jF >= jFmin2 && jF < jFmax2);
            
            if (!is_interior) {
                int cF = F.mesh.cellIndex(iF, jF);
                // Check if scratch (the correction) is nonzero here
                if (std::abs(F.scratch[cF].rho) > 1e-14 || 
                    std::abs(F.scratch[cF].rhou) > 1e-14 ||
                    std::abs(F.scratch[cF].rhov) > 1e-14 ||
                    std::abs(F.scratch[cF].rhoE) > 1e-14) {
                    touches_boundary = true;
                    std::cout << "[BOUNDARY CHECK] WARNING: Prolongation touched boundary cell ("
                              << iF << "," << jF << ")\n";
                }
            }
        }
    }
    
    if (!touches_boundary) {
        std::cout << "[BOUNDARY CHECK] Level " << fineIdx 
                  << ": Prolongation correctly restricted to interior cells only.\n";
    }
    // ============================================================================
    */

    
    for (int jF = jFmin2; jF < jFmax2; ++jF) {
        for (int iF = iFmin2; iF < iFmax2; ++iF) {
            int cF = F.mesh.cellIndex(iF, jF);
            F.U[cF].rho  += F.scratch[cF].rho;
            F.U[cF].rhou += F.scratch[cF].rhou;
            F.U[cF].rhov += F.scratch[cF].rhov;
            F.U[cF].rhoE += F.scratch[cF].rhoE;
            EOS::makePhysical(F.U[cF], F.cfg.gamma);
        }
    }

    F.bc.initializeGhosts(F.U, F.mesh);
    
    // ============================================================================
    // CHECK 1 (continued): Compute residual after correction
    // ============================================================================
    
    F.bc.apply(F.U, F.mesh);
    {
        std::vector<Primitive> W(F.U.size());
        for (std::size_t c = 0; c < F.U.size(); ++c) {
            W[c] = Primitive(F.U[c], F.cfg.gamma);
        }
        F.mesh.computeSpectralRadius(W, F.cfg.gamma);
    }
    std::vector<Conservative> R_after(F.mesh.niTotal * F.mesh.njTotal);
    F.flux.computeResidual(F.U, R_after, F.mesh, F.cfg);
    
    // Add forcing if this is a coarse level
    if (fineIdx > 0) {
        int imin, imax, jmin, jmax;
        F.mesh.getInteriorBounds(imin, imax, jmin, jmax);
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int c = F.mesh.cellIndex(i, j);
                R_after[c].rho  += F.forcing[c].rho;
                R_after[c].rhou += F.forcing[c].rhou;
                R_after[c].rhov += F.forcing[c].rhov;
                R_after[c].rhoE += F.forcing[c].rhoE;
            }
        }
    }
    
    double L2_after = 0.0;
    {
        int imin, imax, jmin, jmax;
        F.mesh.getInteriorBounds(imin, imax, jmin, jmax);
        double acc2 = 0.0, volSum = 0.0;
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int c = F.mesh.cellIndex(i, j);
                double V = F.mesh.cellArea[c];
                const Conservative& R = R_after[c];
                acc2 += V * (R.rho*R.rho + R.rhou*R.rhou + R.rhov*R.rhov + R.rhoE*R.rhoE);
                volSum += V;
            }
        }
        L2_after = std::sqrt(acc2 / volSum);
    }
    
    /**double correction_ratio = L2_after / L2_before;
    std::cout << "[CORRECTION CHECK] Level " << fineIdx 
              << ": L2_before = " << std::scientific << L2_before
              << ", L2_after = " << L2_after
              << ", ratio = " << correction_ratio;*/
    


    // ============================================================================
}
