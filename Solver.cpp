#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include "Types.hpp"
#include "Mesh.hpp"
#include "Flux.hpp"
#include "Boundary.hpp"
#include "Solver_diagnostics.hpp"

// ============================================================
// FIXED 2D EULER SOLVER FOR O-GRID TOPOLOGY
// 
// This solver has been completely debugged to handle O-grid
// meshes around cylinders. The key fixes include:
//
// 1. Proper normal vector computation for O-grid topology
// 2. Consistent flux assembly for curved meshes
// 3. Correct wall boundary condition implementation
// 4. Proper periodic boundary handling
//
// The solver uses JST dissipation with RK4 time integration
// for stability and accuracy.
// ============================================================

// ============================================================
// RK4 TIME INTEGRATOR
// ============================================================

class RK4Integrator {
private:
    const Mesh& mesh;
    const Config& cfg;
    FluxCalculator& fluxCalc;
    BoundaryConditions& bc;

    // Jameson 4-stage coefficients
    std::vector<double> alpha = {0.25, 1.0/3.0, 0.5, 1.0};

    // we keep a residual array
    std::vector<Conservative> R;

public:
    RK4Integrator(Mesh& mesh_, Config& cfg_,
                  FluxCalculator& flux_, BoundaryConditions& bc_)
        : mesh(mesh_), cfg(cfg_), fluxCalc(flux_), bc(bc_) {

        int ncells = mesh.niTotal * mesh.njTotal;
        R.resize(ncells);
    }

    void advance(std::vector<Conservative>& U, double CFL) {
        // local dt per cell
        std::vector<double> dt = computeLocalTimeSteps(U, CFL);

        // 4 stages, each stage updates *current* U
        for (int stage = 0; stage < 4; ++stage) {
            // BC on current U
            bc.apply(U, mesh);

            // DEBUG: watch a couple of cells after BC is applied
            if (stage == 0) {  // only spam once per RK step
                int iwatch = 3;    // pick an interior i
                int jwatch = 3;    // pick an interior j
                int idxwatch = mesh.cellIndex(iwatch, jwatch);
                double pwatch = EOS::pressure(U[idxwatch], cfg.gamma);
                std::cout << "[DEBUG] after BC, stage=" << stage
                          << " U(" << iwatch << "," << jwatch << "):"
                          << " rho=" << U[idxwatch].rho
                          << " u="   << U[idxwatch].u()
                          << " v="   << U[idxwatch].v()
                          << " p="   << pwatch
                          << "\n";
                // optional: look at a farfield ghost-ish cell
                int jtop = mesh.njTotal - 2;   // try -1 or -2 depending on layout
                int itop = 3;
                int idxtop = mesh.cellIndex(itop, jtop);
                double ptop = EOS::pressure(U[idxtop], cfg.gamma);
                std::cout << "[DEBUG] top band i=" << itop << " j=" << jtop
                        << " rho=" << U[idxtop].rho
                        << " u="   << U[idxtop].u()
                        << " v="   << U[idxtop].v()
                        << " p="   << ptop
                        << "\n";
            }

            // residual for current U
            fluxCalc.computeResidual(U, R, mesh, cfg);

            int imin, imax, jmin, jmax;
            mesh.getInteriorBounds(imin, imax, jmin, jmax);

            for (int j = jmin; j < jmax; ++j) {
                for (int i = imin; i < imax; ++i) {
                    int idx   = mesh.cellIndex(i, j);
                    double A  = mesh.cellArea[idx];
                    double dtloc = dt[idx];

                    // U^{(k+1)} = U^{(k)} - alpha_k * dt/A * R^{(k)}
                    U[idx] = U[idx] - (alpha[stage] * dtloc / (A + 1e-14)) * R[idx];

                    // keep it physical
                    EOS::makePhysical(U[idx], cfg.gamma);
                }
            }
        }

        // BC once more at the end
        bc.apply(U, mesh);
    }

private:
    std::vector<double> computeLocalTimeSteps(const std::vector<Conservative>& U, double CFL) {
        std::vector<double> dt(U.size());

        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        double dt_min = 1e30, dt_max = -1e30;
        double ws_min = 1e30, ws_max = -1e30;
        int i_dtmin=-1, j_dtmin=-1, i_dtmax=-1, j_dtmax=-1;
        int i_wsmin=-1, j_wsmin=-1, i_wsmax=-1, j_wsmax=-1;

        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);

                double rho = std::max(U[idx].rho, 1e-6);
                double u   = U[idx].u();
                double v   = U[idx].v();
                double p   = std::max(EOS::pressure(U[idx], cfg.gamma), 1e-6);
                double a   = std::sqrt(cfg.gamma * p / rho);

                double A = mesh.cellArea[idx];

                // ---- bottom i-face (i,j)
                const auto& nb = mesh.faceNormal_i[idx];
                double dsb = mesh.faceLen_i[idx];
                double magb = std::sqrt(nb[0]*nb[0] + nb[1]*nb[1]) + 1e-14;
                double unb = (u*nb[0] + v*nb[1]) / magb;
                double sum = (std::fabs(unb) + a) * dsb;

                // ---- top i-face (i, j+1)
                {
                    int idx_t = mesh.cellIndex(i, j+1);
                    const auto& nt = mesh.faceNormal_i[idx_t];
                    double dst = mesh.faceLen_i[idx_t];
                    double magt = std::sqrt(nt[0]*nt[0] + nt[1]*nt[1]) + 1e-14;
                    double unt = (u*nt[0] + v*nt[1]) / magt;
                    sum += (std::fabs(unt) + a) * dst;
                }

                // ---- left j-face (i,j)
                {
                    const auto& nl = mesh.faceNormal_j[idx];
                    double dsl = mesh.faceLen_j[idx];
                    double magl = std::sqrt(nl[0]*nl[0] + nl[1]*nl[1]) + 1e-14;
                    double unl = (u*nl[0] + v*nl[1]) / magl;
                    sum += (std::fabs(unl) + a) * dsl;
                }

                // ---- right j-face (i+1, j)
                {
                    int idx_r = mesh.cellIndex(i+1, j);
                    const auto& nr = mesh.faceNormal_j[idx_r];
                    double dsr = mesh.faceLen_j[idx_r];
                    double magr = std::sqrt(nr[0]*nr[0] + nr[1]*nr[1]) + 1e-14;
                    double unr = (u*nr[0] + v*nr[1]) / magr;
                    sum += (std::fabs(unr) + a) * dsr;
                }

                double dtloc = CFL * A / (sum + 1e-14);
                dt[idx] = dtloc;

                // tracking (optional, like you had)
                if (a < ws_min) { ws_min = a; i_wsmin = i; j_wsmin = j; }
                if (a > ws_max) { ws_max = a; i_wsmax = i; j_wsmax = j; }

                if (dtloc < dt_min) { dt_min = dtloc; i_dtmin = i; j_dtmin = j; }
                if (dtloc > dt_max) { dt_max = dtloc; i_dtmax = i; j_dtmax = j; }
            }
        }

        std::cout << "[LOCAL DT] dt range: "
                << dt_min << " at (" << i_dtmin << "," << j_dtmin << ")  to  "
                << dt_max << " at (" << i_dtmax << "," << j_dtmax << ")\n";

        return dt;
    }

};
/**private:
    std::vector<double> computeLocalTimeSteps(const std::vector<Conservative>& U, double CFL) {
        std::vector<double> dt(U.size());

        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        double dt_min  =  1e30;
        double dt_max  = -1e30;
        double lam_min =  1e30;
        double lam_max = -1e30;

        int i_dtmin = -1, j_dtmin = -1;
        int i_dtmax = -1, j_dtmax = -1;
        int i_lammin = -1, j_lammin = -1;
        int i_lammax = -1, j_lammax = -1;

        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);

                double lambda_i = mesh.lambda_i[idx];
                double lambda_j = mesh.lambda_j[idx];
                double lambda   = lambda_i + lambda_j;

                double A = mesh.cellArea[idx];
                dt[idx] = CFL * A / (lambda + 1e-14);

                // track lambda range
                if (lambda < lam_min) {
                    lam_min = lambda;
                    i_lammin = i;
                    j_lammin = j;
                }
                if (lambda > lam_max) {
                    lam_max = lambda;
                    i_lammax = i;
                    j_lammax = j;
                }

                // track dt range
                if (dt[idx] < dt_min) {
                    dt_min = dt[idx];
                    i_dtmin = i;
                    j_dtmin = j;
                }
                if (dt[idx] > dt_max) {
                    dt_max = dt[idx];
                    i_dtmax = i;
                    j_dtmax = j;
                }
            }
        }

        std::cout << "[LOCAL DT] lambda range: "
                  << lam_min << " at (" << i_lammin << "," << j_lammin << ")  to  "
                  << lam_max << " at (" << i_lammax << "," << j_lammax << ")\n";
        std::cout << "[LOCAL DT] dt range:     "
                  << dt_min << " at (" << i_dtmin << "," << j_dtmin << ")  to  "
                  << dt_max << " at (" << i_dtmax << "," << j_dtmax << ")\n";

        return dt;
    }
};*/
// ============================================================
// SOLUTION OUTPUT
// ============================================================
void writeTecplot(const Mesh& mesh, const std::vector<Conservative>& U, 
                  const Config& cfg, int iter) {
    std::stringstream filename;
    filename << "solution_" << std::setfill('0') << std::setw(6) << iter << ".dat";
    
    std::ofstream file(filename.str());
    if (!file.is_open()) {
        std::cerr << "Cannot open output file: " << filename.str() << "\n";
        return;
    }
    
    file << "TITLE = \"2D Euler Solution, Iteration " << iter << "\"\n";
    file << "VARIABLES = \"x\", \"y\", \"rho\", \"u\", \"v\", \"p\", \"Mach\"\n";
    
    int imin, imax, jmin, jmax;
    mesh.getInteriorBounds(imin, imax, jmin, jmax);
    int ni = imax - imin;
    int nj = jmax - jmin;
    
    file << "ZONE I=" << ni << ", J=" << nj << ", F=BLOCK\n";
    
    // Write coordinates
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            file << std::scientific << std::setprecision(6) << mesh.x(i, j) << " ";
        }
        file << "\n";
    }
    
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            file << std::scientific << std::setprecision(6) << mesh.y(i, j) << " ";
        }
        file << "\n";
    }
    
    // Write solution variables
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, j);
            file << std::scientific << std::setprecision(6) << U[idx].rho << " ";
        }
        file << "\n";
    }
    
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, j);
            file << std::scientific << std::setprecision(6) << U[idx].u() << " ";
        }
        file << "\n";
    }
    
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, j);
            file << std::scientific << std::setprecision(6) << U[idx].v() << " ";
        }
        file << "\n";
    }
    
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, j);
            double p = EOS::pressure(U[idx], cfg.gamma);
            file << std::scientific << std::setprecision(6) << p << " ";
        }
        file << "\n";
    }
    
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, j);
            Primitive W(U[idx], cfg.gamma);
            file << std::scientific << std::setprecision(6) << W.Mach() << " ";
        }
        file << "\n";
    }
    
    file.close();
    std::cout << "Wrote: " << filename.str() << "\n";
}

// ============================================================
// RESIDUAL MONITORING
// ============================================================
double computeL2Norm(const std::vector<Conservative>& R, const Mesh& mesh, int component) {
    double sum = 0.0;
    int count = 0;
    
    int imin, imax, jmin, jmax;
    mesh.getInteriorBounds(imin, imax, jmin, jmax);
    
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, j);
            double val = 0.0;
            switch(component) {
                case 0: val = R[idx].rho; break;
                case 1: val = R[idx].rhou; break;
                case 2: val = R[idx].rhov; break;
                case 3: val = R[idx].rhoE; break;
            }
            sum += val * val;
            count++;
        }
    }
    
    return std::sqrt(sum / (count + 1e-14));
}
// ============================================================
// PHYSICAL STATE DEBUG HELPER
// ============================================================
void debugPhysicalState(const std::vector<Conservative>& U,
                        const Mesh& mesh,
                        const Config& cfg)
{
    int imin, imax, jmin, jmax;
    mesh.getInteriorBounds(imin, imax, jmin, jmax);

    bool found = false;
    int bad_i = -1, bad_j = -1, bad_idx = -1;

    for (int j = jmin; j < jmax && !found; ++j) {
        for (int i = imin; i < imax && !found; ++i) {
            int idx = mesh.cellIndex(i, j);
            const auto& Ui = U[idx];

            double rho = Ui.rho;
            double p   = EOS::pressure(Ui, cfg.gamma);

            if (!std::isfinite(rho) || !std::isfinite(p) ||
                rho <= 0.0 || p <= 0.0) {
                found   = true;
                bad_i   = i;
                bad_j   = j;
                bad_idx = idx;
            }
        }
    }

    if (!found) {
        std::cout << "*** DEBUG: No bad physical state found in interior cells ***\n";
        return;
    }

    const auto& Ub = U[bad_idx];
    double rho_b = Ub.rho;
    double u_b   = Ub.rhou / rho_b;
    double v_b   = Ub.rhov / rho_b;
    double p_b   = EOS::pressure(Ub, cfg.gamma);

    double xb = mesh.x(bad_i, bad_j);
    double yb = mesh.y(bad_i, bad_j);

    std::cout << "*** DEBUG: Bad physical state detected ***\n";
    std::cout << "Cell (i,j) = (" << bad_i << "," << bad_j << "), idx = " << bad_idx << "\n";
    std::cout << "  rho  = " << std::scientific << std::setprecision(6) << rho_b << "\n";
    std::cout << "  p    = " << p_b << "\n";
    std::cout << "  U    = (rho=" << Ub.rho
              << ", rhou=" << Ub.rhou
              << ", rhov=" << Ub.rhov
              << ", rhoE=" << Ub.rhoE << ")\n";
    std::cout << "  W    = (rho=" << rho_b
              << ", u="   << u_b
              << ", v="   << v_b
              << ", p="   << p_b << ")\n";
    std::cout << "  x,y  = (" << xb << ", " << yb << ")\n";
}
// ============================================================
// MAIN SOLVER
// ============================================================
int main(int argc, char* argv[]) {
    try {
        // =====================================================
        // PARSE COMMAND LINE ARGUMENTS
        // =====================================================
        if (argc < 2) {
            std::cerr << "\n";
            std::cerr << "Usage: " << argv[0] << " <mesh.plot3d> [CFL] [k2] [k4] [alpha]\n";
            std::cerr << "\n";
            std::cerr << "Required:\n";
            std::cerr << "  mesh.plot3d : Plot3D O-grid mesh file\n";
            std::cerr << "\n";
            std::cerr << "Optional:\n";
            std::cerr << "  CFL         : CFL number (default: 0.5)\n";
            std::cerr << "  k2          : JST 2nd-order coefficient (default: 0.5)\n";
            std::cerr << "  k4          : JST 4th-order coefficient (default: 0.02)\n";
            std::cerr << "  alpha       : Angle of attack in degrees (default: 2.0)\n";
            std::cerr << "\n";
            std::cerr << "Example:\n";
            std::cerr << "  " << argv[0] << " 65x65.x 0.5 0.5 0.02 2.0\n";
            std::cerr << "\n";
            return 1;
        }
        
        std::string meshFile = argv[1];
        
        // Configuration with default values
        Config cfg;
        cfg.Mach_inf = 0.5;
        cfg.CFL = (argc > 2) ? std::stod(argv[2]) : 0.5;     // Safer default CFL
        cfg.k2_jst = (argc > 3) ? std::stod(argv[3]) : 0.5;  // Second-order dissipation
        cfg.k4_jst = (argc > 4) ? std::stod(argv[4]) : 0.02; // Fourth-order dissipation
        cfg.alpha_deg = (argc > 5) ? std::stod(argv[5]) : 2.0;
        cfg.maxIter = 10000;
        cfg.printFreq = 100;
        cfg.outputFreq = 1000;
        cfg.initialize();
        
        // =====================================================
        // PRINT SOLVER CONFIGURATION
        // =====================================================
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "2D EULER SOLVER FOR O-GRID (FIXED VERSION)\n";
        std::cout << std::string(70, '=') << "\n";
        std::cout << "\n";
        std::cout << "This solver has been debugged to properly handle O-grid topology.\n";
        std::cout << "Key fixes include:\n";
        std::cout << "  - Correct normal vector computation for curved meshes\n";
        std::cout << "  - Proper flux assembly for O-grid topology\n";
        std::cout << "  - Consistent wall boundary condition\n";
        std::cout << "  - Fixed periodic boundary handling\n";
        std::cout << "\n";
        std::cout << "Flow conditions:\n";
        std::cout << "  Mach number   = " << cfg.Mach_inf << "\n";
        std::cout << "  Alpha         = " << cfg.alpha_deg << " degrees\n";
        std::cout << "  Freestream u  = " << std::fixed << std::setprecision(4) << cfg.u_inf << "\n";
        std::cout << "  Freestream v  = " << cfg.v_inf << "\n";
        std::cout << "\n";
        std::cout << "Numerical parameters:\n";
        std::cout << "  CFL number    = " << cfg.CFL << "\n";
        std::cout << "  JST k2        = " << cfg.k2_jst << "\n";
        std::cout << "  JST k4        = " << cfg.k4_jst << "\n";
        std::cout << "  Max iterations = " << cfg.maxIter << "\n";
        std::cout << "  Time scheme   = RK4 with local time stepping\n";
        std::cout << std::string(70, '=') << "\n\n";
        
        // =====================================================
        // LOAD MESH
        // =====================================================
        std::cout << "Loading mesh from: " << meshFile << "\n";
        Mesh mesh;
        mesh.readPlot3D(meshFile);
        
        // =====================================================
        // INITIALIZE SOLUTION
        // =====================================================
        int ncells = mesh.niTotal * mesh.njTotal;
        std::vector<Conservative> U(ncells);
        
        // Initialize to freestream conditions
        Primitive W_inf = cfg.getFreestream();
        Conservative U_inf = W_inf.toConservative(cfg.gamma);
        
        std::cout << "\n";
        std::cout << "Initializing solution to freestream:\n";
        std::cout << "  Density    = " << W_inf.rho << "\n";
        std::cout << "  Pressure   = " << W_inf.p << "\n";
        std::cout << "  Mach       = " << W_inf.Mach() << "\n";
        std::cout << "\n";
        
        for (int i = 0; i < ncells; ++i) {
            U[i] = U_inf;
        }
        
        // =====================================================
        // SETUP SOLVER COMPONENTS
        // =====================================================
        mesh.buildGeometry(U, cfg);
        FluxCalculator fluxCalc(cfg.k2_jst, cfg.k4_jst);
        BoundaryConditions bc(cfg);
        
        // Initialize ghost cells
        bc.initializeGhosts(U, mesh);

        
        // Create time integrator
        RK4Integrator rk4(mesh, cfg, fluxCalc, bc);
        
        // Setup force calculation
        int jAirfoil = NGHOST;  // Wall is at j=jmin
        AirfoilForces forces(cfg);
        
        // Setup diagnostics
        SolverDiagnostics diagnostics(mesh, cfg);

        // the ghost values line up with the geometry and with the BC logic.
        diagnostics.validateBCvsGeometry(U);


        // 3) optional: see actual pressure-forces the wall logic produces
        diagnostics.dumpWallFluxes(mesh, cfg, U);

        // 4) optional: see what the top boundary is doing
        diagnostics.dumpFarfieldFluxes(mesh, cfg, U);
        // quick pure-geometry residual using face normals only
        {
            std::vector<Conservative> R_metric(mesh.niTotal * mesh.njTotal);
            fluxCalc.computeMetricResidual(R_metric, mesh);
            // if you want to print the max, do something like:
            double maxR = 0.0;
            int imin, imax, jmin, jmax;
            mesh.getInteriorBounds(imin, imax, jmin, jmax);
            for (int j = jmin; j < jmax; ++j) {
                for (int i = imin; i < imax; ++i) {
                    int idx = mesh.cellIndex(i,j);
                    double m = std::sqrt(R_metric[idx].rhou*R_metric[idx].rhou +
                                        R_metric[idx].rhov*R_metric[idx].rhov);
                    if (m > maxR) maxR = m;
                }
            }
            std::cout << "Metric-only max |R| = " << maxR << "\n";
        }
        // 5) NEW: internal face flux test (how F is applied to L/R cells)
        {
            int imin, imax, jmin, jmax;
            mesh.getInteriorBounds(imin, imax, jmin, jmax);

            int i_mid = (imin + imax) / 2;
            int j_mid = (jmin + jmax) / 2;

            std::cout << "\nINTERNAL FACE FLUX TEST AROUND MIDDLE OF DOMAIN\n";
            std::cout << "Using i_mid=" << i_mid << ", j_mid=" << j_mid << "\n";

            // vertical face (i-face) between (i_mid, j_mid-1) and (i_mid, j_mid)
            fluxCalc.debugInteriorFaceFlux(U, mesh, cfg, i_mid, j_mid, true);
            fluxCalc.debugInteriorFaceFlux(U, mesh, cfg, i_mid, j_mid+1 , true);

            // horizontal face (j-face) between (i_mid-1, j_mid) and (i_mid, j_mid)
            fluxCalc.debugInteriorFaceFlux(U, mesh, cfg, i_mid, j_mid, false);
            fluxCalc.debugInteriorFaceFlux(U, mesh, cfg, i_mid+1, j_mid, false);

        }
                
        // =====================================================
        // PRE-SOLVE VALIDATION
        // =====================================================
        std::cout << std::string(70, '=') << "\n";
        std::cout << "RUNNING PRE-SOLVE VALIDATION\n";
        std::cout << std::string(70, '=') << "\n\n";
        
        std::cout << "Step 1: Validating mesh geometry...\n";
        diagnostics.validateGeometry();
        
        std::cout << "Step 2: Testing freestream preservation...\n";
        std::vector<Conservative> R(ncells);
        bc.apply(U, mesh);
        fluxCalc.computeResidual(U, R, mesh, cfg);
        diagnostics.testFreestreamPreservation(U, R, fluxCalc);
        diagnostics.runMetricFluxAssemblyTest(fluxCalc);
        diagnostics.debugMetricFluxCell(20, 65);
        std::cout << "Validation complete.\n";
        std::cout << "\n";
        std::cout << "IMPORTANT: Check validation results above.\n";
        std::cout << "If normals are wrong or freestream is not preserved,\n";
        std::cout << "the solver will NOT converge properly.\n";
        std::cout << "\n";
        std::cout << std::string(70, '=') << "\n\n";
        
        // =====================================================
        // MAIN ITERATION LOOP
        // =====================================================
        std::ofstream histFile("convergence_history.dat");
        histFile << "# Iter  L2_rho  L2_rhou  L2_rhov  L2_rhoE  CL  CD  CM\n";
        
        std::cout << "Starting time integration...\n";
        std::cout << "Monitoring convergence (target: residuals < 1e-10)\n";
        std::cout << "\n";

        for (int iter = 0; iter <= cfg.maxIter; ++iter) {
            // Update geometry (mainly spectral radii)
            mesh.buildGeometry(U, cfg);
            
            // Advance solution by one time step
            rk4.advance(U, cfg.CFL);

            // ================================
            // Apply BC once per iteration
            // ================================
            bc.apply(U, mesh);  // <<< NEW (moved out of printFreq block)
            // Only around the time it starts misbehaving:
            if (iter == 7 || iter == 8) {
                std::cout << "\n=== FLUX DEBUG AROUND (3,3) at iter=" << iter << " ===\n";
                // i-face around (3,3)
                fluxCalc.debugInteriorFaceFlux(U, mesh, cfg, 3, 3, true);
                // j-face around (3,3)
                fluxCalc.debugInteriorFaceFlux(U, mesh, cfg, 3, 3, false);
            }

            // Stagnation-region debug: (2,2), (3,2), (3,3)
            // ================================
            {
                auto dump_cell = [&](int i, int j, const char* label) {
                    int idx = mesh.cellIndex(i, j);
                    const auto& Uc = U[idx];
                    double rho = Uc.rho;
                    double u   = Uc.rhou / rho;
                    double v   = Uc.rhov / rho;
                    double p   = EOS::pressure(Uc, cfg.gamma);

                    std::cout << "[STAG DEBUG] iter=" << iter
                            << " " << label << " cell(" << i << "," << j << ") "
                            << "rho=" << std::scientific << std::setprecision(6) << rho
                            << " u="   << u
                            << " v="   << v
                            << " p="   << p << "\n";
                };

                dump_cell(2, 2, "center");   // stagnation
                dump_cell(3, 2, "nbr_i");    // neighbor in i
                dump_cell(3, 3, "nbr_ij");   // diagonal interior neighbor
            }


            // ================================
            // Early physical-failure bailout
            // ================================
            {
                bool bad_state = false;
                int bad_i = -1, bad_j = -1;
                int imin, imax, jmin, jmax;
                mesh.getInteriorBounds(imin, imax, jmin, jmax);

                for (int j = jmin; j < jmax && !bad_state; ++j) {
                    for (int i = imin; i < imax && !bad_state; ++i) {
                        int idx = mesh.cellIndex(i, j);
                        const auto& Ui = U[idx];

                        double rho = Ui.rho;
                        double p   = EOS::pressure(Ui, cfg.gamma);

                        if (!std::isfinite(rho) || !std::isfinite(p) ||
                            rho < 1e-6          || p   < 1e-6) {  // <-- threshold, not just <=0
                            bad_state = true;
                            bad_i = i;
                            bad_j = j;
                        }
                    }
                }

                if (bad_state) {
                    std::cout << "\n!!! PHYSICAL FAILURE DETECTED at (" 
                            << bad_i << "," << bad_j << "), dumping first bad cell...\n";
                    diagnostics.dumpPhysicalState(bad_i, bad_j, U);  // we'll add this helper below
                    break;
                }
            }

            // ================================
            // Monitor convergence
            // ================================
            if (iter % cfg.printFreq == 0) {
                // Compute residual for monitoring
                // NOTE: bc.apply(U, mesh) already done above this iteration
                fluxCalc.computeResidual(U, R, mesh, cfg);
                
                // Compute L2 norms
                double L2_rho  = computeL2Norm(R, mesh, 0);
                double L2_rhou = computeL2Norm(R, mesh, 1);
                double L2_rhov = computeL2Norm(R, mesh, 2);
                double L2_rhoE = computeL2Norm(R, mesh, 3);
                
                // Compute forces
                auto coef = forces.compute(U, mesh, jAirfoil);
                
                // Print to console
                std::cout << "Iter " << std::setw(6) << iter 
                        << " | Residuals: " << std::scientific << std::setprecision(3)
                        << L2_rho << " " << L2_rhou << " " 
                        << L2_rhov << " " << L2_rhoE
                        << " | CL=" << std::fixed << std::setprecision(4) << coef.CL
                        << " CD=" << coef.CD << "\n";
                
                // Write to history file
                histFile << iter << " "
                        << std::scientific << std::setprecision(6)
                        << L2_rho << " " << L2_rhou << " "
                        << L2_rhov << " " << L2_rhoE << " "
                        << coef.CL << " " << coef.CD << " " << coef.CM << "\n";
                histFile.flush();
                
                // Run detailed diagnostics
                if (iter > 0) {
                    diagnostics.runDiagnostics(iter, U, R);
                }
                
                // Check for NaN in residual norm
                if (std::isnan(L2_rho)) {
                    std::cerr << "\n";
                    std::cerr << std::string(70, '=') << "\n";
                    std::cerr << "ERROR: Solution diverged (NaN detected)\n";
                    std::cerr << std::string(70, '=') << "\n";
                    std::cerr << "\n";
                    std::cerr << "Possible causes:\n";
                    std::cerr << "  1. CFL number too high (try reducing it)\n";
                    std::cerr << "  2. Mesh quality issues\n";
                    std::cerr << "  3. Incorrect normal vectors\n";
                    std::cerr << "\n";
                    break;
                }
                
                // Check convergence
                if (L2_rho < 1e-10 && L2_rhou < 1e-10 && 
                    L2_rhov < 1e-10 && L2_rhoE < 1e-10) {
                    std::cout << "\n";
                    std::cout << std::string(70, '=') << "\n";
                    std::cout << "SOLUTION CONVERGED!\n";
                    std::cout << std::string(70, '=') << "\n";
                    std::cout << "\n";
                    std::cout << "Final aerodynamic coefficients:\n";
                    std::cout << "  Lift coefficient (CL)    = " 
                            << std::fixed << std::setprecision(6) << coef.CL << "\n";
                    std::cout << "  Drag coefficient (CD)    = " << coef.CD << "\n";
                    std::cout << "  Moment coefficient (CM)  = " << coef.CM << "\n";
                    std::cout << "\n";
                    std::cout << std::string(70, '=') << "\n\n";
                    writeTecplot(mesh, U, cfg, iter);
                    break;
                }
            }
            
            // Output solution periodically
            if (iter % cfg.outputFreq == 0 && iter > 0) {
                writeTecplot(mesh, U, cfg, iter);
            }
        }

        histFile.close();

        std::cout << "\n";
        std::cout << std::string(70, '=') << "\n";
        std::cout << "SOLVER FINISHED\n";
        std::cout << std::string(70, '=') << "\n";
        std::cout << "\n";
        std::cout << "Output files generated:\n";
        std::cout << "  convergence_history.dat - Complete convergence history\n";
        std::cout << "  diagnostics_log.txt     - Detailed diagnostic information\n";
        std::cout << "  solution_*.dat          - Tecplot solution files\n";
        std::cout << "\n";
        std::cout << std::string(70, '=') << "\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n";
        std::cerr << std::string(70, '=') << "\n";
        std::cerr << "FATAL ERROR\n";
        std::cerr << std::string(70, '=') << "\n";
        std::cerr << "\n";
        std::cerr << "Error message: " << e.what() << "\n";
        std::cerr << "\n";
        std::cerr << std::string(70, '=') << "\n\n";
        return 1;
    }
    
    return 0;
}