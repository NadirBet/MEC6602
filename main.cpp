// main.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <limits>
#include <chrono>

#include "ImplicitResidualSmoothing.hpp"
#include "mesh.hpp"
#include "Types.hpp"
#include "Flux.hpp"
#include "Boundary.hpp"
#include "FASMultigrid.hpp"  // multigrid

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

// ============================================================================
//                                Small helpers
// ============================================================================

static inline std::string ltrim(std::string s){size_t p=s.find_first_not_of(" \t\r\n");return(p==std::string::npos)?"":s.substr(p);}
static inline std::string rtrim(std::string s){size_t p=s.find_last_not_of(" \t\r\n"); return(p==std::string::npos)?"":s.substr(0,p+1);}
static inline std::string trim(std::string s){return rtrim(ltrim(std::move(s)));}
static inline std::string lower(std::string s){for(char& c:s)c=(char)std::tolower((unsigned char)c);return s;}
static inline std::string strip_dots(std::string s){s.erase(std::remove(s.begin(),s.end(),'.'),s.end());return s;}
static inline bool parse_bool(const std::string& v){std::string t=lower(trim(v));return (t=="1"||t=="true"||t=="yes"||t=="on");}
inline int cidx(int i,int j,int niTotal){return j*niTotal+i;}

// ----- runtime options (from config.txt) -----
struct RunOpts {
    std::string meshFile = "mesh.x";
    bool use_simple_farfield = false;
    bool freestream_everywhere = false;

    // (legacy) list of meshes for multigrid; now unused since MG builds by 2x coarsening
    std::vector<std::string> mg_mesh_files;
    
    // perturbation stuff...
    std::string perturb_mode = "none";
    double perturb_amp = 1e-3;
    int    perturb_waves = 1;
    double perturb_sigma_cells = 6;
};

// --------------------------------------------------------------------
// Write solution on physical cells as VTK UNSTRUCTURED_GRID (QUADS)
// --------------------------------------------------------------------
void writeSolutionVTKPhysical(
    const std::string& filename,
    const Mesh& mesh,
    const std::vector<Conservative>& U,
    const Config& cfg)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[VTK] ERROR: cannot open " << filename << "\n";
        return;
    }

    int imin, imax, jmin, jmax;
    mesh.getInteriorBounds(imin, imax, jmin, jmax);

    const int nnodes     = mesh.niNodes * mesh.njNodes;
    const int nPhysCells = (imax - imin) * (jmax - jmin);

    // header
    out << "# vtk DataFile Version 3.0\n";
    out << "solution dump (physical cells)\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // points
    out << "POINTS " << nnodes << " double\n";
    for (int j = 0; j < mesh.njNodes; ++j) {
        for (int i = 0; i < mesh.niNodes; ++i) {
            int n = mesh.nodeIndex(i, j);
            out << mesh.xNodes[n] << " " << mesh.yNodes[n] << " 0.0\n";
        }
    }

    // cells
    out << "CELLS " << nPhysCells << " " << nPhysCells * 5 << "\n";
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int n0 = mesh.nodeIndex(i,     j);
            int n1 = mesh.nodeIndex(i + 1, j);
            int n2 = mesh.nodeIndex(i + 1, j + 1);
            int n3 = mesh.nodeIndex(i,     j + 1);
            out << 4 << " " << n0 << " " << n1 << " " << n2 << " " << n3 << "\n";
        }
    }

    out << "CELL_TYPES " << nPhysCells << "\n";
    for (int c = 0; c < nPhysCells; ++c) out << 9 << "\n";  // VTK_QUAD

    // cell data
    out << "CELL_DATA " << nPhysCells << "\n";

    // Pressure
    out << "SCALARS p double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int c = mesh.cellIndex(i, j);
            Primitive W(U[c], cfg.gamma);
            out << W.p << "\n";
        }
    }

    // Density
    out << "SCALARS rho double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int c = mesh.cellIndex(i, j);
            out << U[c].rho << "\n";
        }
    }

    // Mach
    out << "SCALARS Mach double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int c = mesh.cellIndex(i, j);
            Primitive W(U[c], cfg.gamma);
            const double a = std::sqrt(std::max(0.0, cfg.gamma * W.p / W.rho));
            const double V = std::sqrt(W.u*W.u + W.v*W.v);
            const double M = (a > 0.0) ? (V / a) : 0.0;
            out << M << "\n";
        }
    }

    // Velocity vector
    out << "VECTORS velocity double\n";
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int c = mesh.cellIndex(i, j);
            Primitive W(U[c], cfg.gamma);
            out << W.u << " " << W.v << " 0.0\n";
        }
    }
    // Pressure coefficient
    out << "SCALARS Cp double 1\n";
    out << "LOOKUP_TABLE default\n";
    const double q_inf = 0.5 * cfg.rho_inf * 
                        (cfg.u_inf * cfg.u_inf + cfg.v_inf * cfg.v_inf);
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int c = mesh.cellIndex(i, j);
            Primitive W(U[c], cfg.gamma);
            const double Cp = (W.p - cfg.p_inf) / std::max(1e-14, q_inf);
            out << Cp << "\n";
        }
    }

    out.close();
    std::cout << "[VTK] wrote " << filename << "\n";
}

// ============================================================================
//                         Diagnostics & numerics helpers
// ============================================================================

struct ResSplit { double interior=0.0, bottom=0.0, top=0.0, periodic=0.0; };

ResSplit residualSplit(const std::vector<Conservative>& R, const Mesh& m)
{
    int imin,imax,jmin,jmax; m.getInteriorBounds(imin,imax,jmin,jmax);
    const int niTot = m.niTotal;
    ResSplit out;
    for(int j=jmin;j<jmax;++j){
        for(int i=imin;i<imax;++i){
            int id = j*niTot + i;
            double rcell = 0.0;
            rcell = std::max(rcell, std::fabs(R[id].rho ));
            rcell = std::max(rcell, std::fabs(R[id].rhou));
            rcell = std::max(rcell, std::fabs(R[id].rhov));
            rcell = std::max(rcell, std::fabs(R[id].rhoE));
            bool isBottom=(j==jmin), isTop=(j==jmax-1), isLeftPer=(i==imin), isRightPer=(i==imax-1);
            if(isLeftPer||isRightPer) out.periodic = std::max(out.periodic, rcell);
            if(isBottom) out.bottom = std::max(out.bottom, rcell);
            if(isTop)    out.top    = std::max(out.top,    rcell);
            if(!isBottom && !isTop && !isLeftPer && !isRightPer) out.interior = std::max(out.interior, rcell);
        }
    }
    return out;
}

static void lambda_stats(const Mesh& m, double& Imin,double& Imax,double& Iavg,
                         double& Jmin,double& Jmax,double& Javg)
{
    int imin,imax,jmin,jmax; m.getInteriorBounds(imin,imax,jmin,jmax);
    Imin=+1e300; Imax=-1e300; Iavg=0.0;
    Jmin=+1e300; Jmax=-1e300; Javg=0.0;
    size_t count=0;
    for(int j=jmin;j<jmax;++j){
        for(int i=imin;i<imax;++i){
            int c=m.cellIndex(i,j);
            Imin = std::min(Imin, m.lambdaI[c]); Imax = std::max(Imax, m.lambdaI[c]); Iavg += m.lambdaI[c];
            Jmin = std::min(Jmin, m.lambdaJ[c]); Jmax = std::max(Jmax, m.lambdaJ[c]); Javg += m.lambdaJ[c];
            ++count;
        }
    }
    if(count){ Iavg/=count; Javg/=count; }
}

// Build per-cell local timesteps: dt_c = CFL * Area / (lambdaI+lambdaJ)
void buildLocalDt(const Mesh& m, const Config& cfg, std::vector<double>& dt)
{
    dt.assign(m.niTotal * m.njTotal, 0.0);
    int imin,imax,jmin,jmax; m.getInteriorBounds(imin,imax,jmin,jmax);
    const double EPS=1e-14;
    for(int j=jmin;j<jmax;++j){
        for(int i=imin;i<imax;++i){
            const int c=m.cellIndex(i,j);
            const double denom = m.lambdaI[c] + m.lambdaJ[c];
            const double Omega = m.cellArea[c];
            dt[c] = cfg.CFL * Omega / std::max(denom, EPS);
        }
    }
}

// 5-stage Jameson RK with LOCAL dt per cell (steady pseudo-time)
void advanceExplicitRK5_LTS(Mesh& m, const Config& cfg, BoundaryConditions& bc,
                            FluxCalculator& flux, std::vector<Conservative>& U,
                            std::vector<Conservative>& R, const std::vector<double>& dt)
{
    static const double alpha[5] = { 0.25, 1.0/6.0, 0.375, 0.5, 1.0 };
    int imin,imax,jmin,jmax; 
    m.getInteriorBounds(imin,imax,jmin,jmax);

    // Freeze base state U0 = W^(0) = U^n for this outer iteration
    std::vector<Conservative> U0 = U;

    // residual smoothing buffer + object
    std::vector<Conservative> R_smooth;
    if (cfg.use_residual_smoothing) {
        R_smooth.resize(U.size());
    }
    ImplicitResidualSmoothing smoother(cfg.smooth_eps_I, cfg.smooth_eps_J);

    for(int s=0; s<5; ++s){
        // Apply BCs to current stage state U
        bc.apply(U, m);

        // Stage spectral radii for JST & diagnostics (based on U^((s)))
        std::vector<Primitive> W(U.size());
        for(size_t c=0; c<U.size();++c)
            W[c] = Primitive(U[c], cfg.gamma);
        m.computeSpectralRadius(W, cfg.gamma);

        // Residual at current stage state U (this is R^(s-1))
        flux.computeResidual(U, R, m, cfg);

        // Optional: apply implicit residual smoothing
        if (cfg.use_residual_smoothing) {
            smoother.smooth(m, R, R_smooth);
        }

        // Stage update: U <- U0 - α_s * (dt/Ω) * R^(s-1)
        for(int j=jmin; j<jmax; ++j){
            for(int i=imin; i<imax; ++i){
                const int c = m.cellIndex(i,j);
                const double factor = alpha[s] * dt[c] / m.cellArea[c];

                // Choose which residual to use
                const Conservative& Rc = cfg.use_residual_smoothing ? R_smooth[c] : R[c];

                U[c].rho  = U0[c].rho  + factor * Rc.rho;
                U[c].rhou = U0[c].rhou + factor * Rc.rhou;
                U[c].rhov = U0[c].rhov + factor * Rc.rhov;
                U[c].rhoE = U0[c].rhoE + factor * Rc.rhoE;
                EOS::makePhysical(U[c], cfg.gamma);
            }
        }
    }

    // After the last stage, ensure R is the residual of the FINAL U
    bc.apply(U, m);

    std::vector<Primitive> W_final(U.size());
    for (size_t c = 0; c < U.size(); ++c)
        W_final[c] = Primitive(U[c], cfg.gamma);
    m.computeSpectralRadius(W_final, cfg.gamma);

    flux.computeResidual(U, R, m, cfg);
}

void advanceExplicitRK5_LTS_MG(
    Mesh& m, 
    const Config& cfg, 
    BoundaryConditions& bc,
    FluxCalculator& flux, 
    std::vector<Conservative>& U,
    std::vector<Conservative>& R, 
    const std::vector<double>& dt,
    const std::vector<Conservative>* forcing = nullptr)  // FAS forcing (nullptr on finest grid)
{
    static const double alpha[5] = { 0.25, 1.0/6.0, 0.375, 0.5, 1.0 };
    int imin, imax, jmin, jmax; 
    m.getInteriorBounds(imin, imax, jmin, jmax);

    // Freeze the base state at the start of this time step
    std::vector<Conservative> U0 = U;

    // Residual smoothing setup (if enabled)
    std::vector<Conservative> R_smooth;
    if (cfg.use_residual_smoothing) {
        R_smooth.resize(U.size());
    }
    ImplicitResidualSmoothing smoother(cfg.smooth_eps_I, cfg.smooth_eps_J);

    // Five-stage Runge-Kutta integration
    for (int s = 0; s < 5; ++s) {
        // Apply boundary conditions to the current stage solution
        bc.apply(U, m);

        // Compute primitive variables and update spectral radii for this stage
        std::vector<Primitive> W(U.size());
        for (size_t c = 0; c < U.size(); ++c)
            W[c] = Primitive(U[c], cfg.gamma);
        m.computeSpectralRadius(W, cfg.gamma);

        // Compute the pure PDE residual for the current stage
        flux.computeResidual(U, R, m, cfg);

        // CRITICAL: Add FAS forcing to residual if we're on a coarse grid
        // This ensures the time integrator sees the complete forced equation
        if (forcing != nullptr) {
            for (int j = jmin; j < jmax; ++j) {
                for (int i = imin; i < imax; ++i) {
                    int c = m.cellIndex(i, j);
                    R[c].rho  += (*forcing)[c].rho;
                    R[c].rhou += (*forcing)[c].rhou;
                    R[c].rhov += (*forcing)[c].rhov;
                    R[c].rhoE += (*forcing)[c].rhoE;
                }
            }
        }

        // Apply implicit residual smoothing if requested
        if (cfg.use_residual_smoothing) {
            smoother.smooth(m, R, R_smooth);
        }

        // Perform the stage update: U^(s) = U^(0) + α_s * (dt/Ω) * R^(s-1)
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                const int c = m.cellIndex(i, j);
                const double factor = alpha[s] * dt[c] / m.cellArea[c];

                // Choose smoothed or unsmoothed residual
                const Conservative& Rc = cfg.use_residual_smoothing ? R_smooth[c] : R[c];

                U[c].rho  = U0[c].rho  + factor * Rc.rho;
                U[c].rhou = U0[c].rhou + factor * Rc.rhou;
                U[c].rhov = U0[c].rhov + factor * Rc.rhov;
                U[c].rhoE = U0[c].rhoE + factor * Rc.rhoE;
                
                // Enforce physical realizability (positive density and pressure)
                EOS::makePhysical(U[c], cfg.gamma);
            }
        }
    }

    // After all stages complete, compute the final residual
    // This gives us R(U^final) which we'll use for convergence monitoring
    bc.apply(U, m);

    std::vector<Primitive> W_final(U.size());
    for (size_t c = 0; c < U.size(); ++c)
        W_final[c] = Primitive(U[c], cfg.gamma);
    m.computeSpectralRadius(W_final, cfg.gamma);

    flux.computeResidual(U, R, m, cfg);
    
    // Add forcing one more time to R for consistency in convergence checking
    if (forcing != nullptr) {
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int c = m.cellIndex(i, j);
                R[c].rho  += (*forcing)[c].rho;
                R[c].rhou += (*forcing)[c].rhou;
                R[c].rhov += (*forcing)[c].rhov;
                R[c].rhoE += (*forcing)[c].rhoE;
            }
        }
    }
}

// global Linf of residual over physical cells
double residualLinf(const std::vector<Conservative>& R, const Mesh& m)
{
    int imin,imax,jmin,jmax; m.getInteriorBounds(imin,imax,jmin,jmax);
    const int niTot=m.niTotal; double maxv=0.0;
    for(int j=jmin;j<jmax;++j)
        for(int i=imin;i<imax;++i){
            int id=j*niTot+i;
            maxv=std::max(maxv,std::fabs(R[id].rho ));
            maxv=std::max(maxv,std::fabs(R[id].rhou));
            maxv=std::max(maxv,std::fabs(R[id].rhov));
            maxv=std::max(maxv,std::fabs(R[id].rhoE));
        }
    return maxv;
}

// ============================================================
// AIRFOIL FORCE / MOMENT COEFFICIENTS (2D O-grid, bottom wall)
// ============================================================

// ============================================================
// AIRFOIL FORCE / MOMENT COEFFICIENTS (2D O-grid, bottom wall)
// ============================================================
struct AeroCoeffs {
    double CL = 0.0;
    double CD = 0.0;
    double CM = 0.0;
};

// Find quarter-chord reference point AND geometric chord c_geom
inline void findQuarterChordRef(
    const Mesh& mesh,
    int jWall,
    double& x_ref,
    double& y_ref,
    double& c_geom   // <-- geometric chord length (x_max - x_min)
)
{
    int imin, imax, jmin, jmax;
    mesh.getInteriorBounds(imin, imax, jmin, jmax);

    auto nidx = [&](int I, int J){ return mesh.nodeIndex(I, J); };

    double x_min =  1e300;
    double x_max = -1e300;

    const int jFaceRow = jWall;

    std::vector<double> x_face, y_face;
    x_face.reserve(imax - imin);
    y_face.reserve(imax - imin);

    for (int i = imin; i < imax; ++i) {
        int nA = nidx(i,   jFaceRow);
        int nB = nidx(i+1, jFaceRow);

        const double xA = mesh.xNodes[nA];
        const double yA = mesh.yNodes[nA];
        const double xB = mesh.xNodes[nB];
        const double yB = mesh.yNodes[nB];

        const double xm = 0.5 * (xA + xB);
        const double ym = 0.5 * (yA + yB);

        x_face.push_back(xm);
        y_face.push_back(ym);

        x_min = std::min(x_min, xm);
        x_max = std::max(x_max, xm);
    }

    // Geometric chord from wall projection in x
    c_geom = x_max - x_min;

    // Quarter-chord x-position in that projection
    const double x_q = x_min + 0.25 * c_geom;

    // Pick the node closest (in x) to x_q as reference
    double bestDist = 1e300;
    x_ref = x_q;
    y_ref = 0.0;

    for (std::size_t k = 0; k < x_face.size(); ++k) {
        const double dx = std::fabs(x_face[k] - x_q);
        if (dx < bestDist) {
            bestDist = dx;
            x_ref    = x_face[k];
            y_ref    = y_face[k];
        }
    }
}

inline AeroCoeffs computeAirfoilCoeffs(
    const std::vector<Conservative>& U,
    const Mesh& mesh,
    const Config& cfg
)
{
    int imin, imax, jmin, jmax;
    mesh.getInteriorBounds(imin, imax, jmin, jmax);

    // Bottom wall is j = jmin in your O-grid convention
    const int jWallCell = jmin;
    const int jWallFace = jWallCell;

    double x_ref = 0.0, y_ref = 0.0;
    double c_geom = 1.0;  // geometric chord length (from mesh)

    // Get quarter-chord ref point + geometric chord
    findQuarterChordRef(mesh, jWallFace, x_ref, y_ref, c_geom);

    // Freestream and reference quantities
    const double u_inf   = cfg.u_inf;
    const double v_inf   = cfg.v_inf;
    const double V_inf2  = u_inf*u_inf + v_inf*v_inf;
    const double q_inf   = 0.5 * cfg.rho_inf * V_inf2;

    // Use geometric chord for both S_ref and c_ref (2D, unit span)
    const double c_ref   = c_geom;
    const double S_ref   = c_geom;   // chord * unit span

    double Fx = 0.0;
    double Fy = 0.0;
    double M_ref = 0.0;

    auto nidx = [&](int I, int J){ return mesh.nodeIndex(I, J); };

    for (int i = imin; i < imax; ++i) {
        // Pressure on wall cell just above the wall face row
        const int cWall   = mesh.cellIndex(i, jWallCell);
        const double p_w  = EOS::pressure(U[cWall], cfg.gamma);

        // Face normal and length along the wall
        const int f       = mesh.iFaceIndex(i, jWallFace);
        const auto& n     = mesh.iFaceNormal[f];
        const double ds   = mesh.iFaceLen[f];

        const double dFx  = -p_w * n[0] * ds;
        const double dFy  = -p_w * n[1] * ds;

        Fx += dFx;
        Fy += dFy;

        // Moment arm: center of the face segment
        const int nA = nidx(i,   jWallFace);
        const int nB = nidx(i+1, jWallFace);

        const double xA = mesh.xNodes[nA];
        const double yA = mesh.yNodes[nA];
        const double xB = mesh.xNodes[nB];
        const double yB = mesh.yNodes[nB];

        const double x_cp = 0.5*(xA + xB);
        const double y_cp = 0.5*(yA + yB);

        const double dx = x_cp - x_ref;
        const double dy = y_cp - y_ref;

        // z-moment (out of plane) about (x_ref,y_ref)
        M_ref -= (dx * dFy - dy * dFx);
    }

    AeroCoeffs C{};

    if (q_inf > 0.0 && S_ref > 0.0 && c_ref > 0.0) {
        const double alpha = cfg.alpha_rad;
        const double ca = std::cos(alpha);
        const double sa = std::sin(alpha);

        // Rotate total force into lift/drag at infinity
        const double Lift = -Fx * sa + Fy * ca;
        const double Drag =  Fx * ca + Fy * sa;

        C.CL = Lift / (q_inf * S_ref);
        C.CD = Drag / (q_inf * S_ref);
        C.CM = M_ref / (q_inf * S_ref * c_ref);
    } else {
        C.CL = C.CD = C.CM = 0.0;
    }

    return C;
}

// ============================================================================
//   Compact convergence + global + dissipation stats
// ============================================================================

struct ResidualStats {
    double L2[4]{}, Linf[4]{}, L2_norm[4]{};
    double L2_total{}, Linf_total{};
};

struct GlobalStats {
    double mass{}, momx{}, momy{}, energy{};
    double massDriftPPM{};
    double minRho{}, minP{}, maxMach{};
    double entropyErr{};
};

struct DissipationStats {
    double D2_L1_total{}, D4_L1_total{}, Dfrac{};
};

static inline void printResidualL2Only(int it, const ResidualStats& RS)
{
    std::printf("[L2] it=%5d  R=[% .6e % .6e % .6e % .6e]  Rtot=% .6e  \n",
        it,
        RS.L2[0], RS.L2[1], RS.L2[2], RS.L2[3],
        RS.L2_total
       );
}

static inline void dt_stats_phys(const Mesh& m, const std::vector<double>& dt, double& dt_min, double& dt_avg){
    int imin,imax,jmin,jmax; m.getInteriorBounds(imin,imax,jmin,jmax);
    dt_min = std::numeric_limits<double>::infinity();
    dt_avg = 0.0; size_t n=0;
    for(int j=jmin;j<jmax;++j)
        for(int i=imin;i<imax;++i){
            int c=m.cellIndex(i,j);
            dt_min = std::min(dt_min, dt[c]);
            dt_avg += dt[c]; ++n;
        }
    if(n) dt_avg /= double(n);
}

static inline ResidualStats computeResidualStatsL2(
    const std::vector<Conservative>& R, const Mesh& m, const double R0_L2[4])
{
    ResidualStats S;
    int imin,imax,jmin,jmax; m.getInteriorBounds(imin,imax,jmin,jmax);
    double wsum = 0.0;
    double acc2[4] = {0,0,0,0};
    double maxAbs[4] = {0,0,0,0};

    for(int j=jmin;j<jmax;++j){
        for(int i=imin;i<imax;++i){
            int id = m.cellIndex(i,j);
            const double V = m.cellArea[id];
            wsum += V;
            const auto &r = R[id];
            const double v[4] = { r.rho, r.rhou, r.rhov, r.rhoE };
            for (int k=0;k<4;++k){
                acc2[k]   += V * v[k]*v[k];
                maxAbs[k]  = std::max(maxAbs[k], std::abs(v[k]));
            }
        }
    }
    const double denom = std::max(1e-300, wsum);
    for (int k=0;k<4;++k){
        S.L2[k]      = std::sqrt(acc2[k]/denom);
        S.Linf[k]    = maxAbs[k];
        S.L2_norm[k] = (R0_L2[k] > 0.0) ? S.L2[k] / R0_L2[k] : 0.0;
        S.L2_total  += S.L2[k]*S.L2[k];
        S.Linf_total = std::max(S.Linf_total, S.Linf[k]);
    }
    S.L2_total = std::sqrt(S.L2_total);
    return S;
}

static inline GlobalStats computeGlobalStats(
    const std::vector<Conservative>& U, const Mesh& m, const Config& cfg, double M0)
{
    GlobalStats G{};
    const double gamma = cfg.gamma;
    const double rho_inf = cfg.rho_inf;
    const double p_inf   = cfg.p_inf;
    const double s_inf   = std::log(p_inf) - gamma*std::log(rho_inf);

    int imin,imax,jmin,jmax; m.getInteriorBounds(imin,imax,jmin,jmax);
    double Vtot=0.0, s_err_sum=0.0;
    double minRho=1e300, minP=1e300, maxM=0.0;

    for(int j=jmin;j<jmax;++j){
        for(int i=imin;i<imax;++i){
            int id = m.cellIndex(i,j);
            const double V = m.cellArea[id];
            Vtot += V;

            const auto &q = U[id];
            const double rho = q.rho;
            const double u = q.rhou / rho;
            const double v = q.rhov / rho;
            const double ke = 0.5 * rho * (u*u + v*v);
            const double p  = (gamma - 1.0) * (q.rhoE - ke);
            const double a2 = std::max(0.0, gamma * p / rho);
            const double a  = std::sqrt(a2);
            const double M  = std::sqrt(u*u + v*v) / std::max(1e-300, a);

            G.mass   += rho * V;
            G.momx   += q.rhou * V;
            G.momy   += q.rhov * V;
            G.energy += q.rhoE * V;

            minRho = std::min(minRho, rho);
            minP   = std::min(minP,   p);
            maxM   = std::max(maxM,   M);

            const double s = std::log(std::max(p,1e-300)) - gamma*std::log(std::max(rho,1e-300));
            s_err_sum += V * std::abs(s - s_inf);
        }
    }
    G.massDriftPPM = ((G.mass - M0) / std::max(1e-300, M0)) * 1.0e6;
    G.minRho = minRho; G.minP = minP; G.maxMach = maxM;
    G.entropyErr = s_err_sum / std::max(1e-300, Vtot);
    return G;
}

static inline DissipationStats extractDiss(const FluxCalculator& flux){
    DissipationStats D{};
    const auto& js = flux.getJSTStats();
    D.D2_L1_total = js.D2_L1_total;
    D.D4_L1_total = js.D4_L1_total;
    const double denom = D.D2_L1_total + D.D4_L1_total;
    D.Dfrac = (denom > 0.0) ? (D.D2_L1_total / denom) : std::numeric_limits<double>::quiet_NaN();
    return D;
}

// ========================== ITERATION TABLE ================================
// Cleaner table: L2 residuals + minRho, maxM, Dfrac + CL,CD,CM
// (R_n/R_0, dM[ppm], and minP/ΔP removed as requested)
// ===========================================================================

static inline void printIterHeader()
{
    std::printf(
        "-----------------------------------------------------------------------------------------------------------------\n");
    std::printf(
        " iter  CFL   dt_min     dt_avg    "
        "|   L2[rho]   L2[rhou]  L2[rhov]  L2[rhoE]  "
        "|   minRho    maxM   Dfrac  "
        "|    CL       CD       CM\n");
    std::printf(
        "-----------------------------------------------------------------------------------------------------------------\n");
}

static inline void printIterLine(
    int it, double CFL, double dt_min, double dt_avg,
    const ResidualStats& RS, const GlobalStats& GS, const DissipationStats& DS,
    const AeroCoeffs& AC)
{
    std::printf(
        "%5d %4.2f %9.3e %9.3e  "
        "| %9.2e %9.2e %9.2e %9.2e  "
        "| %9.2e %7.3f %6.3f  "
        "| %.12f %.12f %.12f\n",
        it, CFL, dt_min, dt_avg,
        RS.L2[0], RS.L2[1], RS.L2[2], RS.L2[3],
        GS.minRho, GS.maxMach, DS.Dfrac,
        AC.CL, AC.CD, AC.CM);
}

// CSV history writer (now includes CL, CD, CM)
struct HistoryCSV {
    std::ofstream f;
    void open(const std::string& path){
        f.open(path, std::ios::out);
        if (f) {
            f << "iter,CFL,dt_min,dt_avg,"
                 "R2_rho,R2_rhou,R2_rhov,R2_rhoE,"
                 "CL,CD,CM\n";
            f.flush();
        }
    }
    void write(int it, double CFL, double dt_min, double dt_avg,
               const ResidualStats& RS, const AeroCoeffs& AC)
    {
        if (!f) return;
        f << it << ","
          << CFL << ","
          << dt_min << ","
          << dt_avg << ","
          << RS.L2[0] << ","
          << RS.L2[1] << ","
          << RS.L2[2] << ","
          << RS.L2[3] << ","
          << AC.CL << ","
          << AC.CD << ","
          << AC.CM << "\n";
        f.flush();
    }
};

// ============================================================================
//                  Pressure perturbations to excite JST sensor
// ============================================================================

inline void setCellPrimitive(std::vector<Conservative>& U, const Mesh& m, int i, int j,
                             double rho, double u, double v, double p, double gamma)
{
    const int c = m.cellIndex(i,j);
    Primitive W(rho,u,v,p,gamma);
    U[c] = W.toConservative(gamma);
}

inline void imprintPressureStripeJ(std::vector<Conservative>& U, const Mesh& m, const Config& cfg,
                                   double A, int waves=1)
{
    int imin,imax,jmin,jmax; m.getInteriorBounds(imin,imax,jmin,jmax);
    const double twoPi = 2.0*M_PI; const int Nj = (jmax - jmin);
    for(int j=jmin;j<jmax;++j){
        const double phi = twoPi * waves * double(j - jmin) / std::max(1, Nj);
        const double p = cfg.p_inf * (1.0 + A * std::sin(phi));
        for(int i=imin;i<imax;++i) setCellPrimitive(U,m,i,j,cfg.rho_inf,cfg.u_inf,cfg.v_inf,p,cfg.gamma);
    }
}

inline void imprintPressureStripeI(std::vector<Conservative>& U, const Mesh& m, const Config& cfg,
                                   double A, int waves=1)
{
    int imin,imax,jmin,jmax; m.getInteriorBounds(imin,imax,jmin,jmax);
    const double twoPi = 2.0*M_PI; const int Ni = (imax - imin);
    for(int i=imin;i<imax;++i){
        const double phi = twoPi * waves * double(i - imin) / std::max(1, Ni);
        const double p = cfg.p_inf * (1.0 + A * std::sin(phi));
        for(int j=jmin;j<jmax;++j) setCellPrimitive(U,m,i,j,cfg.rho_inf,cfg.u_inf,cfg.v_inf,p,cfg.gamma);
    }
}

inline void imprintGaussianBump2D(std::vector<Conservative>& U, const Mesh& m, const Config& cfg,
                                  double A, double sigma_cells=6.0)
{
    int imin,imax,jmin,jmax; m.getInteriorBounds(imin,imax,jmin,jmax);
    const double ic = 0.5*(imin + imax - 1);
    const double jc = 0.5*(jmin + jmax - 1);
    const double s2 = sigma_cells * sigma_cells;

    for(int j=jmin;j<jmax;++j){
        for(int i=imin;i<imax;++i){
            double di=(i-ic), dj=(j-jc);
            double g = std::exp(-(di*di + dj*dj)/(2.0*s2));
            double p = cfg.p_inf * (1.0 + A*g);
            setCellPrimitive(U,m,i,j,cfg.rho_inf,cfg.u_inf,cfg.v_inf,p,cfg.gamma);
        }
    }
}

// ============================================================================
//                       Config parser (config.txt)
// ============================================================================

static bool loadRunFromConfig(const std::string& path, Config& cfg, RunOpts& opts)
{
    std::ifstream in(path);
    if(!in) return false;

    std::string line;
    while(std::getline(in,line)){
        auto pos_hash  = line.find('#');  if(pos_hash  != std::string::npos) line = line.substr(0,pos_hash);
        auto pos_slash = line.find("//"); if(pos_slash != std::string::npos) line = line.substr(0,pos_slash);
        line = trim(line);
        if(line.empty()) continue;

        std::string key,val;
        auto eq=line.find('='); 
        if(eq==std::string::npos){ std::istringstream iss(line); if(!(iss>>key>>val)) continue; }
        else { key=trim(line.substr(0,eq)); val=trim(line.substr(eq+1)); }
        key = strip_dots(lower(key));

        // files / BC toggles
        if(key=="meshfile"||key=="mesh"){ opts.meshFile=val; continue; }
        if(key=="usesimplefarfield"){ opts.use_simple_farfield = parse_bool(val); continue; }

        // flow / gas
        if(key=="mach"||key=="mach_inf"||key=="m"){ cfg.Mach_inf=std::stod(val); continue; }
        if(key=="alpha"||key=="alpha_deg"){ cfg.alpha_deg=std::stod(val); continue; }
        if(key=="gamma"){ cfg.gamma=std::stod(val); continue; }
        if(key=="rho_inf"){ cfg.rho_inf=std::stod(val); continue; }
        if(key=="p_inf"){ cfg.p_inf=std::stod(val); continue; }
        if(key=="freestreameverywhere"){ 
            opts.freestream_everywhere = parse_bool(val); 
            continue; 
        }

        // numerics
        if(key=="cfl"){ cfg.CFL=std::stod(val); continue; }
        if(key=="k2"||key=="k2_jst"){ cfg.k2_jst=std::stod(val); continue; }
        if(key=="k4"||key=="k4_jst"){ cfg.k4_jst=std::stod(val); continue; }
        if(key=="maxiter"||key=="maxiters"){ cfg.maxIter=std::stoi(val); continue; }
        if(key=="printfreq"){ cfg.printFreq=std::stoi(val); continue; }
        if(key=="outputfreq"){ cfg.outputFreq=std::stoi(val); continue; }

        // residual smoothing controls
        if(key=="useresidualsmoothing" || key=="use_residual_smoothing"){
            cfg.use_residual_smoothing = parse_bool(val);
            continue;
        }
        if(key=="smoothepsi" || key=="smooth_eps_i"){
            cfg.smooth_eps_I = std::stod(val);
            continue;
        }
        if(key=="smoothepsj" || key=="smooth_eps_j"){
            cfg.smooth_eps_J = std::stod(val);
            continue;
        }

        // debug flags
        if(key=="debugloglambdas"){ cfg.dbg_log_lambdas = parse_bool(val); continue; }
        if(key=="debuglogjst"){     cfg.dbg_log_jst     = parse_bool(val); continue; }
        if(key=="debuglogsplit"){   cfg.dbg_log_split   = parse_bool(val); continue; }

        // perturbation controls
        if(key=="perturb"||key=="perturbmode"){ opts.perturb_mode=lower(val); continue; }
        if(key=="perturbamp"){ opts.perturb_amp=std::stod(val); continue; }
        if(key=="perturbwaves"){ opts.perturb_waves=std::stoi(val); continue; }
        if(key=="perturbsigma"||key=="perturbsigmacells"){ opts.perturb_sigma_cells=std::stod(val); continue; }

        // ---------- multigrid controls ----------
        if(key=="usemultigrid"){
            cfg.use_multigrid = parse_bool(val);
            continue;
        }
        if(key=="mgcycles" || key=="mg_cycles"){
            cfg.mg_cycles = std::stoi(val);
            continue;
        }
        if(key=="mgrestol" || key=="mg_res_tol"){
            cfg.mg_res_tol = std::stod(val);
            continue;
        }
        if(key=="mgmeshfiles" || key=="mg_mesh_files"){
            // Legacy: comma-separated list; now ignored because MG builds hierarchy by coarsening.
            std::string s = val;
            std::stringstream ss(s);
            std::string token;
            while (std::getline(ss, token, ',')) {
                token = trim(token);
                if (!token.empty())
                    opts.mg_mesh_files.push_back(token);
            }
            continue;
        }
    }
    return true;
}

// ============================================================================
//                                       main
// ============================================================================

int main(int argc, char** argv)
{
    // 1) Config + optional config.txt
    Config cfg;            // defaults in Types.hpp
    cfg.initialize();      // compute derived freestream, etc.
    RunOpts opts;

    const std::string cfgPath = (argc >= 3) ? std::string(argv[2]) : "config.txt";
    if(loadRunFromConfig(cfgPath, cfg, opts)){
        cfg.initialize();  // recompute derived stuff after overrides
        std::cout << "[INFO] loaded " << cfgPath << "\n";
    } else {
        std::cout << "[INFO] using built-in defaults (no config.txt)\n";
    }

    // mesh filename (finest for both single-grid and multigrid)
    std::string meshFile = (argc >= 2) ? argv[1] : opts.meshFile;

    // =====================================================================
    // SINGLE-GRID MODE: explicit RK5-LTS driver (with MG hook)
    // =====================================================================

    Mesh m; 
    m.readPlot3D(meshFile);

    // --- Dump physical mesh once (for ParaView) ---
    m.writeVTKPhysical("mesh_physical.vtk");
    m.writeFaceVTKPhysical("mesh_faces.vtk");
    std::cout << "Mesh VTK written: mesh_physical.vtk and mesh_faces.vtk\n";

    int imin,imax,jmin,jmax; m.getInteriorBounds(imin,imax,jmin,jmax);
    std::cout << "Mesh loaded: " << meshFile << "\n";
    std::cout << "Physical cells: i=" << imin << ".." << (imax-1)
              << ", j=" << jmin << ".." << (jmax-1) << "\n";
    std::cout << "Total cells (with ghosts): " << m.niTotal << " x " << m.njTotal << "\n";

    // 3) BCs & initial state
    BoundaryConditions bc(
        cfg,
        /*use_simple_farfield=*/opts.use_simple_farfield,
        /*freestream_everywhere=*/opts.freestream_everywhere
    );
    Conservative Uinf = cfg.getFreestream().toConservative(cfg.gamma);
    std::vector<Conservative> U(m.niTotal * m.njTotal, Uinf), R(m.niTotal * m.njTotal);

    // Initialize ghosts once from freestream
    bc.initializeGhosts(U, m);

    // Optional perturbation so JST sees something
    if(opts.perturb_mode!="none"){
        if(opts.perturb_mode=="stripej"){
            imprintPressureStripeJ(U, m, cfg, opts.perturb_amp, opts.perturb_waves);
            std::cout << "[PERT] stripeJ A="<<opts.perturb_amp<<" waves="<<opts.perturb_waves<<"\n";
        } else if(opts.perturb_mode=="stripei"){
            imprintPressureStripeI(U, m, cfg, opts.perturb_amp, opts.perturb_waves);
            std::cout << "[PERT] stripeI A="<<opts.perturb_amp<<" waves="<<opts.perturb_waves<<"\n";
        } else if(opts.perturb_mode=="gaussian"){
            imprintGaussianBump2D(U, m, cfg, opts.perturb_amp, opts.perturb_sigma_cells);
            std::cout << "[PERT] gaussian A="<<opts.perturb_amp<<" sigma_cells="<<opts.perturb_sigma_cells<<"\n";
        } else {
            std::cout << "[PERT] unknown mode '"<<opts.perturb_mode<<"' -> none\n";
        }
        bc.initializeGhosts(U, m);
    }

    // initial spectral radii (for first dt)
    {
        std::vector<Primitive> W(U.size());
        for(size_t c=0;c<U.size();++c) W[c]=Primitive(U[c], cfg.gamma);
        m.computeSpectralRadius(W, cfg.gamma);
    }

    FluxCalculator flux(cfg.k2_jst, cfg.k4_jst, opts.freestream_everywhere);

    // 4) Freestream/perturbed residual baseline (iter 0)
    flux.computeResidual(U, R, m, cfg);
    double res0 = residualLinf(R, m);
    ResSplit split0 = residualSplit(R, m);

    std::cout << "Residual Linf (iter 0, all)        = " << std::scientific << res0 << "\n";
    std::cout << "  interior cells max               = " << split0.interior << "\n";
    std::cout << "  bottom-adjacent cells max        = " << split0.bottom   << "\n";
    std::cout << "  top-adjacent cells max           = " << split0.top      << "\n";
    std::cout << "  i-periodic seam cells max        = " << split0.periodic << "\n";

    double R0_L2[4] = {0,0,0,0};
    {
        auto RS0 = computeResidualStatsL2(R, m, R0_L2);
        for (int k=0;k<4;++k) R0_L2[k] = std::max(1e-300, RS0.L2[k]);
    }
    double M0 = 0.0;
    {
        auto GS0 = computeGlobalStats(U, m, cfg, /*M0*/1.0);
        M0 = GS0.mass;
    }

    HistoryCSV hist;
    hist.open("history.csv");
    printIterHeader();

    std::vector<double> dt_local;

    const double RES_TOL = 1e-12;
    cfg.mg_res_tol = RES_TOL;  
    int convergedIter = -1;
    ResidualStats RS_last{};

    auto t_start = std::chrono::high_resolution_clock::now();
    const int MG_START_ITER = 100;        // after 100 explicit sweeps (0..99)
    const int POST_MG_FINE_SWEEPS = 50;   // extra fine-grid RK5 sweeps after MG

    for(int iter=0; iter<cfg.maxIter; ++iter){
        bc.apply(U, m);
        {
            std::vector<Primitive> W(U.size());
            for(size_t c=0;c<U.size();++c) W[c]=Primitive(U[c], cfg.gamma);
            m.computeSpectralRadius(W, cfg.gamma);
        }

        buildLocalDt(m, cfg, dt_local);
        advanceExplicitRK5_LTS(m, cfg, bc, flux, U, R, dt_local);

        // Convergence stats
        auto RS = computeResidualStatsL2(R, m, R0_L2);
        RS_last = RS;

        if (!cfg.use_multigrid && RS_last.L2_total < RES_TOL) {
            convergedIter = iter;
            std::cout << "[STOP] Reached residual tolerance " << RES_TOL
                      << " at iter " << iter
                      << " (L2_total = " << RS_last.L2_total << ")\n";
            break;
        }

        if(iter % cfg.printFreq == 0){
            double dt_min, dt_avg; 
            dt_stats_phys(m, dt_local, dt_min, dt_avg);

            auto GS = computeGlobalStats(U, m, cfg, M0);
            auto DS = extractDiss(flux);
            AeroCoeffs AC = computeAirfoilCoeffs(U, m, cfg);

            printIterLine(iter, cfg.CFL, dt_min, dt_avg, RS, GS, DS, AC);

            if(cfg.dbg_log_split){
                ResSplit sp = residualSplit(R, m);
                std::cout << "  [split] interior="<<sp.interior
                          << " bottom="<<sp.bottom
                          << " top="   <<sp.top
                          << " seam="  <<sp.periodic << "\n";
            }

            if(cfg.dbg_log_lambdas){
                double Imin,Imax,Iavg,Jmin,Jmax,Javg;
                lambda_stats(m, Imin,Imax,Iavg, Jmin,Jmax,Javg);
                std::cout << "  [lambda] I(min/max/avg)= "
                          << Imin << " / " << Imax << " / " << Iavg
                          << " | J(min/max/avg)= "
                          << Jmin << " / " << Jmax << " / " << Javg << "\n";
            }

            if(cfg.dbg_log_jst){
                std::vector<Primitive> W_dbg(U.size());
                for(size_t c=0;c<U.size();++c) W_dbg[c]=Primitive(U[c], cfg.gamma);
                m.computeSpectralRadius(W_dbg, cfg.gamma);

                const auto& js = flux.getJSTStats();
                std::cout << "  [JST] eps2_I max="<< js.max_eps2_I
                          << " eps4_I max="<< js.max_eps4_I
                          << " | cnt(eps2_I>0)=" << js.cnt_eps2_I
                          << " cnt(eps4_I>0)=" << js.cnt_eps4_I << "\n";
                std::cout << "  [JST] eps2_J max="<< js.max_eps2_J
                          << " eps4_J max="<< js.max_eps4_J
                          << " | cnt(eps2_J>0)=" << js.cnt_eps2_J
                          << " cnt(eps4_J>0)=" << js.cnt_eps4_J << "\n";
            }

            hist.write(iter, cfg.CFL, dt_min, dt_avg, RS, AC);
        }

        // ================== SWITCH TO MULTIGRID HERE ==================
        if (cfg.use_multigrid && iter == MG_START_ITER - 1) {
            std::cout << "[MG] Switching to multigrid after "
                      << MG_START_ITER << " explicit iterations.\n";

            // Choose number of MG levels based on mesh; for now, hardcode 3.
            const int NUM_MG_LEVELS = 3;

            FASMultigridSolver::LevelControl defaultCtrl; // details overridden inside
            FASMultigridSolver mg(
                cfg,
                opts.use_simple_farfield,
                opts.freestream_everywhere,
                meshFile,          // finest mesh filename
                NUM_MG_LEVELS,     // number of levels (fine + coarsened)
                defaultCtrl
            );

            // Initialize MG hierarchy (reads same finest mesh, coarsens inside)
            mg.initialize();

            // Inject current single-grid solution into finest level
            mg.finestSolution() = U;

            // Rebuild ghosts/radii/residual/history on finest from this U
            mg.resetFinestFromCurrent();

            const int    nCycles = cfg.mg_cycles;   // e.g. 1 V-cycle
            const double resTol  = cfg.mg_res_tol;  // usually same as RES_TOL

            std::cout << "[MG] Running up to " << nCycles
                      << " V-cycles (resTol = " << resTol << ").\n";

            auto t_mg_start = std::chrono::high_resolution_clock::now();
            mg.runVCycles(nCycles, resTol);
            auto t_mg_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_mg = t_mg_end - t_mg_start;

            auto& levels = mg.getLevels();
            FASMultigridSolver::Level& L0 = levels.front();

            writeSolutionVTKPhysical("solution_mg.vtk", L0.mesh, L0.U, L0.cfg);

            std::cout << "[MG] Finest-level L2 residual after multigrid = "
                      << L0.currentResL2 << "\n";
            std::cout << "[MG] Wall-clock time (multigrid only) = "
                      << elapsed_mg.count() << " s\n";

            // Update U on single-grid side with MG finest solution
            U = L0.U;

            // ----------------- Extra fine-grid sweeps after MG -----------------
            std::cout << "[MG] Running " << POST_MG_FINE_SWEEPS
                      << " extra fine-grid RK5 sweeps after multigrid.\n";
            for (int k = 0; k < POST_MG_FINE_SWEEPS; ++k) {
                bc.apply(U, m);
                std::vector<Primitive> W_tail(U.size());
                for (size_t c = 0; c < U.size(); ++c) {
                    W_tail[c] = Primitive(U[c], cfg.gamma);
                }
                m.computeSpectralRadius(W_tail, cfg.gamma);

                buildLocalDt(m, cfg, dt_local);
                advanceExplicitRK5_LTS_MG(m, cfg, bc, flux, U, R, dt_local, nullptr);
            }

            // Residual / stats AFTER the post-MG sweeps (on finest)
            RS_last = computeResidualStatsL2(R, m, R0_L2);

            if (RS_last.L2_total < RES_TOL) {
                convergedIter = MG_START_ITER - 1 + POST_MG_FINE_SWEEPS;
                std::cout << "[STOP] MG + fine-tail reached global tolerance "
                          << RES_TOL << " (L2_total = " << RS_last.L2_total << ")\n";
            } else {
                convergedIter = -1;
                std::cout << "[WARN] MG + fine-tail did NOT reach global tolerance "
                          << RES_TOL << " (L2_total = " << RS_last.L2_total << ")\n";
            }

            // Exit explicit loop (we're done with iterations either way)
            break;
        }
        // ==============================================================        
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t_end - t_start;

    int finalIter = (convergedIter >= 0) ? convergedIter : (cfg.maxIter - 1);

    std::ostringstream vtkName;
    vtkName << "solution_iter" << finalIter << ".vtk";

    writeSolutionVTKPhysical(vtkName.str(), m, U, cfg);

    std::cout << "[TIME] Wall-clock time = " << elapsed.count() << " s\n";
    if (convergedIter >= 0) {
        std::cout << "[INFO] Converged (by MG/explicit criteria) at approx iteration "
                  << convergedIter
                  << " with L2_total = " << RS_last.L2_total << "\n";
    } else {
        std::cout << "[INFO] Reached maxIter = " << cfg.maxIter
                  << " without hitting tolerance " << RES_TOL << "\n";
    }

    std::cout << "Done.\n";
    return 0;
}
