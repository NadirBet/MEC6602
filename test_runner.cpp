
// Thorough geometry/BC/freestream tests with logging to a single file.
// Build example:
//   g++ -O2 -std=c++17 tests_runner.cpp mesh.cpp Flux.cpp BoundaryConditions.cpp -o tests_runner
//
// Run:
//   ./tests_runner mesh.x
//
// Output:
//   ./tests_output/tests_01_geometry_bc_and_freestream.txt

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <sstream>

#include "mesh.hpp"
#include "Types.hpp"
#include "Flux.hpp"
#include "Boundary.hpp"

static std::string fmt(double v, int w=12, int p=5) {
    std::ostringstream ss;
    ss << std::setw(w) << std::setprecision(p) << std::scientific << v;
    return ss.str();
}

inline int cidx(int i, int j, int niTotal) {
    return j * niTotal + i;
}

// -----------------------------------------------------------------------------
// A minimal central convective flux (NO DISSIPATION) to probe face fluxes
// n is UNIT normal, ds is face length.
// -----------------------------------------------------------------------------

static Conservative centralFluxNoDiss(
    const Conservative& UL,
    const Conservative& UR,
    const std::array<double,2>& n_in, double ds,
    const Config& cfg
) {
    double nn = std::hypot(n_in[0], n_in[1]);
    double nx = (nn > 0.0) ? (n_in[0] / nn) : 0.0;
    double ny = (nn > 0.0) ? (n_in[1] / nn) : 0.0;

    const Primitive WL(UL, cfg.gamma);
    const Primitive WR(UR, cfg.gamma);

    const double unL = WL.u * nx + WL.v * ny;
    const double unR = WR.u * nx + WR.v * ny;

    const double HL = EOS::totalEnthalpy(UL, cfg.gamma);
    const double HR = EOS::totalEnthalpy(UR, cfg.gamma);

    Conservative F;
    F.rho  = 0.5 * (WL.rho * unL + WR.rho * unR) * ds;
    F.rhou = 0.5 * (WL.rho * WL.u * unL + WR.rho * WR.u * unR) * ds
           + 0.5 * (WL.p + WR.p) * nx * ds;
    F.rhov = 0.5 * (WL.rho * WL.v * unL + WR.rho * WR.v * unR) * ds
           + 0.5 * (WL.p + WR.p) * ny * ds;
    F.rhoE = 0.5 * (WL.rho * HL * unL + WR.rho * HR * unR) * ds;
    return F;
}

// -----------------------------------------------------------------------------
// Geometry checks (unit normals, lengths, areas, orientation sanity)
// -----------------------------------------------------------------------------
struct GeoReport {
    double minArea=1e300, maxArea=-1e300;
    double maxUnitErrIF=0.0, maxUnitErrJF=0.0;
    double minLenIF=1e300, maxLenIF=0.0;
    double minLenJF=1e300, maxLenJF=0.0;
};
static void printWallFluxTable(const Mesh& m, const Config& cfg,
                               const std::vector<Conservative>& U,
                               std::ofstream& out)
{
    int imin, imax, jmin, jmax;
    m.getInteriorBounds(imin, imax, jmin, jmax);

    auto dot2 = [](double ax, double ay, double bx, double by){ return ax*bx + ay*by; };
    auto nidx = [&](int I, int J){ return m.nodeIndex(I,J); };

    out << "[WALL TABLE: bottom-row per-cell flux decomposition]\n";
    out << "cols: i | nx ny ds | mx my | cx cy | dot(n,(m-c)) | vn_raw vn_out | "
           "u v un ut p | Fw_rho Fw_mn Fw_mt | Ft_rho Ft_mn Ft_mt | d_mn d_mt | p_face p_face*ds\n";

    for (int i = imin; i < imax; ++i) {
        // --- indices
        const int cGhost = m.cellIndex(i, jmin-1);   // bottom ghost
        const int cUp    = m.cellIndex(i, jmin);     // bottom interior cell

        // --- wall i-face (i, jmin): base normal as stored
        const int fW = m.iFaceIndex(i, jmin);
        const double nxW_raw = m.iFaceNormal[fW][0];
        const double nyW_raw = m.iFaceNormal[fW][1];
        const double dsW     = m.iFaceLen[fW];

        // face endpoints + midpoint
        const double xA = m.xNodes[nidx(i,   jmin)];
        const double yA = m.yNodes[nidx(i,   jmin)];
        const double xB = m.xNodes[nidx(i+1, jmin)];
        const double yB = m.yNodes[nidx(i+1, jmin)];
        const double mx = 0.5*(xA + xB);
        const double my = 0.5*(yA + yB);

        // interior cell center just above the wall
        const double cx = m.x(i, jmin);
        const double cy = m.y(i, jmin);

        // orientation diagnostic and outward normal
        const double dot_nc = nxW_raw*(mx - cx) + nyW_raw*(my - cy);
        const double nx_out = (dot_nc >= 0.0) ? nxW_raw : -nxW_raw;
        const double ny_out = (dot_nc >= 0.0) ? nyW_raw : -nyW_raw;
        const double tx_out = -ny_out, ty_out = nx_out;

        // state & velocities
        const Primitive Wg(U[cGhost], cfg.gamma);
        const Primitive Wu(U[cUp],    cfg.gamma);

        const double vn_raw = Wu.u*nxW_raw + Wu.v*nyW_raw;   // with stored normal
        const double vn_out = Wu.u*nx_out  + Wu.v*ny_out;    // with outward normal
        const double vt_out = Wu.u*tx_out  + Wu.v*ty_out;

        // central face pressure (note: purely diagnostic)
        const double p_face_wall = 0.5*(Wg.p + Wu.p);

        // --- top i-face (between jmin and jmin+1) for balance check
        const int fT   = m.iFaceIndex(i, jmin+1);
        const auto& nT = m.iFaceNormal[fT];
        const double dsT = m.iFaceLen[fT];
        const int cDownTop = m.cellIndex(i, jmin);
        const int cUpTop   = m.cellIndex(i, jmin+1);

        // central fluxes (NO JST), as defined by stored mesh normals
        const Conservative Fw_raw = centralFluxNoDiss(U[cGhost],   U[cUp],    m.iFaceNormal[fW], dsW, cfg);
        const Conservative Ft_raw = centralFluxNoDiss(U[cDownTop], U[cUpTop], nT,                dsT, cfg);

        // project both fluxes on the *outward* (n,t) basis so comparisons are meaningful
        const double Fw_mn = dot2(Fw_raw.rhou, Fw_raw.rhov, nx_out, ny_out);
        const double Fw_mt = dot2(Fw_raw.rhou, Fw_raw.rhov, tx_out, ty_out);
        const double Ft_mn = dot2(Ft_raw.rhou, Ft_raw.rhov, nx_out, ny_out);
        const double Ft_mt = dot2(Ft_raw.rhou, Ft_raw.rhov, tx_out, ty_out);

        // balance in the bottom cell in (n,t) components
        const double d_mn = Fw_mn - Ft_mn;  // should be ~0 in uniform freestream
        const double d_mt = Fw_mt - Ft_mt;  // should be ~0

        out << std::setw(3) << i
            << " | " << fmt(nxW_raw) << " " << fmt(nyW_raw) << " " << fmt(dsW)
            << " | " << fmt(mx) << " " << fmt(my)
            << " | " << fmt(cx) << " " << fmt(cy)
            << " | " << fmt(dot_nc)
            << " | " << fmt(vn_raw) << " " << fmt(vn_out)
            << " | " << fmt(Wu.u) << " " << fmt(Wu.v) << " " << fmt(vn_out) << " " << fmt(vt_out) << " " << fmt(Wu.p)
            << " | " << fmt(Fw_raw.rho) << " " << fmt(Fw_mn) << " " << fmt(Fw_mt)
            << " | " << fmt(Ft_raw.rho) << " " << fmt(Ft_mn) << " " << fmt(Ft_mt)
            << " | " << fmt(d_mn) << " " << fmt(d_mt)
            << " | " << fmt(p_face_wall) << " " << fmt(p_face_wall*dsW)
            << "\n";
    }

    out << "\nNotes:\n"
        << "- dot(n,(m-c)) uses the stored mesh normal; if negative, the code flips it to get the outward normal used for projections.\n"
        << "- vn_raw is with the stored normal; vn_out is with the outward normal. With slip ghosts, vn_out ≈ 0.\n"
        << "- Fw_* and Ft_* are projected on the outward (n,t) basis so d_mn and d_mt should be ~0 in freestream.\n\n";
}
static GeoReport checkGeometry(const Mesh& m, std::ofstream& out)
{
    int imin, imax, jmin, jmax;
    m.getInteriorBounds(imin, imax, jmin, jmax);

    GeoReport rep;

    // Areas: physical cells only
    for (int j = jmin; j < jmax; ++j)
        for (int i = imin; i < imax; ++i) {
            const double A = m.cellArea[m.cellIndex(i,j)];
            rep.minArea = std::min(rep.minArea, A);
            rep.maxArea = std::max(rep.maxArea, A);
        }

    // i-faces: physical strip only (j = jmin..jmax, i = imin..imax-1)
    for (int j = jmin; j <= jmax; ++j)
        for (int i = imin; i < imax; ++i) {
            const int f = m.iFaceIndex(i,j);
            const auto& n = m.iFaceNormal[f];
            const double ds = m.iFaceLen[f];
            const double mag = std::hypot(n[0], n[1]);
            rep.maxUnitErrIF = std::max(rep.maxUnitErrIF, std::fabs(mag - 1.0));
            rep.minLenIF = std::min(rep.minLenIF, ds);
            rep.maxLenIF = std::max(rep.maxLenIF, ds);
        }

    // j-faces: physical strip only (j = jmin..jmax-1, i = imin..imax)
    for (int j = jmin; j < jmax; ++j)
        for (int i = imin; i <= imax; ++i) {
            const int f = m.jFaceIndex(i,j);
            const auto& n = m.jFaceNormal[f];
            const double ds = m.jFaceLen[f];
            const double mag = std::hypot(n[0], n[1]);
            rep.maxUnitErrJF = std::max(rep.maxUnitErrJF, std::fabs(mag - 1.0));
            rep.minLenJF = std::min(rep.minLenJF, ds);
            rep.maxLenJF = std::max(rep.maxLenJF, ds);
        }

    out << "[GEOMETRY]\n";
    out << "cell area min/max = " << fmt(rep.minArea) << "  " << fmt(rep.maxArea) << "\n";
    out << "i-face |n|-1 max  = " << fmt(rep.maxUnitErrIF) << ", len min/max = "
        << fmt(rep.minLenIF) << "  " << fmt(rep.maxLenIF) << "\n";
    out << "j-face |n|-1 max  = " << fmt(rep.maxUnitErrJF) << ", len min/max = "
        << fmt(rep.minLenJF) << "  " << fmt(rep.maxLenJF) << "\n\n";
    return rep;
}

// -----------------------------------------------------------------------------
// BC consistency checks (after bc.apply)
// - Periodic-i: wrap columns equality (interior vs ghosts)
// - Bottom slip wall: normal mass flux at the wall ~ 0 for central flux
// - Top farfield Dirichlet: top ghost equals freestream
// -----------------------------------------------------------------------------
struct BCReport {
    double maxPeriodicDiff=0.0;
    double maxWallMassFlux=0.0;
    double maxFFDiff=0.0;
};

static BCReport checkBCs(const Mesh& m, const Config& cfg,
                         const BoundaryConditions& bcObj,
                         const std::vector<Conservative>& U,
                         std::ofstream& out)
{
    (void)bcObj; // not needed directly
    int imin, imax, jmin, jmax;
    m.getInteriorBounds(imin, imax, jmin, jmax);

    BCReport rep;

    // Periodic wrap (compare ghost columns to sources across all rows)
    for (int j = jmin; j < jmax; ++j) {
        for (int g = 1; g <= NGHOST; ++g) {
            int iL_ghost = imin - g;
            int iL_src   = imax - g;
            int cG = m.cellIndex(iL_ghost, j);
            int cS = m.cellIndex(iL_src,   j);
            const auto d1 = std::fabs(U[cG].rho  - U[cS].rho );
            const auto d2 = std::fabs(U[cG].rhou - U[cS].rhou);
            const auto d3 = std::fabs(U[cG].rhov - U[cS].rhov);
            const auto d4 = std::fabs(U[cG].rhoE - U[cS].rhoE);
            rep.maxPeriodicDiff = std::max({rep.maxPeriodicDiff, d1,d2,d3,d4});
        }
        for (int g = 0; g < NGHOST; ++g) {
            int iR_ghost = imax + g;
            int iR_src   = imin + g;
            int cG = m.cellIndex(iR_ghost, j);
            int cS = m.cellIndex(iR_src,   j);
            const auto d1 = std::fabs(U[cG].rho  - U[cS].rho );
            const auto d2 = std::fabs(U[cG].rhou - U[cS].rhou);
            const auto d3 = std::fabs(U[cG].rhov - U[cS].rhov);
            const auto d4 = std::fabs(U[cG].rhoE - U[cS].rhoE);
            rep.maxPeriodicDiff = std::max({rep.maxPeriodicDiff, d1,d2,d3,d4});
        }
    }

    // Bottom slip wall: check mass flux at wall faces
    for (int i = imin; i < imax; ++i) {
        const int cDown = m.cellIndex(i, jmin-1); // ghost
        const int cUp   = m.cellIndex(i, jmin);   // interior
        const int f  = m.iFaceIndex(i, jmin);
        const auto& n  = m.iFaceNormal[f];
        const double ds = m.iFaceLen[f];
        const Conservative F = centralFluxNoDiss(U[cDown], U[cUp], n, ds, cfg);
        rep.maxWallMassFlux = std::max(rep.maxWallMassFlux, std::fabs(F.rho));
    }

    // Top farfield (Dirichlet assumed for this test): top ghosts equal freestream
    const Conservative Uinf = cfg.getFreestream().toConservative(cfg.gamma);
    for (int i = imin; i < imax; ++i) {
        for (int g = 0; g < NGHOST; ++g) {
            int jg = jmax + g;
            const int c = m.cellIndex(i, jg);
            const auto d1 = std::fabs(U[c].rho  - Uinf.rho );
            const auto d2 = std::fabs(U[c].rhou - Uinf.rhou);
            const auto d3 = std::fabs(U[c].rhov - Uinf.rhov);
            const auto d4 = std::fabs(U[c].rhoE - Uinf.rhoE);
            rep.maxFFDiff = std::max({rep.maxFFDiff, d1,d2,d3,d4});
        }
    }

    out << "[BC CHECKS]\n";
    out << "periodic wrap max |ΔU| = " << fmt(rep.maxPeriodicDiff) << "\n";
    out << "wall mass-flux max |F_rho| (central) = " << fmt(rep.maxWallMassFlux) << "\n";
    out << "farfield ghosts vs Uinf max |ΔU| = " << fmt(rep.maxFFDiff) << "\n\n";
    return rep;
}

// -----------------------------------------------------------------------------
// Residual split (Linf per region)
// -----------------------------------------------------------------------------
struct ResSplit { double interior=0, bottom=0, top=0, periodic=0; };

static ResSplit residualSplit(const std::vector<Conservative>& R, const Mesh& m)
{
    int imin, imax, jmin, jmax;
    m.getInteriorBounds(imin, imax, jmin, jmax);

    ResSplit out;
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            const int id = m.cellIndex(i,j);
            double rcell = 0.0;
            rcell = std::max(rcell, std::fabs(R[id].rho));
            rcell = std::max(rcell, std::fabs(R[id].rhou));
            rcell = std::max(rcell, std::fabs(R[id].rhov));
            rcell = std::max(rcell, std::fabs(R[id].rhoE));

            const bool isBottom   = (j == jmin);
            const bool isTop      = (j == jmax - 1);
            const bool isLeftPer  = (i == imin);
            const bool isRightPer = (i == imax - 1);
            if (isLeftPer || isRightPer) out.periodic = std::max(out.periodic, rcell);
            if (isBottom) out.bottom = std::max(out.bottom, rcell);
            if (isTop)    out.top    = std::max(out.top,    rcell);
            if (!isBottom && !isTop && !isLeftPer && !isRightPer) out.interior = std::max(out.interior, rcell);
        }
    }
    return out;
}

static double residualLinf(const std::vector<Conservative>& R, const Mesh& m)
{
    int imin, imax, jmin, jmax;
    m.getInteriorBounds(imin, imax, jmin, jmax);

    double maxv = 0.0;
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            const int id = m.cellIndex(i,j);
            maxv = std::max(maxv, std::fabs(R[id].rho ));
            maxv = std::max(maxv, std::fabs(R[id].rhou));
            maxv = std::max(maxv, std::fabs(R[id].rhov));
            maxv = std::max(maxv, std::fabs(R[id].rhoE));
        }
    }
    return maxv;
}

// -----------------------------------------------------------------------------
// Main test entry
// -----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    std::string meshFile = (argc >= 2) ? argv[1] : "mesh.x";

    // Mesh
    Mesh m;
    m.readPlot3D(meshFile);

    int imin, imax, jmin, jmax;
    m.getInteriorBounds(imin, imax, jmin, jmax);

    // Output file
    std::system("mkdir -p tests_output");
    std::ofstream out("tests_output/tests_01_geometry_bc_and_freestream.txt");
    out << std::setprecision(6) << std::scientific;

    out << "=== TEST SUITE 01: geometry / BC consistency / freestream residuals ===\n";
    out << "mesh file: " << meshFile << "\n";
    out << "physical cells: i=" << imin << ".." << (imax-1)
        << "  j=" << jmin << ".." << (jmax-1) << "\n\n";

    // Config (nondimensional freestream)
    Config cfg;
    cfg.Mach_inf = 0.5;
    cfg.alpha_deg = 0.0;
    cfg.CFL = 0.5;
    cfg.k2_jst = 0.0;   // TURN OFF DISSIPATION for this test
    cfg.k4_jst = 0.0;
    cfg.initialize();

    out << "[CONFIG]\n";
    out << "M=" << cfg.Mach_inf << "  alpha(deg)=" << cfg.alpha_deg
        << "  gamma=" << cfg.gamma << "  CFL=" << cfg.CFL << "\n";
    out << "Uinf: rho=" << cfg.rho_inf
        << "  u=" << cfg.u_inf
        << "  v=" << cfg.v_inf
        << "  p=" << cfg.p_inf << "\n\n";

    // State arrays
    std::vector<Conservative> U(m.niTotal * m.njTotal);
    std::vector<Conservative> R(m.niTotal * m.njTotal);

    // Initialize freestream everywhere
    const Conservative Uinf = cfg.getFreestream().toConservative(cfg.gamma);
    for (int j = 0; j < m.njTotal; ++j)
        for (int i = 0; i < m.niTotal; ++i)
            U[cidx(i, j, m.niTotal)] = Uinf;

    BoundaryConditions bc(cfg, /*use_simple_farfield=*/true); // use Dirichlet farfield for exact match
    bc.apply(U, m);

    // Geometry checks
    (void)checkGeometry(m, out);

    // BC checks
    (void)checkBCs(m, cfg, bc, U, out);

    // Spectral radii (from current state)
    {
        std::vector<Primitive> W(U.size());
        for (size_t c=0; c<U.size(); ++c) W[c] = Primitive(U[c], cfg.gamma);
        m.computeSpectralRadius(W, cfg.gamma);
    }

    // Flux/residual check (central only, no JST)
    FluxCalculator flux(/*k2*/0.0, /*k4*/0.0);
    flux.computeResidual(U, R, m, cfg);

    const double LinfR = residualLinf(R, m);
    const ResSplit split = residualSplit(R, m);
    out << "[RESIDUALS central-only, no JST]\n";
    out << "Linf(all)   = " << fmt(LinfR) << "\n";
    out << "Linf(inter) = " << fmt(split.interior) << "\n";
    out << "Linf(bottom)= " << fmt(split.bottom)   << "\n";
    out << "Linf(top)   = " << fmt(split.top)      << "\n";
    out << "Linf(period)= " << fmt(split.periodic) << "\n\n";

    // Print sample internal and boundary fluxes (central-only definition)
    out << "[FACE FLUX SAMPLES: central convective + pressure (no JST)]\n";

    // One interior i-face (middle of domain)
    int iMid = std::clamp((imin + imax) / 2, imin + 1, imax - 2);
    int jMid = std::clamp((jmin + jmax) / 2, jmin + 1, jmax - 2);
    {
        int f = m.iFaceIndex(iMid, jMid);
        const auto& n = m.iFaceNormal[f];
        const double ds = m.iFaceLen[f];
        const int cDown = m.cellIndex(iMid, jMid-1);
        const int cUp   = m.cellIndex(iMid, jMid);
        const Conservative F = centralFluxNoDiss(U[cDown], U[cUp], n, ds, cfg);
        out << "interior i-face (i="<<iMid<<", j="<<jMid<<"): "
            << "Frho="<<fmt(F.rho)<<"  Frhou="<<fmt(F.rhou)
            << "  Frhov="<<fmt(F.rhov)<<"  FrhoE="<<fmt(F.rhoE)<<"\n";
    }

    // One interior j-face
    {
        int f = m.jFaceIndex(iMid, jMid);
        const auto& n = m.jFaceNormal[f];
        const double ds = m.jFaceLen[f];
        const int cLeft  = m.cellIndex(iMid-1, jMid);
        const int cRight = m.cellIndex(iMid,   jMid);
        const Conservative F = centralFluxNoDiss(U[cLeft], U[cRight], n, ds, cfg);
        out << "interior j-face (i="<<iMid<<", j="<<jMid<<"): "
            << "Frho="<<fmt(F.rho)<<"  Frhou="<<fmt(F.rhou)
            << "  Frhov="<<fmt(F.rhov)<<"  FrhoE="<<fmt(F.rhoE)<<"\n";
    }

    // Wall face (bottom)
    {
        int iW = iMid;
        int f = m.iFaceIndex(iW, jmin);
        const auto& n = m.iFaceNormal[f];
        const double ds = m.iFaceLen[f];
        const int cDown = m.cellIndex(iW, jmin-1);
        const int cUp   = m.cellIndex(iW, jmin);
        const Conservative F = centralFluxNoDiss(U[cDown], U[cUp], n, ds, cfg);
        out << "bottom wall i-face (i="<<iW<<", j="<<jmin<<"): "
            << "Frho="<<fmt(F.rho)<<"  Frhou="<<fmt(F.rhou)
            << "  Frhov="<<fmt(F.rhov)<<"  FrhoE="<<fmt(F.rhoE)<<"\n";
    }

    // Farfield face (top)
    {
        int iF = iMid;
        int f = m.iFaceIndex(iF, jmax);
        const auto& n = m.iFaceNormal[f];
        const double ds = m.iFaceLen[f];
        const int cDown = m.cellIndex(iF, jmax-1);
        const int cUp   = m.cellIndex(iF, jmax); // ghost
        const Conservative F = centralFluxNoDiss(U[cDown], U[cUp], n, ds, cfg);
        out << "top farfield i-face (i="<<iF<<", j="<<jmax<<"): "
            << "Frho="<<fmt(F.rho)<<"  Frhou="<<fmt(F.rhou)
            << "  Frhov="<<fmt(F.rhov)<<"  FrhoE="<<fmt(F.rhoE)<<"\n";
    }

    out << "\n[NOTES]\n"
        << "These flux samples are computed with a pure central formula (no JST) to\n"
        << "help isolate BC and geometry issues. In uniform freestream, the net\n"
        << "residual should be ~0 up to roundoff.\n";
    printWallFluxTable(m, cfg, U, out);
    out.close();
    std::cout << "Wrote tests_output/tests_01_geometry_bc_and_freestream.txt\n";
    return 0;
}
