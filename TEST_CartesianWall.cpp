// TEST_CartesianWall.cpp
#include <vector>
#include <cmath>
#include <cstdio>

#include "Types.hpp"     // Conservative, Primitive, Config
#include "mesh.hpp"      // Mesh
#include "Boundary.hpp"  // BoundaryConditions
#include "Flux.hpp"      // FluxCalculator


// Simple Linf norm over a Conservative residual vector
static double linfNorm(const std::vector<Conservative>& R,
                       const Mesh& mesh,
                       int imin, int imax,
                       int jmin, int jmax)
{
    double m = 0.0;
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int c = mesh.cellIndex(i,j);
            double v = std::max(
                std::max(std::fabs(R[c].rho),
                         std::fabs(R[c].rhou)),
                std::max(std::fabs(R[c].rhov),
                         std::fabs(R[c].rhoE))
            );
            if (v > m) m = v;
        }
    }
    return m;
}

int main()
{
    // -------------------------------
    // 1) Config: simple freestream
    // -------------------------------
    Config cfg;

    // Flow parameters
    cfg.gamma      = 1.4;
    cfg.Mach_inf   = 0.5;   // M_infinity
    cfg.alpha_deg  = 0.0;   // flow along +x, parallel to wall y=0

    // Reference primitive state
    cfg.rho_inf = 1.0;
    cfg.p_inf   = 1.0;

    // Numerical parameters (if you want to override)
    cfg.k2_jst = 0.5;
    cfg.k4_jst = 0.02;

    // Compute a_inf, u_inf, v_inf, alpha_rad, etc.
    cfg.initialize();

    // -------------------------------
    // 2) Read tiny Cartesian mesh
    // -------------------------------
    Mesh mesh;
    mesh.readPlot3D("cart_3x3.x");   // 1 block, 3x3 nodes → 2x2 physical cells

    const int nCells = mesh.niTotal * mesh.njTotal;
    std::vector<Conservative> U(nCells), R(nCells);

    // -------------------------------
    // 3) Initialize ALL cells to freestream
    // -------------------------------
    Primitive    W_inf = cfg.getFreestream();
    Conservative U_inf = W_inf.toConservative(cfg.gamma);

    for (int c = 0; c < nCells; ++c) {
        U[c] = U_inf;
    }

    // -------------------------------
    // 4) Apply BC to fill ghost cells
    //    - bottom: slip wall (in BoundaryConditions)
    //    - top: SIMPLE Dirichlet farfield (pure constant state)
    // -------------------------------
    BoundaryConditions bc(cfg, /*use_simple_farfield=*/true);
    bc.initializeGhosts(U, mesh);

    // -------------------------------
    // 5) Compute spectral radii (REQUIRED before flux calculation!)
    // -------------------------------
    std::vector<Primitive> W(nCells);
    for (int c = 0; c < nCells; ++c) {
        W[c] = Primitive(U[c], cfg.gamma);
    }
    mesh.computeSpectralRadius(W, cfg.gamma);

    // -------------------------------
    // 6) Compute residual once (Central + JST)
    // -------------------------------
    FluxCalculator flux(cfg.k2_jst, cfg.k4_jst);
    flux.computeResidual(U, R, mesh, cfg);

    int imin, imax, jmin, jmax;
    mesh.getInteriorBounds(imin, imax, jmin, jmax);

    // Linf on ALL interior cells
    double Linf_all    = linfNorm(R, mesh, imin, imax, jmin, jmax);
    // Linf ONLY on bottom interior row j = jmin (wall row)
    double Linf_bottom = linfNorm(R, mesh, imin, imax, jmin, jmin+1);

    // -------------------------------
    // 6) Print
    // -------------------------------
    std::printf("Cartesian constant test (2x2 physical cells)\n");
    std::printf("  ni      = %d, nj      = %d\n", mesh.ni,      mesh.nj);
    std::printf("  niTotal = %d, njTotal = %d\n", mesh.niTotal, mesh.njTotal);
    std::printf("  interior: i=[%d,%d), j=[%d,%d)\n", imin, imax, jmin, jmax);
    std::printf("  Linf(interior)      = %.15e\n", Linf_all);
    std::printf("  Linf(bottom row)    = %.15e\n", Linf_bottom);

    return 0;
}