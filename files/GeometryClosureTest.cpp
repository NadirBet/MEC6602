#include <iostream>
#include <vector>
#include "Types.hpp"
#include "Mesh.hpp"
#include "Solver_diagnostics.hpp"   // your diagnostics header

int main() {
    // 1) Choose a small, clean Cartesian grid
    int ni = 64;
    int nj = 64;

    Mesh mesh(ni, nj);

    // 2) Build a simple unit square [0,1] x [0,1] Cartesian grid
    //    using the helper we just added to Mesh.hpp
    mesh.initCartesianTest(0.0, 1.0, 0.0, 1.0);

    // 3) Minimal Config: just enough to get gamma + freestream
    Config cfg;
    cfg.Mach_inf   = 0.5;
    cfg.alpha_deg  = 0.0;
    cfg.CFL        = 0.0;   // not used here
    cfg.k2_jst     = 0.0;
    cfg.k4_jst     = 0.0;
    cfg.maxIter    = 0;
    cfg.printFreq  = 0;
    cfg.outputFreq = 0;
    cfg.initialize();       // so cfg.gamma, u_inf, v_inf, etc. are consistent

    // 4) Build a uniform state U = freestream everywhere
    int ncells = mesh.niTotal * mesh.njTotal;
    std::vector<Conservative> U(ncells);

    Primitive Winf = cfg.getFreestream();
    Conservative Uinf = Winf.toConservative(cfg.gamma);

    for (int j = 0; j < mesh.njTotal; ++j) {
        for (int i = 0; i < mesh.niTotal; ++i) {
            int idx = mesh.cellIndex(i, j);
            U[idx] = Uinf;
        }
    }

    // 5) Build geometry (areas, normals, etc.) using your usual path
    mesh.buildGeometry(U, cfg);

    // 6) Create diagnostics
    SolverDiagnostics diag(mesh, cfg);

    // 7) Geometry checks (includes Test 6: face-closure)
    diag.validateGeometry();

    // 8) Central-flux residual audit with uniform U
    //    For a perfect Cartesian grid + uniform U, this should be ~0
    diag.centralFluxResidualAudit(mesh, cfg, U);

    return 0;
}