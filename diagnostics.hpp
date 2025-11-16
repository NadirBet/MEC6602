// Diagnostics.hpp
#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include "mesh.hpp"


inline void runMeshDiagnostics(const Mesh& m, bool writeVTK = true) {
    int imin, imax, jmin, jmax;
    m.getInteriorBounds(imin, imax, jmin, jmax);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "mesh diagnostics:\n";
    std::cout << "  physical cells: i=" << imin << ".." << (imax-1)
              << ", j=" << jmin << ".." << (jmax-1) << "\n";
    std::cout << "  total (with ghosts): " << m.niTotal << " x " << m.njTotal << "\n";

    // -----------------------------------------------------
    // 1) PHYSICAL CELLS: dump a few
    // -----------------------------------------------------
    std::cout << "\n=== PHYSICAL CELLS (shared faces) ===\n";
    for (int j = jmin; j < std::min(jmin + 2, jmax); ++j) {
        for (int i = imin; i < std::min(imin + 4, imax); ++i) {
            int c  = m.cellIndex(i,j);
            int fS = m.iFaceIndex(i,   j);
            int fN = m.iFaceIndex(i,   j+1);
            int fW = m.jFaceIndex(i,   j);
            int fE = m.jFaceIndex(i+1, j);

            std::cout << "Cell (" << i << "," << j << ")\n";
            std::cout << "  center : (" << m.x(i,j) << ", " << m.y(i,j) << ")\n";
            std::cout << "  area   : " << m.cellArea[c] << "\n";
            std::cout << "  south i-face: n=(" << m.iFaceNormal[fS][0] << "," << m.iFaceNormal[fS][1]
                      << ") len=" << m.iFaceLen[fS] << "\n";
            std::cout << "  north i-face: n=(" << m.iFaceNormal[fN][0] << "," << m.iFaceNormal[fN][1]
                      << ") len=" << m.iFaceLen[fN] << "\n";
            std::cout << "  west  j-face: n=(" << m.jFaceNormal[fW][0] << "," << m.jFaceNormal[fW][1]
                      << ") len=" << m.jFaceLen[fW] << "\n";
            std::cout << "  east  j-face: n=(" << m.jFaceNormal[fE][0] << "," << m.jFaceNormal[fE][1]
                      << ") len=" << m.jFaceLen[fE] << "\n\n";
        }
    }

    // -----------------------------------------------------
    // 2) LOCAL CELL GEOMETRY TEST (outward face-sum ~ 0)
    // -----------------------------------------------------
    std::cout << "=== LOCAL CELL GEOMETRY TEST ===\n";
    {
        int ic = imin + 1, jc = jmin + 1;
        int fS = m.iFaceIndex(ic,   jc);
        int fN = m.iFaceIndex(ic,   jc+1);
        int fW = m.jFaceIndex(ic,   jc);
        int fE = m.jFaceIndex(ic+1, jc);
        double Rx = 0.0, Ry = 0.0;
        Rx -= m.iFaceNormal[fS][0]; Ry -= m.iFaceNormal[fS][1];
        Rx += m.iFaceNormal[fN][0]; Ry += m.iFaceNormal[fN][1];
        Rx -= m.jFaceNormal[fW][0]; Ry -= m.jFaceNormal[fW][1];
        Rx += m.jFaceNormal[fE][0]; Ry += m.jFaceNormal[fE][1];
        std::cout << "Cell (" << ic << "," << jc << ") outward face-sum: "
                  << "sx=" << Rx << ", sy=" << Ry << "\n";
    }

    // -----------------------------------------------------
    // 3) CONSTANT-FLUX METRIC TEST (Fc = 1)
    // -----------------------------------------------------
    std::cout << "=== CONSTANT-FLUX METRIC TEST (Fc=1) ===\n";
    double maxRes = 0.0;
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int fS = m.iFaceIndex(i,   j);
            int fN = m.iFaceIndex(i,   j+1);
            int fW = m.jFaceIndex(i,   j);
            int fE = m.jFaceIndex(i+1, j);
            double Rx = 0.0, Ry = 0.0;
            Rx -= m.iFaceNormal[fS][0]; Ry -= m.iFaceNormal[fS][1];
            Rx += m.iFaceNormal[fN][0]; Ry += m.iFaceNormal[fN][1];
            Rx -= m.jFaceNormal[fW][0]; Ry -= m.jFaceNormal[fW][1];
            Rx += m.jFaceNormal[fE][0]; Ry += m.jFaceNormal[fE][1];
            double mag = std::sqrt(Rx*Rx + Ry*Ry);
            if (mag > maxRes) maxRes = mag;
        }
    }
    std::cout << "max |residual| over physical cells = " << maxRes << "\n";

    // -----------------------------------------------------
    // 4) OUTER PHYSICAL BOUNDARY WALK
    // -----------------------------------------------------
    std::cout << "=== OUTER PHYSICAL BOUNDARY (one loop) ===\n";
    {
        double bx = 0.0, by = 0.0;
        for (int i = imin; i < imax; ++i) { auto f = m.iFaceIndex(i, jmin); bx -= m.iFaceNormal[f][0]; by -= m.iFaceNormal[f][1]; }
        for (int j = jmin; j < jmax; ++j) { auto f = m.jFaceIndex(imax, j); bx += m.jFaceNormal[f][0]; by += m.jFaceNormal[f][1]; }
        for (int i = imax-1; i >= imin; --i) { auto f = m.iFaceIndex(i, jmax); bx += m.iFaceNormal[f][0]; by += m.iFaceNormal[f][1]; }
        for (int j = jmax-1; j >= jmin; --j) { auto f = m.jFaceIndex(imin, j); bx -= m.jFaceNormal[f][0]; by -= m.jFaceNormal[f][1]; }
        std::cout << "boundary sum: sx=" << bx << ", sy=" << by << "\n";
    }

    // -----------------------------------------------------
    // 5) GLOBAL INTERIOR CLOSURE (sum over all physical cells)
    // -----------------------------------------------------
    std::cout << "=== GLOBAL INTERIOR CLOSURE ===\n";
    {
        double gx = 0.0, gy = 0.0;
        for (int j = jmin; j < jmax; ++j) for (int i = imin; i < imax; ++i) {
            int fS = m.iFaceIndex(i,   j);
            int fN = m.iFaceIndex(i,   j+1);
            int fW = m.jFaceIndex(i,   j);
            int fE = m.jFaceIndex(i+1, j);
            gx -= m.iFaceNormal[fS][0]; gy -= m.iFaceNormal[fS][1];
            gx += m.iFaceNormal[fN][0]; gy += m.iFaceNormal[fN][1];
            gx -= m.jFaceNormal[fW][0]; gy -= m.jFaceNormal[fW][1];
            gx += m.jFaceNormal[fE][0]; gy += m.jFaceNormal[fE][1];
        }
        std::cout << "global interior sum: sx=" << gx << ", sy=" << gy << "\n";
    }

    // -----------------------------------------------------
    // 6) COORDINATE–METRIC TESTS (∮ x n_x ds = A, ∮ y n_y ds = A)
    // -----------------------------------------------------
    std::cout << "=== COORDINATE–METRIC TESTS ===\n";
    {
        auto node = [&](int ii, int jj) {
            int idx = m.nodeIndex(ii, jj);
            return std::pair<double,double>{ m.xNodes[idx], m.yNodes[idx] };
        };

        double maxCoordRes = 0.0;
        for (int j = jmin; j < jmax; ++j) for (int i = imin; i < imax; ++i) {
            int fS = m.iFaceIndex(i,   j);
            int fN = m.iFaceIndex(i,   j+1);
            int fW = m.jFaceIndex(i,   j);
            int fE = m.jFaceIndex(i+1, j);

            auto [xs0, ys0] = node(i,   j);
            auto [xs1, ys1] = node(i+1, j);
            auto [xn0, yn0] = node(i,   j+1);
            auto [xn1, yn1] = node(i+1, j+1);
            auto [xw0, yw0] = node(i,   j);
            auto [xw1, yw1] = node(i,   j+1);
            auto [xe0, ye0] = node(i+1, j);
            auto [xe1, ye1] = node(i+1, j+1);

            double xs = 0.5*(xs0 + xs1), ys = 0.5*(ys0 + ys1);
            double xn = 0.5*(xn0 + xn1), yn = 0.5*(yn0 + yn1);
            double xw = 0.5*(xw0 + xw1), yw = 0.5*(yw0 + yw1);
            double xe = 0.5*(xe0 + xe1), ye = 0.5*(ye0 + ye1);

            double Sx = 0.0, Sy = 0.0;
            Sx -= xs * m.iFaceNormal[fS][0]; Sy -= ys * m.iFaceNormal[fS][1];
            Sx += xn * m.iFaceNormal[fN][0]; Sy += yn * m.iFaceNormal[fN][1];
            Sx -= xw * m.jFaceNormal[fW][0]; Sy -= yw * m.jFaceNormal[fW][1];
            Sx += xe * m.jFaceNormal[fE][0]; Sy += ye * m.jFaceNormal[fE][1];

            double A = m.cellArea[m.cellIndex(i,j)];
            double errX = std::fabs(Sx - A);
            double errY = std::fabs(Sy - A);
            double locMax = (errX > errY) ? errX : errY;
            if (locMax > maxCoordRes) maxCoordRes = locMax;
        }
        std::cout << "max coordinate-metric residual = " << maxCoordRes << "\n";
    }

    // -----------------------------------------------------
    // 7) PERIODIC NODE CHECK (i-direction)
    // -----------------------------------------------------
    std::cout << "=== PERIODIC NODE CHECK ===\n";
    {
        double maxPerNode = 0.0;
        for (int j = 0; j < m.njNodes; ++j) {
            for (int g = 0; g < NGHOST; ++g) {
                int ig = imin - 1 - g, is = imax - g;
                int idxg = m.nodeIndex(ig, j), idxs = m.nodeIndex(is, j);
                double dx = m.xNodes[idxg] - m.xNodes[idxs];
                double dy = m.yNodes[idxg] - m.yNodes[idxs];
                maxPerNode = std::max(maxPerNode, std::hypot(dx,dy));
            }
            for (int g = 0; g < NGHOST; ++g) {
                int ig = imax + 1 + g, is = imin + g;
                int idxg = m.nodeIndex(ig, j), idxs = m.nodeIndex(is, j);
                double dx = m.xNodes[idxg] - m.xNodes[idxs];
                double dy = m.yNodes[idxg] - m.yNodes[idxs];
                maxPerNode = std::max(maxPerNode, std::hypot(dx,dy));
            }
        }
        std::cout << "max periodic node error = " << maxPerNode << "\n";
    }

    // -----------------------------------------------------
    // 8) PERIODIC i-FACE CHECK (ghost faces must match)
    // -----------------------------------------------------
    std::cout << "=== PERIODIC i-FACE CHECK ===\n";
    {
        double maxPerIFaceN = 0.0, maxPerIFaceL = 0.0;
        for (int j = 0; j <= m.njTotal; ++j) {
            int fL1 = m.iFaceIndex(imin - 1, j);
            int fR1 = m.iFaceIndex(imax,     j);
            double dx1 = m.iFaceNormal[fL1][0] - m.iFaceNormal[fR1][0];
            double dy1 = m.iFaceNormal[fL1][1] - m.iFaceNormal[fR1][1];
            maxPerIFaceN = std::max(maxPerIFaceN, std::hypot(dx1,dy1));
            maxPerIFaceL = std::max(maxPerIFaceL, std::fabs(m.iFaceLen[fL1] - m.iFaceLen[fR1]));

            int fL2 = m.iFaceIndex(imin - 2, j);
            int fR2 = m.iFaceIndex(imax - 1, j);
            double dx2 = m.iFaceNormal[fL2][0] - m.iFaceNormal[fR2][0];
            double dy2 = m.iFaceNormal[fL2][1] - m.iFaceNormal[fR2][1];
            maxPerIFaceN = std::max(maxPerIFaceN, std::hypot(dx2,dy2));
            maxPerIFaceL = std::max(maxPerIFaceL, std::fabs(m.iFaceLen[fL2] - m.iFaceLen[fR2]));
        }
        std::cout << "max periodic i-face normal error = " << maxPerIFaceN << "\n";
        std::cout << "max periodic i-face length error = " << maxPerIFaceL << "\n";
    }

    // -----------------------------------------------------
    // 9) GHOST VERTICAL EXTRAPOLATION CHECK (top/bottom)
    // -----------------------------------------------------
    std::cout << "=== GHOST VERTICAL EXTRAPOLATION CHECK ===\n";
    {
        double maxGhostVertErr = 0.0;
        // bottom ghosts
        for (int i = 0; i < m.niNodes; ++i) for (int g = 0; g < NGHOST; ++g) {
            int jg = jmin - 1 - g, j1 = jmin, j2 = jmin + 1 + g;
            int idxg = m.nodeIndex(i, jg), idx1 = m.nodeIndex(i, j1), idx2 = m.nodeIndex(i, j2);
            double xref = 2.0 * m.xNodes[idx1] - m.xNodes[idx2];
            double yref = 2.0 * m.yNodes[idx1] - m.yNodes[idx2];
            maxGhostVertErr = std::max(maxGhostVertErr, std::max(std::fabs(m.xNodes[idxg]-xref), std::fabs(m.yNodes[idxg]-yref)));
        }
        // top ghosts
        for (int i = 0; i < m.niNodes; ++i) for (int g = 0; g < NGHOST; ++g) {
            int jg = jmax + 1 + g, j1 = jmax - g, j2 = jmax - 1 - g;
            int idxg = m.nodeIndex(i, jg), idx1 = m.nodeIndex(i, j1), idx2 = m.nodeIndex(i, j2);
            double xref = 2.0 * m.xNodes[idx1] - m.xNodes[idx2];
            double yref = 2.0 * m.yNodes[idx1] - m.yNodes[idx2];
            maxGhostVertErr = std::max(maxGhostVertErr, std::max(std::fabs(m.xNodes[idxg]-xref), std::fabs(m.yNodes[idxg]-yref)));
        }
        std::cout << "max vertical ghost extrapolation error = " << maxGhostVertErr << "\n";
    }

    // -----------------------------------------------------
    // 10) AREA / FACE LENGTH STATS
    // -----------------------------------------------------
    std::cout << "=== AREA / FACE LENGTH STATS ===\n";
    {
        double minA = 1e30, maxA2 = 0.0;
        for (int j = 0; j < m.njTotal; ++j)
            for (int i = 0; i < m.niTotal; ++i) {
                double A = m.cellArea[m.cellIndex(i,j)];
                minA = std::min(minA, A);
                maxA2 = std::max(maxA2, A);
            }
        double minIF = 1e30, maxIF = 0.0;
        for (double L : m.iFaceLen) { minIF = std::min(minIF, L); maxIF = std::max(maxIF, L); }
        double minJF = 1e30, maxJF = 0.0;
        for (double L : m.jFaceLen) { minJF = std::min(minJF, L); maxJF = std::max(maxJF, L); }
        std::cout << "cell area min=" << minA << " max=" << maxA2 << "\n";
        std::cout << "i-face len min=" << minIF << " max=" << maxIF << "\n";
        std::cout << "j-face len min=" << minJF << " max=" << maxJF << "\n";
    }

    // -----------------------------------------------------
    // 11) VTK dumps (optional)
    // -----------------------------------------------------
    if (writeVTK) {
        m.writeVTKPhysical("mesh.vtk");
        m.writeFaceVTKPhysical("faces.vtk");
        std::cout << "VTK written: mesh.vtk, faces.vtk\n";
    }
}


void testUniformFlow(const Mesh& mesh, const Config& cfg,
                     FluxCalculator& flux,
                     std::vector<Conservative>& U,
                     std::vector<Conservative>& R)
{
    const double rho0=1.0, u0=0.7, v0=0.2, p0=1.0;
    U.assign(U.size(), {});
    for (auto& u : U) {
        u.rho  = rho0;
        u.rhou = rho0*u0;
        u.rhov = rho0*v0;
        u.rhoE = p0/(cfg.gamma-1.0) + 0.5*rho0*(u0*u0+v0*v0);
    }

    // primitives for spectral radius
    std::vector<Primitive> W(U.size());
    for (size_t c=0; c<U.size(); ++c) W[c] = Primitive(U[c], cfg.gamma);
    const_cast<Mesh&>(mesh).computeSpectralRadius(W, cfg.gamma);

#if defined(DIAG_SPECTRAL)
    mesh.debugSpectralRadius();
#endif

    flux.computeResidual(U, R, mesh, cfg);

    double maxabs=0.0;
    for (const auto& r : R) {
        maxabs = std::max(maxabs, std::abs(r.rho));
        maxabs = std::max(maxabs, std::abs(r.rhou));
        maxabs = std::max(maxabs, std::abs(r.rhov));
        maxabs = std::max(maxabs, std::abs(r.rhoE));
    }
    std::printf("Uniform flow max |R| = %.3e  (expect ~1e-12 .. 1e-10)\n", maxabs);
}


// build a small sinusoid in j only (constant in i)
static void set_vary_in_j(const Mesh& mesh, const Config& cfg,
                          std::vector<Conservative>& U)
{
    int imin, imax, jmin, jmax;
    mesh.getInteriorBounds(imin, imax, jmin, jmax);

    const double rho0=1.0, u0=0.3, v0=0.0, p0=1.0;
    for (int j=0; j<mesh.njTotal; ++j) {
        // use *cell* j index mapped to [0,1]
        double t = double(j - jmin) / std::max(1, jmax - jmin);
        double dp = 0.1 * std::sin(2.0*M_PI*t); // small pressure wiggle along j
        for (int i=0; i<mesh.niTotal; ++i) {
            int c = mesh.cellIndex(i,j);
            U[c].rho  = rho0;
            U[c].rhou = rho0*u0;
            U[c].rhov = rho0*v0;
            U[c].rhoE = (p0+dp)/(cfg.gamma-1.0) + 0.5*rho0*(u0*u0+v0*v0);
        }
    }
}

static void set_vary_in_i(const Mesh& mesh, const Config& cfg,
                          std::vector<Conservative>& U)
{
    int imin, imax, jmin, jmax;
    mesh.getInteriorBounds(imin, imax, jmin, jmax);

    const double rho0=1.0, u0=0.3, v0=0.0, p0=1.0;
    for (int j=0; j<mesh.njTotal; ++j) {
        for (int i=0; i<mesh.niTotal; ++i) {
            double s = double(i - imin) / std::max(1, imax - imin);
            double dp = 0.1 * std::sin(2.0*M_PI*s); // wiggle along i
            int c = mesh.cellIndex(i,j);
            U[c].rho  = rho0;
            U[c].rhou = rho0*u0;
            U[c].rhov = rho0*v0;
            U[c].rhoE = (p0+dp)/(cfg.gamma-1.0) + 0.5*rho0*(u0*u0+v0*v0);
        }
    }
}

static double residual_inf_norm(const std::vector<Conservative>& R)
{
    double m=0.0;
    for (auto& r : R) {
        m = std::max(m, std::abs(r.rho));
        m = std::max(m, std::abs(r.rhou));
        m = std::max(m, std::abs(r.rhov));
        m = std::max(m, std::abs(r.rhoE));
    }
    return m;
}

void testDirectionalScaling(Mesh& mesh, const Config& cfg,
                            FluxCalculator& flux,
                            std::vector<Conservative>& U,
                            std::vector<Conservative>& R)
{
    // === Case 1: varies only in j (I-faces active) ===
    set_vary_in_j(mesh, cfg, U);

    // recompute radii
    std::vector<Primitive> W(U.size());
    for (size_t c=0; c<U.size(); ++c) W[c] = Primitive(U[c], cfg.gamma);
    mesh.computeSpectralRadius(W, cfg.gamma);

    flux.computeResidual(U, R, mesh, cfg);
    double base = residual_inf_norm(R);
    std::printf("I-face case: base ||R||_inf = %.3e\n", base);

#if defined(DIAG_LAMBDA_SCALE)
    // scale I by 10, J by 0.1 → only I should matter here
    mesh.scaleLambdas(10.0, 0.1);
    flux.computeResidual(U, R, mesh, cfg);
    double scaled = residual_inf_norm(R);
    std::printf("I-face case: scaled ||R||_inf = %.3e  (ratio=%.2f)\n", scaled, scaled/(base+1e-300));
    // Expect ratio ~10 (± a bit). If you swapped families, ratio will ~0.1.
#endif

    // === Case 2: varies only in i (J-faces active) ===
    set_vary_in_i(mesh, cfg, U);
    for (size_t c=0; c<U.size(); ++c) W[c] = Primitive(U[c], cfg.gamma);
    mesh.computeSpectralRadius(W, cfg.gamma);

    flux.computeResidual(U, R, mesh, cfg);
    base = residual_inf_norm(R);
    std::printf("J-face case: base ||R||_inf = %.3e\n", base);

#if defined(DIAG_LAMBDA_SCALE)
    // scale I by 0.1, J by 10 → only J should matter here
    mesh.scaleLambdas(0.1, 10.0);
    flux.computeResidual(U, R, mesh, cfg);
    scaled = residual_inf_norm(R);
    std::printf("J-face case: scaled ||R||_inf = %.3e  (ratio=%.2f)\n", scaled, scaled/(base+1e-300));
    // Expect ratio ~10 if families are correct.
#endif
}
