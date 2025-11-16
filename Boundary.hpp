// BoundaryConditions.hpp
#pragma once
#include <vector>
#include <cmath>

#include "Types.hpp"
#include "mesh.hpp"    // your mesh class with iFaceIndex, iFaceNormal, ...

// Apply flow BCs to ghost CELLS, consistent with your mesh contract.
// Order: 1) i-periodic  2) bottom wall (slip)  3) top farfield
class BoundaryConditions {
public:
    BoundaryConditions(const Config& cfg_, bool use_simple_farfield_ = false,
                       bool freestream_everywhere_ = false)
        : cfg(cfg_),
          W_inf(cfg_.getFreestream()),
          use_simple_farfield(use_simple_farfield_),
          freestream_everywhere(freestream_everywhere_) {}

    void apply(std::vector<Conservative>& U, const Mesh& mesh) {
        if (freestream_everywhere) {
            applyAllFreestream(U, mesh);
            return;
        }
        applyPeriodicBC(U, mesh);
        applyWallBC(U, mesh);
        if (use_simple_farfield)
            applyFarfieldDirichlet(U, mesh);
        else
            applyFarfieldChar(U, mesh);
    }

    void initializeGhosts(std::vector<Conservative>& U, const Mesh& mesh) {
        apply(U, mesh);
    }


private:
    const Config&   cfg;
    const Primitive W_inf;   // freestream primitive
    bool            use_simple_farfield;
    bool            freestream_everywhere;  





    void applyAllFreestream(std::vector<Conservative>& U, const Mesh& mesh) {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        Conservative Uinf = W_inf.toConservative(cfg.gamma);

        for (int j = 0; j < mesh.njTotal; ++j) {
            for (int i = 0; i < mesh.niTotal; ++i) {
                bool isInterior =
                    (i >= imin && i < imax &&
                     j >= jmin && j < jmax);
                if (!isInterior) {
                    U[ mesh.cellIndex(i,j) ] = Uinf;  // only ghosts
                }
            }
        }
    }
    // ------------------------------------------------------------
    // 1) periodic in i (O-grid)
    //    copy physical columns into left/right ghost columns
    // ------------------------------------------------------------
    void applyPeriodicBC(std::vector<Conservative>& U, const Mesh& mesh) {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        // for *all* rows, including ghost rows
        for (int j = 0; j < mesh.njTotal; ++j) {
            // left ghosts: imin-1, imin-2 <- imax-1, imax-2
            for (int g = 1; g <= NGHOST; ++g) {
                int i_ghost = imin - g;
                int i_src   = imax - g;
                U[ mesh.cellIndex(i_ghost, j) ] = U[ mesh.cellIndex(i_src, j) ];
            }
            // right ghosts: imax, imax+1 <- imin, imin+1
            for (int g = 0; g < NGHOST; ++g) {
                int i_ghost = imax + g;
                int i_src   = imin + g;
                U[ mesh.cellIndex(i_ghost, j) ] = U[ mesh.cellIndex(i_src, j) ];
            }
        }
    }

    // ------------------------------------------------------------
    // 2) bottom wall (slip) at j = jmin
    //    reflect normal velocity using the bottom i-face normal
    // ------------------------------------------------------------
    /**void applyWallBC(std::vector<Conservative>& U, const Mesh& mesh) {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        for (int i = imin; i < imax; ++i) {
            // interior cell touching wall
            int c_int = mesh.cellIndex(i, jmin);
            Primitive Wc(U[c_int], cfg.gamma);

            // bottom face geometry: i-face at (i, jmin)

            // bottom wall
            int f_bottom = mesh.iFaceIndex(i, jmin);
            const auto& n = mesh.iFaceNormal[f_bottom]; // already unit
            double nx = n[0], ny = n[1];
            double vn = Wc.u * nx + Wc.v * ny;
            // slip reflection
            Primitive Wg = Wc;
            Wg.u = Wc.u - 2.0 * vn * nx;
            Wg.v = Wc.v - 2.0 * vn * ny;

            Conservative Ug = Wg.toConservative(cfg.gamma);

            // fill all bottom ghost layers: j = jmin-1, jmin-2, ...
            for (int g = 1; g <= NGHOST; ++g) {
                int jg = jmin - g;
                U[ mesh.cellIndex(i, jg) ] = Ug;
            }
        }
    }*/
    void applyWallBC(std::vector<Conservative>& U, const Mesh& mesh) {

        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);
        
        for (int i = imin; i < imax; ++i) {
            const int c_int = mesh.cellIndex(i, jmin);
            Primitive Wc(U[c_int], cfg.gamma);
            
            // Face geometry and orientation
            const int f = mesh.iFaceIndex(i, jmin);
            double nx = mesh.iFaceNormal[f][0];
            double ny = mesh.iFaceNormal[f][1];
            
            // Flip to outward (same as before)
            const auto nidx = [&](int I, int J){ return mesh.nodeIndex(I,J); };
            const double xA = mesh.xNodes[nidx(i, jmin)];
            const double yA = mesh.yNodes[nidx(i, jmin)];
            const double xB = mesh.xNodes[nidx(i+1, jmin)];
            const double yB = mesh.yNodes[nidx(i+1, jmin)];
            const double mx = 0.5*(xA + xB), my = 0.5*(yA + yB);
            const double cx = mesh.x(i, jmin), cy = mesh.y(i, jmin);
            const double vx = mx - cx, vy = my - cy;
            if (nx*vx + ny*vy < 0.0) { nx = -nx; ny = -ny; }
            
            // NEW: Iterative correction to enforce zero normal velocity at face
            Primitive Wg = Wc;
            for (int iter = 0; iter < 3; ++iter) {  // 2-3 iterations sufficient
                // Average velocity at face
                const double u_face =  Wc.u;
                const double v_face = Wc.v ;
                const double vn_face = u_face*nx + v_face*ny;
                
                // Adjust ghost to zero out face normal velocity
                Wg.u = Wc.u - 2.0 * vn_face * nx;
                Wg.v = Wc.v - 2.0 * vn_face * ny;
            }
            
            Wg.rho = Wc.rho;  // Extrapolate scalars
            Wg.p = Wc.p;
            
            const Conservative Ug = Wg.toConservative(cfg.gamma);
            for (int g = 1; g <= NGHOST; ++g) {
                U[mesh.cellIndex(i, jmin - g)] = Ug;
            }
        }
    }

    // ============================================================
    // 3A) simple top farfield: just put freestream in top ghosts
    //     this is the "make the residual disappear" version
    // ============================================================
    void applyFarfieldDirichlet(std::vector<Conservative>& U, const Mesh& mesh) {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        Conservative Uinf = W_inf.toConservative(cfg.gamma);

        for (int i = imin; i < imax; ++i) {
            for (int g = 0; g < NGHOST; ++g) {
                int jg = jmax + g;
                U[ mesh.cellIndex(i, jg) ] = Uinf;
            }
        }
    }

    // ============================================================
    // 3B) characteristic-style top farfield (your previous code)
    // ============================================================
    void applyFarfieldChar(std::vector<Conservative>& U, const Mesh& mesh) {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        for (int i = imin; i < imax; ++i) {
            // interior cell just below top
            int c_int = mesh.cellIndex(i, jmax - 1);
            Primitive Wd(U[c_int], cfg.gamma);

            // top face geometry: i-face at (i, jmax)
            int f_top = mesh.iFaceIndex(i, jmax);
            const auto& n = mesh.iFaceNormal[f_top];   // points UP, outward

            double nx  = n[0];
            double ny  = n[1];

            // freestream
            const Primitive& Wa = W_inf;

            // interior normal speed and sound
            double un_d = Wd.u * nx + Wd.v * ny;
            double c_d  = std::sqrt(cfg.gamma * Wd.p / Wd.rho);

            bool inflow     = (un_d < 0.0);                 // flow coming from outside in
            bool supersonic = (std::fabs(un_d) >= c_d);

            double rho_b, u_b, v_b, p_b;

            if (inflow && supersonic) {
                // full freestream
                rho_b = Wa.rho;
                u_b   = Wa.u;
                v_b   = Wa.v;
                p_b   = Wa.p;
            } else if (!inflow && supersonic) {
                // fully determined by interior
                rho_b = Wd.rho;
                u_b   = Wd.u;
                v_b   = Wd.v;
                p_b   = Wd.p;
            } else if (inflow && !supersonic) {
                // subsonic inflow: mix interior & freestream
                double pa = Wa.p;
                double pd = Wd.p;
                double corr = nx*(Wa.u - Wd.u) + ny*(Wa.v - Wd.v);
                p_b = 0.5 * (pa + pd - Wd.rho * c_d * corr);

                rho_b = Wa.rho + (p_b - pa)/(c_d*c_d);
                u_b   = Wa.u   - nx*(pa - p_b)/(Wd.rho * c_d);
                v_b   = Wa.v   - ny*(pa - p_b)/(Wd.rho * c_d);
            } else {
                // subsonic outflow: fix pressure to freestream
                double pa = Wa.p;
                double pd = Wd.p;
                p_b   = pa;
                rho_b = Wd.rho + (p_b - pd)/(c_d*c_d);
                u_b   = Wd.u + nx*(pd - p_b)/(Wd.rho * c_d);
                v_b   = Wd.v + ny*(pd - p_b)/(Wd.rho * c_d);
            }

            // build conservative
            double E = p_b / ((cfg.gamma - 1.0) * rho_b) + 0.5*(u_b*u_b + v_b*v_b);
            Conservative Ub;
            Ub.rho  = rho_b;
            Ub.rhou = rho_b * u_b;
            Ub.rhov = rho_b * v_b;
            Ub.rhoE = rho_b * E;

            // write into top ghost rows: j = jmax, jmax+1, ...
            for (int g = 0; g < NGHOST; ++g) {
                int jg = jmax + g;
                U[ mesh.cellIndex(i, jg) ] = Ub;
            }
        }
    }
};
