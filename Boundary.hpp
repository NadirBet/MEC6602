#pragma once
#include <vector>
#include <cmath>
#include "Types.hpp"
#include "Mesh.hpp"

// ============================================================
// CORRECTED BOUNDARY CONDITIONS
// Key points in this version:
// 1. Wall slip condition (normal velocity reflected)
// 2. Ghost cells consistently handled
// 3. Farfield: **freestream Dirichlet BC** (no characteristic BC)
// 4. Proper periodic boundary for O-grid
// ============================================================

class BoundaryConditions {
private:
    const Config& cfg;
    const Primitive U_inf;
    
public:
    BoundaryConditions(const Config& cfg_) 
        : cfg(cfg_), U_inf(cfg_.getFreestream()) {}
    
    // Apply all boundary conditions
    void apply(std::vector<Conservative>& U, const Mesh& mesh) {
        applyWallBC(U, mesh);
        applyFarfieldBC(U, mesh);
        applyPeriodicBC(U, mesh);
    }
    
    // Initialize ghost cells at start
    void initializeGhosts(std::vector<Conservative>& U, const Mesh& mesh) {
        // Apply all BCs to set initial ghost values
        apply(U, mesh);
    }
    
private:
    // =====================================================
    // WALL SLIP BC (bottom boundary at j=jmin)
    // Slip condition: reflect normal velocity component
    // =====================================================
    void applyWallBC(std::vector<Conservative>& U, const Mesh& mesh) {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);
        
        for (int i = imin; i < imax; ++i) {
            // Interior cell state
            int idx_int = mesh.cellIndex(i, jmin);
            Primitive W_int(U[idx_int], cfg.gamma);
            
            // Wall face normal (points into fluid)
            //std::array<double,2> normal = mesh.faceNormal_i[idx_int];
            std::array<double,2> normal = mesh.faceNormal_i_out[idx_int];
            double mag = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
            double nx = normal[0] / (mag + 1e-14);
            double ny = normal[1] / (mag + 1e-14);
            
            // Normal velocity
            double vn = W_int.u * nx + W_int.v * ny;
            
            // Slip reflection: remove normal component
            Primitive W_ghost = W_int;
            W_ghost.u = W_int.u - 2.0 * vn * nx;
            W_ghost.v = W_int.v - 2.0 * vn * ny;
            
            Conservative U_ghost = W_ghost.toConservative(cfg.gamma);
            
            // Same ghost state for all wall ghost layers
            for (int g = 1; g <= NGHOST; ++g) {
                int idx_ghost = mesh.cellIndex(i, jmin - g);
                U[idx_ghost] = U_ghost;
            }
        }
    }


 
    // =====================================================
    // ------------------------------------------------------------
    // 2) FARFIELD (top, j = jmax) — 4 classical cases
    //    based on the slide you sent
    // ------------------------------------------------------------
    void applyFarfieldBC(std::vector<Conservative>& U, const Mesh& mesh) {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        for (int i = imin; i < imax; ++i) {
            // interior state "d"
            int idx_int = mesh.cellIndex(i, jmax - 1);
            Primitive Wd(U[idx_int], cfg.gamma);

            // freestream "a"
            const Primitive& Wa = U_inf;

            // *** use the SAME face the flux uses: face at j = jmax ***
            int faceIdxTop = mesh.cellIndex(i, jmax);
            std::array<double,2> n = mesh.faceNormal_i[faceIdxTop];
            double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]);
            double nx = - n[0] / (mag + 1e-14);
            double ny = - n[1] / (mag + 1e-14);

            // interior sound speed and normal velocity
            double p_d   = Wd.p;
            double rho_d = Wd.rho;
            double c0    = std::sqrt(cfg.gamma * p_d / rho_d);
            double un_d  = Wd.u * nx + Wd.v * ny;

            bool inflow     = (un_d < 0.0);
            bool supersonic = (std::fabs(un_d) >= c0);

            double rho_b, u_b, v_b, p_b;

            if (inflow && supersonic) {
                // supersonic inflow: take freestream
                rho_b = Wa.rho;
                u_b   = Wa.u;
                v_b   = Wa.v;
                p_b   = Wa.p;
            }
            else if (!inflow && supersonic) {
                // supersonic outflow: take interior
                rho_b = rho_d;
                u_b   = Wd.u;
                v_b   = Wd.v;
                p_b   = p_d;
            }
            else if (inflow && !supersonic) {
                // subsonic inflow
                double pa = Wa.p;
                double pd = p_d;
                double corr = nx*(Wa.u - Wd.u) + ny*(Wa.v - Wd.v);
                p_b = 0.5 * (pa + pd - rho_d * c0 * corr);

                rho_b = Wa.rho + (p_b - pa)/(c0*c0);
                u_b   = Wa.u   - nx*(pa - p_b)/(rho_d * c0);
                v_b   = Wa.v   - ny*(pa - p_b)/(rho_d * c0);
            }
            else {
                // subsonic outflow
                double pa = Wa.p;
                double pd = p_d;
                p_b   = pa;
                rho_b = rho_d + (p_b - pd)/(c0*c0);
                u_b   = Wd.u + nx*(pd - p_b)/(rho_d * c0);
                v_b   = Wd.v + ny*(pd - p_b)/(rho_d * c0);
            }

            double E = p_b/((cfg.gamma - 1.0) * rho_b) + 0.5*(u_b*u_b + v_b*v_b);

            Conservative Ub;
            Ub.rho  = rho_b;
            Ub.rhou = rho_b * u_b;
            Ub.rhov = rho_b * v_b;
            Ub.rhoE = rho_b * E;

            for (int g = 0; g < NGHOST; ++g) {
                int jg = jmax + g;
                U[ mesh.cellIndex(i, jg) ] = Ub;
            }
        }
    }

    
    // =====================================================
    // PERIODIC BC (I-direction for O-grid)
    // =====================================================
    void applyPeriodicBC(std::vector<Conservative>& U, const Mesh& mesh) {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);
        
        // For all j levels (including ghosts)
        for (int j = 0; j < mesh.njTotal; ++j) {
            // Left ghosts receive from right interior
            for (int g = 1; g <= NGHOST; ++g) {
                int i_ghost = imin - g;
                int i_donor = imax - g;  // Wrap from right side
                int idx_ghost = mesh.cellIndex(i_ghost, j);
                int idx_donor = mesh.cellIndex(i_donor, j);
                U[idx_ghost] = U[idx_donor];
            }
            
            // Right ghosts receive from left interior
            for (int g = 0; g < NGHOST; ++g) {
                int i_ghost = imax + g;
                int i_donor = imin + g;  // Wrap from left side
                int idx_ghost = mesh.cellIndex(i_ghost, j);
                int idx_donor = mesh.cellIndex(i_donor, j);
                U[idx_ghost] = U[idx_donor];
            }
        }
    }
};

// ============================================================
// AIRFOIL FORCE CALCULATION
// ============================================================
class AirfoilForces {
private:
    const Config& cfg;
    double x_ref, y_ref;  // Reference point for moments (quarter-chord)
    
public:
    struct Coefficients {
        double CL = 0.0;  // Lift coefficient
        double CD = 0.0;  // Drag coefficient
        double CM = 0.0;  // Moment coefficient
    };
    
    AirfoilForces(const Config& cfg_) : cfg(cfg_), x_ref(0.0), y_ref(0.0) {}
    
    Coefficients compute(
        const std::vector<Conservative>& U,
        const Mesh& mesh,
        int jAirfoil  // j-index of airfoil surface
    ) {
        // Find reference point (quarter-chord)
        findQuarterChord(mesh, jAirfoil);
        
        Coefficients coef;
        double Fx = 0.0, Fy = 0.0, M_ref = 0.0;
        
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);
        
        // Reference quantities
        const double q_inf = 0.5 * cfg.rho_inf * cfg.Mach_inf * cfg.Mach_inf 
                           * cfg.a_inf * cfg.a_inf;
        const double S_ref = cfg.chord;
        
        // Integrate pressure forces along airfoil
        for (int i = imin; i < imax; ++i) {
            // Wall pressure (use interior cell value)
            int idx_wall = mesh.cellIndex(i, jAirfoil);
            double p_wall = EOS::pressure(U[idx_wall], cfg.gamma);
            
            // Face geometry (normal into fluid)
            //std::array<double,2> normal = mesh.faceNormal_i[idx_wall];
            std::array<double,2> normal = mesh.faceNormal_i_out[idx_wall];
            double ds = mesh.faceLen_i[idx_wall];
            
            // Normalize
            double mag = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
            double nx = normal[0] / (mag + 1e-14);
            double ny = normal[1] / (mag + 1e-14);
            
            // Pressure forces (wall pushes on fluid)
            double dFx = p_wall * nx * ds;
            double dFy = p_wall * ny * ds;
            
            // Accumulate forces
            Fx += dFx;
            Fy += dFy;
            
            // Moment about reference point
            double x_cp = mesh.x(i, jAirfoil);
            double y_cp = mesh.y(i, jAirfoil);
            double dx = x_cp - x_ref;
            double dy = y_cp - y_ref;
            M_ref += dx * dFy - dy * dFx;
        }
        
        // Transform to lift/drag axes
        const double alpha = cfg.alpha_rad;
        const double cos_a = std::cos(alpha);
        const double sin_a = std::sin(alpha);
        
        double Lift = -Fx * sin_a + Fy * cos_a;
        double Drag =  Fx * cos_a + Fy * sin_a;
        
        // Non-dimensionalize
        coef.CL = Lift / (q_inf * S_ref);
        coef.CD = Drag / (q_inf * S_ref);
        coef.CM = M_ref / (q_inf * S_ref * cfg.chord);
        
        return coef;
    }
    
private:
    void findQuarterChord(const Mesh& mesh, int jAirfoil) {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);
        
        // Find leading and trailing edge x-coordinates
        double x_min = 1e10, x_max = -1e10;
        for (int i = imin; i < imax; ++i) {
            double x = mesh.x(i, jAirfoil);
            x_min = std::min(x_min, x);
            x_max = std::max(x_max, x);
        }
        
        // Quarter-chord location
        x_ref = x_min + 0.25 * (x_max - x_min);
        y_ref = 0.0;
    }
};
