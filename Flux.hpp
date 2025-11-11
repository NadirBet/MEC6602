#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>
#include <iostream>

#include "Types.hpp"
#include "Mesh.hpp"

// ============================================================
// FIXED FLUX CALCULATOR FOR O-GRID
// ============================================================
class FluxCalculator {
private:
    double k2, k4;

public:
    FluxCalculator(double k2_ = 0.5, double k4_ = 0.02)
        : k2(k2_), k4(k4_) {}

    // ------------------------------------------------------------
    // main Euler residual
    // ------------------------------------------------------------
    void computeResidual(
        const std::vector<Conservative>& U,
        std::vector<Conservative>& R,
        const Mesh& mesh,
        const Config& cfg
    ) {
        // Zero residual
        for (auto& r : R) {
            r.rho = 0.0; r.rhou = 0.0; r.rhov = 0.0; r.rhoE = 0.0;
        }

        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);


        // ================= I-FACES =================
        // ================= I-FACES =================
        for (int j = jmin; j <= jmax; ++j) {
            for (int i = imin; i < imax; ++i) {

                if (j == jmin) {
                    int idInt   = mesh.cellIndex(i, jmin);     // interior
                    int idGhost = mesh.cellIndex(i, jmin - 1); // wall ghost you already set in BC

                    // face geometry: bottom i-face of the interior cell
                    int faceIdx = mesh.cellIndex(i, jmin);
                    const auto& n  = mesh.faceNormal_i_out[faceIdx]; // outward normal
                    double      ds = mesh.faceLen_i[faceIdx];

                    Primitive Wint(U[idInt],   cfg.gamma);
                    Primitive Wgst(U[idGhost], cfg.gamma);

                    // slide’s W̄_01 = 1/2 (W0 + W1)
                    Primitive Wbar;
                    Wbar.rho = 0.5 * (Wint.rho + Wgst.rho);
                    Wbar.u   = 0.5 * (Wint.u   + Wgst.u);
                    Wbar.v   = 0.5 * (Wint.v   + Wgst.v);
                    Wbar.p   = 0.5 * (Wint.p   + Wgst.p);

                    // Euler flux with averaged state
                    double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]) + 1e-14;
                    double nx  = n[0] / mag;
                    double ny  = n[1] / mag;

                    double un = Wbar.u * nx + Wbar.v * ny;
                    double H  = cfg.gamma/(cfg.gamma-1.0) * (Wbar.p / Wbar.rho) + 0.5*(Wbar.u*Wbar.u + Wbar.v*Wbar.v);

                    Conservative F;
                    F.rho  = (Wbar.rho * un) * ds;
                    F.rhou = (Wbar.rho * Wbar.u * un + Wbar.p * nx) * ds;
                    F.rhov = (Wbar.rho * Wbar.v * un + Wbar.p * ny) * ds;
                    F.rhoE = (Wbar.rho * H * un) * ds;

                    // add to interior cell
                    R[idInt].rho  += F.rho;
                    R[idInt].rhou += F.rhou;
                    R[idInt].rhov += F.rhov;
                    R[idInt].rhoE += F.rhoE;
                }
                else if (j < jmax) {
                    // interior face between (i,j-1) and (i,j)
                    int idL = mesh.cellIndex(i, j-1);
                    int idR = mesh.cellIndex(i, j);
                    int faceIdx = mesh.cellIndex(i, j);

                    Conservative F = computeInteriorFlux(
                        U[idL], U[idR],
                        mesh.faceNormal_i[faceIdx],
                        mesh.faceLen_i[faceIdx],
                        cfg,
                        /*alongI=*/true,
                        i, j, U, mesh
                    );

                    R[idL].rho  -= F.rho;
                    R[idL].rhou -= F.rhou;
                    R[idL].rhov -= F.rhov;
                    R[idL].rhoE -= F.rhoE;

                    R[idR].rho  += F.rho;
                    R[idR].rhou += F.rhou;
                    R[idR].rhov += F.rhov;
                    R[idR].rhoE += F.rhoE;
                }
                else { // j == jmax : farfield, but ghost already has the char BC
                    int idInt   = mesh.cellIndex(i, jmax-1); // interior
                    int idGhost = mesh.cellIndex(i, jmax);   // first ghost
                    int faceIdx = mesh.cellIndex(i, jmax);   // top face geom

                      // normale géométrique
                    //const auto& n0 = mesh.faceNormal_i[faceIdx];
                    //std::array<double,2> n = { -n0[0], -n0[1] };  // ← on la retourne pour matcher le BC
                    // use the geometric face normal as is
                    const auto& n = mesh.faceNormal_i[faceIdx];
                    double ds     = mesh.faceLen_i[faceIdx];
                    Conservative F = computeInteriorFlux(
                        U[idInt], U[idGhost],
                        n,
                        mesh.faceLen_i[faceIdx],
                        cfg,
                        /*alongI=*/true,
                        i, jmax, U, mesh
                    );

                    // flux leaves the interior
                    R[idInt].rho  -= F.rho;
                    R[idInt].rhou -= F.rhou;
                    R[idInt].rhov -= F.rhov;
                    R[idInt].rhoE -= F.rhoE;
                }
            }
        }

        // ================= J-FACES =================


        
        // ================= J-FACES =================
        for (int j = jmin; j < jmax; ++j) {
            // interior faces between (i-1,j) and (i,j)
            for (int i = imin + 1; i < imax; ++i) {
                int idL = mesh.cellIndex(i - 1, j);  // left cell
                int idR = mesh.cellIndex(i,     j);  // right cell

                int faceIdx = mesh.cellIndex(i, j);

                Conservative F = computeInteriorFlux(
                    U[idL], U[idR],
                    mesh.faceNormal_j[faceIdx],
                    mesh.faceLen_j[faceIdx],
                    cfg,
                    /*alongI=*/false,
                    i, j, U, mesh
                );

                // J-faces interior
                R[idL].rho  -= F.rho;
                R[idL].rhou -= F.rhou;
                R[idL].rhov -= F.rhov;
                R[idL].rhoE -= F.rhoE;

                R[idR].rho  += F.rho;
                R[idR].rhou += F.rhou;
                R[idR].rhov += F.rhov;
                R[idR].rhoE += F.rhoE;
            }

            // periodic closure between (imax-1,j) [left] and (imin,j) [right]
            {
                int idL = mesh.cellIndex(imax - 1, j);
                int idR = mesh.cellIndex(imin,     j);

                int faceIdx = mesh.cellIndex(imin, j);

                Conservative F = computeInteriorFlux(
                    U[idL], U[idR],
                    mesh.faceNormal_j[faceIdx],
                    mesh.faceLen_j[faceIdx],
                    cfg,
                    /*alongI=*/false,
                    imax, j, U, mesh
                );

                // J-faces interior
                R[idL].rho  -= F.rho;
                R[idL].rhou -= F.rhou;
                R[idL].rhov -= F.rhov;
                R[idL].rhoE -= F.rhoE;

                R[idR].rho  += F.rho;
                R[idR].rhou += F.rhou;
                R[idR].rhov += F.rhov;
                R[idR].rhoE += F.rhoE;
            }
        }
    }

private:


    // ------------------------------------------------------------
    // interior flux + JST
    // ------------------------------------------------------------
    Conservative computeInteriorFlux(
        const Conservative& UL,
        const Conservative& UR,
        const std::array<double,2>& normal,
        double ds,
        const Config& cfg,
        bool alongI,
        int i, int j,
        const std::vector<Conservative>& U,
        const Mesh& mesh
    ) const {
        double mag = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
        double nx = normal[0] / (mag + 1e-14);
        double ny = normal[1] / (mag + 1e-14);

        Primitive WL(UL, cfg.gamma);
        Primitive WR(UR, cfg.gamma);

        double unL = WL.u * nx + WL.v * ny;
        double unR = WR.u * nx + WR.v * ny;

        double HL = EOS::totalEnthalpy(UL, cfg.gamma);
        double HR = EOS::totalEnthalpy(UR, cfg.gamma);

        Conservative F;
        F.rho  = 0.5 * (WL.rho * unL + WR.rho * unR) * ds;
        F.rhou = 0.5 * (WL.rho*WL.u*unL + WL.p*nx +
                        WR.rho*WR.u*unR + WR.p*nx) * ds;
        F.rhov = 0.5 * (WL.rho*WL.v*unL + WL.p*ny +
                        WR.rho*WR.v*unR + WR.p*ny) * ds;
        F.rhoE = 0.5 * (WL.rho*HL*unL + WR.rho*HR*unR) * ds;

        addJSTDissipation(F, UL, UR, WL, WR, nx, ny, ds,
                          alongI, i, j, U, mesh, cfg);

        return F;
    }


// ------------------------------------------------------------
// JST with pressure sensor and boundary-safe fallback
// ------------------------------------------------------------
    void addJSTDissipation(
        Conservative& F,
        const Conservative& UL,
        const Conservative& UR,
        const Primitive& WL,
        const Primitive& WR,
        double nx, double ny, double ds,
        bool alongI, int i, int j,
        const std::vector<Conservative>& U,
        const Mesh& mesh,
        const Config& cfg
    ) const
    {
        // 1) face spectral radius
        double unL = WL.u * nx + WL.v * ny;
        double unR = WR.u * nx + WR.v * ny;
        double lambda = std::max(
            std::abs(unL) + WL.a,
            std::abs(unR) + WR.a
        );

        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        // 2) pressure sensor around this face
        auto cellPressure = [&](int ii, int jj) -> double {
            int idx = mesh.cellIndex(ii, jj);
            return EOS::pressure(U[idx], cfg.gamma);
        };

        auto sensor3 = [&](double pm1, double p0, double pp1) -> double {
            double num = std::fabs(pp1 - 2.0 * p0 + pm1);
            double den = (pp1 + 2.0 * p0 + pm1) + 1e-14;
            return num / den;
        };

        double Y_face = 0.0;
        if (alongI) {
            // face between (i, j-1) and (i, j)
            double Y_lower = 0.0;
            double Y_upper = 0.0;

            // lower cell sensor: needs j-2, j-1, j
            if (j >= jmin + 2 && j <= jmax) {
                Y_lower = sensor3(
                    cellPressure(i, j-2),
                    cellPressure(i, j-1),
                    cellPressure(i, j)
                );
            }
            // upper cell sensor: needs j-1, j, j+1
            if (j >= jmin + 1 && j <= jmax - 1) {
                Y_upper = sensor3(
                    cellPressure(i, j-1),
                    cellPressure(i, j),
                    cellPressure(i, j+1)
                );
            }
            Y_face = std::max(Y_lower, Y_upper);
        } else {
            // j-face between (i-1, j) and (i, j) with periodic in i
            int iL   = (i == imin) ? (imax - 1) : (i - 1);
            int iR   = (i == imax) ? imin : i;
            int iLm1 = (iL == imin) ? (imax - 1) : (iL - 1);
            int iRp1 = (iR == imax - 1) ? imin : (iR + 1);

            double Y_left = sensor3(
                cellPressure(iLm1, j),
                cellPressure(iL,  j),
                cellPressure(iR,  j)
            );
            double Y_right = sensor3(
                cellPressure(iL,  j),
                cellPressure(iR,  j),
                cellPressure(iRp1,j)
            );
            Y_face = std::max(Y_left, Y_right);
        }

        // 3) base variable coefficients from the slide
        double eps2 = k2 * Y_face;
        double eps4 = std::max(0.0, k4 - eps2);

        // 4) check if 4th-order stencil is actually available
        bool hasD4 = false;
        if (alongI) {
            // need j-2 and j+1
            hasD4 = (j > jmin + 1) && (j < jmax - 1);
        } else {
            // in j-direction we have periodic wrap → ok
            hasD4 = true;
        }

        // 5) boundary-safe fallback:
        // if we can't do 4th order, force a plain 2nd-order JST
        if (!hasD4) {
            eps2 = k2;
            eps4 = 0.0;
        }

        // 6) 2nd-order dissipation
        if (eps2 > 0.0) {
            Conservative D2;
            D2.rho  = eps2 * lambda * (UR.rho  - UL.rho)  * ds;
            D2.rhou = eps2 * lambda * (UR.rhou - UL.rhou) * ds;
            D2.rhov = eps2 * lambda * (UR.rhov - UL.rhov) * ds;
            D2.rhoE = eps2 * lambda * (UR.rhoE - UL.rhoE) * ds;

            F.rho  -= D2.rho;
            F.rhou -= D2.rhou;
            F.rhov -= D2.rhov;
            F.rhoE -= D2.rhoE;
        }

        // 7) 4th-order dissipation
        if (eps4 > 0.0 && hasD4) {
            Conservative D4;
            D4.rho = D4.rhou = D4.rhov = D4.rhoE = 0.0;

            if (alongI) {
                // i fixed, vary j
                int idLm1 = mesh.cellIndex(i, j-2);
                int idRp1 = mesh.cellIndex(i, j+1);
                const Conservative& ULm1 = U[idLm1];
                const Conservative& URp1 = U[idRp1];

                D4.rho  = eps4 * lambda *
                    (URp1.rho  - 3.0*UR.rho  + 3.0*UL.rho  - ULm1.rho) * ds;
                D4.rhou = eps4 * lambda *
                    (URp1.rhou - 3.0*UR.rhou + 3.0*UL.rhou - ULm1.rhou) * ds;
                D4.rhov = eps4 * lambda *
                    (URp1.rhov - 3.0*UR.rhov + 3.0*UL.rhov - ULm1.rhov) * ds;
                D4.rhoE = eps4 * lambda *
                    (URp1.rhoE - 3.0*UR.rhoE + 3.0*UL.rhoE - ULm1.rhoE) * ds;
            } else {
                // j fixed, vary i (periodic)
                int iL   = (i == imin) ? (imax - 1) : (i - 1);
                int iLm1 = (iL == imin) ? (imax - 1) : (iL - 1);
                int iR   = (i == imax) ? imin : i;
                int iRp1 = (iR == imax - 1) ? imin : (iR + 1);

                int idLm1 = mesh.cellIndex(iLm1, j);
                int idRp1 = mesh.cellIndex(iRp1, j);

                const Conservative& ULm1 = U[idLm1];
                const Conservative& URp1 = U[idRp1];

                D4.rho  = eps4 * lambda *
                    (URp1.rho  - 3.0*UR.rho  + 3.0*UL.rho  - ULm1.rho) * ds;
                D4.rhou = eps4 * lambda *
                    (URp1.rhou - 3.0*UR.rhou + 3.0*UL.rhou - ULm1.rhou) * ds;
                D4.rhov = eps4 * lambda *
                    (URp1.rhov - 3.0*UR.rhov + 3.0*UL.rhov - ULm1.rhov) * ds;
                D4.rhoE = eps4 * lambda *
                    (URp1.rhoE - 3.0*UR.rhoE + 3.0*UL.rhoE - ULm1.rhoE) * ds;
            }

            F.rho  += D4.rho;
            F.rhou += D4.rhou;
            F.rhov += D4.rhov;
            F.rhoE += D4.rhoE;
        }
    }



















public:
    // Debug helper: print flux and residual contributions for one internal face
    void debugInteriorFaceFlux(
        const std::vector<Conservative>& U,
        const Mesh& mesh,
        const Config& cfg,
        int i, int j,
        bool alongI   // true = i-face (vertical), false = j-face (horizontal)
    ) const {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        int idL, idR;
        std::array<double,2> n;
        double ds;

        int i_face = i;
        int j_face = j;

        if (alongI) {
            // I-face between (i, j-1) [L] and (i, j) [R]
            j_face = std::min(std::max(j, jmin + 1), jmax - 1);
            i_face = std::min(std::max(i, imin),     imax - 1);

            idL = mesh.cellIndex(i_face,     j_face - 1);  // lower
            idR = mesh.cellIndex(i_face,     j_face);      // upper

            int faceIdx = mesh.cellIndex(i_face, j_face);  // bottom of upper cell
            n  = mesh.faceNormal_i[faceIdx];
            ds = mesh.faceLen_i[faceIdx];
        } else {
            // J-face between (i-1, j) [L] and (i, j) [R]
            i_face = std::min(std::max(i, imin + 1), imax - 1);
            j_face = std::min(std::max(j, jmin),     jmax - 1);

            idL = mesh.cellIndex(i_face - 1, j_face);     // left
            idR = mesh.cellIndex(i_face,     j_face);     // right

            int faceIdx = mesh.cellIndex(i_face, j_face); // left face of right cell
            n  = mesh.faceNormal_j[faceIdx];
            ds = mesh.faceLen_j[faceIdx];
        }

        const Conservative& UL = U[idL];
        const Conservative& UR = U[idR];

        double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]);
        double nx  = n[0] / (mag + 1e-14);
        double ny  = n[1] / (mag + 1e-14);

        Primitive WL(UL, cfg.gamma);
        Primitive WR(UR, cfg.gamma);

        double unL = WL.u * nx + WL.v * ny;
        double unR = WR.u * nx + WR.v * ny;

        // Use the actual scheme (central + JST)
        Conservative F = computeInteriorFlux(
            UL, UR, n, ds, cfg,
            alongI,
            i_face, j_face,
            U, mesh
        );

        Conservative dR_L, dR_R;
        dR_L.rho  = -F.rho;   dR_R.rho  = +F.rho;
        dR_L.rhou = -F.rhou;  dR_R.rhou = +F.rhou;
        dR_L.rhov = -F.rhov;  dR_R.rhov = +F.rhov;
        dR_L.rhoE = -F.rhoE;  dR_R.rhoE = +F.rhoE;

        std::cout << "\n=== INTERNAL FACE FLUX DEBUG ===\n";
        std::cout << "alongI = " << (alongI ? "true (i-face)" : "false (j-face)") << "\n";
        std::cout << "Face indices (i_face,j_face) = (" << i_face << "," << j_face << ")\n";

        std::cout << "Left cell id=" << idL << "  Right cell id=" << idR << "\n";

        std::cout << "Normal n = (" << nx << ", " << ny << "), ds = " << ds << "\n";
        std::cout << "unL = " << unL << ", unR = " << unR << "\n";

        std::cout << "UL = (rho="  << UL.rho
                << ", rhou="     << UL.rhou
                << ", rhov="     << UL.rhov
                << ", rhoE="     << UL.rhoE << ")\n";
        std::cout << "UR = (rho="  << UR.rho
                << ", rhou="     << UR.rhou
                << ", rhov="     << UR.rhov
                << ", rhoE="     << UR.rhoE << ")\n";

        std::cout << "Flux F = ("
                << "rho="  << F.rho  << ", "
                << "rhou=" << F.rhou << ", "
                << "rhov=" << F.rhov << ", "
                << "rhoE=" << F.rhoE << ")\n";

        std::cout << "dR_L = (-F) = ("
                << dR_L.rho  << ", "
                << dR_L.rhou << ", "
                << dR_L.rhov << ", "
                << dR_L.rhoE << ")\n";
        std::cout << "dR_R = (+F) = ("
                << dR_R.rho  << ", "
                << dR_R.rhou << ", "
                << dR_R.rhov << ", "
                << dR_R.rhoE << ")\n";

        std::cout << "dR_L + dR_R (should be ~0) = ("
                << dR_L.rho  + dR_R.rho  << ", "
                << dR_L.rhou + dR_R.rhou << ", "
                << dR_L.rhov + dR_R.rhov << ", "
                << dR_L.rhoE + dR_R.rhoE << ")\n";
        std::cout << "================================\n";
    }

    // ============================================================
    // METRIC / GEOMETRY CHECK IMPLEMENTATION
    // ============================================================

    void computeMetricResidual(
        std::vector<Conservative>& R,
        const Mesh& mesh
    ) const {
        // Zero out all residuals first
        for (auto& r : R) {
            r.rho = r.rhou = r.rhov = r.rhoE = 0.0;
        }

        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        // Now we need to match EXACTLY what the mesh geometry does:
        // Each cell owns its bottom i-face and left j-face
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int c = mesh.cellIndex(i, j);
                
                // The key insight: we need to gather the four faces that bound this cell
                // using the EXACT ownership pattern from computeCellAreasAndFaces()
                
                // Bottom i-face: owned by this cell
                const auto& n_bottom = mesh.faceNormal_i[c];
                double ds_bottom = mesh.faceLen_i[c];
                
                // Top i-face: owned by the cell above
                int c_above = mesh.cellIndex(i, j+1);
                const auto& n_top = mesh.faceNormal_i[c_above];
                double ds_top = mesh.faceLen_i[c_above];
                
                // Left j-face: owned by this cell
                const auto& n_left = mesh.faceNormal_j[c];
                double ds_left = mesh.faceLen_j[c];
                
                // Right j-face: owned by the cell to the right
                int c_right = mesh.cellIndex(i+1, j);
                const auto& n_right = mesh.faceNormal_j[c_right];
                double ds_right = mesh.faceLen_j[c_right];
                
                // Now sum them exactly as the geometry module does
                R[c].rhou = n_bottom[0]*ds_bottom + n_top[0]*ds_top + 
                        n_left[0]*ds_left + n_right[0]*ds_right;
                R[c].rhov = n_bottom[1]*ds_bottom + n_top[1]*ds_top + 
                        n_left[1]*ds_left + n_right[1]*ds_right;
            }
        }
    }


// ------------------------------------------------------------
// METRIC FLUX TEST: use *flux loops* with F = n * ds
//
// This MUST mirror computeResidual's face loops, but instead of
// physical fluxes we inject a purely geometric flux:
//
//   F.rho  = 0
//   F.rhou = n_x * ds
//   F.rhov = n_y * ds
//   F.rhoE = 0
//
// If geometry + residual assembly are consistent, then for a closed
// surface the per–cell sum of F over its faces must be ≈ 0, so
// R.rhou, R.rhov should be ~1e-12 everywhere in the interior.
// ------------------------------------------------------------
    void computeMetricResidualFluxAssembly(
        std::vector<Conservative>& R,
        const Mesh& mesh
    ) const {
        // zero out
        for (auto& r : R) {
            r.rho = r.rhou = r.rhov = r.rhoE = 0.0;
        }

        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        Conservative F;

        // ================= I-FACES (vertical, along i) =================
        // ================= I-FACES (vertical, along i) =================
        for (int j = jmin; j <= jmax; ++j) {
            for (int i = imin; i < imax; ++i) {

                if (j == jmin) {
                    // bottom boundary: interface ghost ↔ interior
                    int c = mesh.cellIndex(i, jmin);

                    // use the SAME normal as the wall flux uses
                    const auto& n_out = mesh.faceNormal_i_out[c];
                    double ds = mesh.faceLen_i[c];

                    double mag = std::sqrt(n_out[0]*n_out[0] + n_out[1]*n_out[1]) + 1e-14;
                    double nx  = n_out[0] / mag;
                    double ny  = n_out[1] / mag;

                    F.rho  = 0.0;
                    F.rhou = nx * ds;
                    F.rhov = ny * ds;
                    F.rhoE = 0.0;

                    // same assembly as computeResidual's wall branch
                    R[c].rhou += F.rhou;
                    R[c].rhov += F.rhov;
                }
                else if (j < jmax) {
                    // interior face between (i,j-1) [lower] and (i,j) [upper]
                    int idL = mesh.cellIndex(i, j-1);
                    int idR = mesh.cellIndex(i, j);

                    // geometry of this face: bottom i-face of the *upper* cell (i,j)
                    int faceIdx   = mesh.cellIndex(i, j);
                    const auto& n = mesh.faceNormal_i[faceIdx];
                    double ds     = mesh.faceLen_i[faceIdx];

                    double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]) + 1e-14;
                    double nx  = n[0] / mag;
                    double ny  = n[1] / mag;

                    F.rho  = 0.0;
                    F.rhou = nx * ds;
                    F.rhov = ny * ds;
                    F.rhoE = 0.0;

                    // same pattern as computeResidual:
                    R[idL].rhou -= F.rhou;
                    R[idL].rhov -= F.rhov;

                    R[idR].rhou += F.rhou;
                    R[idR].rhov += F.rhov;
                }
                else { // j == jmax : farfield boundary band
                    int idInt   = mesh.cellIndex(i, jmax-1);
                    int faceIdx = mesh.cellIndex(i, jmax);

                    // géométrique
                    const auto& n0 = mesh.faceNormal_i[faceIdx];
                    // mais on la retourne comme dans le flux
                    double mag = std::sqrt(n0[0]*n0[0] + n0[1]*n0[1]) + 1e-14;
                    double nx  = -n0[0] / mag;
                    double ny  = -n0[1] / mag;
                    double ds  = mesh.faceLen_i[faceIdx];

                    F.rho  = 0.0;
                    F.rhou = nx * ds;
                    F.rhov = ny * ds;
                    F.rhoE = 0.0;

                    R[idInt].rhou += F.rhou;
                    R[idInt].rhov += F.rhov;
                }
            }
        }

        // ================= J-FACES (horizontal, along j) =================
        for (int j = jmin; j < jmax; ++j) {
            // interior faces i = imin+1..imax-1
            for (int i = imin + 1; i < imax; ++i) {
                int idL = mesh.cellIndex(i - 1, j);  // left cell
                int idR = mesh.cellIndex(i,     j);  // right cell

                int faceIdx   = mesh.cellIndex(i, j);
                const auto& n = mesh.faceNormal_j[faceIdx];
                double      ds = mesh.faceLen_j[faceIdx];

                double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]) + 1e-14;
                double nx  = n[0] / mag;
                double ny  = n[1] / mag;

                F.rho  = 0.0;
                F.rhou = nx * ds;
                F.rhov = ny * ds;
                F.rhoE = 0.0;

                // *** flip signs to match computeResidual ***

                R[idL].rhou -= F.rhou;
                R[idL].rhov -= F.rhov;
                R[idR].rhou += F.rhou;
                R[idR].rhov += F.rhov;
            }

            // periodic closure ...
            {
                int idL = mesh.cellIndex(imax - 1, j);
                int idR = mesh.cellIndex(imin,     j);

                int faceIdx   = mesh.cellIndex(imin, j);
                const auto& n = mesh.faceNormal_j[faceIdx];
                double      ds = mesh.faceLen_j[faceIdx];

                double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]) + 1e-14;
                double nx  = n[0] / mag;
                double ny  = n[1] / mag;

                F.rho  = 0.0;
                F.rhou = nx * ds;
                F.rhov = ny * ds;
                F.rhoE = 0.0;
                //  match computeResidual
                R[idL].rhou -= F.rhou;
                R[idL].rhov -= F.rhov;
                R[idR].rhou += F.rhou;
                R[idR].rhov += F.rhov;
            }
        }
    }

};
