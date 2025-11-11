#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include "Types.hpp"
#include "Mesh.hpp"
#include "Flux.hpp"

// ============================================================================
// COMPREHENSIVE SOLVER DIAGNOSTICS (O-grid friendly)
// ============================================================================

class SolverDiagnostics {
private:
    const Mesh& mesh;
    const Config& cfg;
    std::ofstream log_file;
    bool geometry_validated;
    
public:
    SolverDiagnostics(const Mesh& mesh_, const Config& cfg_) 
        : mesh(mesh_), cfg(cfg_), geometry_validated(false) {
        log_file.open("diagnostics_log.txt");
        log_file << "# Iteration,TotalMass,SumR_rho,SumR_rhou,SumR_rhov,SumR_rhoE,"
                 << "MaxWallVn,MinP,MaxP,MinRho,MaxRho,L2_rho,WallFluxY\n";
    }
    
    ~SolverDiagnostics() {
        if (log_file.is_open()) log_file.close();
    }
    
    // ========================================================================
    // VALIDATE GEOMETRY (run once after mesh setup)
    // ========================================================================
    void validateGeometry() {
        if (geometry_validated) return;
        
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "GEOMETRIC CONSISTENCY VALIDATION\n";
        std::cout << std::string(70, '=') << "\n\n";
        
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        // ------------------------------------------------------------
        // 0) infer center from inner row
        // ------------------------------------------------------------
        double xc = 0.0, yc = 0.0;
        {
            int cnt = 0;
            for (int i = imin; i < imax; ++i) {
                xc += mesh.x(i, jmin);
                yc += mesh.y(i, jmin);
                ++cnt;
            }
            xc /= std::max(1, cnt);
            yc /= std::max(1, cnt);
        }
        std::cout << "Inferred center (from j=jmin):  xc=" << std::setprecision(6) << xc
                  << ", yc=" << yc << "\n\n";
    










        // ------------------------------------------------------------
        // Test 1: i-face normals outward
        // ------------------------------------------------------------
        int i_fail = 0;
        std::cout << "Test 1: I-face normal directions (should be outward w.r.t. center)\n";
        
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);
                double rx = mesh.x(i, j) - xc;
                double ry = mesh.y(i, j) - yc;
                const auto& n = mesh.faceNormal_i_out[idx];
                double dot = n[0]*rx + n[1]*ry;
                if (dot < 0.0) {
                    ++i_fail;
                    if (i_fail < 5) {
                        std::cout << "    ✗ Cell("<<i<<","<<j<<") n·r="<<dot<<" (not outward)\n";
                    }
                }
            }
        }
        if (i_fail == 0) {
            std::cout << "  ✓ All i-face normals are outward (radial)\n";
        } else {
            std::cout << "  ✗ " << i_fail << " i-face normals point inward\n";
        }

        // ------------------------------------------------------------
        // Test 2: j-face normals not degenerate
        // ------------------------------------------------------------
        int j_fail = 0;
        std::cout << "\nTest 2: J-face normal directions (O-grid sanity)\n";
        int i_test = imin + (imax - imin) / 2;
        for (int j = jmin; j < jmax; ++j) {
            int idx = mesh.cellIndex(i_test, j);
            double nx = mesh.faceNormal_j[idx][0];
            double ny = mesh.faceNormal_j[idx][1];
            double mag = std::sqrt(nx*nx + ny*ny);
            if (mag < 1e-12) {
                ++j_fail;
                if (j_fail < 5) {
                    std::cout << "    ✗ Cell("<<i_test<<","<<j<<") |n_j|=0 (degenerate)\n";
                }
            }
        }
        if (j_fail == 0) {
            std::cout << "  ✓ All checked j-face normals are non-degenerate\n";
        } else {
            std::cout << "  ✗ " << j_fail << " j-face normals have issues\n";
        }

        // ------------------------------------------------------------
        // Test 3: wall normals outward
        // ------------------------------------------------------------
        std::cout << "\nTest 3: Wall boundary normals (should be outward from center)\n";
        int wall_fail = 0;
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, jmin);
            double rx = mesh.x(i, jmin) - xc;
            double ry = mesh.y(i, jmin) - yc;
            const auto& n = mesh.faceNormal_i_out[idx];
            double dot = n[0]*rx + n[1]*ry;
            if (dot < 0.0) {
                ++wall_fail;
                if (wall_fail < 5) {
                    std::cout << "    ✗ Wall i="<<i<<": n·r="<<dot<<" (not outward)\n";
                }
            }
        }
        if (wall_fail == 0) {
            std::cout << "  ✓ All wall normals are outward\n";
        } else {
            std::cout << "  ✗ " << wall_fail << " wall normals are inward\n";
            std::cout << "  THIS WILL CAUSE WRONG SIGN WALL FORCES!\n";
        }

        // ------------------------------------------------------------
        // Test 4: cell areas + richest print
        // ------------------------------------------------------------
        std::cout << "\nTest 4: Cell areas\n";
        int area_fail = 0;
        double min_area = 1e10, max_area = -1e10;
        int min_i = imin, min_j = jmin, max_i = imin, max_j = jmin;
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);
                double area = mesh.cellArea[idx];
                if (area <= 0) {
                    ++area_fail;
                    if (area_fail < 5) {
                        std::cout << "    ✗ Cell("<<i<<","<<j<<") area="<<area<<"\n";
                    }
                }
                if (area < min_area) { min_area = area; min_i = i; min_j = j; }
                if (area > max_area) { max_area = area; max_i = i; max_j = j; }
            }
        }
        if (area_fail == 0) {
            std::cout << "  ✓ All cell areas are positive\n";
            std::cout << "  Area range: [" << std::scientific << std::setprecision(3) 
                      << min_area << " at ("<<min_i<<","<<min_j<<"), "
                      << max_area << " at ("<<max_i<<","<<max_j<<")]\n";
        } else {
            std::cout << "  ✗ " << area_fail << " cells have non-positive area\n";
        }

        // ------------------------------------------------------------
        // Test 5: periodic consistency + sample pairs
        // ------------------------------------------------------------
        std::cout << "\nTest 5: Periodic boundary consistency\n";
        bool periodic_ok = true;
        for (int j = jmin; j < jmax; ++j) {
            double xL = mesh.x(imin,   j);
            double yL = mesh.y(imin,   j);
            double xR = mesh.x(imax-1, j);
            double yR = mesh.y(imax-1, j);
            double dist = std::sqrt((xR-xL)*(xR-xL) + (yR-yL)*(yR-yL));
            double dx = mesh.x(imin+1, j) - mesh.x(imin, j);
            double dy = mesh.y(imin+1, j) - mesh.y(imin, j);
            double h  = std::sqrt(dx*dx + dy*dy);
            if (dist > 2.0*h) {
                periodic_ok = false;
                if (j == jmin) {
                    std::cout << "    ✗ Large gap at j="<<j<<": dist="<<dist<<", h="<<h<<"\n";
                }
            }
        }
        if (periodic_ok) {
            std::cout << "  ✓ Periodic boundaries properly connected\n";
        } else {
            std::cout << "  ✗ Periodic boundaries have gaps\n";
        }

        // ------------------------------------------------------------
        // EXTRA: dump a few wall and farfield entries so we can eyeball them
        // ------------------------------------------------------------
        std::cout << "\nSample wall entries (j=jmin):\n";
        {
            int printed = 0;
            for (int i = imin; i < imax && printed < 5; ++i, ++printed) {
                int idx = mesh.cellIndex(i, jmin);
                const auto& n = mesh.faceNormal_i[idx];
                double ds = mesh.faceLen_i[idx];
                double rx = mesh.x(i, jmin) - xc;
                double ry = mesh.y(i, jmin) - yc;
                double dot = n[0]*rx + n[1]*ry;
                std::cout << "  i="<<i
                          << "  x="<<mesh.x(i,jmin)
                          << "  y="<<mesh.y(i,jmin)
                          << "  n=("<<n[0]<<","<<n[1]<<")"
                          << "  |n|="<<std::sqrt(n[0]*n[0]+n[1]*n[1])
                          << "  ds="<<ds
                          << "  n·r="<<dot << "\n";
            }
        }


        std::cout << "\nSample farfield entries (j=jmax):\n";
        {
            int printed = 0;
            for (int i = imin; i < imax && printed < 5; ++i, ++printed) {
                int idx = mesh.cellIndex(i, jmax);  // <-- face at ghost row
                const auto& n = mesh.faceNormal_i[idx];
                double ds = mesh.faceLen_i[idx];
                std::cout << "  i="<<i
                        << "  x="<<mesh.x(i,jmax)
                        << "  y="<<mesh.y(i,jmax)
                        << "  n=("<<n[0]<<","<<n[1]<<")"
                        << "  |n|="<<std::sqrt(n[0]*n[0]+n[1]*n[1])
                        << "  ds="<<ds << "\n";
            }
        }

        std::cout << "\nSample periodic pairs (left/right):\n";
        {
            int printed = 0;
            for (int j = jmin; j < jmax && printed < 3; ++j, ++printed) {
                double xL = mesh.x(imin,   j);
                double yL = mesh.y(imin,   j);
                double xR = mesh.x(imax-1, j);
                double yR = mesh.y(imax-1, j);
                double dist = std::sqrt((xR-xL)*(xR-xL) + (yR-yL)*(yR-yL));
                std::cout << "  j="<<j<<"  L=("<<xL<<","<<yL<<")  R=("<<xR<<","<<yR
                          << ")  dist="<<dist << "\n";
            }
        }

        // ------------------------------------------------------------
        // NEW: angles + ds ranges
        // ------------------------------------------------------------
        std::cout << "\nWall normal angles (deg):\n";
        {
            int printed = 0;
            for (int i = imin; i < imax && printed < 8; ++i, ++printed) {
                int idx = mesh.cellIndex(i, jmin);
                const auto& n = mesh.faceNormal_i[idx];
                double ang = std::atan2(n[1], n[0]) * 180.0 / M_PI;
                std::cout << "  i=" << i << "  angle=" << ang << "°\n";
            }
        }

        std::cout << "\nFarfield normal angles (deg):\n";
        {
            int printed = 0;
            for (int i = imin; i < imax && printed < 8; ++i, ++printed) {
                int idx = mesh.cellIndex(i, jmax);
                const auto& n = mesh.faceNormal_i[idx];
                double ang = std::atan2(n[1], n[0]) * 180.0 / M_PI;
                std::cout << "  i=" << i << "  angle=" << ang << "°\n";
            }
        }

        // Face-length summary
        double min_wall_ds = 1e30, max_wall_ds = -1e30;
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, jmin);
            min_wall_ds = std::min(min_wall_ds, mesh.faceLen_i[idx]);
            max_wall_ds = std::max(max_wall_ds, mesh.faceLen_i[idx]);
        }
        double min_ff_ds = 1e30, max_ff_ds = -1e30;
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, jmax-1);
            min_ff_ds = std::min(min_ff_ds, mesh.faceLen_i[idx]);
            max_ff_ds = std::max(max_ff_ds, mesh.faceLen_i[idx]);
        }
        std::cout << "\nWall ds range:    [" << min_wall_ds << ", " << max_wall_ds << "]\n";
        std::cout << "Farfield ds range: [" << min_ff_ds << ", " << max_ff_ds << "]\n";


        checkCellFaceClosure();
        std::cout << "\n" << std::string(70, '=') << "\n\n";
        
        geometry_validated = true;
    }






    // ========================================================================
    // TEST 6: CELL FACE CLOSURE (METRIC IDENTITY)
    //
    // For each interior cell, we form the sum of its four face area vectors
    // using the stored face normals + lengths. If the geometry is metric-
    // consistent (i.e. faces are built from a closed polygon), we expect
    //   Σ (n̂ * ds) ≈ 0
    // for each cell. Large values indicate a geometric inconsistency.
    //
    // NOTE:
    //  - We skip the wall band (jmin) and farfield band (jmax-1)
    //  - We also skip the first/last interior i to avoid periodic index
    //    gymnastics; this still samples the whole "bulk" of the domain.
    // ========================================================================



    void checkCellFaceClosure() const {
            int imin, imax, jmin, jmax;
            mesh.getInteriorBounds(imin, imax, jmin, jmax);

            double max_mag = 0.0;
            int i_max = -1, j_max = -1;

            // Need (i+1,j) and (i,j+1) to exist as interior cells:
            int i_start = imin;
            int i_end   = imax - 2;   // because right face uses i+1
            int j_start = jmin;
            int j_end   = jmax - 2;   // because top face uses j+1

            for (int j = j_start; j <= j_end; ++j) {
                for (int i = i_start; i <= i_end; ++i) {

                    double Sx = 0.0;
                    double Sy = 0.0;

                    auto addFace = [&](const std::array<double,2>& n,
                                    double ds,
                                    double sign) {
                        double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]) + 1e-14;
                        double nx  = n[0] / mag;
                        double ny  = n[1] / mag;
                        Sx += sign * nx * ds;
                        Sy += sign * ny * ds;
                    };

                    // Bottom face: stored at (i, j), outward = +n
                    {
                        int idxB = mesh.cellIndex(i, j);
                        addFace(mesh.faceNormal_i[idxB], mesh.faceLen_i[idxB], +1.0);
                    }

                    // Top face: stored at (i, j+1), outward for cell (i,j) is opposite of that
                    {
                        int idxT = mesh.cellIndex(i, j+1);
                        addFace(mesh.faceNormal_i[idxT], mesh.faceLen_i[idxT], -1.0);
                    }

                    // Left face: stored at (i, j), outward for cell (i,j) is opposite of stored
                    {
                        int idxL = mesh.cellIndex(i, j);
                        addFace(mesh.faceNormal_j[idxL], mesh.faceLen_j[idxL], -1.0);
                    }

                    // Right face: stored at (i+1, j), outward for cell (i,j) is same as stored
                    {
                        int idxR = mesh.cellIndex(i+1, j);
                        addFace(mesh.faceNormal_j[idxR], mesh.faceLen_j[idxR], +1.0);
                    }

                    double magS = std::sqrt(Sx*Sx + Sy*Sy);
                    if (magS > max_mag) {
                        max_mag = magS;
                        i_max = i;
                        j_max = j;
                    }
                }
            }

            std::cout << "\nTest 6: Cell face-closure (metric identity)\n";
            if (i_max >= 0) {
                std::cout << "  max |Σ n̂ ds| = " << std::scientific << std::setprecision(4)
                        << max_mag << "  at cell (i,j) = ("
                        << i_max << "," << j_max << ")\n";
                std::cout << "  If this is O(1) instead of ~1e-10, face geometry is NOT\n"
                        << "  consistent with a closed polygon per cell and freestream\n"
                        << "  cannot be preserved.\n";

                // Also dump detailed breakdown for the worst cell
                debugCellClosure(i_max, j_max);
            } else {
                std::cout << "  (no interior cells tested)\n";
            }
        }


        

    
    
    void debugCellClosure(int i0, int j0) const {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        std::cout << "\nDEBUG: Cell face-closure details at ("
                  << i0 << "," << j0 << ")\n";

        double Sx = 0.0, Sy = 0.0;

        auto printFace = [&](const char* label,
                             const std::array<double,2>& n,
                             double ds,
                             double sign,
                             int idx) {
            double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]) + 1e-14;
            double nx  = n[0] / mag;
            double ny  = n[1] / mag;
            double fx  = sign * nx * ds;
            double fy  = sign * ny * ds;
            Sx += fx;
            Sy += fy;

            std::cout << "  " << label
                      << "  idx=" << idx
                      << "  n=(" << std::scientific << nx << "," << ny << ")"
                      << "  ds=" << ds
                      << "  contrib=(" << fx << "," << fy << ")\n";
        };

        // Bottom: i-face at (i0, j0), outward = +n
        {
            int idxB = mesh.cellIndex(i0, j0);
            printFace("Bottom i-face", mesh.faceNormal_i[idxB],
                      mesh.faceLen_i[idxB], +1.0, idxB);
        }

        // Top: i-face at (i0, j0+1), outward for cell (i0,j0) = -n
        {
            int idxT = mesh.cellIndex(i0, j0+1);
            printFace("Top    i-face", mesh.faceNormal_i[idxT],
                      mesh.faceLen_i[idxT], -1.0, idxT);
        }

        // Left: j-face at (i0, j0), outward = -n
        {
            int idxL = mesh.cellIndex(i0, j0);
            printFace("Left   j-face", mesh.faceNormal_j[idxL],
                      mesh.faceLen_j[idxL], -1.0, idxL);
        }

        // Right: j-face at (i0+1, j0), outward = +n
        {
            int idxR = mesh.cellIndex(i0+1, j0);
            printFace("Right  j-face", mesh.faceNormal_j[idxR],
                      mesh.faceLen_j[idxR], +1.0, idxR);
        }

        double magS = std::sqrt(Sx*Sx + Sy*Sy);
        std::cout << "  ==> Sum S = (" << Sx << "," << Sy
                  << "), |S| = " << magS << "\n";
    }
        
        
        // ========================================================================
    // METRIC FLUX-ASSEMBLY TEST:
    // Compare direct ∑ n ds (Test 6) with residual-assembly version
    // ========================================================================
    void runMetricFluxAssemblyTest(const FluxCalculator& fluxCalc) const {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        std::vector<Conservative> R_metric(mesh.niTotal * mesh.njTotal);
        fluxCalc.computeMetricResidualFluxAssembly(R_metric, mesh);

        double max_mag      = 0.0;
        double sum_rhou_all = 0.0;
        double sum_rhov_all = 0.0;
        int    i_max = -1, j_max_idx = -1;

        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);

                double sx = R_metric[idx].rhou;
                double sy = R_metric[idx].rhov;
                double mag = std::sqrt(sx*sx + sy*sy);

                sum_rhou_all += sx;
                sum_rhov_all += sy;

                if (mag > max_mag) {
                    max_mag   = mag;
                    i_max     = i;
                    j_max_idx = j;
                }
            }
        }

        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "METRIC FLUX-ASSEMBLY TEST (F = n * ds)\n";
        std::cout << std::string(70, '=') << "\n";

        std::cout << "  Global sum of metric residuals:\n"
                  << "    sum R_rhou = " << sum_rhou_all << "\n"
                  << "    sum R_rhov = " << sum_rhov_all << "\n";

        std::cout << "  Expectation for a consistent closed surface:\n"
                  << "    max |R|      ~ O(1e-12)\n"
                  << "    sum R_rhou   ~ 0\n"
                  << "    sum R_rhov   ~ 0\n";

        if (max_mag > 1e-8) {
            std::cout << "  ✗ Metric-assembly test FAILED: residual assembly and\n"
                      << "    geometry (face normals / ds) are inconsistent.\n";
            std::cout << "    Use debugCellClosure() at the reported (i,j) and\n"
                      << "    also inspect the flux loops around that cell.\n";
        } else {
            std::cout << "  ✓ Metric-assembly test PASSED: residual assembly uses\n"
                      << "    the same geometry as the direct ∑ n ds test.\n";
        }

        std::cout << std::string(70, '=') << "\n\n";
    }
    

    // ========================================================================
    // DETAILED METRIC FLUX-ASSEMBLY DEBUG FOR ONE CELL
    //
    // This walks the *same loops* as computeMetricResidualFluxAssembly,
    // but only accumulates / prints contributions to the chosen cell (i0,j0).
    //
    // Then it prints the sum and calls debugCellClosure(i0,j0) so you can
    // compare:
    //
    //   - direct ∑ n ds   (debugCellClosure)      <-- pure geometry
    //   - flux-based sum  via residual assembly   <-- what computeResidual does
    // ========================================================================
    void debugMetricFluxCell(int i0, int j0) const {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        int target = mesh.cellIndex(i0, j0);

        double Rx = 0.0, Ry = 0.0;

        auto makeF = [&](const std::array<double,2>& n, double ds) {
            Conservative F;
            double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]) + 1e-14;
            double nx  = n[0] / mag;
            double ny  = n[1] / mag;
            F.rho  = 0.0;
            F.rhou = nx * ds;
            F.rhov = ny * ds;
            F.rhoE = 0.0;
            return F;
        };

        std::cout << "\n------------------------------------------------------------\n";
        std::cout << "DEBUG METRIC FLUX-ASSEMBLY at cell (i,j) = (" 
                << i0 << "," << j0 << ")\n";
        std::cout << "Cell index = " << target << "\n";

        // ================= I-FACES =================
        for (int j = jmin; j <= jmax; ++j) {
            for (int i = imin; i < imax; ++i) {

                if (j == jmin) {
                    // bottom boundary: bottom i-face of (i,jmin)
                    int c = mesh.cellIndex(i, jmin);
                    if (c != target) continue;

                    const auto& n_out = mesh.faceNormal_i_out[c];
                    double      ds = mesh.faceLen_i[c];
                    Conservative F = makeF(n_out, ds);

                    Rx += F.rhou;
                    Ry += F.rhov;

                    std::cout << "  I-face: bottom wall face at (i="<<i<<", j="<<jmin<<")\n"
                            << "    contrib = (+F) = (" << F.rhou << ", " << F.rhov << ")\n";
                }
                else if (j < jmax) {
                    // interior i-face between (i,j-1) [L] and (i,j) [R]
                    int idL = mesh.cellIndex(i, j-1);
                    int idR = mesh.cellIndex(i, j);
                    int faceIdx = mesh.cellIndex(i, j);
                    const auto& n  = mesh.faceNormal_i[faceIdx];
                    double      ds = mesh.faceLen_i[faceIdx];
                    Conservative F = makeF(n, ds);

                    if (idL == target) {
                        // target is LEFT cell -> this is its RIGHT face, outward = -F
                        Rx -= F.rhou;
                        Ry -= F.rhov;
                        std::cout << "  J-face: RIGHT of target ... contrib = (-F) = ("
                                << -F.rhou << ", " << -F.rhov << ")\n";
                    }
                    if (idR == target) {
                        // target is RIGHT cell -> this is its LEFT face, outward = +(-F) i.e. also -F
                        Rx -= F.rhou;
                        Ry -= F.rhov;
                        std::cout << "  J-face: LEFT of target ... contrib = (-F) = ("
                                << -F.rhou << ", " << -F.rhov << ")\n";
                    }
                }
                else { // j == jmax (farfield band)
                    int idInt   = mesh.cellIndex(i, jmax-1);
                    if (idInt != target) continue;

                    int faceIdx = mesh.cellIndex(i, jmax);
                    const auto& n  = mesh.faceNormal_i_out[faceIdx];
                    double      ds = mesh.faceLen_i[faceIdx];
                    Conservative F = makeF(n, ds);

                    Rx -= F.rhou;
                    Ry -= F.rhov;

                    std::cout << "  I-face: TOP farfield face for target (i="<<i<<", j="<<jmax-1<<")\n"
                            << "    contrib = (-F) = (" << -F.rhou << ", " << -F.rhov << ")\n";
                }
            }
        }

        // ================= J-FACES =================
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin + 1; i < imax; ++i) {
                int idL = mesh.cellIndex(i - 1, j);
                int idR = mesh.cellIndex(i,     j);
                int faceIdx   = mesh.cellIndex(i, j);
                const auto& n = mesh.faceNormal_j[faceIdx];
                double      ds = mesh.faceLen_j[faceIdx];
                Conservative F = makeF(n, ds);

                if (idL == target) {
                    // target is LEFT cell -> this is its RIGHT face, contrib = -F
                    Rx -= F.rhou;
                    Ry -= F.rhov;
                    std::cout << "  J-face: RIGHT of target (between ("<<i-1<<","<<j<<") and ("
                            << i << "," << j << "))\n"
                            << "    target is LEFT cell => contrib = (-F) = ("
                            << -F.rhou << ", " << -F.rhov << ")\n";
                }
                if (idR == target) {
                    // target is RIGHT cell -> this is its LEFT face, contrib = +F
                    Rx += F.rhou;
                    Ry += F.rhov;
                    std::cout << "  J-face: LEFT of target (between ("<<i-1<<","<<j<<") and ("
                            << i << "," << j << "))\n"
                            << "    target is RIGHT cell => contrib = (+F) = ("
                            << F.rhou << ", " << F.rhov << ")\n";
                }
            }

            // periodic closure: (imax-1,j) [L] -> (imin,j) [R]
            {
                int idL = mesh.cellIndex(imax - 1, j);
                int idR = mesh.cellIndex(imin,     j);
                int faceIdx   = mesh.cellIndex(imin, j);
                const auto& n = mesh.faceNormal_j[faceIdx];
                double      ds = mesh.faceLen_j[faceIdx];
                Conservative F = makeF(n, ds);
                if (idL == target) {
                    // target is LEFT cell -> this is its RIGHT face, outward = -F
                    Rx -= F.rhou;
                    Ry -= F.rhov;
                    std::cout << "  J-face: RIGHT of target ... contrib = (-F) = ("
                            << -F.rhou << ", " << -F.rhov << ")\n";
                }
                if (idR == target) {
                    // target is RIGHT cell -> this is its LEFT face, outward = +(-F) i.e. also -F
                    Rx -= F.rhou;
                    Ry -= F.rhov;
                    std::cout << "  J-face: LEFT of target ... contrib = (-F) = ("
                            << -F.rhou << ", " << -F.rhov << ")\n";
                }
            }
        }

        double magR = std::sqrt(Rx*Rx + Ry*Ry);
        std::cout << "\n  ==> Sum of metric flux contributions at cell ("
                << i0 << "," << j0 << ") via residual assembly:\n"
                << "      R_metric = (" << Rx << ", " << Ry << "), |R| = " << magR << "\n";

        std::cout << "\n  Now direct geometry closure for the same cell:\n";
        debugCellClosure(i0, j0);

        std::cout << "------------------------------------------------------------\n\n";
    }
    // Find and dump the first cell with NaN / negative p or rho
    void debugPhysicalState(
        const std::vector<Conservative>& U,
        const Mesh& mesh,
        const Config& cfg
    ) const {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);
                const auto& Ui = U[idx];

                double rho = Ui.rho;
                double p   = EOS::pressure(Ui, cfg.gamma);

                bool bad = std::isnan(rho) || std::isnan(p) ||
                        (rho <= 0.0) || (p <= 0.0);

                if (bad) {
                    Primitive W(Ui, cfg.gamma);

                    std::cout << "\n*** DEBUG: Bad physical state detected ***\n";
                    std::cout << "Cell (i,j) = (" << i << "," << j << "), idx = " << idx << "\n";
                    std::cout << "  rho  = " << rho << "\n";
                    std::cout << "  p    = " << p   << "\n";
                    std::cout << "  U    = (rho=" << Ui.rho
                            << ", rhou=" << Ui.rhou
                            << ", rhov=" << Ui.rhov
                            << ", rhoE=" << Ui.rhoE << ")\n";
                    std::cout << "  W    = (rho=" << W.rho
                            << ", u="   << W.u
                            << ", v="   << W.v
                            << ", p="   << W.p << ")\n";
                    std::cout << "  x,y  = (" << mesh.x(i,j) << ", " << mesh.y(i,j) << ")\n";

                    return; // stop at first bad cell
                }
            }
        }

        std::cout << "\nDEBUG: No bad cells (NaN/negative p or rho) found.\n";
    }
    // ========================================================================
    // BC vs GEOMETRY CROSS-CHECK
    // Call this right AFTER you apply boundary conditions
    // ========================================================================
    void validateBCvsGeometry(const std::vector<Conservative>& U) {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        // re-use the inferred center to define "outward"
        double xc = 0.0, yc = 0.0;
        {
            int cnt = 0;
            for (int i = imin; i < imax; ++i) {
                xc += mesh.x(i, jmin);
                yc += mesh.y(i, jmin);
                ++cnt;
            }
            xc /= std::max(1, cnt);
            yc /= std::max(1, cnt);
        }

        std::cout << "\n--- BC vs GEOMETRY CHECK ---------------------------------\n";

        // ------------------------------------------------------------
        // 1) WALL GHOSTS: slip = reflected normal velocity
        //    interior at (i, jmin), ghosts at (i, jmin-1), (i, jmin-2)
        // ------------------------------------------------------------
        std::cout << "Wall ghosts (slip reflection):\n";
        for (int i = imin; i < std::min(imax, imin+6); ++i) {
            int idx_wall = mesh.cellIndex(i, jmin);
            Primitive Ww(U[idx_wall], cfg.gamma);

            // wall normal from interior
            std::array<double,2> n = mesh.faceNormal_i_out[idx_wall];
            double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]);
            double nx = n[0] / (mag + 1e-14);
            double ny = n[1] / (mag + 1e-14);

            // interior vn
            double vn_int = Ww.u * nx + Ww.v * ny;

            // first ghost
            int idx_g1 = mesh.cellIndex(i, jmin-1);
            Primitive Wg1(U[idx_g1], cfg.gamma);
            double vn_g1 = Wg1.u * nx + Wg1.v * ny;

            // second ghost
            int idx_g2 = mesh.cellIndex(i, jmin-2);
            Primitive Wg2(U[idx_g2], cfg.gamma);
            double vn_g2 = Wg2.u * nx + Wg2.v * ny;

            // outwardness of ghost geometry
            double rx = mesh.x(i, jmin-1) - xc;
            double ry = mesh.y(i, jmin-1) - yc;
            const auto& ng = mesh.faceNormal_i[idx_g1]; // geometry for ghost row
            double dot_out = ng[0]*rx + ng[1]*ry;

            std::cout << "  i=" << i
                      << "  vn_int=" << vn_int
                      << "  vn_g1=" << vn_g1
                      << "  vn_g2=" << vn_g2
                      << "  p_int=" << Ww.p
                      << "  p_g1=" << Wg1.p
                      << "  geom n·r(ghost)=" << dot_out
                      << "\n";
        }

        // ------------------------------------------------------------
        // 2) FARFIELD GHOSTS: should look like freestream-ish
        //    interior at (i, jmax-1), ghosts at (i, jmax), (i, jmax+1)
        // ------------------------------------------------------------
        std::cout << "\nFarfield ghosts:\n";
        Primitive Winf = cfg.getFreestream();
        for (int i = imin; i < std::min(imax, imin+6); ++i) {
            int idx_int = mesh.cellIndex(i, jmax-1);
            Primitive Wi(U[idx_int], cfg.gamma);

            int idx_g1 = mesh.cellIndex(i, jmax);
            int idx_g2 = mesh.cellIndex(i, jmax+1);
            Primitive Wg1(U[idx_g1], cfg.gamma);
            Primitive Wg2(U[idx_g2], cfg.gamma);

            // geometry outward check for ghost row
            double rx = mesh.x(i, jmax) - xc;
            double ry = mesh.y(i, jmax) - yc;
            const auto& ng = mesh.faceNormal_i[ mesh.cellIndex(i, jmax) ]; 
            double dot_out = ng[0]*rx + ng[1]*ry;

            std::cout << "  i=" << i
                      << "  INT(p=" << Wi.p << ",u=" << Wi.u << ",v=" << Wi.v << ")"
                      << "  G1(p=" << Wg1.p << ",u=" << Wg1.u << ",v=" << Wg1.v << ")"
                      << "  G2(p=" << Wg2.p << ",u=" << Wg2.u << ",v=" << Wg2.v << ")"
                      << "  INF(p=" << Winf.p << ",u=" << Winf.u << ",v=" << Winf.v << ")"
                      << "  geom n·r(top)=" << dot_out
                      << "\n";
        }

        // ------------------------------------------------------------
        // 3) PERIODIC GHOSTS: left/right wrap
        //    left ghosts: i = imin-1, imin-2  <-  from right interior
        //    right ghosts: i = imax, imax+1   <-  from left interior
        // We check a few j’s, including jmin and jmax-1
        // ------------------------------------------------------------
        std::cout << "\nPeriodic ghosts:\n";
        for (int j = jmin-1; j <= jmax; j += (jmax - jmin + 1)) { // bottom ghosts and top interior
            if (j < 0 || j >= mesh.njTotal) continue;
            // left side
            for (int g = 1; g <= NGHOST; ++g) {
                int i_ghost = imin - g;
                int i_donor = imax - g;
                int idx_g = mesh.cellIndex(i_ghost, j);
                int idx_d = mesh.cellIndex(i_donor, j);
                Primitive Pg(U[idx_g], cfg.gamma);
                Primitive Pd(U[idx_d], cfg.gamma);
                std::cout << "  j=" << j << " Lg" << g
                          << "  ghost(rho=" << Pg.rho << ",u=" << Pg.u << ",v=" << Pg.v << ",p=" << Pg.p << ")"
                          << "  donor(rho=" << Pd.rho << ",u=" << Pd.u << ",v=" << Pd.v << ",p=" << Pd.p << ")\n";
            }
            // right side
            for (int g = 0; g < NGHOST; ++g) {
                int i_ghost = imax + g;
                int i_donor = imin + g;
                int idx_g = mesh.cellIndex(i_ghost, j);
                int idx_d = mesh.cellIndex(i_donor, j);
                Primitive Pg(U[idx_g], cfg.gamma);
                Primitive Pd(U[idx_d], cfg.gamma);
                std::cout << "  j=" << j << " Rg" << g
                          << "  ghost(rho=" << Pg.rho << ",u=" << Pg.u << ",v=" << Pg.v << ",p=" << Pg.p << ")"
                          << "  donor(rho=" << Pd.rho << ",u=" << Pd.u << ",v=" << Pd.v << ",p=" << Pd.p << ")\n";
            }
        }

        std::cout << "-----------------------------------------------------------\n\n";
    }
    
    // ========================================================================
    // MAIN DIAGNOSTICS (run periodically during solve)
    // ========================================================================
    void runDiagnostics(
        int iter,
        const std::vector<Conservative>& U,
        const std::vector<Conservative>& R
    ) {
        if (iter % 100 != 0) return;
        
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "DIAGNOSTICS: Iteration " << iter << "\n";
        std::cout << std::string(70, '=') << "\n\n";
        
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        // recompute center for wall checks
        double xc = 0.0, yc = 0.0;
        {
            int cnt = 0;
            for (int i = imin; i < imax; ++i) {
                xc += mesh.x(i, jmin);
                yc += mesh.y(i, jmin);
                ++cnt;
            }
            xc /= std::max(1, cnt);
            yc /= std::max(1, cnt);
        }
        
        // ====================================================================
        // CONSERVATION LAWS
        // ====================================================================
        std::cout << "CONSERVATION LAWS\n";
        std::cout << std::string(50, '-') << "\n";
        
        double total_mass = 0.0;
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);
                total_mass += U[idx].rho * mesh.cellArea[idx];
            }
        }
        std::cout << "Total mass:     " << std::scientific << std::setprecision(6) 
                  << total_mass << "\n";
        
        double sum_rho = 0.0, sum_rhou = 0.0, sum_rhov = 0.0, sum_rhoE = 0.0;
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);
                sum_rho  += R[idx].rho;
                sum_rhou += R[idx].rhou;
                sum_rhov += R[idx].rhov;
                sum_rhoE += R[idx].rhoE;
            }
        }
        
        std::cout << "Sum(R_rho):     " << sum_rho << "\n";
        std::cout << "Sum(R_rhou):    " << sum_rhou << "\n";
        std::cout << "Sum(R_rhov):    " << sum_rhov << " ";
        if (std::abs(cfg.alpha_deg) < 0.1 && std::abs(sum_rhov) > 1e-10) {
            std::cout << "✗ Should be ~0 for alpha=0°\n";
        } else {
            std::cout << "\n";
        }
        std::cout << "Sum(R_rhoE):    " << sum_rhoE << "\n";
        
        // ====================================================================
        // WALL BOUNDARY CHECKS
        // ====================================================================
        std::cout << "\nWALL BOUNDARY\n";
        std::cout << std::string(50, '-') << "\n";
        
        double max_vn = 0.0;
        double total_wall_flux_y = 0.0;
        
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, jmin);
            Primitive W(U[idx], cfg.gamma);
            
            std::array<double,2> normal =mesh.faceNormal_i_out[idx];;
            double ds  = mesh.faceLen_i[idx];
            double mag = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
            double nx = normal[0] / (mag + 1e-14);
            double ny = normal[1] / (mag + 1e-14);
            
            // radial outward (for info)
            double rx = mesh.x(i, jmin) - xc;
            double ry = mesh.y(i, jmin) - yc;
            double dot_out = nx*rx + ny*ry;
            
            // Normal velocity
            double vn = W.u * nx + W.v * ny;
            max_vn = std::max(max_vn, std::abs(vn));
            
            double p_wall = W.p;
            double Fy_wall = p_wall * ny * ds;
            total_wall_flux_y += Fy_wall;

            if (dot_out < 0.0) {
                std::cout << "  ✗ Wall normal at i=" << i << " is not outward (n·r=" 
                          << dot_out << ")\n";
            }
        }
        
        std::cout << "Max |v·n| at wall: " << std::scientific << max_vn;
        if (max_vn < 1e-10) {
            std::cout << " ✓ (slip condition satisfied)\n";
        } else {
            std::cout << " ✗ (slip condition violated)\n";
        }
        
        // we still print y-force, but we don't enforce sign on O-grid
        std::cout << "Total wall y-force: " << std::fixed << std::setprecision(4) 
                  << total_wall_flux_y << " (info)\n";
        
        // ====================================================================
        // PHYSICAL STATE CHECKS
        // ====================================================================
        std::cout << "\nPHYSICAL STATE\n";
        std::cout << std::string(50, '-') << "\n";
        
        double min_p = 1e10, max_p = -1e10;
        double min_rho = 1e10, max_rho = -1e10;
        int nan_count = 0, neg_p = 0, neg_rho = 0;
        
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);
                double p = EOS::pressure(U[idx], cfg.gamma);
                double rho = U[idx].rho;
                
                if (std::isnan(p) || std::isnan(rho)) nan_count++;
                if (p < 0) neg_p++;
                if (rho < 0) neg_rho++;
                
                min_p = std::min(min_p, p);
                max_p = std::max(max_p, p);
                min_rho = std::min(min_rho, rho);
                max_rho = std::max(max_rho, rho);
            }
        }
        
        std::cout << "Pressure range:  [" << std::fixed << std::setprecision(4) 
                  << min_p << ", " << max_p << "]";
        if (neg_p > 0) {
            std::cout << " ✗ " << neg_p << " negative!\n";
        } else {
            std::cout << " ✓\n";
        }
        
        std::cout << "Density range:   [" << min_rho << ", " << max_rho << "]";
        if (neg_rho > 0) {
            std::cout << " ✗ " << neg_rho << " negative!\n";
        } else {
            std::cout << " ✓\n";
        }
        
        if (nan_count > 0) {
            std::cout << "✗ " << nan_count << " NaN values detected!\n";
        }
        // If anything is physically broken, dump the first bad cell
        if (std::isnan(total_mass) || nan_count > 0 || neg_p > 0 || neg_rho > 0) {
            std::cout << "\n!!! PHYSICAL FAILURE DETECTED, dumping first bad cell...\n";
            debugPhysicalState(U, mesh, cfg);
        }


        // ====================================================================
        // RESIDUAL NORMS
        // ====================================================================
        std::cout << "\nRESIDUAL NORMS\n";
        std::cout << std::string(50, '-') << "\n";
        
        double L2_rho = 0.0, L2_rhou = 0.0, L2_rhov = 0.0, L2_rhoE = 0.0;
        int count = 0;
        
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);
                L2_rho  += R[idx].rho  * R[idx].rho;
                L2_rhou += R[idx].rhou * R[idx].rhou;
                L2_rhov += R[idx].rhov * R[idx].rhov;
                L2_rhoE += R[idx].rhoE * R[idx].rhoE;
                count++;
            }
        }
        
        L2_rho  = std::sqrt(L2_rho  / (count + 1e-14));
        L2_rhou = std::sqrt(L2_rhou / (count + 1e-14));
        L2_rhov = std::sqrt(L2_rhov / (count + 1e-14));
        L2_rhoE = std::sqrt(L2_rhoE / (count + 1e-14));
        
        std::cout << "L2(R_rho):  " << std::scientific << std::setprecision(3) << L2_rho << "\n";
        std::cout << "L2(R_rhou): " << L2_rhou << "\n";
        std::cout << "L2(R_rhov): " << L2_rhov << "\n";
        std::cout << "L2(R_rhoE): " << L2_rhoE << "\n";
        
        double max_R_mag = 0.0;
        int i_max = imin, j_max_idx = jmin;
        
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);
                double mag = std::sqrt(
                    R[idx].rho*R[idx].rho + R[idx].rhou*R[idx].rhou +
                    R[idx].rhov*R[idx].rhov + R[idx].rhoE*R[idx].rhoE
                );
                if (mag > max_R_mag) {
                    max_R_mag = mag;
                    i_max = i;
                    j_max_idx = j;
                }
            }
        }
        
        std::cout << "Max |R| at cell (" << i_max << "," << j_max_idx 
                  << "): " << max_R_mag << "\n";
        
        // ====================================================================
        // PHYSICAL CHECKS (for airfoil)
        // ====================================================================
        if (cfg.alpha_deg < 5.0) {
            std::cout << "\nPHYSICAL CHECKS\n";
            std::cout << std::string(50, '-') << "\n";
            
            int i_stag = imin;
            double min_vmag = 1e10;
            
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, jmin);
                Primitive W(U[idx], cfg.gamma);
                double vmag = std::sqrt(W.u*W.u + W.v*W.v);
                if (vmag < min_vmag) {
                    min_vmag = vmag;
                    i_stag = i;
                }
            }
            
            int idx_stag = mesh.cellIndex(i_stag, jmin);
            double p_stag = EOS::pressure(U[idx_stag], cfg.gamma);
            double q_inf = 0.5 * cfg.rho_inf * cfg.Mach_inf * cfg.Mach_inf 
                         * cfg.a_inf * cfg.a_inf;
            double Cp_stag = (p_stag - cfg.p_inf) / (q_inf + 1e-14);
            
            std::cout << "Stagnation point at i=" << i_stag 
                      << ", x=" << std::fixed << std::setprecision(4) 
                      << mesh.x(i_stag, jmin) << "\n";
            std::cout << "Cp at stagnation: " << std::setprecision(3) << Cp_stag;
            
            double expected_Cp = 1.0 - 0.25 * cfg.Mach_inf * cfg.Mach_inf;
            if (std::abs(Cp_stag - expected_Cp) < 0.2) {
                std::cout << " ✓ (expected ~" << expected_Cp << ")\n";
            } else {
                std::cout << " (expected ~" << expected_Cp << ")\n";
            }
        }
        
        std::cout << "\n" << std::string(70, '=') << "\n\n";
        
        // ====================================================================
        // LOG TO FILE
        // ====================================================================
        log_file << iter << ","
                 << std::scientific << std::setprecision(6)
                 << total_mass << ","
                 << sum_rho << ","
                 << sum_rhou << ","
                 << sum_rhov << ","
                 << sum_rhoE << ","
                 << max_vn << ","
                 << min_p << ","
                 << max_p << ","
                 << min_rho << ","
                 << max_rho << ","
                 << L2_rho << ","
                 << total_wall_flux_y << "\n";
        log_file.flush();
    }

    // ========================================================================
    // FREESTREAM AND GHOST CELL TESTING 
    // ========================================================================
    void centralFluxResidualAudit(
        const Mesh& mesh,
        const Config& cfg,
        const std::vector<Conservative>& U
    ) {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        // residual-style accumulator
        std::vector<Conservative> Rc(mesh.niTotal * mesh.njTotal);
        for (auto& r : Rc) r = {0,0,0,0};

        auto centralFaceFlux = [&](const Conservative& UL,
                                const Conservative& UR,
                                const std::array<double,2>& n,
                                double ds) {
            // normalize
            double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]);
            double nx = n[0] / (mag + 1e-14);
            double ny = n[1] / (mag + 1e-14);

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
            return F;
        };

        // ------------------------------
        // I-FACES
        // ------------------------------
        for (int j = jmin; j <= jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                if (j == jmin) {
                    // wall face: use interior cell twice (like wall BC),
                    // or, better, use the existing interior cell vs itself
                    int idx = mesh.cellIndex(i, jmin);
                    Primitive W(U[idx], cfg.gamma);

                    // use the SAME geometry as fluxCalculator::computeResidual for j==jmin
                    int faceIdx = mesh.cellIndex(i, jmin);
                    const auto& n = mesh.faceNormal_i_out[faceIdx];
                    double ds = mesh.faceLen_i[faceIdx];


                    double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]);
                    double nx = n[0] / (mag + 1e-14);
                    double ny = n[1] / (mag + 1e-14);

                    Conservative Fw;
                    Fw.rho  = 0.0;
                    Fw.rhou = W.p * nx * ds;
                    Fw.rhov = W.p * ny * ds;
                    Fw.rhoE = 0.0;

                    Rc[idx].rho  += Fw.rho;
                    Rc[idx].rhou += Fw.rhou;
                    Rc[idx].rhov += Fw.rhov;
                    Rc[idx].rhoE += Fw.rhoE;
                }
                else if (j < jmax) {
                    int idL = mesh.cellIndex(i, j-1);
                    int idR = mesh.cellIndex(i, j);
                    const auto& n = mesh.faceNormal_i[idL];
                    double ds = mesh.faceLen_i[idL];

                    Conservative F = centralFaceFlux(U[idL], U[idR], n, ds);

                    Rc[idL].rho  -= F.rho;
                    Rc[idL].rhou -= F.rhou;
                    Rc[idL].rhov -= F.rhov;
                    Rc[idL].rhoE -= F.rhoE;

                    Rc[idR].rho  += F.rho;
                    Rc[idR].rhou += F.rhou;
                    Rc[idR].rhov += F.rhov;
                    Rc[idR].rhoE += F.rhoE;
                } else {
                    // top: interior vs top ghost
                    int idL = mesh.cellIndex(i, jmax-1);
                    int idR = mesh.cellIndex(i, jmax);
                    const auto& n0 = mesh.faceNormal_i[idL];

                    double ds = mesh.faceLen_i[idL];
                    std::array<double,2> n = { -n0[0], -n0[1] };

                    Conservative F = centralFaceFlux(U[idL], U[idR], n, ds);

                    Rc[idL].rho  -= F.rho;
                    Rc[idL].rhou -= F.rhou;
                    Rc[idL].rhov -= F.rhov;
                    Rc[idL].rhoE -= F.rhoE;
                }
            }
        }

        // ------------------------------
        // J-FACES (single periodic closure per j)
        // ------------------------------
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin + 1; i < imax; ++i) {
                int idL = mesh.cellIndex(i-1, j);
                int idR = mesh.cellIndex(i,   j);
                const auto& n = mesh.faceNormal_j[idL];
                double ds = mesh.faceLen_j[idL];

                Conservative F = centralFaceFlux(U[idL], U[idR], n, ds);

                Rc[idL].rho  -= F.rho;
                Rc[idL].rhou -= F.rhou;
                Rc[idL].rhov -= F.rhov;
                Rc[idL].rhoE -= F.rhoE;

                Rc[idR].rho  += F.rho;
                Rc[idR].rhou += F.rhou;
                Rc[idR].rhov += F.rhov;
                Rc[idR].rhoE += F.rhoE;
            }

            // periodic closure (imax-1 -> imin)
            {
                int idL = mesh.cellIndex(imax - 1, j);
                int idR = mesh.cellIndex(imin,     j);
                const auto& n = mesh.faceNormal_j[idL];
                double ds = mesh.faceLen_j[idL];

                Conservative F = centralFaceFlux(U[idL], U[idR], n, ds);

                Rc[idL].rho  -= F.rho;
                Rc[idL].rhou -= F.rhou;
                Rc[idL].rhov -= F.rhov;
                Rc[idL].rhoE -= F.rhoE;

                Rc[idR].rho  += F.rho;
                Rc[idR].rhou += F.rhou;
                Rc[idR].rhov += F.rhov;
                Rc[idR].rhoE += F.rhoE;
            }
        }

        // ------------ report ------------
        double maxWall = 0.0, maxInterior = 0.0, maxFar = 0.0;
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);
                double m = std::sqrt(
                    Rc[idx].rho*Rc[idx].rho +
                    Rc[idx].rhou*Rc[idx].rhou +
                    Rc[idx].rhov*Rc[idx].rhov +
                    Rc[idx].rhoE*Rc[idx].rhoE
                );
                if (j == jmin)        maxWall     = std::max(maxWall,     m);
                else if (j == jmax-1) maxFar      = std::max(maxFar,      m);
                else                  maxInterior = std::max(maxInterior, m);
            }
        }

        std::cout << "\nCENTRAL-FLUX GEOMETRY AUDIT\n";
        std::cout << "  max |R| wall:     " << maxWall << "\n";
        std::cout << "  max |R| interior: " << maxInterior << "\n";
        std::cout << "  max |R| farfield: " << maxFar << "\n";
        std::cout << "If these are big for uniform U, the issue is still in face setup/BC, not JST.\n\n";
    }



    void dumpFarfieldFluxes(const Mesh& mesh,
                        const Config& cfg,
                        const std::vector<Conservative>& U)
    {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        std::cout << "\nFARFIELD FLUX DUMP:\n";
        for (int i = imin; i < imax; ++i) {
            int idL = mesh.cellIndex(i, jmax-1);
            int idR = mesh.cellIndex(i, jmax);   // ghost
            Primitive WL(U[idL], cfg.gamma);
            Primitive WR(U[idR], cfg.gamma);

            int idx = mesh.cellIndex(i, jmax);
            const auto& n0 = mesh.faceNormal_i[idx];
            double mag = std::sqrt(n0[0]*n0[0] + n0[1]*n0[1]);
            double ds = mesh.faceLen_i[idx];
            double nx  = -n0[0] / (mag + 1e-14);
            double ny  = -n0[1] / (mag + 1e-14);

            double unL = WL.u * nx + WL.v * ny;
            double unR = WR.u * nx + WR.v * ny;

            // central mass flux
            double Fm = 0.5 * (WL.rho * unL + WR.rho * unR) * ds;

            std::cout << "  i=" << i
                    << " unL=" << unL
                    << " unR=" << unR
                    << " Fm="  << Fm
                    << " ds="  << ds
                    << "\n";
        }
    }
    void dumpWallFluxes(const Mesh& mesh,
                        const Config& cfg,
                        const std::vector<Conservative>& U)
    {
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        std::cout << "\nWALL FLUX DUMP:\n";
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, jmin);
            Primitive W(U[idx], cfg.gamma);
            const auto& n = mesh.faceNormal_i_out[idx];
            double ds = mesh.faceLen_i[idx];
            double mag = std::sqrt(n[0]*n[0] + n[1]*n[1]);
            double nx = n[0] / (mag + 1e-14);
            double ny = n[1] / (mag + 1e-14);

            double Fx = W.p * nx * ds;
            double Fy = W.p * ny * ds;

            std::cout << "  i=" << i
                    << " p=" << W.p
                    << " n=(" << nx << "," << ny << ")"
                    << " ds=" << ds
                    << " Fx=" << Fx
                    << " Fy=" << Fy
                    << "\n";
        }
    }
    
    // ========================================================================
    // FREESTREAM PRESERVATION TEST + FACE-BY-FACE DECOMPOSITION
    // ========================================================================
    void testFreestreamPreservation(
        const std::vector<Conservative>& U,
        const std::vector<Conservative>& R,
        const FluxCalculator& fluxCalc
    ) {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "FREESTREAM PRESERVATION TEST\n";
        std::cout << std::string(70, '=') << "\n";
        std::cout << "Testing if uniform freestream produces zero residual\n";
        std::cout << "(This validates flux discretization and geometry)\n\n";
        
        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);
        
        double max_R_wall = 0.0, max_R_interior = 0.0, max_R_farfield = 0.0;
        
        // --- wall band (j = jmin) ---
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, jmin);
            double R_mag = std::sqrt(
                R[idx].rho*R[idx].rho + R[idx].rhou*R[idx].rhou +
                R[idx].rhov*R[idx].rhov + R[idx].rhoE*R[idx].rhoE
            );
            max_R_wall = std::max(max_R_wall, R_mag);
        }
        
        // --- interior band ---
        for (int j = jmin+1; j < jmax-1; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);
                double R_mag = std::sqrt(
                    R[idx].rho*R[idx].rho + R[idx].rhou*R[idx].rhou +
                    R[idx].rhov*R[idx].rhov + R[idx].rhoE*R[idx].rhoE
                );
                max_R_interior = std::max(max_R_interior, R_mag);
            }
        }
        
        // --- farfield band (j = jmax-1) ---
        for (int i = imin; i < imax; ++i) {
            int idx = mesh.cellIndex(i, jmax-1);
            double R_mag = std::sqrt(
                R[idx].rho*R[idx].rho + R[idx].rhou*R[idx].rhou +
                R[idx].rhov*R[idx].rhov + R[idx].rhoE*R[idx].rhoE
            );
            max_R_farfield = std::max(max_R_farfield, R_mag);
        }
        
        std::cout << "Max |R| at wall:     " << std::scientific << std::setprecision(3) 
                << max_R_wall;
        if (max_R_wall < 1e-10) {
            std::cout << " ✓\n";
        } else {
            std::cout << " ✗ (wall flux error)\n";
        }
        
        std::cout << "Max |R| in interior: " << max_R_interior;
        if (max_R_interior < 1e-10) {
            std::cout << " ✓\n";
        } else {
            std::cout << " ✗ (flux discretization error)\n";
        }
        
        std::cout << "Max |R| at farfield: " << max_R_farfield;
        if (max_R_farfield < 1e-10) {
            std::cout << " ✓\n";
        } else {
            std::cout << " ✗ (farfield BC error)\n";
        }
        
        double max_R_overall = std::max({max_R_wall, max_R_interior, max_R_farfield});
        
        std::cout << "\nOverall: ";
        if (max_R_overall < 1e-10) {
            std::cout << "✓ Freestream preserved (flux discretization correct)\n";
        } else {
            std::cout << "✗ Freestream NOT preserved (geometric/flux error)\n";
            std::cout << "  This indicates fundamental discretization problems!\n";
        }
        
        // ====================================================================
        // NEW: pick an interior cell with largest residual and decompose
        // its residual into contributions from its four faces
        // ====================================================================
        double maxR_int = -1.0;
        int    i_max_int = -1, j_max_int = -1;

        for (int j = jmin+1; j < jmax-1; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = mesh.cellIndex(i, j);
                double R_mag = std::sqrt(
                    R[idx].rho*R[idx].rho + R[idx].rhou*R[idx].rhou +
                    R[idx].rhov*R[idx].rhov + R[idx].rhoE*R[idx].rhoE
                );
                if (R_mag > maxR_int) {
                    maxR_int = R_mag;
                    i_max_int = i;
                    j_max_int = j;
                }
            }
        }

        if (i_max_int >= 0) {
            std::cout << "\n======================================================================\n";
            std::cout << "FREESTREAM RESIDUAL DECOMPOSITION (INTERIOR CELL)\n";
            std::cout << "======================================================================\n";
            std::cout << "Max interior residual cell: (i,j) = ("
                    << i_max_int << "," << j_max_int << "), |R| = "
                    << std::scientific << maxR_int << "\n";

            int i0 = i_max_int;
            int j0 = j_max_int;

            std::cout << "\nDEBUGGING FACE FLUXES AROUND CELL ("
                    << i0 << "," << j0 << ")\n";

            // -----------------------------------------------------------------
            // Bottom i-face: between (i0, j0-1) [L] and (i0, j0) [R]
            // -----------------------------------------------------------------
            if (j0 > jmin) {
                std::cout << "\n--- Bottom i-face (between (" 
                        << i0 << "," << j0-1 << ") [L] and ("
                        << i0 << "," << j0   << ") [R]) ---\n";
                fluxCalc.debugInteriorFaceFlux(U, mesh, cfg, i0, j0, true);
            } else {
                std::cout << "\n[Bottom face touches wall band; skipping interior-face debug]\n";
            }

            // -----------------------------------------------------------------
            // Top i-face: between (i0, j0) [L] and (i0, j0+1) [R]
            // -----------------------------------------------------------------
            if (j0+1 < jmax) {
                std::cout << "\n--- Top i-face (between (" 
                        << i0 << "," << j0   << ") [L] and ("
                        << i0 << "," << j0+1 << ") [R]) ---\n";
                fluxCalc.debugInteriorFaceFlux(U, mesh, cfg, i0, j0+1, true);
            } else {
                std::cout << "\n[Top face touches farfield band; skipping interior-face debug]\n";
            }

            // -----------------------------------------------------------------
            // Left j-face: between (i0-1, j0) [L] and (i0, j0) [R]
            // -----------------------------------------------------------------
            if (i0 > imin) {
                std::cout << "\n--- Left j-face (between (" 
                        << i0-1 << "," << j0 << ") [L] and ("
                        << i0   << "," << j0 << ") [R]) ---\n";
                fluxCalc.debugInteriorFaceFlux(U, mesh, cfg, i0, j0, false);
            } else {
                std::cout << "\n[Left face hits periodic closure; interior-face debug not used here]\n";
            }

            // -----------------------------------------------------------------
            // Right j-face: between (i0, j0) [L] and (i0+1, j0) [R]
            // -----------------------------------------------------------------
            if (i0+1 < imax) {
                std::cout << "\n--- Right j-face (between (" 
                        << i0   << "," << j0 << ") [L] and ("
                        << i0+1 << "," << j0 << ") [R]) ---\n";
                fluxCalc.debugInteriorFaceFlux(U, mesh, cfg, i0+1, j0, false);
            } else {
                std::cout << "\n[Right face hits periodic closure; interior-face debug not used here]\n";
            }

            std::cout << "======================================================================\n\n";
        }

        std::cout << "\n" << std::string(70, '=') << "\n\n";
    }
    void dumpPhysicalState(int i, int j,
                                          const std::vector<Conservative>& U) const
    {
        int idx = mesh.cellIndex(i, j);
        const auto& Ui = U[idx];
        double p = EOS::pressure(Ui, cfg.gamma);
        double u = Ui.rhou / Ui.rho;
        double v = Ui.rhov / Ui.rho;

        std::cout << "*** DEBUG: Bad physical state detected ***\n";
        std::cout << "Cell (i,j) = (" << i << "," << j << "), idx = " << idx << "\n";
        std::cout << "  rho  = " << Ui.rho << "\n";
        std::cout << "  p    = " << p << "\n";
        std::cout << "  U    = (rho=" << Ui.rho
                << ", rhou=" << Ui.rhou
                << ", rhov=" << Ui.rhov
                << ", rhoE=" << Ui.rhoE << ")\n";
        std::cout << "  W    = (rho=" << Ui.rho
                << ", u="   << u
                << ", v="   << v
                << ", p="   << p << ")\n";
        std::cout << "  x,y  = (" << mesh.x(i,j) << ", " << mesh.y(i,j) << ")\n";
    }
};