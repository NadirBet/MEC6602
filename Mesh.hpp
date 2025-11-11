#pragma once
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "Types.hpp"

constexpr int NGHOST = 2;

// ============================================================
// CORRECTED MESH FOR O-GRID TOPOLOGY
// Key fixes:
// 1. Proper quadrilateral area calculation using shoelace formula
// 2. Consistent normal vector orientations
// 3. Ensures discrete geometric conservation laws
// ============================================================
class Mesh {
public:
    int ni, nj;               // Interior cell counts
    int niTotal, njTotal;     // Total including ghosts
    int niNodes, njNodes;     // Node dimensions
    
    // Node coordinates
    std::vector<double> xNodes, yNodes;
    
    // Cell centers
    std::vector<double> xCells, yCells;
    
    // Cell area
    std::vector<double> cellArea;
    
    // Face normals and lengths (GEOMETRIC)
    std::vector<std::array<double,2>> faceNormal_i;
    std::vector<std::array<double,2>> faceNormal_j;
    std::vector<double> faceLen_i, faceLen_j;

    // "Averaged" quantities (GEOMETRIC)
    std::vector<std::array<double,2>> faceNormal_i_avg;
    std::vector<std::array<double,2>> faceNormal_j_avg;
    std::vector<double> faceLen_i_avg, faceLen_j_avg;

    // NEW: “oriented” normals (for BCs / forces etc.)
    // These are the ones we are allowed to flip for “outward” conventions.
    std::vector<std::array<double,2>> faceNormal_i_out;
    std::vector<std::array<double,2>> faceNormal_j_out;
    
    // Spectral radii
    std::vector<double> lambda_i, lambda_j;
    
    Mesh() = default;
    
    Mesh(int ni_, int nj_) : ni(ni_), nj(nj_) {
        niTotal = ni + 2*NGHOST;
        njTotal = nj + 2*NGHOST;
        niNodes = niTotal + 1;
        njNodes = njTotal + 1;
        
        int ncells = niTotal * njTotal;
        int nnodes = niNodes * njNodes;
        
        xNodes.resize(nnodes);
        yNodes.resize(nnodes);
        xCells.resize(ncells);
        yCells.resize(ncells);
        
        cellArea.resize(ncells);
        faceNormal_i.resize(ncells);
        faceNormal_j.resize(ncells);
        faceLen_i.resize(ncells);
        faceLen_j.resize(ncells);

        faceNormal_i_avg.resize(ncells);
        faceNormal_j_avg.resize(ncells);
        faceLen_i_avg.resize(ncells);
        faceLen_j_avg.resize(ncells);

        // NEW: oriented normals
        faceNormal_i_out.resize(ncells);
        faceNormal_j_out.resize(ncells);

        lambda_i.resize(ncells);
        lambda_j.resize(ncells);
    }
    
    // Index helpers
    int nodeIndex(int i, int j) const { return j * niNodes + i; }
    int cellIndex(int i, int j) const { return j * niTotal + i; }
    
    // Cell center access
    double x(int i, int j) const { return xCells[cellIndex(i,j)]; }
    double y(int i, int j) const { return yCells[cellIndex(i,j)]; }
    



    // ------------------------------------------------------------
    // Simple Cartesian test grid initializer
    // Builds a uniform [x0,x1] x [y0,y1] grid in *cell indices*,
    // including ghosts, so interior cells are a clean Cartesian mesh.
    // ------------------------------------------------------------
    void initCartesianTest(double x0, double x1, double y0, double y1) {
        // interior cell counts are ni x nj
        // so there are (ni+1) x (nj+1) physical nodes
        double dx = (x1 - x0) / static_cast<double>(ni);
        double dy = (y1 - y0) / static_cast<double>(nj);

        // Fill all nodes, including ghosts, by extending the uniform grid
        for (int j = 0; j < njNodes; ++j) {
            // physical index j_phys = j - NGHOST  (0..nj are interior nodes)
            double y = y0 + static_cast<double>(j - NGHOST) * dy;
            for (int i = 0; i < niNodes; ++i) {
                // physical index i_phys = i - NGHOST  (0..ni are interior nodes)
                double x = x0 + static_cast<double>(i - NGHOST) * dx;

                int idx = nodeIndex(i, j);
                xNodes[idx] = x;
                yNodes[idx] = y;
            }
        }

        // Compute cell centers from these nodes
        computeCellCenters();
    }
        // ------------------------------------------------------------
    // Read Plot3D mesh (single block, ASCII), O-grid layout
    // ------------------------------------------------------------
    void readPlot3D(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open mesh file: " + filename);
        }
        
        int nBlocks;
        file >> nBlocks;
        if (nBlocks != 1) {
            throw std::runtime_error("Only single-block meshes supported");
        }
        
        int ni_nodes, nj_nodes;
        file >> ni_nodes >> nj_nodes;
        
        // Rebuild with the right physical size
        *this = Mesh(ni_nodes - 1, nj_nodes - 1);
        
        // Read x-coordinates
        for (int j = 0; j < nj_nodes; ++j) {
            for (int i = 0; i < ni_nodes; ++i) {
                double x;
                file >> x;
                int idx = nodeIndex(i + NGHOST, j + NGHOST);
                xNodes[idx] = x;
            }
        }
        
        // Read y-coordinates
        for (int j = 0; j < nj_nodes; ++j) {
            for (int i = 0; i < ni_nodes; ++i) {
                double y;
                file >> y;
                int idx = nodeIndex(i + NGHOST, j + NGHOST);
                yNodes[idx] = y;
            }
        }
        
        file.close();
        
        // Fill ghost nodes
        fillGhostNodes();
        
        // Compute cell centers
        computeCellCenters();
    }
    
    // Interior bounds
    void getInteriorBounds(int& imin, int& imax, int& jmin, int& jmax) const {
        imin = NGHOST;
        imax = NGHOST + ni;
        jmin = NGHOST;
        jmax = NGHOST + nj;
    }
    
    // ------------------------------------------------------------
    // Fill ghost nodes for boundary conditions
    // ------------------------------------------------------------
    void fillGhostNodes() {
        int imin = NGHOST;
        int imax = NGHOST + ni;
        int jmin = NGHOST;
        int jmax = NGHOST + nj;
        
        // 1) Bottom (wall) -> mirror in y
        for (int i = 0; i < niNodes; ++i) {
            for (int g = 0; g < NGHOST; ++g) {
                int j_ghost  = jmin - 1 - g;
                int j_wall   = jmin;
                int j_mirror = jmin + 1 + g;
                
                int idx_g  = nodeIndex(i, j_ghost);
                int idx_w  = nodeIndex(i, j_wall);
                int idx_m  = nodeIndex(i, j_mirror);
                
                xNodes[idx_g] = xNodes[idx_w];
                yNodes[idx_g] = 2.0 * yNodes[idx_w] - yNodes[idx_m];
            }
        }
        
        // 2) Top (farfield) -> linear extrapolation
        for (int i = 0; i < niNodes; ++i) {
            for (int g = 0; g < NGHOST; ++g) {
                int j_ghost = jmax + 1 + g;
                int j1      = jmax - g;
                int j2      = jmax - 1 - g;
                
                int idx_g = nodeIndex(i, j_ghost);
                int idx1  = nodeIndex(i, j1);
                int idx2  = nodeIndex(i, j2);
                
                xNodes[idx_g] = 2.0 * xNodes[idx1] - xNodes[idx2];
                yNodes[idx_g] = 2.0 * yNodes[idx1] - yNodes[idx2];
            }
        }
        
        // 3) i-periodic for O-grid, skipping the duplicated TE
        for (int j = 0; j < njNodes; ++j) {
            // Left ghosts from right
            for (int g = 0; g < NGHOST; ++g) {
                int i_ghost  = imin - 1 - g;
                int i_source = imax - 1 - g;
                int idx_g = nodeIndex(i_ghost, j);
                int idx_s = nodeIndex(i_source, j);
                xNodes[idx_g] = xNodes[idx_s];
                yNodes[idx_g] = yNodes[idx_s];
            }
            // Right ghosts from left
            for (int g = 0; g < NGHOST; ++g) {
                int i_ghost  = imax + 1 + g;
                int i_source = imin + 1 + g;
                int idx_g = nodeIndex(i_ghost, j);
                int idx_s = nodeIndex(i_source, j);
                xNodes[idx_g] = xNodes[idx_s];
                yNodes[idx_g] = yNodes[idx_s];
            }
        }
    }
    
    // ------------------------------------------------------------
    // Compute cell centers from nodes
    // ------------------------------------------------------------
    void computeCellCenters() {
        for (int j = 0; j < njTotal; ++j) {
            for (int i = 0; i < niTotal; ++i) {
                int n00 = nodeIndex(i,     j);
                int n10 = nodeIndex(i + 1, j);
                int n01 = nodeIndex(i,     j + 1);
                int n11 = nodeIndex(i + 1, j + 1);
                
                int cidx = cellIndex(i, j);
                xCells[cidx] = 0.25 * (xNodes[n00] + xNodes[n10] + xNodes[n01] + xNodes[n11]);
                yCells[cidx] = 0.25 * (yNodes[n00] + yNodes[n10] + yNodes[n01] + yNodes[n11]);
            }
        }
    }
    
    // ------------------------------------------------------------
    // Full geometry build
    // ------------------------------------------------------------
    void buildGeometry(const std::vector<Conservative>& U, const Config& cfg) {
        computeCellAreasAndFaces();
        correctNormalOrientations();
        fillGhostGeometry();
        applyPeriodicBoundaries();
        computeAveragedQuantities();
        computeSpectralRadii(U, cfg);
    }
    
    // ------------------------------------------------------------
    // CORRECTED: Area + face normals with proper formulas
    // ------------------------------------------------------------

    void computeCellAreasAndFaces() {
        for (int j = 0; j < njTotal; ++j) {
            for (int i = 0; i < niTotal; ++i) {
                // Get 4 corner nodes of the cell
                int n00 = nodeIndex(i,     j);      // SW
                int n10 = nodeIndex(i + 1, j);      // SE 
                int n01 = nodeIndex(i,     j + 1);  // NW
                int n11 = nodeIndex(i + 1, j + 1);  // NE
                
                double x00 = xNodes[n00], y00 = yNodes[n00];
                double x10 = xNodes[n10], y10 = yNodes[n10];
                double x01 = xNodes[n01], y01 = yNodes[n01];
                double x11 = xNodes[n11], y11 = yNodes[n11];
                
                // Area via shoelace formula
                double A = 0.5 * std::abs(
                    x00*(y10 - y01) + 
                    x10*(y11 - y00) + 
                    x11*(y01 - y10) + 
                    x01*(y00 - y11)
                );
                
                int idx = cellIndex(i, j);
                cellArea[idx] = A;
                
                // Bottom i-face: edge from n00 -> n10
                double dx_i = x10 - x00;
                double dy_i = y10 - y00;
                double nx_i =  dy_i;
                double ny_i = -dx_i;
                faceNormal_i[idx] = {nx_i, ny_i};
                faceLen_i[idx]    = std::sqrt(dx_i*dx_i + dy_i*dy_i);
                
                // Left j-face: edge from n00 -> n01
                double dx_j = x01 - x00;
                double dy_j = y01 - y00;
                double nx_j =  dy_j;
                double ny_j = -dx_j;
                faceNormal_j[idx] = {nx_j, ny_j};
                faceLen_j[idx]    = std::sqrt(dx_j*dx_j + dy_j*dy_j);

                // DEFAULT: oriented normals = geometric normals
                faceNormal_i_out[idx] = faceNormal_i[idx];
                faceNormal_j_out[idx] = faceNormal_j[idx];
            }
        }
    }
        
    // ------------------------------------------------------------
    // Make i-faces radial and j-faces continuous & periodic
    // ------------------------------------------------------------



    void correctNormalOrientations() {
        int imin, imax, jmin, jmax;
        getInteriorBounds(imin, imax, jmin, jmax);

        // 1) Estimate center from the wall row
        double xc0 = 0.0, yc0 = 0.0;
        int count = 0;
        for (int i = imin; i < imax; ++i) {
            int idx = cellIndex(i, jmin);
            xc0 += xCells[idx];
            yc0 += yCells[idx];
            ++count;
        }
        if (count > 0) {
            xc0 /= count;
            yc0 /= count;
        }
        
        // 2) Force i-faces (in oriented version) to be radial outward
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = cellIndex(i, j);

                double rx  = xCells[idx] - xc0;
                double ry  = yCells[idx] - yc0;

                // start from geometric i-normal
                double nx = faceNormal_i[idx][0];
                double ny = faceNormal_i[idx][1];

                double dot = nx * rx + ny * ry;
                if (dot < 0.0) {
                    nx = -nx;
                    ny = -ny;
                }
                faceNormal_i_out[idx][0] = nx;
                faceNormal_i_out[idx][1] = ny;
            }
        }

        // 3) j-face oriented normals just copy geometric everywhere
        for (int j = 0; j < njTotal; ++j) {
            for (int i = 0; i < niTotal; ++i) {
                int idx = cellIndex(i, j);
                faceNormal_j_out[idx] = faceNormal_j[idx];
            }
        }
    }
        
    // ------------------------------------------------------------
    // Fill ghost cell geometry
    // ------------------------------------------------------------
    void fillGhostGeometry() {
        int imin, imax, jmin, jmax;
        getInteriorBounds(imin, imax, jmin, jmax);
        
    // Bottom ghosts (wall)
    for (int i = imin; i < imax; ++i) {
        int idx_i = cellIndex(i, jmin);
        for (int g = 1; g <= NGHOST; ++g) {
            int idx_g = cellIndex(i, jmin - g);

            cellArea[idx_g]     = cellArea[idx_i];
            faceNormal_i[idx_g] = faceNormal_i[idx_i];
            faceNormal_j[idx_g] = faceNormal_j[idx_i];
            faceLen_i[idx_g]    = faceLen_i[idx_i];
            faceLen_j[idx_g]    = faceLen_j[idx_i];

            // NEW: oriented normals for ghosts
            faceNormal_i_out[idx_g] = faceNormal_i_out[idx_i];
            faceNormal_j_out[idx_g] = faceNormal_j_out[idx_i];
        }
    }

        // Top ghosts (farfield)
        // IMPORTANT: keep the first ghost ring (j = jmax) as computed from node extrapolation.
        // Only copy its geometry to deeper ghost cells.
        // Top ghosts (farfield)
        for (int i = imin; i < imax; ++i) {
            int idx_ref = cellIndex(i, jmax); // first ghost
            for (int g = 1; g < NGHOST; ++g) {
                int idx_g = cellIndex(i, jmax + g);
                cellArea[idx_g]     = cellArea[idx_ref];
                faceNormal_i[idx_g] = faceNormal_i[idx_ref];
                faceNormal_j[idx_g] = faceNormal_j[idx_ref];
                faceLen_i[idx_g]    = faceLen_i[idx_ref];
                faceLen_j[idx_g]    = faceLen_j[idx_ref];

                // NEW: oriented normals
                faceNormal_i_out[idx_g] = faceNormal_i_out[idx_ref];
                faceNormal_j_out[idx_g] = faceNormal_j_out[idx_ref];
            }
        }
    }
    
    // ------------------------------------------------------------
    // Apply periodic boundaries
    // ------------------------------------------------------------
    void applyPeriodicBoundaries() {
        int imin, imax, jmin, jmax;
        getInteriorBounds(imin, imax, jmin, jmax);
        
        // Copy geometry periodically in i-direction for all j
        for (int j = 0; j < njTotal; ++j) {
            // Left ghosts from right
            // Left ghosts from right
            for (int g = 1; g <= NGHOST; ++g) {
                int i_ghost = imin - g;
                int i_donor = imax - g;
                int idx_g = cellIndex(i_ghost, j);
                int idx_d = cellIndex(i_donor, j);
                
                cellArea[idx_g]     = cellArea[idx_d];
                faceNormal_i[idx_g] = faceNormal_i[idx_d];
                faceNormal_j[idx_g] = faceNormal_j[idx_d];
                faceLen_i[idx_g]    = faceLen_i[idx_d];
                faceLen_j[idx_g]    = faceLen_j[idx_d];

                // NEW
                faceNormal_i_out[idx_g] = faceNormal_i_out[idx_d];
                faceNormal_j_out[idx_g] = faceNormal_j_out[idx_d];
            }

            // Right ghosts from left
            for (int g = 0; g < NGHOST; ++g) {
                int i_ghost = imax + g;
                int i_donor = imin + g;
                int idx_g = cellIndex(i_ghost, j);
                int idx_d = cellIndex(i_donor, j);
                
                cellArea[idx_g]     = cellArea[idx_d];
                faceNormal_i[idx_g] = faceNormal_i[idx_d];
                faceNormal_j[idx_g] = faceNormal_j[idx_d];
                faceLen_i[idx_g]    = faceLen_i[idx_d];
                faceLen_j[idx_g]    = faceLen_j[idx_d];

                // NEW
                faceNormal_i_out[idx_g] = faceNormal_i_out[idx_d];
                faceNormal_j_out[idx_g] = faceNormal_j_out[idx_d];
            }
        }
    }
    
    // ------------------------------------------------------------
    // Compute averaged quantities for wave speed
    // ------------------------------------------------------------
    void computeAveragedQuantities() {
        int imin, imax, jmin, jmax;
        getInteriorBounds(imin, imax, jmin, jmax);

        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = cellIndex(i, j);

                // Average i-direction metrics
                if (i > imin && i < imax - 1) {
                    int idxm1 = cellIndex(i - 1, j);
                    faceLen_i_avg[idx] = 0.5 * (faceLen_i[idx] + faceLen_i[idxm1]);
                    faceNormal_i_avg[idx][0] = 0.5 * (faceNormal_i[idx][0] + faceNormal_i[idxm1][0]);
                    faceNormal_i_avg[idx][1] = 0.5 * (faceNormal_i[idx][1] + faceNormal_i[idxm1][1]);
                } else {
                    faceLen_i_avg[idx]    = faceLen_i[idx];
                    faceNormal_i_avg[idx] = faceNormal_i[idx];
                }

                // Average j-direction metrics
                if (j > jmin && j < jmax - 1) {
                    int idxm1 = cellIndex(i, j - 1);
                    faceLen_j_avg[idx] = 0.5 * (faceLen_j[idx] + faceLen_j[idxm1]);
                    faceNormal_j_avg[idx][0] = 0.5 * (faceNormal_j[idx][0] + faceNormal_j[idxm1][0]);
                    faceNormal_j_avg[idx][1] = 0.5 * (faceNormal_j[idx][1] + faceNormal_j[idxm1][1]);
                } else {
                    faceLen_j_avg[idx]    = faceLen_j[idx];
                    faceNormal_j_avg[idx] = faceNormal_j[idx];
                }
            }
        }
    }    
    // ------------------------------------------------------------
    // Compute spectral radii for local time stepping
    // ------------------------------------------------------------
    void computeSpectralRadii(const std::vector<Conservative>& U, const Config& cfg) {
        int imin, imax, jmin, jmax;
        getInteriorBounds(imin, imax, jmin, jmax);

        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int idx = cellIndex(i, j);
                Primitive W(U[idx], cfg.gamma);

                // i-direction
                const auto& n_i = faceNormal_i_avg[idx];
                double nx_i = n_i[0] / (faceLen_i_avg[idx] + 1e-14);
                double ny_i = n_i[1] / (faceLen_i_avg[idx] + 1e-14);
                double un_i = W.u * nx_i + W.v * ny_i;
                lambda_i[idx] = (std::abs(un_i) + W.a) * faceLen_i_avg[idx];

                // j-direction
                const auto& n_j = faceNormal_j_avg[idx];
                double nx_j = n_j[0] / (faceLen_j_avg[idx] + 1e-14);
                double ny_j = n_j[1] / (faceLen_j_avg[idx] + 1e-14);
                double un_j = W.u * nx_j + W.v * ny_j;
                lambda_j[idx] = (std::abs(un_j) + W.a) * faceLen_j_avg[idx];
            }
        }
    }
};