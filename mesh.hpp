#pragma once
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include "Types.hpp"
constexpr int NGHOST = 2;

// ============================================================
// Structured O-grid-like mesh with ghost layers
// - i-direction periodic
// - 2 ghost layers around
// - shared-face storage (i-faces and j-faces)
// - can write VTK to visualize
// ============================================================
class Mesh {
public:
    // physical cells
    int ni = 0, nj = 0;

    // total cells (with ghosts)
    int niTotal = 0, njTotal = 0;

    // nodes (cells+1)
    int niNodes = 0, njNodes = 0;

    // node coordinates
    std::vector<double> xNodes;
    std::vector<double> yNodes;

    // cell centers + areas (per CELL)
    std::vector<double> xCells;
    std::vector<double> yCells;
    std::vector<double> cellArea;

    // ---------- shared faces ----------
    // i-faces: horizontal, from node(i,j) -> node(i+1,j), j = 0..njTotal
    // count = niTotal * (njTotal+1)
    std::vector<std::array<double,2>> iFaceNormal;
    std::vector<double>              iFaceLen;

    // j-faces: vertical, from node(i,j) -> node(i,j+1), i = 0..niNodes-1
    // count = (niTotal+1) * njTotal = niNodes * njTotal
    std::vector<std::array<double,2>> jFaceNormal;
    std::vector<double>              jFaceLen;
    std::vector<double> lambdaI;  // spectral radius in i-direction
    std::vector<double> lambdaJ;  // spectral radius in j-direction
    std::vector<double> lambda;   // total spectral radius (lambdaI + lambdaJ)

    Mesh() = default;
    explicit Mesh(int ni_, int nj_) { allocate(ni_, nj_); }

    // --------------------------------------------------------
    // index helpers
    // --------------------------------------------------------
    inline int nodeIndex(int i, int j) const {
        return j * niNodes + i;
    }
    inline int cellIndex(int i, int j) const {
        return j * niTotal + i;
    }
    // shared-face indices
    // i-face: i = 0..niTotal-1, j = 0..njTotal
    inline int iFaceIndex(int i, int j) const {
        return j * niTotal + i;
    }
    // j-face: i = 0..niNodes-1, j = 0..njTotal-1
    inline int jFaceIndex(int i, int j) const {
        return j * niNodes + i;
    }

    // --------------------------------------------------------
    // allocate
    // --------------------------------------------------------
    void allocate(int ni_, int nj_) {
        ni = ni_;
        nj = nj_;
        niTotal = ni + 2*NGHOST;
        njTotal = nj + 2*NGHOST;
        niNodes = niTotal + 1;
        njNodes = njTotal + 1;

        int ncells = niTotal * njTotal;
        int nnodes = niNodes * njNodes;

        xNodes.assign(nnodes, 0.0);
        yNodes.assign(nnodes, 0.0);

        xCells.assign(ncells, 0.0);
        yCells.assign(ncells, 0.0);
        cellArea.assign(ncells, 0.0);

        // faces
        int nIFaces = niTotal * (njTotal + 1);
        int nJFaces = niNodes * njTotal;

        iFaceNormal.assign(nIFaces, {0.0, 0.0});
        iFaceLen.assign(nIFaces, 0.0);

        jFaceNormal.assign(nJFaces, {0.0, 0.0});
        jFaceLen.assign(nJFaces, 0.0);
    }

    // --------------------------------------------------------
    // interior bounds in CELL indices
    // --------------------------------------------------------
    void getInteriorBounds(int& imin, int& imax,
                           int& jmin, int& jmax) const {
        imin = NGHOST;
        imax = NGHOST + ni;
        jmin = NGHOST;
        jmax = NGHOST + nj;
    }

    // --------------------------------------------------------
    // read PLOT3D: 1 block, ASCII
    // --------------------------------------------------------

    void normalizePhysicalNodes(double Lref) {
        if (Lref <= 0.0) return;
        // scale *all* nodes we currently have
        for (double &x : xNodes) x /= Lref;
        for (double &y : yNodes) y /= Lref;
    }

    void readPlot3D(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open mesh file: " + filename);
        }

        int nBlocks;
        file >> nBlocks;
        if (nBlocks != 1) {
            throw std::runtime_error("Only single-block Plot3D supported");
        }

        int ni_nodes_phys, nj_nodes_phys;
        file >> ni_nodes_phys >> nj_nodes_phys;

        // physical cells = nodes - 1
        allocate(ni_nodes_phys - 1, nj_nodes_phys - 1);

        // read X, Y into the shifted physical region
        for (int j = 0; j < nj_nodes_phys; ++j) {
            for (int i = 0; i < ni_nodes_phys; ++i) {
                double x; file >> x;
                xNodes[nodeIndex(i + NGHOST, j + NGHOST)] = x;
            }
        }
        for (int j = 0; j < nj_nodes_phys; ++j) {
            for (int i = 0; i < ni_nodes_phys; ++i) {
                double y; file >> y;
                yNodes[nodeIndex(i + NGHOST, j + NGHOST)] = y;
            }
        }
        file.close();

        // ---- NEW: build a length scale from physical nodes ----
        double xmin = 1e30, xmax = -1e30;
        double ymin = 1e30, ymax = -1e30;
        for (int j = 0; j < nj_nodes_phys; ++j) {
            for (int i = 0; i < ni_nodes_phys; ++i) {
                int ii = i + NGHOST;
                int jj = j + NGHOST;
                double x = xNodes[nodeIndex(ii, jj)];
                double y = yNodes[nodeIndex(ii, jj)];
                if (x < xmin) xmin = x;
                if (x > xmax) xmax = x;
                if (y < ymin) ymin = y;
                if (y > ymax) ymax = y;
            }
        }
        double Lx = xmax - xmin;
        double Ly = ymax - ymin;
        double Lref = (Lx > Ly) ? Lx : Ly;
        if (Lref <= 0.0) Lref = 1.0;  // safety

        // scale all nodes we read
        normalizePhysicalNodes(Lref);

        // ghosts -> centers -> areas -> faces (all in normalized units)
        fillGhostNodes();
        computeCellCenters();
        computeCellAreas();
        computeFaceGeometry();
    }

    // --------------------------------------------------------
    // cartesian init (optional)
    // --------------------------------------------------------
    void initCartesianTest(double x0, double x1, double y0, double y1) {
        // assumes ni, nj already set
        allocate(ni, nj);

        double dx = (x1 - x0) / static_cast<double>(ni);
        double dy = (y1 - y0) / static_cast<double>(nj);

        for (int j = 0; j < njNodes; ++j) {
            double y = y0 + (j - NGHOST) * dy;
            for (int i = 0; i < niNodes; ++i) {
                double x = x0 + (i - NGHOST) * dx;
                xNodes[nodeIndex(i,j)] = x;
                yNodes[nodeIndex(i,j)] = y;
            }
        }

        computeCellCenters();
        computeCellAreas();
        computeFaceGeometry();
    }

    // --------------------------------------------------------
    // fill ghost nodes (fixed + clearer order)
    // --------------------------------------------------------
    void fillGhostNodes() {
        int imin = NGHOST;
        int imax = NGHOST + ni;
        int jmin = NGHOST;
        int jmax = NGHOST + nj;

        // 1) i-periodic first: so left/right ghost columns are valid
        for (int j = 0; j < njNodes; ++j) {
            // left ghosts: copy from rightmost physical nodes
            for (int g = 0; g < NGHOST; ++g) {
                int ig = imin - 1 - g;     // ghost cols: ..., 1, 0
                int is = imax - g;         // physical: imax, imax-1
                xNodes[nodeIndex(ig, j)] = xNodes[nodeIndex(is, j)];
                yNodes[nodeIndex(ig, j)] = yNodes[nodeIndex(is, j)];
            }
            // right ghosts: copy from leftmost physical nodes
            for (int g = 0; g < NGHOST; ++g) {
                int ig = imax + 1 + g;     // ghost cols: imax+1, imax+2...
                int is = imin + g;         // physical: imin, imin+1
                xNodes[nodeIndex(ig, j)] = xNodes[nodeIndex(is, j)];
                yNodes[nodeIndex(ig, j)] = yNodes[nodeIndex(is, j)];
            }
        }

        // 2) bottom ghosts (symmetric extrapolation, like top)
        for (int i = 0; i < niNodes; ++i) {
            for (int g = 0; g < NGHOST; ++g) {
                int jg = jmin - 1 - g;       // bottom ghost rows
                int j1 = jmin;               // first physical row
                int j2 = jmin + 1 + g;       // row above
                xNodes[nodeIndex(i, jg)] =
                    2.0 * xNodes[nodeIndex(i, j1)] - xNodes[nodeIndex(i, j2)];
                yNodes[nodeIndex(i, jg)] =
                    2.0 * yNodes[nodeIndex(i, j1)] - yNodes[nodeIndex(i, j2)];
            }
        }

        // 3) top ghosts (as before)
        for (int i = 0; i < niNodes; ++i) {
            for (int g = 0; g < NGHOST; ++g) {
                int jg = jmax + 1 + g;
                int j1 = jmax - g;
                int j2 = jmax - 1 - g;
                xNodes[nodeIndex(i, jg)] =
                    2.0 * xNodes[nodeIndex(i, j1)] - xNodes[nodeIndex(i, j2)];
                yNodes[nodeIndex(i, jg)] =
                    2.0 * yNodes[nodeIndex(i, j1)] - yNodes[nodeIndex(i, j2)];
            }
        }
    }

    // --------------------------------------------------------
    // cell centers
    // --------------------------------------------------------
    void computeCellCenters() {
        for (int j = 0; j < njTotal; ++j) {
            for (int i = 0; i < niTotal; ++i) {
                int n00 = nodeIndex(i,     j);
                int n10 = nodeIndex(i + 1, j);
                int n01 = nodeIndex(i,     j + 1);
                int n11 = nodeIndex(i + 1, j + 1);

                int c = cellIndex(i, j);
                xCells[c] = 0.25 * (xNodes[n00] + xNodes[n10] +
                                    xNodes[n01] + xNodes[n11]);
                yCells[c] = 0.25 * (yNodes[n00] + yNodes[n10] +
                                    yNodes[n01] + yNodes[n11]);
            }
        }
    }

    // --------------------------------------------------------
    // cell areas (shoelace)
    // --------------------------------------------------------
    void computeCellAreas() {
        for (int j = 0; j < njTotal; ++j) {
            for (int i = 0; i < niTotal; ++i) {
                int n00 = nodeIndex(i,     j);
                int n10 = nodeIndex(i + 1, j);
                int n01 = nodeIndex(i,     j + 1);
                int n11 = nodeIndex(i + 1, j + 1);

                double x00 = xNodes[n00], y00 = yNodes[n00];
                double x10 = xNodes[n10], y10 = yNodes[n10];
                double x01 = xNodes[n01], y01 = yNodes[n01];
                double x11 = xNodes[n11], y11 = yNodes[n11];

                double A = 0.5 * std::abs(
                    x00*(y10 - y01) +
                    x10*(y11 - y00) +
                    x11*(y01 - y10) +
                    x01*(y00 - y11)
                );
                cellArea[cellIndex(i,j)] = A;
            }
        }
    }

    // --------------------------------------------------------
    // shared-face geometry from nodes
    // --------------------------------------------------------
    void computeFaceGeometry() {
        // i-faces
        for (int j = 0; j <= njTotal; ++j) {
            for (int i = 0; i < niTotal; ++i) {
                int f = iFaceIndex(i,j);
                int n0 = nodeIndex(i,     j);
                int n1 = nodeIndex(i + 1, j);

                // i-faces
                double dx = xNodes[n1] - xNodes[n0];
                double dy = yNodes[n1] - yNodes[n0];
                double L  = std::sqrt(dx*dx + dy*dy);
                iFaceLen[f] = L;
                if (L > 0.0) {
                    iFaceNormal[f] = { -dy / L,  dx / L };  // unit normal
                } else {
                    iFaceNormal[f] = { 0.0, 0.0 };
                }
                
            }
        }

        // j-faces
        for (int j = 0; j < njTotal; ++j) {
            for (int i = 0; i < niNodes; ++i) {
                int f = jFaceIndex(i,j);
                int n0 = nodeIndex(i, j);
                int n1 = nodeIndex(i, j + 1);

                
                double dx = xNodes[n1] - xNodes[n0];
                double dy = yNodes[n1] - yNodes[n0];
                double L  = std::sqrt(dx*dx + dy*dy);
                jFaceLen[f] = L;
                if (L > 0.0) {
                    jFaceNormal[f] = {  dy / L, -dx / L };  // unit normal
                } else {
                    jFaceNormal[f] = { 0.0, 0.0 };
                }
            }
        }
    }



    void computeSpectralRadius(const std::vector<Primitive>& W, double gamma) {
        lambdaI.assign(niTotal * njTotal, 0.0);
        lambdaJ.assign(niTotal * njTotal, 0.0);
        lambda.assign(niTotal * njTotal, 0.0);

        int imin, imax, jmin, jmax;
        getInteriorBounds(imin, imax, jmin, jmax);

        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int c = cellIndex(i, j);

  
                // ----- I direction -----
                int fiL = iFaceIndex(i, j);
                int fiR = iFaceIndex(i, j + 1);

                double nxL = iFaceNormal[fiL][0], nyL = iFaceNormal[fiL][1];
                double nxR = iFaceNormal[fiR][0], nyR = iFaceNormal[fiR][1];
                double LiL = iFaceLen[fiL],        LiR = iFaceLen[fiR];

                double unL = std::abs(W[c].u * nxL + W[c].v * nyL);
                double unR = std::abs(W[c].u * nxR + W[c].v * nyR);
                double a   = W[c].a;

                lambdaI[c] = 0.5 * ( (unL + a) * LiL + (unR + a) * LiR );

                // ----- J direction -----
                int fjB = jFaceIndex(i,     j);
                int fjT = jFaceIndex(i + 1, j);

                double nxB = jFaceNormal[fjB][0], nyB = jFaceNormal[fjB][1];
                double nxT = jFaceNormal[fjT][0], nyT = jFaceNormal[fjT][1];
                double LjB = jFaceLen[fjB],        LjT = jFaceLen[fjT];

                double unB = std::abs(W[c].u * nxB + W[c].v * nyB);
                double unT = std::abs(W[c].u * nxT + W[c].v * nyT);

                lambdaJ[c] = 0.5 * ( (unB + a) * LjB + (unT + a) * LjT );

                lambda[c] = lambdaI[c] + lambdaJ[c];
            }
        }
    }




    // --------------------------------------------------------
    // write VTK legacy unstructured grid (QUADS)
    // --------------------------------------------------------
    void writeVTKPhysical(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot write VTK: " + filename);
        }

        int imin = NGHOST;
        int imax = NGHOST + ni;
        int jmin = NGHOST;
        int jmax = NGHOST + nj;

        // we can still write ALL nodes
        int nnodes = niNodes * njNodes;

        out << "# vtk DataFile Version 3.0\n";
        out << "mesh dump (physical)\n";
        out << "ASCII\n";
        out << "DATASET UNSTRUCTURED_GRID\n";

        out << "POINTS " << nnodes << " double\n";
        for (int j = 0; j < njNodes; ++j) {
            for (int i = 0; i < niNodes; ++i) {
                int n = nodeIndex(i,j);
                out << xNodes[n] << " " << yNodes[n] << " 0.0\n";
            }
        }

        int nPhysCells = (imax - imin) * (jmax - jmin);
        out << "CELLS " << nPhysCells << " " << nPhysCells * 5 << "\n";
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int n0 = nodeIndex(i,     j);
                int n1 = nodeIndex(i + 1, j);
                int n2 = nodeIndex(i + 1, j + 1);
                int n3 = nodeIndex(i,     j + 1);
                out << 4 << " " << n0 << " " << n1 << " " << n2 << " " << n3 << "\n";
            }
        }

        out << "CELL_TYPES " << nPhysCells << "\n";
        for (int c = 0; c < nPhysCells; ++c) out << 9 << "\n";

        out << "CELL_DATA " << nPhysCells << "\n";
        out << "SCALARS area double 1\n";
        out << "LOOKUP_TABLE default\n";
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                out << cellArea[cellIndex(i,j)] << "\n";
            }
        }

        out.close();
    }
    // --------------------------------------------------------
    // write all shared faces as POLYDATA with a vector = face normal
    // so ParaView can Glyph them

    void writeFaceVTKPhysical(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot write face VTK: " + filename);
        }

        int imin = NGHOST;
        int imax = NGHOST + ni;
        int jmin = NGHOST;
        int jmax = NGHOST + nj;

        // still write all nodes
        int nnodes = niNodes * njNodes;

        // number of physical i-faces: i = imin..imax-1, j = jmin..jmax (note jmax included)
        int nIFaces = (imax - imin) * (jmax - jmin + 1);
        // number of physical j-faces: i = imin..imax, j = jmin..jmax-1
        int nJFaces = (imax - imin + 1) * (jmax - jmin);
        int nFaces  = nIFaces + nJFaces;

        out << "# vtk DataFile Version 3.0\n";
        out << "physical faces\n";
        out << "ASCII\n";
        out << "DATASET POLYDATA\n";

        // points
        out << "POINTS " << nnodes << " double\n";
        for (int j = 0; j < njNodes; ++j)
            for (int i = 0; i < niNodes; ++i) {
                int n = nodeIndex(i,j);
                out << xNodes[n] << " " << yNodes[n] << " 0.0\n";
            }

        // lines
        out << "LINES " << nFaces << " " << nFaces * 3 << "\n";

        // physical i-faces
        for (int j = jmin; j <= jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int p0 = nodeIndex(i,     j);
                int p1 = nodeIndex(i + 1, j);
                out << 2 << " " << p0 << " " << p1 << "\n";
            }
        }
        // physical j-faces
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i <= imax; ++i) {
                int p0 = nodeIndex(i, j);
                int p1 = nodeIndex(i, j + 1);
                out << 2 << " " << p0 << " " << p1 << "\n";
            }
        }

        // cell data: normals
        out << "CELL_DATA " << nFaces << "\n";
        out << "VECTORS face_normal double\n";

        // same order as lines:
        // i-faces
        for (int j = jmin; j <= jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                int f = iFaceIndex(i, j);
                out << iFaceNormal[f][0] << " " << iFaceNormal[f][1] << " 0.0\n";
            }
        }
        // j-faces
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i <= imax; ++i) {
                int f = jFaceIndex(i, j);
                out << jFaceNormal[f][0] << " " << jFaceNormal[f][1] << " 0.0\n";
            }
        }

        out.close();
    }


    // --------------------------------------------------------
    // centers
    // --------------------------------------------------------
    double x(int i, int j) const { return xCells[cellIndex(i,j)]; }
    double y(int i, int j) const { return yCells[cellIndex(i,j)]; }
};
