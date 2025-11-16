#pragma once

#include <vector>
#include <string>
#include <deque>
#include <limits>
#include <cmath>

#include "mesh.hpp"
#include "Types.hpp"
#include "Flux.hpp"
#include "Boundary.hpp"

// -----------------------------------------------------------------------------
// FASMultigridSolver
//
// - Builds a hierarchy of O-grid levels by 2x coarsening from a *single* finest
//   mesh file (geometric multigrid).
// - Uses Explicit RK5 + local time stepping as smoother on *all* levels.
// - Full Approximation Storage (FAS):
//     * coarse initial guess from restricted fine solution,
//     * coarse forcing Q_F built from restricted fine residual,
//     * coarse residual used for diagnostics is R + Q_F (RK5 still uses PDE
//       residual internally via FluxCalculator).
// -----------------------------------------------------------------------------
class FASMultigridSolver {
public:
    // Per-level control parameters
    struct LevelControl {
        int    minSweeps       = 3;    // Minimum smoother iterations per visit
        int    maxSweeps       = 30;   // Hard cap per visit
        int    stallWindow     = 5;    // # of last residuals for stall detection
        double stallTol        = 0.02; // Stalled if residual drop < 2% over window
        double targetReduction = 0.5;  // R_new / R0 <= targetReduction
    };

    // One grid level in the multigrid hierarchy
    struct Level {
        int   index = 0;          // 0 = finest, L-1 = coarsest

        Mesh               mesh;
        Config             cfg;
        BoundaryConditions bc;
        FluxCalculator     flux;

        // State and residuals
        std::vector<Conservative> U;     // State (conservative variables)
        std::vector<Conservative> R;     // Residual (includes Q_F on coarse levels)
        std::vector<double>        dt;   // Local dt per cell

        // FAS forcing term Q_F on this level (zero on finest)
        std::vector<Conservative> forcing;

        // For prolongation / temporary buffers
        std::vector<Conservative> scratch;

        // Control + residual tracking
        LevelControl       ctrl;
        std::deque<double> resHistory;   // last residuals
        double             currentResL2 = std::numeric_limits<double>::infinity();
        double             bestResL2    = std::numeric_limits<double>::infinity();

        Level(const Config&       baseCfg,
              bool                useSimpleFarfield,
              bool                freestreamEverywhere,
              int                 levelIndex,
              const LevelControl& ctrlIn)
            : index(levelIndex)
            , mesh()
            , cfg(baseCfg)
            , bc(cfg, useSimpleFarfield, freestreamEverywhere)
            , flux(cfg.k2_jst, cfg.k4_jst, freestreamEverywhere)
            , U()
            , R()
            , dt()
            , forcing()
            , scratch()
            , ctrl(ctrlIn)
            , resHistory()
            , currentResL2(std::numeric_limits<double>::infinity())
            , bestResL2(std::numeric_limits<double>::infinity())
        {
        }
    };

    // -------------------------------------------------------------------------
    // Constructor:
    //
    //  baseCfg        : your Config (Mach, alpha, gamma, CFL, k2, k4, etc.)
    //  use_simple_ff  : same flag you pass to BoundaryConditions
    //  freestream_all : same as opts.freestream_everywhere
    //  finestMeshFile : ONE Plot3D file for the finest grid (e.g. "257x257.x")
    //  numLevels      : total number of multigrid levels (1=single grid)
    //  defaultCtrl    : base LevelControl, can be overridden per level
    // -------------------------------------------------------------------------
    FASMultigridSolver(const Config&      baseCfg,
                       bool               use_simple_ff,
                       bool               freestream_all,
                       const std::string& finestMeshFile,
                       int                numLevels,
                       const LevelControl& defaultCtrl);

    // Access to levels to tune ctrl per level, etc.
    std::vector<Level>&       getLevels()       { return levels; }
    const std::vector<Level>& getLevels() const { return levels; }

    // Initialize all levels:
    // - read finest mesh
    // - build coarser meshes by 2× coarsening
    // - initialize freestream state on each level
    // - compute initial spectral radii and forced residual
    void initialize();

    // Run up to n V-cycles, stop early if finest L2 < resTol
    void runVCycles(int nCycles, double resTol = 1e-12);

    // Rebuild residual / history on finest after an external modification of U
    void resetFinestFromCurrent();

    // Convenience: solution on finest grid
    std::vector<Conservative>&       finestSolution();
    const std::vector<Conservative>& finestSolution() const;

private:
    std::vector<Level> levels;
    bool               useSimpleFarfield;
    bool               freestreamEverywhere;
    int                nLevels = 0;
    std::string        finestMeshFile_;

    // ---------- core multigrid internals ----------

    // One V-cycle starting at level ell (0 = finest)
    void vCycle(int ell);

    // Smoother with min/max sweeps + stall detection
    void smoothWithControl(Level& L);

    // One smoothing sweep on a level:
    //  - ALL levels: RK5-LTS smoother using PDE residual.
    //    (FAS forcing enters via computeResidual() only.)
    void smootherSweep(Level& L);

    // Build local dt from lambdaI + lambdaJ
    void buildLocalDt(Level& L);

    // Compute forced residual on a level:
    //  - applies BC
    //  - recomputes spectral radii
    //  - flux residual R
    //  - adds Q_F on coarse levels
    //  - updates L.currentResL2
    void computeResidual(Level& L);

    // Volume-weighted L2 residual over physical cells
    double computeL2Residual(const Level& L) const;

    // Stall detection based on resHistory and ctrl.stallTol
    bool isStalled(const Level& L) const;

    // Restrict solution & residual fine → coarse (FAS construction of coarse level)
    void restrictSolution(int fineIdx, int coarseIdx);

    // Prolongate correction coarse → fine
    void prolongateCorrection(int coarseIdx, int fineIdx);

    // Shortcuts
    Level&       finestLevel()         { return levels.front(); }
    const Level& finestLevel() const   { return levels.front(); }
    Level&       coarsestLevel()       { return levels.back(); }
    const Level& coarsestLevel() const { return levels.back(); }
};
