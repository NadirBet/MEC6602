// Flux.hpp
#pragma once
#include <vector>
#include <array>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstddef>   // size_t
#include "Types.hpp"
#include "mesh.hpp"

class FluxCalculator {
public:
    // -----------------------------------------------------------
    // JST diagnostics (captured over the last computeResidual())
    // - max_eps*: running maxima of actual eps2/eps4 used
    // - cnt_eps*: number of faces where eps > COUNT_EPS_TOL
    // -----------------------------------------------------------

    struct JSTStats {
        double max_eps2_I = 0.0, max_eps4_I = 0.0;
        double max_eps2_J = 0.0, max_eps4_J = 0.0;
        std::size_t cnt_eps2_I = 0, cnt_eps4_I = 0;
        std::size_t cnt_eps2_J = 0, cnt_eps4_J = 0;

        // NEW: L1 totals of JST dissipation (face-integrated)
        double D2_L1_total = 0.0;   // sum over faces of ds * ||D2||_1
        double D4_L1_total = 0.0;   // sum over faces of ds * ||D4||_1

        void reset() {
            max_eps2_I = max_eps4_I = max_eps2_J = max_eps4_J = 0.0;
            cnt_eps2_I = cnt_eps4_I = cnt_eps2_J = cnt_eps4_J = 0;
            // NEW:
            D2_L1_total = 0.0;
            D4_L1_total = 0.0;
        }

        std::size_t cnt_eps2_total() const { return cnt_eps2_I + cnt_eps2_J; }
        std::size_t cnt_eps4_total() const { return cnt_eps4_I + cnt_eps4_J; }
    };

    // Threshold for counting eps>0 (ignore FP dust)
    static constexpr double COUNT_EPS_TOL = 1e-12;


    FluxCalculator(double k2_ = 0.5, double k4_ = 0.02, bool skip_wall_flux_ = false) 
        : skip_wall_pressure_flux(skip_wall_flux_) {
        (void)k2_; (void)k4_; // cfg carries k2/k4; keep ctor for compatibility
    }

    const JSTStats& getJSTStats() const { return jstStats; }
    void resetJSTStats() { jstStats.reset(); }


    // Central (Rusanov-free) Euler flux, per unit length (NO ds here).
    // WL, WR are Primitive states; n is a UNIT normal; gamma is EOS γ.
    static inline Conservative centralFlux(
        const Primitive& WL,
        const Primitive& WR,
        const std::array<double,2>& n,
        double gamma
    ) {
        const double nx = n[0], ny = n[1];

        const double unL = WL.u * nx + WL.v * ny;
        const double unR = WR.u * nx + WR.v * ny;

        // total enthalpy H = h + 0.5|u|^2, with h = γ/(γ−1) * p/ρ
        const double qL2 = WL.u*WL.u + WL.v*WL.v;
        const double qR2 = WR.u*WR.u + WR.v*WR.v;
        const double hL  = (gamma/(gamma-1.0)) * (WL.p / WL.rho);
        const double hR  = (gamma/(gamma-1.0)) * (WR.p / WR.rho);
        const double HL  = hL + 0.5*qL2;
        const double HR  = hR + 0.5*qR2;

        Conservative F;
        // mass
        F.rho  = 0.5 * (WL.rho * unL + WR.rho * unR);
        // momentum
        F.rhou = 0.5 * (WL.rho*WL.u*unL + WR.rho*WR.u*unR) + 0.5*(WL.p + WR.p)*nx;
        F.rhov = 0.5 * (WL.rho*WL.v*unL + WR.rho*WR.v*unR) + 0.5*(WL.p + WR.p)*ny;
        // energy
        F.rhoE = 0.5 * (WL.rho * HL * unL + WR.rho * HR * unR);

        return F; // NOTE: caller multiplies by ds
    }



    // --------------------------------------------------------------
    // Assemble residual R from shared faces (Central + JST)
    // Assumes: mesh.iFaceNormal/jFaceNormal are UNIT normals;
    //          mesh.iFaceLen/jFaceLen are face lengths (ds).
    // --------------------------------------------------------------
    void computeResidual(
        const std::vector<Conservative>& U,
        std::vector<Conservative>& R,
        const Mesh& mesh,
        const Config& cfg
    ) {
        // zero residuals and reset JST stats
        for (auto& r : R) r = Conservative();
        jstStats.reset();

        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);

        auto isInterior = [&](int i, int j) {
            return (i >= imin && i < imax && j >= jmin && j < jmax);
        };

        // =========================
        // 1) I-FACES (horizontal, normal points UP)
        // faces (i,j), j = jmin..jmax
        // cDown = (i, j-1), cUp = (i, j)
        // =========================
        for (int j = jmin; j <= jmax; ++j) {
            for (int i = imin; i < imax; ++i) {


                 // ---- Bottom slip wall (I-face at j=jmin): pressure-only traction, no convective flux ----
 
                /**   if (!skip_wall_pressure_flux && j == jmin) {
                    // 1) Extrapolate p at the wall from interior cells
                    const int c2 = mesh.cellIndex(i, jmin    );
                    const int c3 = mesh.cellIndex(i, jmin + 1);
                    const int c4 = mesh.cellIndex(i, jmin + 2);

                    const double p2 = Primitive(U[c2], cfg.gamma).p;
                    const double p3 = Primitive(U[c3], cfg.gamma).p;
                    const double p4 = Primitive(U[c4], cfg.gamma).p;

                    const bool have3   = (jmin + 2 < jmax);
                    /**   const double p_wall = have3
                        ? 0.125*(15.0*p2 - 10.0*p3 + 3.0*p4)
                        : p2;
                    */ /**
                   const double p_wall = p2;
                    // 2) Wall face geometry
                    const int f_wall = mesh.iFaceIndex(i, jmin);
                    double nx = mesh.iFaceNormal[f_wall][0];
                    double ny = mesh.iFaceNormal[f_wall][1];
                    const double ds = mesh.iFaceLen[f_wall];

                    // 3) Force normal to be outward from interior cell
                    const auto nidx = [&](int I, int J){ return mesh.nodeIndex(I,J); };
                    const double xA = mesh.xNodes[nidx(i,   jmin)];
                    const double yA = mesh.yNodes[nidx(i,   jmin)];
                    const double xB = mesh.xNodes[nidx(i+1, jmin)];
                    const double yB = mesh.yNodes[nidx(i+1, jmin)];
                    const double mx = 0.5*(xA + xB);
                    const double my = 0.5*(yA + yB);

                    const double cx = mesh.x(i, jmin);
                    const double cy = mesh.y(i, jmin);

                    const double vx = mx - cx;
                    const double vy = my - cy;

                    if (nx*vx + ny*vy < 0.0) { nx = -nx; ny = -ny; }

                    // 4) Pressure-only traction on wall cell
                    Conservative Fw{};
                    Fw.rho  = 0.0;
                    Fw.rhou = +p_wall * nx * ds;
                    Fw.rhov = +p_wall * ny * ds;
                    Fw.rhoE =  0.0;  // stationary adiabatic wall

                    const int cUp = mesh.cellIndex(i, jmin);  // interior cell
                    R[cUp].rho  += Fw.rho;
                    R[cUp].rhou += Fw.rhou;
                    R[cUp].rhov += Fw.rhov;
                    R[cUp].rhoE += Fw.rhoE;

                    // Important: skip the “interior” flux handling for this face
                    continue;
                }*/ 

                // ---- END bottom wall override ----

                const int cDown = (j == jmin) ? mesh.cellIndex(i, jmin - 1)
                                              : mesh.cellIndex(i, j - 1);
                const int cUp   = (j == jmax) ? mesh.cellIndex(i, jmax)
                                              : mesh.cellIndex(i, j);

                const int f  = mesh.iFaceIndex(i, j);
                const auto& n = mesh.iFaceNormal[f];   // unit
                const double ds = mesh.iFaceLen[f];

                const Conservative F = computeFaceFlux(
                    U[cDown], U[cUp], n, ds, cfg,
                    /*alongI=*/true,  i, j, U, mesh
                );

                const int jDown = (j == jmin) ? (jmin - 1) : (j - 1);
                const bool downInterior = isInterior(i, jDown);
                const bool upInterior   = isInterior(i, (j == jmax) ? jmax : j);

                if (downInterior) {
                    R[cDown].rho  -= F.rho;
                    R[cDown].rhou -= F.rhou;
                    R[cDown].rhov -= F.rhov;
                    R[cDown].rhoE -= F.rhoE;
                }
                if (upInterior) {
                    R[cUp].rho  += F.rho;
                    R[cUp].rhou += F.rhou;
                    R[cUp].rhov += F.rhov;
                    R[cUp].rhoE += F.rhoE;
                }
            }
        }

        // =========================
        // 2) J-FACES (vertical, normal points RIGHT)
        // faces (i,j), i = imin..imax-1  (seam handled at i==imin)
        // cLeft = (iLphys, j), cRight = (iRphys, j) with periodic wrap
        // =========================
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {

                const int iLphys = (i == imin) ? (imax - 1) : (i - 1);
                const int iRphys = i;

                const int cLeft  = mesh.cellIndex(iLphys, j);
                const int cRight = mesh.cellIndex(iRphys, j);

                const int f  = mesh.jFaceIndex(i, j);
                const auto& n = mesh.jFaceNormal[f];   // unit
                const double ds = mesh.jFaceLen[f];

                const Conservative F = computeFaceFlux(
                    U[cLeft], U[cRight], n, ds, cfg,
                    /*alongI=*/false, i, j, U, mesh
                );

                // conservative two-sided update on PHYSICAL cells
                R[cLeft].rho  -= F.rho;   R[cRight].rho  += F.rho;
                R[cLeft].rhou -= F.rhou;  R[cRight].rhou += F.rhou;
                R[cLeft].rhov -= F.rhov;  R[cRight].rhov += F.rhov;
                R[cLeft].rhoE -= F.rhoE;  R[cRight].rhoE += F.rhoE;

                /** // seam diagnostics
                if (i == imin) {
                    // per-face conservation check
                    const double s_rho  = (-F.rho)  + (F.rho);
                    const double s_rhou = (-F.rhou) + (F.rhou);
                    const double s_rhov = (-F.rhov) + (F.rhov);
                    const double s_rhoE = (-F.rhoE) + (F.rhoE);

                    static bool seam_cons_header = false;
                    if (!seam_cons_header) {
                        std::printf("---- i-periodic seam per-face conservation (should be ~0) ----\n");
                        std::printf("%4s %12s %12s %12s %12s\n","j","sum_Frho","sum_Frhou","sum_Frhov","sum_FrhoE");
                        seam_cons_header = true;
                    }
                    std::printf("%4d % .3e % .3e % .3e % .3e\n",
                                j, s_rho, s_rhou, s_rhov, s_rhoE);

                    // seam state diffs and flux
                    static bool seam_header = false;
                    if (!seam_header) {
                        std::printf("---- i-periodic seam diagnostics (face at i=imin=%d) ----\n", imin);
                        std::printf("%4s %9s %9s | %9s %9s | %9s %9s %9s %9s | %12s %12s %12s %12s\n",
                                    "j","nx","ny","ds","|dU|","d_rho","d_u","d_v","d_p",
                                    "Frho","Frhou","Frhov","FrhoE");
                        seam_header = true;
                    }
                    const Primitive WL(U[cLeft],  cfg.gamma);
                    const Primitive WR(U[cRight], cfg.gamma);

                    const double d_rho = WR.rho - WL.rho;
                    const double d_u   = WR.u   - WL.u;
                    const double d_v   = WR.v   - WL.v;
                    const double d_p   = WR.p   - WL.p;
                    const double dUmag = std::sqrt(d_rho*d_rho + d_u*d_u + d_v*d_v + d_p*d_p);

                    std::printf("%4d % .3e % .3e | % .3e % .3e | % .3e % .3e % .3e % .3e | % .6e % .6e % .6e % .6e\n",
                                j, n[0], n[1], ds, dUmag, d_rho, d_u, d_v, d_p,
                                F.rho, F.rhou, F.rhov, F.rhoE);
                }*/
            }
        }

        // ---- DEBUG: bottom row residual Linf per component ----
        /**{
            int imin_, imax_, jmin_, jmax_;
            mesh.getInteriorBounds(imin_, imax_, jmin_, jmax_);

            double linf_rho  = 0.0;
            double linf_rhou = 0.0;
            double linf_rhov = 0.0;
            double linf_rhoE = 0.0;

            for (int i = imin_; i < imax_; ++i) {
                const int c = mesh.cellIndex(i, jmin_);
                linf_rho  = std::max(linf_rho,  std::abs(R[c].rho));
                linf_rhou = std::max(linf_rhou, std::abs(R[c].rhou));
                linf_rhov = std::max(linf_rhov, std::abs(R[c].rhov));
                linf_rhoE = std::max(linf_rhoE, std::abs(R[c].rhoE));
            }
            for (int i = imin_; i < imax_; ++i) {
                const int c = mesh.cellIndex(i, jmin_);
                std::printf("[DEBUG] bottom R(i=%d): rho=% .3e rhou=% .3e rhov=% .3e rhoE=% .3e\n",
                            i, R[c].rho, R[c].rhou, R[c].rhov, R[c].rhoE);
            }
            std::printf("[DEBUG] bottom-row Linf |R|: rho=%.3e rhou=%.3e rhov=%.3e rhoE=%.3e\n",
                        linf_rho, linf_rhou, linf_rhov, linf_rhoE);
        }

        // ---- DEBUG: seam conservation check (sum of residuals across seam) ----
        {
            int imin_, imax_, jmin_, jmax_;
            mesh.getInteriorBounds(imin_, imax_, jmin_, jmax_);

            double linf_sum_rho=0, linf_sum_rhou=0, linf_sum_rhov=0, linf_sum_rhoE=0;

            for (int j = jmin_; j < jmax_; ++j) {
                const int cL = mesh.cellIndex(imax_-1, j);  // left of seam
                const int cR = mesh.cellIndex(imin_,   j);  // right of seam

                const double s_rho  = R[cL].rho  + R[cR].rho;
                const double s_rhou = R[cL].rhou + R[cR].rhou;
                const double s_rhov = R[cL].rhov + R[cR].rhov;
                const double s_rhoE = R[cL].rhoE + R[cR].rhoE;

                linf_sum_rho  = std::max(linf_sum_rho,  std::abs(s_rho));
                linf_sum_rhou = std::max(linf_sum_rhou, std::abs(s_rhou));
                linf_sum_rhov = std::max(linf_sum_rhov, std::abs(s_rhov));
                linf_sum_rhoE = std::max(linf_sum_rhoE, std::abs(s_rhoE));
            }
            std::printf("[DEBUG] seam sum(L+R) Linf |R|: rho=%.3e rhou=%.3e rhov=%.3e rhoE=%.3e\n",
                        linf_sum_rho, linf_sum_rhou, linf_sum_rhov, linf_sum_rhoE);
        }*/
    }

private:
    // -------------------------------------------------------
    // Central flux + JST scalar dissipation on a shared face.
    // n: UNIT normal; ds: face length.
    // alongI: true  => I-face (stencil varies j)
    //         false => J-face (stencil varies i)
    // i_face,j_face: face indices in their respective face grids.
    // -------------------------------------------------------
    Conservative computeFaceFlux(
        const Conservative& UL,
        const Conservative& UR,
        const std::array<double,2>& n,  // unit
        double ds,
        const Config& cfg,
        bool alongI,
        int i_face, int j_face,
        const std::vector<Conservative>& U,
        const Mesh& mesh
    ) const
    {
        const double nx = n[0];
        const double ny = n[1];

        const Primitive WL(UL, cfg.gamma);
        const Primitive WR(UR, cfg.gamma);

        const double unL = WL.u * nx + WL.v * ny;
        const double unR = WR.u * nx + WR.v * ny;

        const double HL = EOS::totalEnthalpy(UL, cfg.gamma);
        const double HR = EOS::totalEnthalpy(UR, cfg.gamma);

        Conservative F;

        // Mass
        F.rho  = 0.5 * (WL.rho * unL + WR.rho * unR) * ds;

        // Momentum-x
        F.rhou = 0.5 * (WL.rho * WL.u * unL + WR.rho * WR.u * unR) * ds
               + 0.5 * (WL.p + WR.p) * nx * ds;

        // Momentum-y
        F.rhov = 0.5 * (WL.rho * WL.v * unL + WR.rho * WR.v * unR) * ds
               + 0.5 * (WL.p + WR.p) * ny * ds;

        // Energy
        F.rhoE = 0.5 * (WL.rho * HL * unL + WR.rho * HR * unR) * ds;

        // JST dissipation (scalar)

        const Conservative D = computeJSTDissipation(
            U, mesh, cfg, alongI, i_face, j_face, ds
        );

        F.rho  -= D.rho;
        F.rhou -= D.rhou;
        F.rhov -= D.rhov;
        F.rhoE -= D.rhoE;

        return F;
    }


    // -----------------------------------------------------------------
    // JST scalar dissipation at a face:
    // D = λ̂ [ ε2 (U_{+1} - U_0) - ε4 (U_{+2} - 3U_{+1} + 3U_0 - U_{-1}) ]
    // Direction per alongI (I-face varies j), alongJ (J-face varies i).
    // λ̂ is a face-averaged scalar spectral radius from mesh.lambdaI/J.
    //
    // NEW: we split D into D2 and D4 parts to accumulate L1 totals:
    //   D2_L1_total += ds * (|D2.rho| + |D2.rhou| + |D2.rhov| + |D2.rhoE|)
    //   D4_L1_total += ds * (|D4.rho| + |D4.rhou| + |D4.rhov| + |D4.rhoE|)
    // -----------------------------------------------------------------
    Conservative computeJSTDissipation(
        const std::vector<Conservative>& U,
        const Mesh& mesh,
        const Config& cfg,
        bool alongI,
        int i, int j,
        double ds          // NEW: needed to make L1 totals face-integrated
    ) const
    {
        const int imin = NGHOST;
        const int imax = NGHOST + mesh.ni;
        const int jmin = NGHOST;
        const int jmax = NGHOST + mesh.nj;

        auto wrapI = [&](int ii)->int {
            if (ii < 0) return ii + mesh.niTotal;
            if (ii >= mesh.niTotal) return ii - mesh.niTotal;
            return ii;
        };
        auto mirrorJ = [&](int jj)->int {
            if (jj <  jmin) return 2*jmin - 1 - jj;
            if (jj >= jmax) return 2*jmax - 1 - jj;
            return jj;
        };
        auto idxPhys = [&](int ii, int jj)->int {
            return mesh.cellIndex(wrapI(ii), mirrorJ(jj));
        };

        const int di = alongI ? 0 : 1;
        const int dj = alongI ? 1 : 0;

        const int iLL = i - 2*di, jLL = j - 2*dj;
        const int iL  = i - 1*di, jL  = j - 1*dj;
        const int iR  = i,        jR  = j;
        const int iRR = i + 1*di, jRR = j + 1*dj;

        const Conservative& ULL = U[idxPhys(iLL, jLL)];
        const Conservative& UL  = U[idxPhys(iL , jL )];
        const Conservative& UR  = U[idxPhys(iR , jR )];
        const Conservative& URR = U[idxPhys(iRR, jRR)];

        auto p_of = [&](const Conservative& C)->double {
            return Primitive(C, cfg.gamma).p;
        };
        auto sens_j = [&](int ii, int jj)->double {
            const double pm1 = p_of( U[idxPhys(ii, jj-1)] );
            const double p0  = p_of( U[idxPhys(ii, jj  )] );
            const double pp1 = p_of( U[idxPhys(ii, jj+1)] );
            const double num = std::fabs(pp1 - 2.0*p0 + pm1);
            const double den = std::fabs(pp1) + 2.0*std::fabs(p0) + std::fabs(pm1) + 1e-14;
            return num / den;
        };
        auto sens_i = [&](int ii, int jj)->double {
            const double pm1 = p_of( U[idxPhys(ii-1, jj)] );
            const double p0  = p_of( U[idxPhys(ii,   jj)] );
            const double pp1 = p_of( U[idxPhys(ii+1, jj)] );
            const double num = std::fabs(pp1 - 2.0*p0 + pm1);
            const double den = std::fabs(pp1) + 2.0*std::fabs(p0) + std::fabs(pm1) + 1e-14;
            return num / den;
        };

        double eps2 = 0.0, eps4 = 0.0;
        if (alongI) {
            const double YL = sens_j(i, j-1);
            const double YR = sens_j(i, j  );
            eps2 = cfg.k2_jst * std::max(YL, YR);
        } else {
            const double YL = sens_i(i-1, j);
            const double YR = sens_i(i,   j);
            eps2 = cfg.k2_jst * std::max(YL, YR);
        }
        eps4 = std::max(0.0, cfg.k4_jst - eps2);

        // Disable D4 if 4-pt stencil not available along j (near top/bottom)
        if (alongI) {
            const bool hasD4 = ((j - 2) >= jmin) && ((j + 1) < jmax);
            if (!hasD4) eps4 = 0.0;
        }

        // Face-averaged scalar spectral radius
        double lambda_hat = 0.0;
        if (alongI) {
            const int cDown = mesh.cellIndex(i, mirrorJ(j-1));
            const int cUp   = mesh.cellIndex(i, mirrorJ(j  ));
            lambda_hat = 0.5 * (mesh.lambdaI[cDown] + mesh.lambdaI[cUp]);
        } else {
            const int iLcell = (i == imin) ? (imax - 1) : (i - 1);
            const int iRcell = (i == imax) ? (imin)     : (i);
            const int cLeft  = mesh.cellIndex(wrapI(iLcell), j);
            const int cRight = mesh.cellIndex(wrapI(iRcell), j);
            lambda_hat = 0.5 * (mesh.lambdaJ[cLeft] + mesh.lambdaJ[cRight]);
        }

        // --- JST instrumentation (eps maxima + counts) ---
        if (alongI) {
            if (eps2 > jstStats.max_eps2_I) jstStats.max_eps2_I = eps2;
            if (eps4 > jstStats.max_eps4_I) jstStats.max_eps4_I = eps4;
            if (eps2 > COUNT_EPS_TOL) ++jstStats.cnt_eps2_I;
            if (eps4 > COUNT_EPS_TOL) ++jstStats.cnt_eps4_I;
        } else {
            if (eps2 > jstStats.max_eps2_J) jstStats.max_eps2_J = eps2;
            if (eps4 > jstStats.max_eps4_J) jstStats.max_eps4_J = eps4;
            if (eps2 > COUNT_EPS_TOL) ++jstStats.cnt_eps2_J;
            if (eps4 > COUNT_EPS_TOL) ++jstStats.cnt_eps4_J;
        }

        // --- Build D2 and D4 separately, then combine ---
        Conservative D2; // zero
        Conservative D4; // zero

        if (eps2 > 0.0) {
            const double a2 = eps2 * lambda_hat;
            D2.rho  = a2 * (UR.rho  - UL.rho );
            D2.rhou = a2 * (UR.rhou - UL.rhou);
            D2.rhov = a2 * (UR.rhov - UL.rhov);
            D2.rhoE = a2 * (UR.rhoE - UL.rhoE);
        }
        if (eps4 > 0.0) {
            const double a4 = eps4 * lambda_hat;
            // NOTE the minus sign per JST formula
            D4.rho  = -a4 * (URR.rho  - 3.0*UR.rho  + 3.0*UL.rho  - ULL.rho );
            D4.rhou = -a4 * (URR.rhou - 3.0*UR.rhou + 3.0*UL.rhou - ULL.rhou);
            D4.rhov = -a4 * (URR.rhov - 3.0*UR.rhov + 3.0*UL.rhov - ULL.rhov);
            D4.rhoE = -a4 * (URR.rhoE - 3.0*UR.rhoE + 3.0*UL.rhoE - ULL.rhoE);
        }

        // NEW: L1 accumulation (face-integrated with ds)
        auto L1 = [](const Conservative& C)->double {
            return std::fabs(C.rho) + std::fabs(C.rhou) + std::fabs(C.rhov) + std::fabs(C.rhoE);
        };
        jstStats.D2_L1_total += ds * L1(D2);
        jstStats.D4_L1_total += ds * L1(D4);

        // Return combined dissipation (what subtracts from central flux)
        Conservative D;
        D.rho  = D2.rho  + D4.rho;
        D.rhou = D2.rhou + D4.rhou;
        D.rhov = D2.rhov + D4.rhov;
        D.rhoE = D2.rhoE + D4.rhoE;
        return D;
    }

private:
    mutable JSTStats jstStats;
    bool skip_wall_pressure_flux; 
};