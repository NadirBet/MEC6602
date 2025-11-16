#pragma once

#include <vector>
#include <Eigen/Dense>

#include "mesh.hpp"      // Mesh::getInteriorBounds, cellIndex(...)
#include "Types.hpp"    // struct Conservative { double rho,rhou,rhov,rhoE; }

// -----------------------------------------------------------------------------
// 1. Scalar Thomas solver for a tridiagonal system
//    a(i) * x(i-1) + b(i) * x(i) + c(i) * x(i+1) = d(i)
// -----------------------------------------------------------------------------
inline Eigen::VectorXd thomasAlgorithm(const Eigen::VectorXd& a,
                                       const Eigen::VectorXd& b,
                                       const Eigen::VectorXd& c,
                                       const Eigen::VectorXd& d)
{
    const int n = static_cast<int>(b.size());
    Eigen::VectorXd c_prime(n);
    Eigen::VectorXd d_prime(n);
    Eigen::VectorXd x(n);

    // First row
    c_prime(0) = c(0) / b(0);
    d_prime(0) = d(0) / b(0);

    // Forward sweep
    for (int i = 1; i < n; ++i) {
        const double m = b(i) - a(i) * c_prime(i - 1);
        c_prime(i) = c(i) / m;
        d_prime(i) = (d(i) - a(i) * d_prime(i - 1)) / m;
    }

    // Back substitution
    x(n - 1) = d_prime(n - 1);
    for (int i = n - 2; i >= 0; --i) {
        x(i) = d_prime(i) - c_prime(i) * x(i + 1);
    }

    return x;
}

// -----------------------------------------------------------------------------
// 2. Helper: access components of a Conservative state by index
//    comp = 0: rho, 1: rhou, 2: rhov, 3: rhoE
// -----------------------------------------------------------------------------
inline double getComponent(const Conservative& Q, int comp)
{
    switch (comp) {
    case 0: return Q.rho;
    case 1: return Q.rhou;
    case 2: return Q.rhov;
    default: return Q.rhoE;
    }
}

inline void setComponent(Conservative& Q, int comp, double value)
{
    switch (comp) {
    case 0: Q.rho  = value; break;
    case 1: Q.rhou = value; break;
    case 2: Q.rhov = value; break;
    default: Q.rhoE = value; break;
    }
}

// -----------------------------------------------------------------------------
// 3. Implicit residual smoothing module (structured I-then-J line solves)
// -----------------------------------------------------------------------------
class ImplicitResidualSmoothing
{
public:
    ImplicitResidualSmoothing(double epsI_, double epsJ_)
        : epsI(epsI_), epsJ(epsJ_) {}

    // R_in  : residuals on physical cells, size = mesh.niTotal * mesh.njTotal
    // R_out : smoothed residuals (output). Can alias R_in if desired.
    void smooth(const Mesh& mesh,
                const std::vector<Conservative>& R_in,
                std::vector<Conservative>&       R_out) const
    {
        // If no smoothing, just copy and return
        if (epsI <= 0.0 && epsJ <= 0.0) {
            R_out = R_in;
            return;
        }

        // Make sure output has correct size
        if (R_out.size() != R_in.size()) {
            R_out.resize(R_in.size());
        }

        // Intermediate storage after I-direction smoothing
        std::vector<Conservative> R_star(R_in.size());
        R_star = R_in; // initialize all components

        int imin, imax, jmin, jmax;
        mesh.getInteriorBounds(imin, imax, jmin, jmax);
        const int Ni = imax - imin;
        const int Nj = jmax - jmin;

        // ------------------------------
        // 3.1 I-direction line smoothing
        // ------------------------------
        if (epsI > 0.0) {
            Eigen::VectorXd a(Ni), b(Ni), c(Ni), d(Ni);

            // For each scalar component
            for (int comp = 0; comp < 4; ++comp) {

                // For each j-line (row)
                for (int j = jmin; j < jmax; ++j) {

                    // Build RHS d from R_in along this line (fixed j)
                    for (int li = 0; li < Ni; ++li) {
                        int i   = imin + li;
                        int cid = mesh.cellIndex(i, j);
                        d(li)   = getComponent(R_in[cid], comp);
                    }

                    // Build tridiagonal coefficients with Neumann BC at ends
                    if (Ni == 1) {
                        // Degenerate case: just copy RHS
                        a(0) = 0.0; b(0) = 1.0; c(0) = 0.0;
                    } else {
                        // Left boundary
                        a(0) = 0.0;
                        b(0) = 1.0 + epsI;
                        c(0) = -epsI;

                        // Interior
                        for (int li = 1; li < Ni - 1; ++li) {
                            a(li) = -epsI;
                            b(li) = 1.0 + 2.0 * epsI;
                            c(li) = -epsI;
                        }

                        // Right boundary
                        a(Ni - 1) = -epsI;
                        b(Ni - 1) = 1.0 + epsI;
                        c(Ni - 1) = 0.0;
                    }

                    // Solve line system
                    Eigen::VectorXd u = thomasAlgorithm(a, b, c, d);

                    // Write back into R_star
                    for (int li = 0; li < Ni; ++li) {
                        int i   = imin + li;
                        int cid = mesh.cellIndex(i, j);
                        setComponent(R_star[cid], comp, u(li));
                    }
                } // end loop j
            }     // end loop comp
        } // end I-smoothing

        // ------------------------------
        // 3.2 J-direction line smoothing
        // ------------------------------
        if (epsJ > 0.0) {
            Eigen::VectorXd a(Nj), b(Nj), c(Nj), d(Nj);

            for (int comp = 0; comp < 4; ++comp) {

                // For each i-column (fixed i)
                for (int i = imin; i < imax; ++i) {

                    // Build RHS from R_star along this column
                    for (int lj = 0; lj < Nj; ++lj) {
                        int j   = jmin + lj;
                        int cid = mesh.cellIndex(i, j);
                        d(lj)   = getComponent(R_star[cid], comp);
                    }

                    // Build tridiagonal coefficients with Neumann BC at ends
                    if (Nj == 1) {
                        a(0) = 0.0; b(0) = 1.0; c(0) = 0.0;
                    } else {
                        // Bottom boundary (wall)
                        a(0) = 0.0;
                        b(0) = 1.0 + epsJ;
                        c(0) = -epsJ;

                        // Interior
                        for (int lj = 1; lj < Nj - 1; ++lj) {
                            a(lj) = -epsJ;
                            b(lj) = 1.0 + 2.0 * epsJ;
                            c(lj) = -epsJ;
                        }

                        // Top boundary (farfield)
                        a(Nj - 1) = -epsJ;
                        b(Nj - 1) = 1.0 + epsJ;
                        c(Nj - 1) = 0.0;
                    }

                    // Solve line system
                    Eigen::VectorXd u = thomasAlgorithm(a, b, c, d);

                    // Write back into R_out
                    for (int lj = 0; lj < Nj; ++lj) {
                        int j   = jmin + lj;
                        int cid = mesh.cellIndex(i, j);
                        setComponent(R_out[cid], comp, u(lj));
                    }
                } // i
            }     // comp
        } else {
            // No J-smoothing: just copy R_star into R_out
            R_out = R_star;
        }
    }

private:
    double epsI;  // smoothing strength in I-direction
    double epsJ;  // smoothing strength in J-direction
};
