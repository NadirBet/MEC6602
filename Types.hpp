#pragma once
#include <cmath>
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================
// Conservative state: U = [rho, rho*u, rho*v, rho*E]
// ============================================================
struct Conservative {
    double rho;   // density
    double rhou;  // x-momentum
    double rhov;  // y-momentum
    double rhoE;  // total energy

    Conservative() : rho(0.0), rhou(0.0), rhov(0.0), rhoE(0.0) {}
    
    Conservative(double rho_, double rhou_, double rhov_, double rhoE_)
        : rho(rho_), rhou(rhou_), rhov(rhov_), rhoE(rhoE_) {}
    
    // Init from primitives (rho, u, v, p) + gamma
    Conservative(double rho_, double u_, double v_, double p_, double gamma)
        : rho(rho_), rhou(rho_ * u_), rhov(rho_ * v_) {
        double V2 = u_*u_ + v_*v_;
        double e_int = p_ / ((gamma - 1.0) * rho_);
        rhoE = rho_ * (e_int + 0.5 * V2);
    }

    // Velocities (safe division)
    double u() const { return rhou / (rho + 1e-14); }
    double v() const { return rhov / (rho + 1e-14); }

    // Operators for RK updates
    Conservative& operator+=(const Conservative& rhs) {
        rho += rhs.rho; rhou += rhs.rhou; rhov += rhs.rhov; rhoE += rhs.rhoE;
        return *this;
    }
    
    Conservative& operator-=(const Conservative& rhs) {
        rho -= rhs.rho; rhou -= rhs.rhou; rhov -= rhs.rhov; rhoE -= rhs.rhoE;
        return *this;
    }
};

// Free operators for RK stages
inline Conservative operator+(const Conservative& a, const Conservative& b) {
    return Conservative(a.rho + b.rho, a.rhou + b.rhou, 
                       a.rhov + b.rhov, a.rhoE + b.rhoE);
}

inline Conservative operator-(const Conservative& a, const Conservative& b) {
    return Conservative(a.rho - b.rho, a.rhou - b.rhou, 
                       a.rhov - b.rhov, a.rhoE - b.rhoE);
}

inline Conservative operator*(const Conservative& a, double s) {
    return Conservative(a.rho * s, a.rhou * s, a.rhov * s, a.rhoE * s);
}

inline Conservative operator*(double s, const Conservative& a) {
    return a * s;
}

// ============================================================
// Primitive state: W = [rho, u, v, p, a]
// ============================================================
struct Primitive {
    double rho, u, v, p, a;  // a = sound speed

    Primitive() : rho(0.0), u(0.0), v(0.0), p(0.0), a(0.0) {}
    
    Primitive(double rho_, double u_, double v_, double p_, double gamma)
        : rho(rho_), u(u_), v(v_), p(p_) {
        a = std::sqrt(gamma * p / rho);
    }
    
    Primitive(const Conservative& U, double gamma) {
        rho = std::max(U.rho, 1e-14);
        u = U.rhou / rho;
        v = U.rhov / rho;
        double ke = 0.5 * (u*u + v*v);
        p = std::max((gamma - 1.0) * (U.rhoE - rho * ke), 1e-14);
        a = std::sqrt(gamma * p / rho);
    }

    Conservative toConservative(double gamma) const {
        double V2 = u*u + v*v;
        double rhoE_ = p/(gamma - 1.0) + 0.5 * rho * V2;
        return Conservative(rho, rho*u, rho*v, rhoE_);
    }
    
    double Mach() const {
        double V = std::sqrt(u*u + v*v);
        return V / (a + 1e-14);
    }
};

// ============================================================
// EOS utilities (perfect gas)
// ============================================================
namespace EOS {
    inline double pressure(const Conservative& U, double gamma) {
        double rho = std::max(U.rho, 1e-14);
        double u = U.rhou / rho;
        double v = U.rhov / rho;
        double ke = 0.5 * rho * (u*u + v*v);
        return std::max((gamma - 1.0) * (U.rhoE - ke), 1e-14);
    }

    inline double soundSpeed(const Conservative& U, double gamma) {
        double p = pressure(U, gamma);
        double rho = std::max(U.rho, 1e-14);
        return std::sqrt(gamma * p / rho);
    }

    inline double totalEnthalpy(const Conservative& U, double gamma) {
        double p = pressure(U, gamma);
        double rho = std::max(U.rho, 1e-14);
        return (U.rhoE + p) / rho;
    }

    inline void makePhysical(Conservative& U, double gamma) {
        U.rho = std::max(U.rho, 1e-14);
        double p = pressure(U, gamma);
        if (p < 1e-14) {
            double rho = U.rho;
            double u = U.rhou / rho;
            double v = U.rhov / rho;
            double ke = 0.5 * rho * (u*u + v*v);
            U.rhoE = (1e-14 / (gamma - 1.0)) + ke;
        }
    }
}

// ============================================================
// Configuration
// ============================================================
struct Config {
    // Flow conditions
    double Mach_inf = 0.5;
    double alpha_deg = 0.0;
    double alpha_rad;
    double gamma = 1.4;
    
    // Reference values
    double rho_inf = 1.0;
    double p_inf = 1.0;
    double a_inf;
    double u_inf, v_inf;
    double chord = 1.0;
    
    // Numerical parameters
    double CFL    = 0.1;
    double k2_jst = 0.5;     // Second-order dissipation
    double k4_jst = 0.02;    // Fourth-order dissipation
    int    maxIter    = 1000;
    int    printFreq  = 100;
    int    outputFreq = 500;
    int    dbg_log_jst      = 1;
    int    dbg_log_lambdas  = 1;
    int    dbg_log_split    = 1;

    // Residual smoothing controls
    bool   use_residual_smoothing = false;
    double smooth_eps_I           = 0.0;   // I-direction IRS eps
    double smooth_eps_J           = 0.0;   // J-direction IRS eps

    // -------- Multigrid controls (NEW) --------
    bool use_multigrid = false;   // if true: use FASMultigridSolver path
    int  mg_cycles     = 10;      // number of V-cycles to run
    double mg_res_tol  = 1e-12; 

    
    void initialize() {
        alpha_rad = alpha_deg * M_PI / 180.0;
        a_inf = std::sqrt(gamma * p_inf / rho_inf);
        u_inf = Mach_inf * a_inf * std::cos(alpha_rad);
        v_inf = Mach_inf * a_inf * std::sin(alpha_rad);
    }
    
    Primitive getFreestream() const {
        return Primitive(rho_inf, u_inf, v_inf, p_inf, gamma);
    }
};
