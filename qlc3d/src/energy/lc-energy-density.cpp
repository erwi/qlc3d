#include <energy/lc-energy-density.h>
#include <geom/vec3.h>
#include <cmath>

// Physical constants
static constexpr double e0 = 8.8541878176e-12; // vacuum permittivity [F/m]

// Tensor basis constants
static const double rt2 = std::sqrt(2.0);
static const double rt3 = std::sqrt(3.0);
static const double rt6 = std::sqrt(6.0);

double elasticEnergyDensity(const GaussPointData &p, const EnergyMaterialParams &m) {
    const double S0 = m.S0;
    const double q0 = m.q0;

    const double q1 = p.q1, q2 = p.q2, q3 = p.q3, q4 = p.q4, q5 = p.q5;
    const double q1x = p.q1x, q1y = p.q1y, q1z = p.q1z;
    const double q2x = p.q2x, q2y = p.q2y, q2z = p.q2z;
    const double q3x = p.q3x, q3y = p.q3y, q3z = p.q3z;
    const double q4x = p.q4x, q4y = p.q4y, q4z = p.q4z;
    const double q5x = p.q5x, q5y = p.q5y, q5z = p.q5z;

    // G4: twist invariant: (9S²/4) * n·curl(n) in TTensor form
    double G4 = (q2*q4x - q4*q2x - q3*q5x + q5*q3x + q2*q5y + q3*q4y - q4*q3y - q5*q2y
                 - 2*q2*q3z + 2*q3*q2z + q4*q5z - q5*q4z) * 0.5
              + (3*q1*q4x - 3*q4*q1x - 3*q1*q5y + 3*q5*q1y) / (rt2 * rt6);

    double R = q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5;
    // aa = (3/2)*R is the chirality scale factor ≈ (9S²/4) at ground state
    double aa = 1.5 * R;
    double F_twist = (G4 - aa * q0) / (9.0 * S0 * S0) * 4.0;
    F_twist *= F_twist; // squared: (twist - chirality_offset)²

    // G1: sum of all squared TTensor spatial derivatives (isotropic bulk term)
    double G1 = (6.0 / (rt6*rt6)) * (q1x*q1x + q1y*q1y + q1z*q1z)
              + (2.0 / (rt2*rt2)) * (q2x*q2x + q3x*q3x + q4x*q4x + q5x*q5x
                                   + q2y*q2y + q3y*q3y + q4y*q4y + q5y*q5y
                                   + q2z*q2z + q3z*q3z + q4z*q4z + q5z*q5z);

    // G2: splay-related invariant in TTensor form
    double G2 = (1.0/(rt6*rt6)) * (q1x*q1x + q1y*q1y + 4.0*q1z*q1z)
              + (1.0/(rt2*rt2)) * (2*q2x*q3y - 2*q3x*q2y + 2*q5x*q4y + 2*q2x*q5z + 2*q3x*q4z
                                  - 2*q2y*q4z + 2*q3y*q5z
                                  + q2x*q2x + q3x*q3x + q5x*q5x + q2y*q2y + q3y*q3y
                                  + q4y*q4y + q4z*q4z + q5z*q5z)
              - (2*q1x*q2x + 2*q1x*q3y + 2*q3x*q1y - 2*q1y*q2y + 2*q1x*q5z
                 - 4*q5x*q1z + 2*q1y*q4z - 4*q4y*q1z) / (rt2*rt6);

    // G6: higher-order elastic invariant in TTensor form (contributes to splay + bend split)
    double G6 = (1.0/(rt2*rt2*rt2)) * (
          2*q2*q2x*q2x + 2*q2*q3x*q3x + 2*q2*q4x*q4x + 2*q2*q5x*q5x
        - 2*q2*q2y*q2y - 2*q2*q3y*q3y - 2*q2*q4y*q4y - 2*q2*q5y*q5y
        + 4*q3*q2x*q2y + 4*q3*q3x*q3y + 4*q3*q4x*q4y + 4*q3*q5x*q5y
        + 4*q5*q2x*q2z + 4*q5*q3x*q3z + 4*q5*q4x*q4z + 4*q5*q5x*q5z
        + 4*q4*q2y*q2z + 4*q4*q3y*q3z + 4*q4*q4y*q4z + 4*q4*q5y*q5z)
      + (6.0/(rt6*rt6)) * (q2*(q1x*q1x) - q2*(q1y*q1y) + 2*q3*q1x*q1y + 2*q5*q1x*q1z + 2*q4*q1y*q1z) / rt2
      - q1 * (1.0/(rt6*rt6*rt6)) * 6.0 * (q1x*q1x + q1y*q1y - 2.0*q1z*q1z)
      - (q1 * (1.0/(rt2*rt2)) * 2.0 * (q2x*q2x + q3x*q3x + q4x*q4x + q5x*q5x
                                       + q2y*q2y + q3y*q3y + q4y*q4y + q5y*q5y
                                       - 2*q2z*q2z - 2*q3z*q3z - 2*q4z*q4z - 2*q5z*q5z)) / rt6;

    // Decompose into splay (K11), twist (K22), and bend (K33) contributions
    double F_splay = 4.0 * G2  / (9.0  * S0*S0)
                   - 2.0 * G1  / (27.0 * S0*S0)
                   - 4.0 * G6  / (27.0 * S0*S0*S0);
    double F_bend  = 2.0 * G1  / (27.0 * S0*S0)
                   + 4.0 * G6  / (27.0 * S0*S0*S0);

    return 0.5 * m.K11 * F_splay
         + 0.5 * m.K22 * F_twist
         + 0.5 * m.K33 * F_bend;
}

double thermotropicEnergyDensity(const GaussPointData &p, const EnergyMaterialParams &m) {
    const double A = m.A, B = m.B, C = m.C;

    const double q1 = p.q1, q2 = p.q2, q3 = p.q3, q4 = p.q4, q5 = p.q5;
    double R = q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5;

    // Bulk thermotropic free energy in TTensor form (raw, no ground-state offset).
    // The cubic (B) term represents (B/3)*tr(Q³) expressed via TTensor DOFs.
    // At the equilibrium state this is negative (not zero).
    double Fth = A * R / 2.0
               + B * (  q5*q5*q1*rt6/4.0
                      - q1*rt6*q2*q2/2.0
                      - q3*q3*q1*rt6/2.0
                      + 3.0/4.0*q5*q5*q2*rt2
                      + 3.0/2.0*q3*rt2*q5*q4
                      + q4*q4*q1*rt6/4.0
                      - 3.0/4.0*q4*q4*q2*rt2
                      + q1*q1*q1*rt6/6.0) / 3.0
               + C * (R*R) / 4.0;

    return Fth;
}

double electricEnergyDensity(const GaussPointData &p, const EnergyMaterialParams &m) {
    const double S0 = m.S0;

    // Effective permittivity coefficients in TTensor convention
    double epsav  = m.eps_per / S0;
    double deleps = (m.eps_par - m.eps_per) / S0;
    // Flexoelectric coupling constants in TTensor convention
    double efe  = 2.0 / (9.0 * S0)      * (m.e1 + 2.0 * m.e3);
    double efe2 = 4.0 / (9.0 * S0 * S0) * (m.e1 - m.e3);

    // Electric field is stored as E_i = -∂φ/∂x_i.
    // The original energy.cpp used Vx=∂φ/∂x (potential gradient); convert here.
    double Vx = -p.Ex, Vy = -p.Ey, Vz = -p.Ez;

    const double q1 = p.q1, q2 = p.q2, q3 = p.q3, q4 = p.q4, q5 = p.q5;

    // Dielectric contribution: -½ε₀ε·E² in TTensor form
    double Fdiel = e0 * (-Vx*Vx - Vy*Vy - Vz*Vz) * epsav * 0.5
                 + e0 * deleps * (  Vx*Vx * q1*rt6/12.0
                                  - Vx*Vx * q2*rt2/4.0
                                  - Vx*Vy * q3*rt2/2.0
                                  - Vx*Vz * q5*rt2/2.0
                                  - Vy*Vz * q4*rt2/2.0
                                  - Vz*Vz * q1*rt6/6.0
                                  + Vy*Vy * q2/4.0);

    double Fflx = 0.0;
    if (efe != 0.0) {
        // First flexoelectric term: linear in (e1+2e3), linear in V·∇q
        const double q1x = p.q1x, q1y = p.q1y, q1z = p.q1z;
        const double q2x = p.q2x, q2y = p.q2y;
        const double q3x = p.q3x, q3y = p.q3y;
        const double q4y = p.q4y, q4z = p.q4z;
        const double q5x = p.q5x, q5z = p.q5z;
        Fflx += efe * ( (2*Vz*q1z - Vx*q1x - Vy*q1y) / rt6
                      + (Vx*(q2x+q3y+q5z) + Vy*(-q2y+q3x+q4z) + Vz*(q4y+q5x)) / rt2);
    }
    if (efe2 != 0.0) {
        // Second flexoelectric term: linear in (e1-e3), bilinear in q·V·∇q
        const double q1x = p.q1x, q1y = p.q1y, q1z = p.q1z;
        const double q2x = p.q2x, q2y = p.q2y;
        const double q3x = p.q3x, q3y = p.q3y;
        const double q4y = p.q4y, q4z = p.q4z;
        const double q5x = p.q5x, q5y = p.q5y, q5z = p.q5z;
        Fflx += efe2 * (
              ( Vx*(-q1*(q2x+q3y+q5z) - q2*q1x - q3*q1y + 2*q5*q1z)
              + Vy*(q1*(q2y-q3x-q4z) + q2*q1y - q3*q1x + 2*q4*q1z)
              + Vz*(2*q1*(q4y+q5x) - q4*q1y - q5*q1x) ) * rt3/6.0
              + Vx*(q2*(q2x+q3y+q5z) + q3*(-q2y+q3x+q4z) + q5*(q4y+q5x))/2.0
              + (Vy*q2-Vz*q4)*(q2y-q3x-q4z)/2.0
              + (Vy*q3+Vz*q5)*(q2x+q3y+q5z)/2.0
              + Vy*q4*(q4y+q5x)/2.0
              + q1*(Vx*q1x + 4*Vz*q1z + Vy*q1y)/6.0);
    }
    return Fdiel + Fflx;
}

EnergyMaterialParams lcToEnergyMaterialParams(const LC &lc) {
    constexpr double pi = 3.14159265358979323846;
    double q0 = 0.0;
    if (lc.p0() > 0.0) {
        q0 = 2.0 * pi / lc.p0();
    }
    return EnergyMaterialParams{
        .S0      = lc.S0(),
        .K11     = lc.K11(),
        .K22     = lc.K22(),
        .K33     = lc.K33(),
        .A       = lc.A(),
        .B       = lc.B(),
        .C       = lc.C(),
        .eps_par = lc.eps_par(),
        .eps_per = lc.eps_per(),
        .e1      = lc.e1(),
        .e3      = lc.e3(),
        .q0      = q0
    };
}

double surfaceEnergyDensity(double q1, double q2, double q3, double q4, double q5,
                            const Vec3 &v1, const Vec3 &v2,
                            double W, double K1, double K2, double S0) {
    const double A = (K1 + K2) / (S0 * 6.0);

    // Isotropic orientation-independent term
    double fIso = W * A * (q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5);

    // v̂·Q·v̂ in qlc3d TTensor basis form:
    //   q1/rt6 * (-vx²-vy²+2vz²) + q2/rt2 * (vx²-vy²) + q3*rt2*vx*vy + q4*rt2*vy*vz + q5*rt2*vx*vz
    // v̂₁·Q·v̂₁
    double vQv1 = q1 / rt6 * (-v1.x()*v1.x() - v1.y()*v1.y() + 2*v1.z()*v1.z())
                + q2 / rt2 * (v1.x()*v1.x() - v1.y()*v1.y())
                + q3 * rt2 * v1.x() * v1.y()
                + q4 * rt2 * v1.y() * v1.z()
                + q5 * rt2 * v1.x() * v1.z();

    // v̂₂·Q·v̂₂
    double vQv2 = q1 / rt6 * (-v2.x()*v2.x() - v2.y()*v2.y() + 2*v2.z()*v2.z())
                + q2 / rt2 * (v2.x()*v2.x() - v2.y()*v2.y())
                + q3 * rt2 * v2.x() * v2.y()
                + q4 * rt2 * v2.y() * v2.z()
                + q5 * rt2 * v2.x() * v2.z();

    return fIso + W * K1 * vQv1 + W * K2 * vQv2;
}


