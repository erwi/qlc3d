#include <catch.h>
#include <energy/lc-energy-density.h>
#include <lc.h>
#include <lc-representation.h>
#include <geom/vec3.h>
#include <cmath>
#include <memory>

// Ground state director along z: n=(0,0,1), S=S0
// TTensor from QTensor conversion:
//   Q11 = S/2*(3*0-1) = -S/2,  Q22 = -S/2,  Q33 = S,  Q12=Q23=Q13=0
//   t1 = (Q11+Q22)*(-sqrt(6)/2) = -S*(-sqrt(6)/2) = S*sqrt(6)/2
//   t2 = (Q11+(Q11+Q22)/(-2))*sqrt(2) = (-S/2 + S/2)*sqrt(2) = 0
//   t3 = Q12*sqrt(2) = 0, t4=t5=0
static double groundStateQ1(double S0) { return S0 * std::sqrt(6.0) / 2.0; }

// Default test LC material (typical nematic parameters)
// Returns the default LCBuilder parameters as an LC object.
static std::unique_ptr<LC> makeTestLC() {
    return std::unique_ptr<LC>(LCBuilder{}
        .K11(10e-12)
        .K22(10e-12)
        .K33(10e-12)
        .A(-1.2e5)
        .B(-2.1333e6)
        .C(1.7333e6)
        .eps_par(18.5)
        .eps_per(7.0)
        .e1(0)
        .e3(0)
        .build());
}

TEST_CASE("elasticEnergyDensity returns zero for uniform ground state (no gradients)", "[energy-density]") {
    // ARRANGE
    auto lc = makeTestLC();
    auto mat = lcToEnergyMaterialParams(*lc);
    // Ground state along z: q1 = S0*sqrt(6)/2, all others and all derivatives zero
    GaussPointData p;
    p.q1 = groundStateQ1(mat.S0);

    // ACT
    double result = elasticEnergyDensity(p, mat);

    // ASSERT: no spatial gradients → zero elastic energy
    REQUIRE(result == Approx(0.0).margin(1e-30));
}

TEST_CASE("thermotropicEnergyDensity: ground state is a local minimum w.r.t. order parameter perturbations", "[energy-density]") {
    // GIVEN: typical nematic LC material
    // WHEN:  director is at the ground state (uniform uniaxial S0, along z)
    // THEN:  any perturbation away from ground state raises the energy:
    //        - reducing the scalar order parameter (less ordered)
    //        - increasing the scalar order parameter (over-ordered)
    //        - adding a biaxial component q2 (off-axis deviation at constant order)
    auto lc = makeTestLC();
    auto mat = lcToEnergyMaterialParams(*lc);

    GaussPointData p_gs;
    p_gs.q1 = groundStateQ1(mat.S0);
    double f_gs = thermotropicEnergyDensity(p_gs, mat);

    // Perturbation 1: reduce order parameter by 10%
    GaussPointData p_reduced;
    p_reduced.q1 = groundStateQ1(mat.S0) * 0.9;
    REQUIRE(thermotropicEnergyDensity(p_reduced, mat) > f_gs);

    // Perturbation 2: increase order parameter by 10%
    GaussPointData p_increased;
    p_increased.q1 = groundStateQ1(mat.S0) * 1.1;
    REQUIRE(thermotropicEnergyDensity(p_increased, mat) > f_gs);

    // Perturbation 3: add a biaxial component q2 (keeps |q| roughly same but breaks uniaxial symmetry)
    GaussPointData p_biaxial;
    p_biaxial.q1 = groundStateQ1(mat.S0);
    p_biaxial.q2 = 0.1 * mat.S0;
    REQUIRE(thermotropicEnergyDensity(p_biaxial, mat) > f_gs);

    // Sanity: ground state energy is negative for these typical nematic parameters
    REQUIRE(f_gs < 0.0);
}

TEST_CASE("elasticEnergyDensity is consistent with single non-zero gradient q1z", "[energy-density]") {
    // ARRANGE: only q1z ≠ 0, all other values and derivatives zero
    // Check against the manually simplified contribution from the L1-like term.
    //
    // With only q1z ≠ 0:
    //   G1 = (6/(rt6²))*(q1z²) = q1z²
    //   G2 = (1/(rt6²))*(4*q1z²) - (−4*q4y*q1z)/(rt2*rt6) = (2/3)*q1z²  [q4y=0]
    //        Wait: G2 numerically = (1/6)*4*q1z² = (2/3)*q1z²
    //   G6 = 0  (all q=0)
    //   G4 = 0  (all values and x/y derivatives zero)
    //   F_twist = 0
    //   F_splay = 4*G2/(9*S0²) - 2*G1/(27*S0²) = (8/3)/(9*S0²)*q1z² - 2/(27*S0²)*q1z²
    //           = (8/3 - 2/3 (wait) = let me compute numerically below
    EnergyMaterialParams mat;
    mat.S0  = 0.6;
    mat.K11 = 10e-12;
    mat.K22 = 10e-12;
    mat.K33 = 10e-12;
    mat.A = -1.2e5;
    mat.B = -2.1333e6;
    mat.C = 1.7333e6;
    mat.eps_par = 18.5;
    mat.eps_per = 7.0;
    mat.e1 = 0; mat.e3 = 0; mat.q0 = 0;

    const double dq1z = 1e3; // arbitrary non-zero gradient (units: S0/m in SI)
    const double S0 = mat.S0;

    GaussPointData p;
    p.q1z = dq1z;
    // all other values zero

    // ACT
    double result = elasticEnergyDensity(p, mat);

    // ASSERT: manually compute expected value
    // With only q1z != 0, all qi = 0:
    //   G4 = 0 (no twist), aa=0 → F_twist = 0 (q=0 so aa=0)
    //   G1 = 6/(rt6²)*q1z² = 6/6*q1z² = q1z²
    //   G2 = 1/(rt6²)*(4*q1z²) = 4/6*q1z² = (2/3)*q1z²
    //   G6 = 0  (all qi=0)
    //   F_splay = 4*G2/(9*S0²) - 2*G1/(27*S0²)
    //           = 4*(2/3)*q1z²/(9*S0²) - 2*q1z²/(27*S0²)
    //           = (8/3)/(9*S0²)*q1z² - 2/(27*S0²)*q1z²
    //           = (8/(27) - 2/(27)) * q1z²/S0²
    //           = (6/27)*q1z²/S0²
    //           = (2/9)*q1z²/S0²
    //   F_bend  = 2*G1/(27*S0²) = 2*q1z²/(27*S0²)
    //   result = 0.5*K11*(2/9)*q1z²/S0² + 0.5*K33*2/(27)*q1z²/S0²
    //          = q1z²/S0² * (K11/9 + K33/27)
    double expected = (dq1z * dq1z) / (S0 * S0) * (mat.K11 / 9.0 + mat.K33 / 27.0);
    REQUIRE(result == Approx(expected).epsilon(1e-10));
}

TEST_CASE("electricEnergyDensity with Ez only and ground state along z matches analytic value", "[energy-density]") {
    // ...existing test body...
    // ARRANGE: ground state q along z, E field along z only, no flexoelectric
    const double S0 = 0.6;
    EnergyMaterialParams mat;
    mat.S0      = S0;
    mat.K11     = 10e-12;
    mat.K22     = 10e-12;
    mat.K33     = 10e-12;
    mat.A       = -1.2e5;
    mat.B       = -2.1333e6;
    mat.C       = 1.7333e6;
    mat.eps_par = 18.5;
    mat.eps_per = 7.0;
    mat.e1      = 0;
    mat.e3      = 0;
    mat.q0      = 0;

    const double Ez0 = 1e6; // 1 MV/m field along z
    GaussPointData p;
    p.q1 = groundStateQ1(S0);
    p.Ez = Ez0;
    // Ex=Ey=0, no gradients needed for dielectric term (no flexo with e1=e3=0)

    // ACT
    double result = electricEnergyDensity(p, mat);

    // ASSERT: manually compute expected value
    // Vx=-Ex=0, Vy=-Ey=0, Vz=-Ez=-Ez0
    // epsav = eps_per/S0, deleps = (eps_par-eps_per)/S0
    // Fdiel = e0*(-Vz²)*epsav*0.5 + e0*deleps*(-Vz²*q1*rt6/6)
    //       = -e0*Ez0²*(eps_per/S0)/2 + e0*(eps_par-eps_per)/S0*(-Ez0²*(S0*rt6/2)*rt6/6)
    //       = -e0*Ez0²*eps_per/(2*S0) - e0*(eps_par-eps_per)/S0 * Ez0² * S0/2
    //       = -e0*Ez0²*eps_per/(2*S0) - e0*(eps_par-eps_per)*Ez0²/2
    constexpr double e0_val = 8.8541878176e-12;
    const double rt6_val = std::sqrt(6.0);
    double expected = -e0_val * Ez0 * Ez0 * mat.eps_per / (2.0 * S0)
                    - e0_val * (mat.eps_par - mat.eps_per) / S0
                      * (Ez0 * Ez0 * groundStateQ1(S0) * rt6_val / 6.0);
    REQUIRE(result == Approx(expected).epsilon(1e-10));
}

// ============================================================================
// surfaceEnergyDensity tests
//
// Geometry: v1=(0,1,0), v2=(0,0,1), easy axis e = v1×v2 = (1,0,0).
//
// The Rapini-Papoular energy in the qlc3d TTensor basis at uniaxial order S0:
//   f_S(θ) = f_min + (3/2)·W·K1·S0·sin²(θ₁) + (3/2)·W·K2·S0·sin²(θ₂)
// where θ₁ = angle of tilt towards v1, θ₂ = angle of tilt towards v2 and
// f_min = -W·S0·(K1+K2)/4.
//
// This follows from v̂·Q·v̂ = S0/2·(3·(v̂·n)² − 1) for uniaxial Q = S0/2·(3n⊗n − I):
//   rotating towards v1 by θ: v1·n = sin(θ), v2·n = 0
//     f(θ) - f(0) = W·K1·[v1·Q·v1(θ) - v1·Q·v1(0)]
//                 = W·K1·S0/2·[3·sin²(θ) − 0] = (3/2)·W·K1·S0·sin²(θ)
//   rotating towards v2 by θ: v1·n = 0, v2·n = sin(θ)  (analogous, gives K2)
// ============================================================================


TEST_CASE("surfaceEnergyDensity: energy varies as Rapini-Papoular sin² when rotating director towards v1", "[energy-density]") {
    // GIVEN: planar WeakAnchoring with K1=1, K2=0  (only v1 term active)
    //        W=1, S0=0.6, v1=(0,1,0), v2=(0,0,1), easy axis e=(1,0,0)
    const double W  = 1.0;
    const double K1 = 1.0;
    const double K2 = 0.0;
    const double S0 = 0.6;
    const Vec3 v1(0, 1, 0);
    const Vec3 v2(0, 0, 1);

    // WHEN: director rotates from e=(1,0,0) towards v1=(0,1,0) by twist angle θ:
    //       n(θ) = (cos θ, sin θ, 0) = fromRadianAngles(tilt=0, twist=θ).
    //       At each θ, compute f_S(θ) and compare to f_S(0) + (3/2)·W·K1·S0·sin²(θ).
    const qlc3d::TTensor t0 = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(0, 0, S0));
    double f_easy = surfaceEnergyDensity(t0[0], t0[1], t0[2], t0[3], t0[4], v1, v2, W, K1, K2, S0);

    // THEN: sin² shape holds at 10 angles from 0 to π/2
    const int N = 10;
    for (int i = 1; i <= N; i++) {
        double theta = i * M_PI / (2.0 * N);
        // tilt=0, twist=θ gives n=(cos θ, sin θ, 0): rotation in the x-y plane towards v1
        const qlc3d::TTensor t = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(0, theta, S0));
        double f = surfaceEnergyDensity(t[0], t[1], t[2], t[3], t[4], v1, v2, W, K1, K2, S0);

        double expectedVariation = 1.5 * W * K1 * S0 * std::sin(theta) * std::sin(theta);
        REQUIRE((f - f_easy) == Approx(expectedVariation).epsilon(1e-10));
    }
    // THEN: the easy axis is the minimum (all variations positive)
    REQUIRE(f_easy < surfaceEnergyDensity(t0[0], t0[1], t0[2], t0[3], t0[4], v1, v2, W, K1, K2, S0) + 1e-15);
}

TEST_CASE("surfaceEnergyDensity: energy varies as Rapini-Papoular sin² when rotating director towards v2", "[energy-density]") {
    // GIVEN: planar WeakAnchoring with K1=0, K2=1  (only v2 term active)
    //        W=1, S0=0.6, v1=(0,1,0), v2=(0,0,1), easy axis e=(1,0,0)
    const double W  = 1.0;
    const double K1 = 0.0;
    const double K2 = 1.0;
    const double S0 = 0.6;
    const Vec3 v1(0, 1, 0);
    const Vec3 v2(0, 0, 1);

    // WHEN: director rotates from e=(1,0,0) towards v2=(0,0,1) by tilt angle θ:
    //       n(θ) = (cos θ, 0, sin θ) = fromRadianAngles(tilt=θ, twist=0).
    //       At each θ, energy variation must follow (3/2)·W·K2·S0·sin²(θ).
    const qlc3d::TTensor t0 = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(0, 0, S0));
    double f_easy = surfaceEnergyDensity(t0[0], t0[1], t0[2], t0[3], t0[4], v1, v2, W, K1, K2, S0);

    // THEN: sin² shape holds at 10 angles from 0 to π/2
    const int N = 10;
    for (int i = 1; i <= N; i++) {
        double theta = i * M_PI / (2.0 * N);
        // tilt=θ, twist=0 gives n=(cos θ, 0, sin θ): rotation in the x-z plane towards v2
        const qlc3d::TTensor t = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(theta, 0, S0));
        double f = surfaceEnergyDensity(t[0], t[1], t[2], t[3], t[4], v1, v2, W, K1, K2, S0);

        double expectedVariation = 1.5 * W * K2 * S0 * std::sin(theta) * std::sin(theta);
        REQUIRE((f - f_easy) == Approx(expectedVariation).epsilon(1e-10));
    }
}

TEST_CASE("surfaceEnergyDensity: K1 and K2 independently control energy variation along v1 and v2", "[energy-density]") {
    // GIVEN: planar WeakAnchoring, W=1, S0=0.6, v1=(0,1,0), v2=(0,0,1), θ=45°
    //        easy axis e=(1,0,0)
    // WHEN:  K1=1, K2=0: rotation towards v1 costs energy; rotation towards v2 does not.
    // WHEN:  K1=0, K2=1: rotation towards v2 costs energy; rotation towards v1 does not.
    // WHEN:  K1=2, K2=1: rotation towards v1 costs twice as much as towards v2.
    const double W  = 1.0;
    const double S0 = 0.6;
    const Vec3 v1(0, 1, 0);
    const Vec3 v2(0, 0, 1);
    const double theta = M_PI / 4.0;  // 45°

    const qlc3d::TTensor t_e  = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(0, 0, S0));
    // Director tilted 45° towards v1: twist=45°, tilt=0 → n=(cos θ, sin θ, 0)
    const qlc3d::TTensor t_v1 = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(0, theta, S0));
    // Director tilted 45° towards v2: tilt=45°, twist=0 → n=(cos θ, 0, sin θ)
    const qlc3d::TTensor t_v2 = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(theta, 0, S0));

    {
        // K1=1, K2=0: only v1-tilt has energy cost
        const double K1 = 1.0, K2 = 0.0;
        double f_e  = surfaceEnergyDensity(t_e[0],  t_e[1],  t_e[2],  t_e[3],  t_e[4],  v1, v2, W, K1, K2, S0);
        double f_v1 = surfaceEnergyDensity(t_v1[0], t_v1[1], t_v1[2], t_v1[3], t_v1[4], v1, v2, W, K1, K2, S0);
        double f_v2 = surfaceEnergyDensity(t_v2[0], t_v2[1], t_v2[2], t_v2[3], t_v2[4], v1, v2, W, K1, K2, S0);

        REQUIRE(f_v1 > f_e);                                                          // v1-tilt raises energy
        REQUIRE(f_v2 == Approx(f_e).epsilon(1e-10));                                  // v2-tilt has zero cost (K2=0)
        REQUIRE((f_v1 - f_e) == Approx(1.5 * W * K1 * S0 * 0.5).epsilon(1e-10));    // sin²(45°) = 0.5
    }
    {
        // K1=0, K2=1: only v2-tilt has energy cost
        const double K1 = 0.0, K2 = 1.0;
        double f_e  = surfaceEnergyDensity(t_e[0],  t_e[1],  t_e[2],  t_e[3],  t_e[4],  v1, v2, W, K1, K2, S0);
        double f_v1 = surfaceEnergyDensity(t_v1[0], t_v1[1], t_v1[2], t_v1[3], t_v1[4], v1, v2, W, K1, K2, S0);
        double f_v2 = surfaceEnergyDensity(t_v2[0], t_v2[1], t_v2[2], t_v2[3], t_v2[4], v1, v2, W, K1, K2, S0);

        REQUIRE(f_v1 == Approx(f_e).epsilon(1e-10));                                  // v1-tilt has zero cost (K1=0)
        REQUIRE(f_v2 > f_e);                                                          // v2-tilt raises energy
        REQUIRE((f_v2 - f_e) == Approx(1.5 * W * K2 * S0 * 0.5).epsilon(1e-10));    // sin²(45°) = 0.5
    }
    {
        // K1=2, K2=1: v1-tilt costs twice as much as v2-tilt
        const double K1 = 2.0, K2 = 1.0;
        double f_e  = surfaceEnergyDensity(t_e[0],  t_e[1],  t_e[2],  t_e[3],  t_e[4],  v1, v2, W, K1, K2, S0);
        double f_v1 = surfaceEnergyDensity(t_v1[0], t_v1[1], t_v1[2], t_v1[3], t_v1[4], v1, v2, W, K1, K2, S0);
        double f_v2 = surfaceEnergyDensity(t_v2[0], t_v2[1], t_v2[2], t_v2[3], t_v2[4], v1, v2, W, K1, K2, S0);

        REQUIRE((f_v1 - f_e) == Approx(2.0 * (f_v2 - f_e)).epsilon(1e-10));  // K1/K2 = 2 ratio
    }
}

// ============================================================================
// Homeotropic anchoring tests (WeakHomeotropic convention: K1=0, K2=1, W<0)
//
// The preferred direction is the surface normal v2.
// Energy formula:
//   f_S(θ) = W·S0·[−1/4 + (3/2)·cos²(θ)]
// where θ = angle between n and v2.
// Variation from the minimum at θ=0:
//   f_S(θ) − f_S(0) = −(3/2)·W·K2·S0·sin²(θ)
// Since W<0 this is positive, i.e. energy rises as the director tilts away.
// ============================================================================


TEST_CASE("surfaceEnergyDensity: homeotropic anchoring minimum is at surface normal (n=v2)", "[energy-density]") {
    // GIVEN: WeakHomeotropic conventions: K1=0, K2=1, W<0, v2 = surface normal = (0,0,1)
    const double W  = -1.0;  // negative per WeakHomeotropic convention
    const double K1 = 0.0;
    const double K2 = 1.0;
    const double S0 = 0.6;
    const Vec3 v1(0, 0, 0); // K1=0, so v1 is unused  (zero vector, contributes nothing)
    const Vec3 v2(0, 0, 1); // surface normal

    // WHEN: director is aligned with v2 (homeotropic, preferred orientation)
    // tilt=90°, twist=0 gives n=(0,0,1), i.e. along v2
    const qlc3d::TTensor t_normal = qlc3d::TTensor::fromDirector(qlc3d::Director::fromDegreeAngles(90, 0, S0));
    double f_normal = surfaceEnergyDensity(t_normal[0], t_normal[1], t_normal[2], t_normal[3], t_normal[4],
                                           v1, v2, W, K1, K2, S0);

    // THEN: tilting the director away from v2 always increases the energy
    const int N = 8;
    for (int i = 1; i <= N; i++) {
        double theta = i * M_PI / (2.0 * N);
        // tilt = π/2 − θ, twist=0 gives n=(sin θ, 0, cos θ): rotating away from v2 towards x-axis
        const qlc3d::TTensor t = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(M_PI / 2.0 - theta, 0, S0));
        double f = surfaceEnergyDensity(t[0], t[1], t[2], t[3], t[4], v1, v2, W, K1, K2, S0);
        REQUIRE(f > f_normal); // energy increases as director tilts away from normal
    }
    // Sanity: at n=v2 energy is negative (the most negative value for W<0, K2=1)
    REQUIRE(f_normal < 0.0);
}

TEST_CASE("surfaceEnergyDensity: homeotropic anchoring follows Rapini-Papoular sin² when tilting away from normal", "[energy-density]") {
    // GIVEN: WeakHomeotropic conventions: K1=0, K2=1, W<0, v2=(0,0,1)
    // WHEN:  director rotates from v2=(0,0,1) towards (1,0,0) by angle θ:
    //        n(θ) = (sin θ, 0, cos θ) = fromRadianAngles(tilt=π/2−θ, twist=0)
    //        v2·n = cos(θ),  v2·Q·v2 = S0/2·(3·cos²(θ)−1)
    // THEN:  f(θ) − f(0) = −(3/2)·W·K2·S0·sin²(θ)  (positive, since W<0)
    const double W  = -1.0;
    const double K1 = 0.0;
    const double K2 = 1.0;
    const double S0 = 0.6;
    const Vec3 v1(0, 0, 0);
    const Vec3 v2(0, 0, 1);

    // tilt=90°, twist=0 gives n=(0,0,1), aligned with v2
    const qlc3d::TTensor t0 = qlc3d::TTensor::fromDirector(qlc3d::Director::fromDegreeAngles(90, 0, S0));
    double f_normal = surfaceEnergyDensity(t0[0], t0[1], t0[2], t0[3], t0[4], v1, v2, W, K1, K2, S0);

    const int N = 10;
    for (int i = 1; i <= N; i++) {
        double theta = i * M_PI / (2.0 * N);
        // tilt = π/2 − θ, twist=0 gives n=(sin θ, 0, cos θ)
        const qlc3d::TTensor t = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(M_PI / 2.0 - theta, 0, S0));
        double f = surfaceEnergyDensity(t[0], t[1], t[2], t[3], t[4], v1, v2, W, K1, K2, S0);

        // Expected variation: −(3/2)·W·K2·S0·sin²(θ)  (positive since W<0)
        double expectedVariation = -1.5 * W * K2 * S0 * std::sin(theta) * std::sin(theta);
        REQUIRE((f - f_normal) == Approx(expectedVariation).epsilon(1e-10));
    }
}

// ============================================================================
// Planar degenerate anchoring tests (Degenerate convention: K1=0, K2=1, W>0)
//
// Any director in the plane perpendicular to the surface normal v2 is preferred
// (the energy is degenerate in that plane). The normal direction v2 is the maximum.
//
// Parametrising by γ = angle of tilt away from the degenerate plane (γ=0 → in plane):
//   n(γ) = cos(γ)·p + sin(γ)·v2,  where p is any unit vector ⊥ v2
//   v2·n = sin(γ),  v2·Q·v2 = S0/2·(3·sin²(γ)−1)
//   f_S(γ) − f_S(0) = (3/2)·W·K2·S0·sin²(γ)   (positive, energy rises toward normal)
// ============================================================================


TEST_CASE("surfaceEnergyDensity: planar degenerate anchoring is energy-degenerate for all in-plane directors", "[energy-density]") {
    // GIVEN: Degenerate anchoring: K1=0, K2=1, W>0, v2=(0,0,1) (surface normal to be avoided)
    // WHEN:  director is constrained to the x-y plane (any direction ⊥ v2)
    // THEN:  all such directors give the same energy (true degeneracy)
    const double W  = 1.0;
    const double K1 = 0.0;
    const double K2 = 1.0;
    const double S0 = 0.6;
    const Vec3 v1(0, 0, 0); // K1=0, v1 unused
    const Vec3 v2(0, 0, 1);

    // Reference: n = (1, 0, 0) = fromRadianAngles(tilt=0, twist=0)
    const qlc3d::TTensor t_ref = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(0, 0, S0));
    double f_ref = surfaceEnergyDensity(t_ref[0], t_ref[1], t_ref[2], t_ref[3], t_ref[4],
                                         v1, v2, W, K1, K2, S0);

    // Test 7 other in-plane directions at different azimuthal (twist) angles
    const int N = 7;
    for (int i = 1; i <= N; i++) {
        double phi = i * M_PI / N;  // azimuthal angle in x-y plane
        // tilt=0, twist=φ gives n=(cos φ, sin φ, 0): pure in-plane direction
        const qlc3d::TTensor t = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(0, phi, S0));
        double f = surfaceEnergyDensity(t[0], t[1], t[2], t[3], t[4], v1, v2, W, K1, K2, S0);
        REQUIRE(f == Approx(f_ref).epsilon(1e-10)); // all in-plane directors are degenerate
    }
    // Sanity: in-plane energy is the minimum (negative for W>0, K2=1)
    REQUIRE(f_ref < 0.0);
}

TEST_CASE("surfaceEnergyDensity: planar degenerate anchoring follows Rapini-Papoular sin² when tilting out of plane", "[energy-density]") {
    // GIVEN: Degenerate anchoring: K1=0, K2=1, W>0, v2=(0,0,1)
    // WHEN:  director tilts from in-plane (x-axis) towards v2 by angle γ:
    //        n(γ) = (cos γ, 0, sin γ)
    //        v2·n = sin(γ),  v2·Q·v2 = S0/2·(3·sin²(γ)−1)
    // THEN:  f(γ) − f(0) = (3/2)·W·K2·S0·sin²(γ)   (always positive: energy rises)
    const double W  = 1.0;
    const double K1 = 0.0;
    const double K2 = 1.0;
    const double S0 = 0.6;
    const Vec3 v1(0, 0, 0);
    const Vec3 v2(0, 0, 1);

    // Reference at γ=0: director in-plane along x = fromRadianAngles(tilt=0, twist=0)
    const qlc3d::TTensor t0 = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(0, 0, S0));
    double f_planar = surfaceEnergyDensity(t0[0], t0[1], t0[2], t0[3], t0[4], v1, v2, W, K1, K2, S0);

    const int N = 10;
    for (int i = 1; i <= N; i++) {
        double gamma = i * M_PI / (2.0 * N);
        // tilt=γ, twist=0 gives n=(cos γ, 0, sin γ): tilting out of the x-y plane towards v2
        const qlc3d::TTensor t = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(gamma, 0, S0));
        double f = surfaceEnergyDensity(t[0], t[1], t[2], t[3], t[4], v1, v2, W, K1, K2, S0);

        // Expected variation: (3/2)·W·K2·S0·sin²(γ)
        double expectedVariation = 1.5 * W * K2 * S0 * std::sin(gamma) * std::sin(gamma);
        REQUIRE((f - f_planar) == Approx(expectedVariation).epsilon(1e-10));
    }
}

TEST_CASE("surfaceEnergyDensity: planar degenerate anchoring: tilt-out-of-plane cost is independent of in-plane reference direction", "[energy-density]") {
    // GIVEN: Degenerate anchoring: K1=0, K2=1, W>0, v2=(0,0,1)
    // WHEN:  tilting by the same out-of-plane angle γ but starting from different azimuthal directions
    // THEN:  the energy variation f(γ) − f(0) is the same for all starting in-plane directions
    const double W  = 1.0;
    const double K1 = 0.0;
    const double K2 = 1.0;
    const double S0 = 0.6;
    const Vec3 v1(0, 0, 0);
    const Vec3 v2(0, 0, 1);

    const double gamma = M_PI / 4.0;  // 45° out of plane
    const double expectedVariation = 1.5 * W * K2 * S0 * std::sin(gamma) * std::sin(gamma);

    // Test 5 different azimuthal starting directions p = (cos φ, sin φ, 0) = fromRadianAngles(0, φ, S0)
    const int N = 5;
    for (int i = 0; i < N; i++) {
        double phi = i * 2.0 * M_PI / N;
        // In-plane reference: tilt=0, twist=φ
        const qlc3d::TTensor t_ref = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(0, phi, S0));
        double f_ref = surfaceEnergyDensity(t_ref[0], t_ref[1], t_ref[2], t_ref[3], t_ref[4],
                                             v1, v2, W, K1, K2, S0);
        // Tilted: n = cos(γ)·p + sin(γ)·v2 = fromRadianAngles(tilt=γ, twist=φ)
        const qlc3d::TTensor t = qlc3d::TTensor::fromDirector(qlc3d::Director::fromRadianAngles(gamma, phi, S0));
        double f = surfaceEnergyDensity(t[0], t[1], t[2], t[3], t[4], v1, v2, W, K1, K2, S0);
        REQUIRE((f - f_ref) == Approx(expectedVariation).epsilon(1e-10));
    }
}
