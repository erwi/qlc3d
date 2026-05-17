#include <catch.h>
#include <energy/lc-energy-integrator.h>
#include <energy/lc-energy-density.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <mesh/mesh.h>
#include <solutionvector.h>
#include <lc.h>
#include <lc-representation.h>
#include <alignment.h>
#include <material_numbers.h>
#include <globals.h>
#include <cmath>
#include <memory>

// ============================================================================
// Helper: build a unit tetrahedron geometry suitable for energy integration tests.
//
// Nodes at (0,0,0), (1,0,0), (0,1,0), (0,0,1) in µm.
// Volume = 1/6 µm³ = 1/6 × 10⁻¹⁸ m³.
//
// @param elementOrder 1 for TET4, 2 for TET10.
// ============================================================================
static Geometry makeSingleTetGeometry(unsigned int elementOrder) {
    // Corner nodes in µm (qlc3d coordinate unit)
    std::vector<Vec3> pts = {
        Vec3(0.0, 0.0, 0.0),  // n0
        Vec3(1.0, 0.0, 0.0),  // n1
        Vec3(0.0, 1.0, 0.0),  // n2
        Vec3(0.0, 0.0, 1.0),  // n3
    };
    std::vector<unsigned int> tetNodes, triNodes;
    if (elementOrder == 1) {
        // TET4: 4 corner nodes
        tetNodes = {0, 1, 2, 3};
        // A dummy linear triangle face (required by setMeshData)
        triNodes = {0, 1, 2};
    } else {
        // TET10: 4 corners + 6 mid-edge nodes
        pts.push_back(Vec3(0.5, 0.0, 0.0)); // n4 = mid(n0,n1)
        pts.push_back(Vec3(0.5, 0.5, 0.0)); // n5 = mid(n1,n2)
        pts.push_back(Vec3(0.0, 0.5, 0.0)); // n6 = mid(n0,n2)
        pts.push_back(Vec3(0.0, 0.0, 0.5)); // n7 = mid(n0,n3)
        pts.push_back(Vec3(0.5, 0.0, 0.5)); // n8 = mid(n1,n3)  (Gmsh [8])
        pts.push_back(Vec3(0.0, 0.5, 0.5)); // n9 = mid(n2,n3)  (Gmsh [9])
        tetNodes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        // A dummy quadratic triangle face (TRI6): corners n0,n1,n2 + midpoints n4,n5,n6
        triNodes = {0, 1, 2, 4, 5, 6};
    }
    auto coords = std::make_shared<Coordinates>(std::move(pts));
    Geometry geom;
    geom.setMeshData(elementOrder, coords,
                     std::move(tetNodes), {MAT_DOMAIN1},
                     std::move(triNodes), {MAT_FIXLC1});
    return geom;
}

// Returns an LC with typical nematic parameters and no flexoelectric effect
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

// Set all Q-tensor DOFs to a uniform ground state aligned along z
// q1 = S0*sqrt(6)/2, q2=q3=q4=q5=0
static void setGroundStateQ(SolutionVector &q, double S0) {
    const qlc3d::Director d(0, 0, 1, S0);
    const qlc3d::TTensor t = qlc3d::TTensor::fromDirector(d);
    for (idx i = 0; i < q.getnDoF(); i++) {
        q.setValue(i, t);
    }
}

// ============================================================================
// Test: TET4 — uniform ground state, no field
// All energy contributions should be ~0.
// ============================================================================
TEST_CASE("integrateVolumeEnergy: TET4 uniform ground state gives ~zero energy", "[volume-integrator]") {
    // ARRANGE
    auto lc = makeTestLC();
    Geometry geom = makeSingleTetGeometry(1);
    const unsigned int np = geom.getnp();
    const unsigned int npLC = geom.getnpLC();

    // Zero potential
    SolutionVector v(np, 1);
    v.setValuesTo(0.0);

    // Ground state Q-tensor
    SolutionVector q(npLC, 5);
    setGroundStateQ(q, lc->S0());

    // ACT
    EnergyResult result = integrateVolumeEnergy(*lc, geom, v, q);

    // ASSERT: uniform state → zero elastic (no gradients), zero electric (no field).
    // Thermotropic is the raw ground-state value (negative, not zero).
    const double volume_m3 = (1.0 / 6.0) * 1e-18;
    EnergyMaterialParams mat = lcToEnergyMaterialParams(*lc);
    const double S0 = mat.S0;
    double f0 = (3.0*mat.A/4.0)*(S0*S0) + (mat.B/4.0)*(S0*S0*S0) + (9.0*mat.C/16.0)*(S0*S0*S0*S0);
    REQUIRE(result.elastic      == Approx(0.0).margin(1e-30));
    REQUIRE(result.thermotropic == Approx(f0 * volume_m3).epsilon(1e-10));
    REQUIRE(result.electric     == Approx(0.0).margin(1e-30));
}

// ============================================================================
// Test: TET4 — uniform linear potential (uniform E field) with ground-state Q along z
// Electric energy should match analytical value: electricEnergyDensity * volume.
// ============================================================================
TEST_CASE("integrateVolumeEnergy: TET4 uniform E-field electric energy matches analytic", "[volume-integrator]") {
    // ARRANGE
    auto lc = makeTestLC();
    Geometry geom = makeSingleTetGeometry(1);
    const unsigned int np = geom.getnp();
    const unsigned int npLC = geom.getnpLC();

    // Linear potential φ = E0*z applied at nodes (z in µm → convert)
    // Coordinates are (0,0,0),(1,0,0),(0,1,0),(0,0,1) in µm
    // We set φ = E0_SI * z_µm * 1e-6 at each node to create uniform E_z = -E0_SI V/m
    const double E0_SI = 1e6; // 1 MV/m in z direction
    SolutionVector v(np, 1);
    // Node z-coordinates (in µm): n0=0, n1=0, n2=0, n3=1µm
    // φ_i = E0 * z_i in V (z in m), so φ_3 = E0 * 1e-6
    v.setValue(0, 0, 0.0);
    v.setValue(1, 0, 0.0);
    v.setValue(2, 0, 0.0);
    v.setValue(3, 0, E0_SI * 1e-6); // E_z = -∂φ/∂z → ∂φ/∂z = +E0 → φ increases with z

    // Ground state Q-tensor along z
    SolutionVector q(npLC, 5);
    setGroundStateQ(q, lc->S0());

    // ACT
    EnergyResult result = integrateVolumeEnergy(*lc, geom, v, q);

    // ASSERT: electric energy should match density * volume (in m³)
    // Tetrahedron volume = 1/6 µm³ = 1/6 × 10⁻¹⁸ m³
    // Gradient Vz = ∂φ/∂z_m = E0_SI V/m → Ex=0, Ey=0, Ez = -Vz = -E0_SI
    const double volume_m3 = (1.0 / 6.0) * 1e-18;
    EnergyMaterialParams mat = lcToEnergyMaterialParams(*lc);
    GaussPointData gp;
    gp.q1 = std::sqrt(6.0) / 2.0 * lc->S0();
    gp.Ez = -E0_SI; // E = -∇φ, ∇φ points in +z, so E points in -z
    // Note: the source sets φ increasing with z → E_z = -∂φ/∂z < 0 if ∂φ/∂z > 0
    // But in the test we set φ_3 = E0*1e-6 > 0, so ∂φ/∂z > 0 → E_z < 0... but the sign
    // of E² doesn't matter for energy — just check magnitude consistency.
    double expectedDensity = electricEnergyDensity(gp, mat);
    double expected = expectedDensity * volume_m3;

    REQUIRE(result.electric == Approx(expected).epsilon(1e-6));
}

// ============================================================================
// Test: TET10 — same geometry as TET4 (same corners, midpoints at midpoints),
// thermotropic and elastic energies should be equal to the TET4 result.
// This confirms element-order independence for geometry-equivalent meshes.
// ============================================================================
TEST_CASE("integrateVolumeEnergy: TET10 gives same energy as TET4 for identical geometry", "[volume-integrator]") {
    // ARRANGE
    auto lc = makeTestLC();
    Geometry geom4  = makeSingleTetGeometry(1);  // TET4
    Geometry geom10 = makeSingleTetGeometry(2);  // TET10

    // Zero potential for both
    SolutionVector v4(geom4.getnp(), 1);
    v4.setValuesTo(0.0);
    SolutionVector v10(geom10.getnp(), 1);
    v10.setValuesTo(0.0);

    // Ground state Q-tensor for both (all nodes, including midpoints)
    SolutionVector q4(geom4.getnpLC(), 5);
    setGroundStateQ(q4, lc->S0());
    SolutionVector q10(geom10.getnpLC(), 5);
    setGroundStateQ(q10, lc->S0());

    // ACT
    EnergyResult r4  = integrateVolumeEnergy(*lc, geom4,  v4,  q4);
    EnergyResult r10 = integrateVolumeEnergy(*lc, geom10, v10, q10);

    // ASSERT: both should give the same thermotropic energy (ground state, no field)
    // and zero elastic and electric
    REQUIRE(r4.elastic       == Approx(0.0).margin(1e-30));
    REQUIRE(r10.elastic      == Approx(0.0).margin(1e-30));
    REQUIRE(r4.thermotropic  == Approx(r10.thermotropic).epsilon(1e-10));
    REQUIRE(r4.electric      == Approx(r10.electric).margin(1e-30));
}

// ============================================================================
// Test: Surface energy integration with a single weak-anchoring triangle face.
//
// Build a unit tetrahedron; face {n0,n1,n2} is in the z=0 plane, area = 0.5 µm².
// Use tilt=0, twist=0 so that v1=(0,1,0), v2=(0,0,1) (from Surface::calculateV1/V2).
// Easy axis e = v1×v2 = (1,0,0).
// Use K1=0, K2=1 to isolate the v2·Q·v2 contribution.
//
// Ground state Q along z: v̂₂=(0,0,1) so v̂₂·Q·v̂₂ = S0 (Q33 in qlc3d convention).
//
// Raw energy density at ground state (no ground-state offset applied):
//   fIso    = W*(K1+K2)/(S0*6) * R   = W*K2/(S0*6) * (3*S0²/2) = W*K2*S0/4
//   fOrient = W*K2*(v̂₂·Q·v̂₂)        = W*K2*S0
//   f_raw   = W*5*K2*S0/4
//
// Expected total = f_raw * area_m2
// ============================================================================
TEST_CASE("integrateSurfaceEnergy: single weak-anchoring face with ground-state Q matches analytic", "[surface-integrator]") {
    // ARRANGE
    const double S0 = 0.6;
    const double W = 1e-3;  // anchoring strength [J/m²]
    const double K1 = 0.0;  // no penalty along v1
    const double K2 = 1.0;  // penalise deviation from v2=(0,0,1)

    // Build a single tetrahedron; the triangle {0,1,2} lies in the z=0 plane.
    Geometry geom = makeSingleTetGeometry(1); // TET4

    // Ground state Q along z
    SolutionVector q(geom.getnpLC(), 5);
    setGroundStateQ(q, S0);

    // Create a weak surface on FIXLC1. tilt=0, twist=0 gives:
    //   v1=(0,1,0), v2=(0,0,1)  [from Surface::calculateV1/V2 with a=0,b=0,g=0]
    // With K1=0 only the K2*(v̂₂·Q·v̂₂) term contributes.
    Alignment alignment;
    alignment.addSurface(Surface::ofWeakAnchoring(1, 0.0, 0.0, W, K1, K2));

    // ACT
    double surfaceEnergy = integrateSurfaceEnergy(alignment, geom, q, S0);

    // ASSERT: raw analytical expectation (no ground-state offset subtracted).
    // Physical area of triangle = 0.5 µm² = 0.5e-12 m²
    //
    // Raw density at n=z (v2 direction):
    //   fIso    = W*K2/(S0*6) * (3*S0²/2) = W*K2*S0/4
    //   fOrient = W*K2*S0  (v̂₂·Q·v̂₂ = S0 when n∥v̂₂)
    //   f_raw   = W*(S0/4 + S0) = W*5*S0/4
    const double fRaw    = W * (K1 + K2) / (S0 * 6.0) * (3.0 * S0 * S0 / 2.0)  // fIso
                         + W * K2 * S0;                                            // fOrient (v2·Q·v2=S0)
    const double area_m2 = 0.5e-12; // 0.5 µm² in m²
    const double expected = fRaw * area_m2;

    REQUIRE(surfaceEnergy == Approx(expected).margin(std::abs(expected) * 1e-4));
}

TEST_CASE("integrateSurfaceEnergy: surface energy is minimum (negative) when LC is aligned with easy axis", "[surface-integrator]") {
    // GIVEN: same geometry but Q aligned along the easy axis e=(1,0,0) = v1×v2
    //        WeakAnchoring with K1=0, K2=1, W=1e-3, tilt=0, twist=0.
    //        At n=(1,0,0): v1·Q·v1 = v2·Q·v2 = -S0/2 (both axes ⊥ easy axis).
    //        f_raw = W*K2/(S0*6)*(3S0²/2)  +  W*K2*(-S0/2)
    //              = W*K2*S0/4 - W*K2*S0/2 = -W*K2*S0/4
    // WHEN:  integrateSurfaceEnergy is called
    // THEN:  the result is negative and equals -W*K2*S0/4 * area
    const double S0 = 0.6;
    const double W  = 1e-3;
    const double K1 = 0.0;
    const double K2 = 1.0;

    Geometry geom = makeSingleTetGeometry(1);

    // Set Q along x = easy axis e = v1×v2 = (0,1,0)×(0,0,1) = (1,0,0)
    SolutionVector q(geom.getnpLC(), 5);
    const qlc3d::Director dx(1, 0, 0, S0);
    const qlc3d::TTensor tx = qlc3d::TTensor::fromDirector(dx);
    for (idx i = 0; i < q.getnDoF(); i++) {
        q.setValue(i, tx);
    }

    Alignment alignment;
    alignment.addSurface(Surface::ofWeakAnchoring(1, 0.0, 0.0, W, K1, K2));

    double surfaceEnergy = integrateSurfaceEnergy(alignment, geom, q, S0);

    // Expected: -W*K2*S0/4 * area_m2
    const double area_m2 = 0.5e-12;
    const double expected = -W * K2 * S0 / 4.0 * area_m2;
    REQUIRE(surfaceEnergy == Approx(expected).margin(std::abs(expected) * 1e-4));
    REQUIRE(surfaceEnergy < 0.0); // minimum is negative, not zero
}
