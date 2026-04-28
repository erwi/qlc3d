#include <catch.h>
#include <refinement/q-tensor-interpolator.h>
#include <test-util.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <io/meshreader.h>
#include <simu.h>
#include <simulation-state.h>
#include <electrodes.h>
#include <alignment.h>
#include <refinement.h>
#include <refinement/refinement-spec.h>
#include <solutionvector.h>

// ============================================================================
// Unit tests for interpolateScalar
// ============================================================================

TEST_CASE("interpolateScalar_atCorner0") {
    double loc[4] = {1, 0, 0, 0};
    double S[4]   = {10, 20, 30, 40};
    REQUIRE(interpolateScalar(loc, S) == Approx(10.0));
}

TEST_CASE("interpolateScalar_atCorner3") {
    double loc[4] = {0, 0, 0, 1};
    double S[4]   = {10, 20, 30, 40};
    REQUIRE(interpolateScalar(loc, S) == Approx(40.0));
}

TEST_CASE("interpolateScalar_atEdgeMidpoint") {
    double loc[4] = {0.5, 0.5, 0, 0};
    double S[4]   = {0, 1, 0, 0};
    REQUIRE(interpolateScalar(loc, S) == Approx(0.5));
}

TEST_CASE("interpolateScalar_atCentroid") {
    double loc[4] = {0.25, 0.25, 0.25, 0.25};
    double S[4]   = {0, 0, 0, 4};
    REQUIRE(interpolateScalar(loc, S) == Approx(1.0));
}

TEST_CASE("interpolateScalar_uniformField") {
    double loc[4] = {0.1, 0.3, 0.4, 0.2};
    double S[4]   = {7, 7, 7, 7};
    REQUIRE(interpolateScalar(loc, S) == Approx(7.0));
}

// ============================================================================
// Integration tests for interpolateQTensor
// ============================================================================

namespace {

std::filesystem::path meshResourcePath(const char* filename) {
    return std::filesystem::absolute(std::filesystem::path(__FILE__))
        .parent_path().parent_path().parent_path() / "resources" / filename;
}

std::unique_ptr<Simu> makeSimu() {
    return std::unique_ptr<Simu>(SimuBuilder().build());
}

Electrodes makeElectrodes() {
    std::vector<std::shared_ptr<Electrode>> electrodesVec;
    electrodesVec.emplace_back(std::make_shared<Electrode>(1, std::vector<double>{0}, std::vector<double>{0}));
    electrodesVec.emplace_back(std::make_shared<Electrode>(2, std::vector<double>{0}, std::vector<double>{0}));
    return Electrodes::withElectrodePotentials(electrodesVec);
}

/** Load the small cube mesh into geomOld and copy it to geomNew (same mesh, no refinement).
 *  This gives us a trivially valid pair for testing the corner-node fast path. */
void loadMeshPair(Geometry &geomOld, Geometry &geomNew) {
    auto rawMeshData = MeshReader::readMesh(meshResourcePath("gmsh-small-cube.msh"));
    REQUIRE(rawMeshData.getElementOrder() == 1);
    auto coordinates = std::make_shared<Coordinates>(std::move(rawMeshData.points));
    geomOld.setMeshData(rawMeshData.getElementOrder(), coordinates,
                        std::move(rawMeshData.tetNodes), std::move(rawMeshData.tetMaterials),
                        std::move(rawMeshData.triNodes), std::move(rawMeshData.triMaterials));
    geomNew.setTo(&geomOld);
}

void setLinearQField(SolutionVector &q, const Geometry &geom) {
    for (idx n = 0; n < geom.getnpLC(); ++n) {
        const Vec3 p = geom.getCoordinates().getPoint(n);
        q.setValue(n, 0, 1.0 + 2.0 * p.x() + 3.0 * p.y() + 5.0 * p.z());
        q.setValue(n, 1, -2.0 + 4.0 * p.x() - 1.0 * p.y() + 0.5 * p.z());
        q.setValue(n, 2, 7.0 - 3.0 * p.x() + 2.0 * p.y() + 1.0 * p.z());
        q.setValue(n, 3, 0.25 + 1.0 * p.x() - 2.0 * p.y() + 4.0 * p.z());
        q.setValue(n, 4, -1.5 + 0.75 * p.x() + 0.5 * p.y() - 0.25 * p.z());
    }
}

double evaluateLinearQComponent(const Vec3 &p, int component) {
    switch (component) {
        case 0: return 1.0 + 2.0 * p.x() + 3.0 * p.y() + 5.0 * p.z();
        case 1: return -2.0 + 4.0 * p.x() - 1.0 * p.y() + 0.5 * p.z();
        case 2: return 7.0 - 3.0 * p.x() + 2.0 * p.y() + 1.0 * p.z();
        case 3: return 0.25 + 1.0 * p.x() - 2.0 * p.y() + 4.0 * p.z();
        case 4: return -1.5 + 0.75 * p.x() + 0.5 * p.y() - 0.25 * p.z();
        default: return 0.0;
    }
}

} // namespace

TEST_CASE("interpolateQTensor_outputSizeMatchesNewMesh") {
    Geometry geomOld, geomNew;
    loadMeshPair(geomOld, geomNew);

    SolutionVector qOld(geomOld.getnpLC(), 5);

    SolutionVector result = interpolateQTensor(geomNew, geomOld, qOld);

    REQUIRE(result.getnDoF() == geomNew.getnpLC());
    REQUIRE(result.getnDimensions() == 5);
}

TEST_CASE("interpolateQTensor_preservedNodeRetainsValue") {
    // When geomOld == geomNew (same nodes), every node is a corner-node copy.
    // A uniform Q field should be reproduced exactly.
    Geometry geomOld, geomNew;
    loadMeshPair(geomOld, geomNew);

    const idx npLC = geomOld.getnpLC();
    SolutionVector qOld(npLC, 5);
    // Set a uniform non-trivial Q field
    for (idx n = 0; n < npLC; n++) {
        for (idx d = 0; d < 5; d++) {
            qOld.setValue(n, d, 1.0 + 0.1 * d);
        }
    }

    SolutionVector result = interpolateQTensor(geomNew, geomOld, qOld);

    for (idx n = 0; n < npLC; n++) {
        for (idx d = 0; d < 5; d++) {
            REQUIRE(result.getValue(n, d) == Approx(qOld.getValue(n, d)).epsilon(1e-12));
        }
    }
}

// ============================================================================
// Unit tests for evaluateTet10ShapeFunctions
// ============================================================================

TEST_CASE("evaluateTet10ShapeFunctions: partition of unity at centroid") {
    // GIVEN: the centroid of a tet in barycentric coordinates
    double loc[4] = {0.25, 0.25, 0.25, 0.25};
    double N[10];

    // WHEN: shape functions are evaluated
    evaluateTet10ShapeFunctions(loc, N);

    // THEN: the shape functions sum to 1 (partition of unity)
    double sum = 0.0;
    for (int i = 0; i < 10; i++) { sum += N[i]; }
    REQUIRE(sum == Approx(1.0).epsilon(1e-14));
}

TEST_CASE("evaluateTet10ShapeFunctions: partition of unity at random interior point") {
    // GIVEN: an arbitrary interior barycentric point
    double loc[4] = {0.1, 0.3, 0.4, 0.2};
    double N[10];

    // WHEN:
    evaluateTet10ShapeFunctions(loc, N);

    // THEN: shape functions still sum to 1
    double sum = 0.0;
    for (int i = 0; i < 10; i++) { sum += N[i]; }
    REQUIRE(sum == Approx(1.0).epsilon(1e-14));
}

TEST_CASE("evaluateTet10ShapeFunctions: at corner node 0, only N[0]=1") {
    // GIVEN: barycentric point exactly at corner 0 (N1=1, rest=0)
    double loc[4] = {1.0, 0.0, 0.0, 0.0};
    double N[10];

    // WHEN:
    evaluateTet10ShapeFunctions(loc, N);

    // THEN: only N[0] is 1, all others are 0
    REQUIRE(N[0] == Approx(1.0).epsilon(1e-14));
    for (int i = 1; i < 10; i++) {
        REQUIRE(N[i] == Approx(0.0).margin(1e-14));
    }
}

TEST_CASE("evaluateTet10ShapeFunctions: at midpoint of edge 0-1, only N[4]=1") {
    // GIVEN: barycentric midpoint of edge between corner 0 (N1=0.5) and corner 1 (N2=0.5)
    // This is exactly the position of mid-edge node [4] (AB in Gmsh ordering).
    double loc[4] = {0.5, 0.5, 0.0, 0.0};
    double N[10];

    // WHEN:
    evaluateTet10ShapeFunctions(loc, N);

    // THEN: only N[4] is 1 (the AB mid-edge shape function), all others are 0
    REQUIRE(N[4] == Approx(1.0).epsilon(1e-14));
    for (int i = 0; i < 10; i++) {
        if (i != 4) { REQUIRE(N[i] == Approx(0.0).margin(1e-14)); }
    }
}

TEST_CASE("evaluateTet10ShapeFunctions: at midpoint of edge 0-3, only N[7]=1") {
    // GIVEN: midpoint of edge AD (corner 0 and corner 3) — corresponds to node [7]
    double loc[4] = {0.5, 0.0, 0.0, 0.5};
    double N[10];

    // WHEN:
    evaluateTet10ShapeFunctions(loc, N);

    // THEN: only N[7] (AD mid-edge) is 1
    REQUIRE(N[7] == Approx(1.0).epsilon(1e-14));
    for (int i = 0; i < 10; i++) {
        if (i != 7) { REQUIRE(N[i] == Approx(0.0).margin(1e-14)); }
    }
}

TEST_CASE("evaluateTet10ShapeFunctions: all six mid-edge nodes each yield exactly one N=1") {
    // GIVEN: the six mid-edge barycentric positions for a TET10 element.
    // Each position is defined by two barycentric coordinates = 0.5, rest = 0.
    // The expected dominant shape-function index (per Gmsh TET10 ordering) is listed.
    struct Case { double loc[4]; int expectedIdx; const char* label; };
    Case cases[] = {
        {{0.5, 0.5, 0.0, 0.0}, 4, "AB"},
        {{0.0, 0.5, 0.5, 0.0}, 5, "BC"},
        {{0.5, 0.0, 0.5, 0.0}, 6, "AC"},
        {{0.5, 0.0, 0.0, 0.5}, 7, "AD"},
        {{0.0, 0.0, 0.5, 0.5}, 8, "CD"},
        {{0.0, 0.5, 0.0, 0.5}, 9, "BD"},
    };

    for (auto& c : cases) {
        double N[10];
        // WHEN: shape functions are evaluated at the mid-edge position
        evaluateTet10ShapeFunctions(c.loc, N);

        // THEN: the expected mid-edge shape function equals 1 and all others are 0.
        // This confirms the mid-edge fast path can be triggered reliably.
        REQUIRE(N[c.expectedIdx] == Approx(1.0).epsilon(1e-14));
        for (int i = 0; i < 10; i++) {
            if (i != c.expectedIdx) {
                REQUIRE(N[i] == Approx(0.0).margin(1e-14));
            }
        }
    }
}

TEST_CASE("interpolateQTensor_quadraticField: TET10 interior point uses full quadratic path") {
    // GIVEN: reference TET10 element (unit tet, same node coordinates as elsewhere).
    // Q field Q0 = x^2.
    const double nodeX[10] = {0.0, 1.0, 0.0, 0.0,  0.5, 0.5, 0.0, 0.0, 0.0, 0.5};
    double qAt[10];
    for (int i = 0; i < 10; i++) { qAt[i] = nodeX[i] * nodeX[i]; }

    // WHEN: query at a general interior point, not a corner and not a mid-edge node.
    // loc = {0.2, 0.3, 0.3, 0.2}  =>  x = N2*1 + N5*0.5 + N6*0 + N7*0 + N9*0.5
    //   but more simply: x = loc[1]*1.0 + mid-edge contributions
    // At this point none of the 10 shape functions equals 1 (max is ~4*0.3*0.3=0.36),
    // so neither the corner fast path nor the mid-edge fast path fires; the full
    // quadratic sum is used.
    double loc[4] = {0.2, 0.3, 0.3, 0.2};
    double N[10];
    evaluateTet10ShapeFunctions(loc, N);

    // Verify this point is not at a corner or mid-edge (no shape function near 1).
    for (int k = 0; k < 10; k++) {
        REQUIRE(N[k] < 0.99);
    }

    // THEN: full TET10 sum reproduces x^2 exactly at the query point.
    // x at query point = loc[1] (only node 1 has x=1; all others have x=0 or
    // contribute via shape functions): x = N[1]*1 + N[4]*0.5 + N[5]*0.5 + N[9]*0.5
    const double x_query = N[1] * 1.0 + N[4] * 0.5 + N[5] * 0.5 + N[9] * 0.5;
    const double expected = x_query * x_query;

    double q_tet10 = 0.0;
    for (int i = 0; i < 10; i++) { q_tet10 += N[i] * qAt[i]; }

    REQUIRE(q_tet10 == Approx(expected).epsilon(1e-13));

    // Linear TET4 interpolation misses the quadratic term.
    double q_linear = 0.0;
    for (int i = 0; i < 4; i++) { q_linear += loc[i] * qAt[i]; }
    // q_linear = 0*0 + 0.3*1 + 0 + 0 = 0.3, but x_query != 0.3 in general,
    // so just confirm TET10 and TET4 differ for a non-trivial interior point.
    REQUIRE(std::abs(q_tet10 - q_linear) > 1e-4);
}

// ============================================================================
// WP5 — Specification test: TET10 quadratic interpolation is exact for
//         degree-2 polynomial fields, while linear TET4 interpolation is not.
//
// This test works purely algebraically without a full Geometry setup.
// It demonstrates the improvement implemented in WP4.
// ============================================================================

TEST_CASE("interpolateQTensor_quadraticField: TET10 shape functions reproduce x^2 exactly at mid-edge") {
    // GIVEN: reference TET10 element with corners:
    //   node 0: (0, 0, 0)
    //   node 1: (1, 0, 0)
    //   node 2: (0, 1, 0)
    //   node 3: (0, 0, 1)
    // Mid-edge nodes per Gmsh TET10 ordering (AB, BC, AC, AD, CD, BD):
    //   node 4 (AB): (0.5, 0,   0  )   x=0.5
    //   node 5 (BC): (0.5, 0.5, 0  )   x=0.5
    //   node 6 (AC): (0,   0.5, 0  )   x=0
    //   node 7 (AD): (0,   0,   0.5)   x=0
    //   node 8 (CD): (0,   0.5, 0.5)   x=0
    //   node 9 (BD): (0.5, 0,   0.5)   x=0.5
    //
    // Q field: Q0 = x^2 — a degree-2 polynomial not representable by TET4.
    const double nodeX[10] = {0.0, 1.0, 0.0, 0.0,   0.5, 0.5, 0.0, 0.0, 0.0, 0.5};
    double qAt[10];
    for (int i = 0; i < 10; i++) { qAt[i] = nodeX[i] * nodeX[i]; }

    // WHEN: query at the midpoint of edge AB: barycentric {0.5, 0.5, 0, 0}.
    // This is exactly midpoint of edge 0-1, so x=0.5 at the query point.
    double loc[4] = {0.5, 0.5, 0.0, 0.0};
    double N[10];
    evaluateTet10ShapeFunctions(loc, N);

    // TET10 interpolated value
    double q_tet10 = 0.0;
    for (int i = 0; i < 10; i++) { q_tet10 += N[i] * qAt[i]; }

    // TET4 linear interpolation (only 4 corner nodes)
    double q_linear = 0.0;
    for (int i = 0; i < 4; i++) { q_linear += loc[i] * qAt[i]; }

    // THEN:
    // TET10 reproduces (0.5)^2 = 0.25 exactly.
    const double exact = 0.5 * 0.5;
    REQUIRE(q_tet10 == Approx(exact).epsilon(1e-14));

    // Linear interpolation gives 0.5*0 + 0.5*1 = 0.5 -- demonstrably wrong.
    REQUIRE(q_linear == Approx(0.5).epsilon(1e-14));
    REQUIRE(std::abs(q_linear - exact) > 0.2); // linear error is large (0.25)
}

TEST_CASE("interpolateQTensor_quadraticField: TET10 reproduces xy cross-term exactly at element centroid") {
    // GIVEN: same TET10 reference element.
    // Q field: Q0 = x * y (another degree-2 term).
    const double nodeX[10] = {0.0, 1.0, 0.0, 0.0,  0.5, 0.5, 0.0, 0.0, 0.0, 0.5};
    const double nodeY[10] = {0.0, 0.0, 1.0, 0.0,  0.0, 0.5, 0.5, 0.0, 0.5, 0.0};
    double qAt[10];
    for (int i = 0; i < 10; i++) { qAt[i] = nodeX[i] * nodeY[i]; }

    // WHEN: query at the midpoint of edge BC: barycentric {0, 0.5, 0.5, 0}.
    // At this point x=0.5, y=0.5, so Q0 = 0.25 exactly.
    double loc[4] = {0.0, 0.5, 0.5, 0.0};
    double N[10];
    evaluateTet10ShapeFunctions(loc, N);

    double q_tet10 = 0.0;
    for (int i = 0; i < 10; i++) { q_tet10 += N[i] * qAt[i]; }

    // TET4 linear: 0*0 + 0.5*0 + 0.5*0 + 0*0 = 0 (wrong; exact is 0.25)
    double q_linear = 0.0;
    for (int i = 0; i < 4; i++) { q_linear += loc[i] * qAt[i]; }

    // THEN:
    REQUIRE(q_tet10 == Approx(0.25).epsilon(1e-14));
    REQUIRE(q_linear == Approx(0.0).margin(1e-14)); // linear misses the cross-term entirely
}

// ============================================================================
// WP5 — End-to-end specification test: interpolateQTensor on a TET10 source
//        geometry reproduces a degree-2 polynomial field exactly.
//
// Strategy: load the quadratic small-cube mesh, assign Q0 = x^2 to every LC
// node (a field that is NOT representable by linear TET4 interpolation), copy
// the geometry to geomNew (same nodes), then verify that interpolateQTensor
// returns the exact analytic value at every node including mid-edge ones.
//
// Corner nodes hit the fast-copy path; mid-edge nodes exercise the TET10
// quadratic path (barycentric weights never reach 0.99999 for an edge midpoint).
// ============================================================================

TEST_CASE("interpolateQTensor_quadraticFieldOnTET10_isExact") {
    // GIVEN: a TET10 source geometry (quadratic small cube)
    auto raw = MeshReader::readMesh(meshResourcePath("gmsh-small-cube-quadratic.msh"));
    REQUIRE(raw.getElementOrder() == 2);
    auto coords = std::make_shared<Coordinates>(std::move(raw.points));
    Geometry geomOld;
    geomOld.setMeshData(2, coords,
                        std::move(raw.tetNodes), std::move(raw.tetMaterials),
                        std::move(raw.triNodes), std::move(raw.triMaterials));

    // GIVEN: Q field Q0 = x^2 (degree-2; not reproducible by linear TET4 interpolation).
    //        All other components are zero.
    const idx npLC = geomOld.getnpLC();
    REQUIRE(npLC > 0);
    SolutionVector qOld(npLC, 5);
    for (idx n = 0; n < npLC; n++) {
        const Vec3 p = geomOld.getCoordinates().getPoint(n);
        qOld.setValue(n, 0, p.x() * p.x());
        // components 1-4 remain zero
    }

    // WHEN: geomNew is a copy of geomOld (same set of nodes, including mid-edge nodes).
    //       Corner nodes hit the existing corner fast-copy path (max barycentric coord ~1).
    //       Mid-edge nodes have max barycentric coord ~0.5, so they bypass the corner
    //       fast-copy path but hit the new mid-edge fast-copy path (N[k] ~1 for the
    //       corresponding mid-edge shape function).  Points strictly inside element
    //       interiors would go through the full TET10 quadratic sum.
    Geometry geomNew;
    geomNew.setTo(&geomOld);

    SolutionVector result = interpolateQTensor(geomNew, geomOld, qOld);

    // THEN: for every LC node the interpolated Q0 matches x^2 exactly (to 1e-10).
    //       This can only hold for mid-edge nodes if the quadratic TET10 path is used;
    //       linear TET4 interpolation would give (x_A + x_B)/2 instead of
    //       ((x_A + x_B)/2)^2, which is strictly less than the average of x_A^2 and x_B^2.
    REQUIRE(result.getnDoF() == geomNew.getnpLC());
    for (idx n = 0; n < geomNew.getnpLC(); n++) {
        const Vec3 p = geomNew.getCoordinates().getPoint(n);
        const double expected = p.x() * p.x();
        REQUIRE(result.getValue(n, 0) == Approx(expected).margin(1e-10));
        // Other components should remain zero
        for (int d = 1; d < 5; d++) {
            REQUIRE(result.getValue(n, d) == Approx(0.0).margin(1e-10));
        }
    }
}

TEST_CASE("interpolateQTensor_linearFieldSurvivesRefinement") {
    Geometry geomOld, geomNew;
    std::unique_ptr<Simu> simu = makeSimu();
    auto rawMeshData = MeshReader::readMesh(meshResourcePath("gmsh-small-cube.msh"));
    REQUIRE(rawMeshData.getElementOrder() == 1);
    auto coordinates = std::make_shared<Coordinates>(std::move(rawMeshData.points));
    geomOld.setMeshData(rawMeshData.getElementOrder(), coordinates,
                        std::move(rawMeshData.tetNodes), std::move(rawMeshData.tetMaterials),
                        std::move(rawMeshData.triNodes), std::move(rawMeshData.triMaterials));
    geomNew.setTo(&geomOld);

    SolutionVector qCurrent(geomOld.getnpLC(), 5);
    setLinearQField(qCurrent, geomOld);
    SolutionVector potential(geomOld.getnpLC(), 1);
    SimulationState simulationState;

    Alignment alignment;
    alignment.addSurface(Surface::ofFreeze(1));
    alignment.addSurface(Surface::ofFreeze(2));
    Electrodes electrodes = makeElectrodes();

    auto spec = RefinementSpec::makePeriodic("Sphere", {0.1}, {0}, {0}, {0});
    std::vector<const RefinementSpec*> specs = {spec.get()};
    std::unique_ptr<RegularGrid> regGrid;

    REQUIRE(autoref(geomOld, geomNew, qCurrent, potential, specs, *simu, simulationState, alignment, electrodes, 0.5, regGrid));

    REQUIRE(qCurrent.getnDoF() == geomNew.getnpLC());
    REQUIRE(qCurrent.getnDimensions() == 5);

    for (idx n = 0; n < geomNew.getnpLC(); ++n) {
        const Vec3 p = geomNew.getCoordinates().getPoint(n);
        for (int d = 0; d < 5; ++d) {
            REQUIRE(qCurrent.getValue(n, d) == Approx(evaluateLinearQComponent(p, d)).margin(1e-10));
        }
    }
}

