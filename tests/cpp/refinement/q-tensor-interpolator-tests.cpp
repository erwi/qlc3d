#include <catch.h>
#include <refinement/q-tensor-interpolator.h>
#include <test-util.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <inits.h>
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

Alignment makeAlignment() {
    Alignment alignment;
    alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
    alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
    return alignment;
}

Alignment makeAlignmentWithoutOverride() {
    Alignment alignment;
    alignment.addSurface(Surface::ofFreeze(1));
    alignment.addSurface(Surface::ofFreeze(2));
    return alignment;
}

std::unique_ptr<Simu> makeSimu() {
    return std::unique_ptr<Simu>(SimuBuilder().build());
}

Electrodes makeElectrodes() {
    std::vector<std::shared_ptr<Electrode>> electrodesVec;
    electrodesVec.emplace_back(std::shared_ptr<Electrode>(new Electrode(1, {0}, {0})));
    electrodesVec.emplace_back(std::shared_ptr<Electrode>(new Electrode(2, {0}, {0})));
    return Electrodes::withElectrodePotentials(electrodesVec);
}

/** Load the small cube mesh into geomOld and copy it to geomNew (same mesh, no refinement).
 *  This gives us a trivially valid pair for testing the corner-node fast path. */
void loadMeshPair(Geometry &geomOld, Geometry &geomNew) {
    Alignment alignment = makeAlignment();
    Electrodes electrodes = makeElectrodes();
    prepareGeometry(geomOld, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, electrodes, alignment, {1, 1, 1});
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

    const int npLC = geomOld.getnpLC();
    SolutionVector qOld(npLC, 5);
    // Set a uniform non-trivial Q field
    for (int n = 0; n < npLC; n++) {
        for (int d = 0; d < 5; d++) {
            qOld.setValue(n, d, 1.0 + 0.1 * d);
        }
    }

    SolutionVector result = interpolateQTensor(geomNew, geomOld, qOld);

    for (int n = 0; n < npLC; n++) {
        for (int d = 0; d < 5; d++) {
            REQUIRE(result.getValue(n, d) == Approx(qOld.getValue(n, d)).epsilon(1e-12));
        }
    }
}

TEST_CASE("interpolateQTensor_linearFieldSurvivesRefinement") {
    Geometry geomOld, geomNew;
    std::unique_ptr<Simu> simu = makeSimu();
    Alignment alignment = makeAlignmentWithoutOverride();
    Electrodes electrodes = makeElectrodes();
    prepareGeometry(geomOld, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, electrodes, alignment, {1, 1, 1});
    geomNew.setTo(&geomOld);

    SolutionVector qCurrent(geomOld.getnpLC(), 5);
    setLinearQField(qCurrent, geomOld);
    SolutionVector potential(geomOld.getnpLC(), 1);
    SimulationState simulationState;

    auto spec = RefinementSpec::makePeriodic("Sphere", {0.1}, {0}, {0}, {0});
    std::vector<const RefinementSpec*> specs = {spec.get()};

    REQUIRE(autoref(geomOld, geomNew, qCurrent, potential, specs, *simu, simulationState, alignment, electrodes, 0.5));

    REQUIRE(qCurrent.getnDoF() == geomNew.getnpLC());
    REQUIRE(qCurrent.getnDimensions() == 5);

    for (idx n = 0; n < geomNew.getnpLC(); ++n) {
        const Vec3 p = geomNew.getCoordinates().getPoint(n);
        for (int d = 0; d < 5; ++d) {
            REQUIRE(qCurrent.getValue(n, d) == Approx(evaluateLinearQComponent(p, d)).margin(1e-10));
        }
    }
}

