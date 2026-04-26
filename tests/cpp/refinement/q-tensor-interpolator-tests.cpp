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

