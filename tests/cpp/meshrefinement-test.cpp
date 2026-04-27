#include <catch.h>
#include <test-util.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <io/meshreader.h>
#include <simu.h>
#include <electrodes.h>
#include <refinement.h>
#include <refinement/refinement-spec.h>
#include <eventhandler.h>
#include <simulation-state.h>
#include <regulargrid.h>
#include <memory>
#include <list>

namespace {

template <typename MeshType>
double sumDeterminants(const MeshType &mesh) {
    double sum = 0.0;
    for (int i = 0; i < mesh.getnElements(); i++) {
        sum += mesh.getDeterminant(i);
    }
    return sum;
}

std::unique_ptr<Simu> makeSimu() {
    return std::unique_ptr<Simu>(SimuBuilder().build());
}

Alignment makeAlignment() {
    Alignment alignment;
    alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
    alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
    return alignment;
}

Electrodes makeElectrodes() {
    std::vector<std::shared_ptr<Electrode>> electrodesVec;
    electrodesVec.emplace_back(std::make_shared<Electrode>(1, std::vector<double>{0}, std::vector<double>{0}));
    electrodesVec.emplace_back(std::make_shared<Electrode>(2, std::vector<double>{0}, std::vector<double>{0}));
    return Electrodes::withElectrodePotentials(electrodesVec);
}

struct RefinementTestSetup {
    Geometry originalGeometry;
    Geometry workingGeometry;
    std::unique_ptr<Simu> simu;
    Alignment alignment;
    Electrodes electrodes;
    SimulationState simulationState;
    std::unique_ptr<RegularGrid> regGrid;
};

void prepareGeometryPair(Geometry &originalGeometry, Geometry &workingGeometry, const std::string &meshPath,
                         unsigned int expectedElementOrder) {
    auto rawMeshData = MeshReader::readMesh(meshPath);
    REQUIRE(rawMeshData.getElementOrder() == expectedElementOrder);
    const auto coordinates = std::make_shared<Coordinates>(std::move(rawMeshData.points));
    originalGeometry.setMeshData(rawMeshData.getElementOrder(), coordinates,
                                 std::move(rawMeshData.tetNodes), std::move(rawMeshData.tetMaterials),
                                 std::move(rawMeshData.triNodes), std::move(rawMeshData.triMaterials));
    workingGeometry.setTo(&originalGeometry);
}

RefinementTestSetup makeRefinementTestSetup(const std::string &meshPath, unsigned int expectedElementOrder) {
    RefinementTestSetup setup;
    setup.simu = makeSimu();
    setup.alignment = makeAlignment();
    setup.electrodes = makeElectrodes();
    prepareGeometryPair(setup.originalGeometry, setup.workingGeometry, meshPath, expectedElementOrder);
    return setup;
}

void populateChangeField(SolutionVector &qCurrent) {
    for (idx n = 0; n < qCurrent.getnDoF(); ++n) {
        for (idx dim = 0; dim < qCurrent.getnDimensions(); ++dim) {
            qCurrent.setValue(n, dim, static_cast<double>(n) + 0.1 * static_cast<double>(dim));
        }
    }
}

} // namespace

TEST_CASE("mesh refinement") {
    // ARRANGE - create required two geometries that are identical
    auto setup = makeRefinementTestSetup(TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, 1);

    unsigned int npLC = setup.originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    SolutionVector potential(npLC, 1);

    // refinement region is a single spherical volume
    auto spec = RefinementSpec::makePeriodic("Sphere", {0.1}, {0}, {0}, {0});
    std::vector<const RefinementSpec*> specs = {spec.get()};

    double S0 = 0.5;

    // ACT - run mesh refinement
    REQUIRE(autoref(setup.originalGeometry, setup.workingGeometry, qCurrent, potential, specs, *setup.simu,
                    setup.simulationState, setup.alignment, setup.electrodes, S0, setup.regGrid));

    // ASSERT
    // check new refined mesh size
    REQUIRE(setup.workingGeometry.getnpLC() == 50);
    REQUIRE(setup.workingGeometry.getTetrahedra().getnElements() == 140);
    REQUIRE(setup.workingGeometry.getTriangles().getnElements() == 80);

    // check initialisations. Sum of element determinants should not change, as total mesh volume/surface areas do not change
    auto &originalTets = setup.originalGeometry.getTetrahedra();
    auto &originalTris = setup.originalGeometry.getTriangles();

    double sumDeterminantsOriginal = sumDeterminants(originalTets);
    double sumDeterminantsOriginalTris = sumDeterminants(originalTris);

    auto &workingTets = setup.workingGeometry.getTetrahedra();
    auto &workingTris = setup.workingGeometry.getTriangles();

    double sumDeterminantsWorking = sumDeterminants(workingTets);
    double sumDeterminantsWorkingTris = sumDeterminants(workingTris);

    REQUIRE(sumDeterminantsOriginal > 0);
    REQUIRE(sumDeterminantsOriginal == Approx(sumDeterminantsWorking).margin(1e-18));

    REQUIRE(sumDeterminantsOriginalTris > 0);
    REQUIRE(sumDeterminantsOriginalTris == Approx(sumDeterminantsWorkingTris).margin(1e-18));

    // time step should be restricted to minimum after fine registration
    REQUIRE(setup.simulationState.restrictedTimeStep());
    REQUIRE(setup.simulationState.dt() == setup.simu->getMindt());
    REQUIRE(setup.simulationState.meshModified());
}

TEST_CASE("quadratic mesh refinement") {
    // ARRANGE - create matching quadratic geometries and refinement inputs
    auto setup = makeRefinementTestSetup(TestUtil::RESOURCE_SMALL_CUBE_QUADRATIC_GMSH_MESH, 2);

    unsigned int npLC = setup.originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    SolutionVector potential(npLC, 1);

    auto spec = RefinementSpec::makePeriodic("Sphere", {0.1}, {0}, {0}, {0});
    std::vector<const RefinementSpec*> specs = {spec.get()};

    double S0 = 0.5;

    // ACT - run mesh refinement on quadratic input
    REQUIRE(autoref(setup.originalGeometry, setup.workingGeometry, qCurrent, potential, specs, *setup.simu,
                    setup.simulationState, setup.alignment, setup.electrodes, S0, setup.regGrid));

    // ASSERT - the refined volume mesh must be quadratic and conserve total measures
    REQUIRE(setup.workingGeometry.getTetrahedra().getElementType() == ElementType::QUADRATIC_TETRAHEDRON);
    REQUIRE(setup.workingGeometry.getnpLC() > setup.originalGeometry.getnpLC());
    REQUIRE(setup.workingGeometry.getTetrahedra().getnElements() > setup.originalGeometry.getTetrahedra().getnElements());
    REQUIRE(setup.workingGeometry.getTriangles().getnElements() > setup.originalGeometry.getTriangles().getnElements());

    REQUIRE(sumDeterminants(setup.originalGeometry.getTetrahedra()) > 0);
    REQUIRE(sumDeterminants(setup.originalGeometry.getTetrahedra()) ==
            Approx(sumDeterminants(setup.workingGeometry.getTetrahedra())).margin(1e-14));
    REQUIRE(sumDeterminants(setup.originalGeometry.getTriangles()) > 0);
    REQUIRE(sumDeterminants(setup.originalGeometry.getTriangles()) ==
            Approx(sumDeterminants(setup.workingGeometry.getTriangles())).margin(1e-14));

    REQUIRE(setup.simulationState.restrictedTimeStep());
    REQUIRE(setup.simulationState.dt() == setup.simu->getMindt());
    REQUIRE(setup.simulationState.meshModified());
}

TEST_CASE("box refinement end-to-end") {
    auto setup = makeRefinementTestSetup(TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, 1);

    unsigned int npLC = setup.originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    SolutionVector potential(npLC, 1);

    auto spec = RefinementSpec::makePeriodic("Box", {}, {0.0, 0.5}, {0.0, 1.0}, {0.0, 1.0});
    std::vector<const RefinementSpec*> specs = {spec.get()};

    double S0 = 0.5;

    REQUIRE(autoref(setup.originalGeometry, setup.workingGeometry, qCurrent, potential, specs, *setup.simu,
                    setup.simulationState, setup.alignment, setup.electrodes, S0, setup.regGrid));

    REQUIRE(setup.workingGeometry.getnpLC() > setup.originalGeometry.getnpLC());
    REQUIRE(setup.workingGeometry.getTetrahedra().getnElements() > setup.originalGeometry.getTetrahedra().getnElements());
    REQUIRE(setup.workingGeometry.getTriangles().getnElements() > setup.originalGeometry.getTriangles().getnElements());

    REQUIRE(sumDeterminants(setup.originalGeometry.getTetrahedra()) > 0);
    REQUIRE(sumDeterminants(setup.originalGeometry.getTetrahedra()) ==
            Approx(sumDeterminants(setup.workingGeometry.getTetrahedra())).margin(1e-18));
    REQUIRE(sumDeterminants(setup.originalGeometry.getTriangles()) > 0);
    REQUIRE(sumDeterminants(setup.originalGeometry.getTriangles()) ==
            Approx(sumDeterminants(setup.workingGeometry.getTriangles())).margin(1e-18));
    REQUIRE(setup.simulationState.meshModified());
}

TEST_CASE("change refinement end-to-end") {
    auto setup = makeRefinementTestSetup(TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, 1);

    unsigned int npLC = setup.originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    populateChangeField(qCurrent);
    SolutionVector potential(npLC, 1);

    auto spec = RefinementSpec::makePeriodic("Change", {0.0}, {}, {}, {});
    std::vector<const RefinementSpec*> specs = {spec.get()};

    double S0 = 0.5;

    REQUIRE(autoref(setup.originalGeometry, setup.workingGeometry, qCurrent, potential, specs, *setup.simu,
                    setup.simulationState, setup.alignment, setup.electrodes, S0, setup.regGrid));

    REQUIRE(setup.workingGeometry.getnpLC() > setup.originalGeometry.getnpLC());
    REQUIRE(setup.workingGeometry.getTetrahedra().getnElements() > setup.originalGeometry.getTetrahedra().getnElements());
    REQUIRE(setup.workingGeometry.getTriangles().getnElements() > setup.originalGeometry.getTriangles().getnElements());
    REQUIRE(sumDeterminants(setup.originalGeometry.getTetrahedra()) ==
            Approx(sumDeterminants(setup.workingGeometry.getTetrahedra())).margin(1e-18));
    REQUIRE(sumDeterminants(setup.originalGeometry.getTriangles()) ==
            Approx(sumDeterminants(setup.workingGeometry.getTriangles())).margin(1e-18));
    REQUIRE(setup.simulationState.meshModified());
}

TEST_CASE("quadratic box refinement") {
    // GIVEN a quadratic TET10 geometry and a box refinement spec covering the lower-x half
    auto setup = makeRefinementTestSetup(TestUtil::RESOURCE_SMALL_CUBE_QUADRATIC_GMSH_MESH, 2);

    unsigned int npLC = setup.originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    SolutionVector potential(npLC, 1);

    // Box covers lower half in x: [0, 0.5] × [0, 1] × [0, 1]
    auto spec = RefinementSpec::makePeriodic("Box", {}, {0.0, 0.5}, {0.0, 1.0}, {0.0, 1.0});
    std::vector<const RefinementSpec*> specs = {spec.get()};

    double S0 = 0.5;

    // WHEN autoref is called on a quadratic mesh
    REQUIRE(autoref(setup.originalGeometry, setup.workingGeometry, qCurrent, potential, specs, *setup.simu,
                    setup.simulationState, setup.alignment, setup.electrodes, S0, setup.regGrid));

    // THEN the refined mesh is still quadratic, has more elements and nodes, and conserves volume
    REQUIRE(setup.workingGeometry.getTetrahedra().getElementType() == ElementType::QUADRATIC_TETRAHEDRON);
    REQUIRE(setup.workingGeometry.getnpLC() > setup.originalGeometry.getnpLC());
    REQUIRE(setup.workingGeometry.getTetrahedra().getnElements() > setup.originalGeometry.getTetrahedra().getnElements());
    REQUIRE(setup.workingGeometry.getTriangles().getnElements() > setup.originalGeometry.getTriangles().getnElements());

    REQUIRE(sumDeterminants(setup.originalGeometry.getTetrahedra()) > 0);
    REQUIRE(sumDeterminants(setup.originalGeometry.getTetrahedra()) ==
            Approx(sumDeterminants(setup.workingGeometry.getTetrahedra())).margin(1e-14));
    REQUIRE(sumDeterminants(setup.originalGeometry.getTriangles()) > 0);
    REQUIRE(sumDeterminants(setup.originalGeometry.getTriangles()) ==
            Approx(sumDeterminants(setup.workingGeometry.getTriangles())).margin(1e-14));

    REQUIRE(setup.simulationState.meshModified());
}

TEST_CASE("quadratic change refinement") {
    // GIVEN a quadratic TET10 geometry and a change-based refinement spec with a non-trivial Q-tensor field
    auto setup = makeRefinementTestSetup(TestUtil::RESOURCE_SMALL_CUBE_QUADRATIC_GMSH_MESH, 2);

    unsigned int npLC = setup.originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    populateChangeField(qCurrent);   // give a non-uniform Q tensor so Change detection has data to work with
    SolutionVector potential(npLC, 1);

    auto spec = RefinementSpec::makePeriodic("Change", {0.0}, {}, {}, {});
    std::vector<const RefinementSpec*> specs = {spec.get()};

    double S0 = 0.5;

    // WHEN autoref is called on a quadratic mesh
    REQUIRE(autoref(setup.originalGeometry, setup.workingGeometry, qCurrent, potential, specs, *setup.simu,
                    setup.simulationState, setup.alignment, setup.electrodes, S0, setup.regGrid));

    // THEN the refined mesh is still quadratic, has more elements and nodes, and conserves volume
    REQUIRE(setup.workingGeometry.getTetrahedra().getElementType() == ElementType::QUADRATIC_TETRAHEDRON);
    REQUIRE(setup.workingGeometry.getnpLC() > setup.originalGeometry.getnpLC());
    REQUIRE(setup.workingGeometry.getTetrahedra().getnElements() > setup.originalGeometry.getTetrahedra().getnElements());
    REQUIRE(setup.workingGeometry.getTriangles().getnElements() > setup.originalGeometry.getTriangles().getnElements());

    REQUIRE(sumDeterminants(setup.originalGeometry.getTetrahedra()) > 0);
    REQUIRE(sumDeterminants(setup.originalGeometry.getTetrahedra()) ==
            Approx(sumDeterminants(setup.workingGeometry.getTetrahedra())).margin(1e-14));
    REQUIRE(sumDeterminants(setup.originalGeometry.getTriangles()) > 0);
    REQUIRE(sumDeterminants(setup.originalGeometry.getTriangles()) ==
            Approx(sumDeterminants(setup.workingGeometry.getTriangles())).margin(1e-14));

    REQUIRE(setup.simulationState.meshModified());
}

TEST_CASE("handleMeshRefinement smoke test with RefinementSpec") {
    // WP-5 smoke test: verify handleMeshRefinement works end-to-end using the new
    // typed RefinementSpec event path (no RefInfo, no void* cast).
    auto setup = makeRefinementTestSetup(TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, 1);

    unsigned int npLC = setup.originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    SolutionVector potential(npLC, 1);

    // Build a typed refinement event using makeRefinement (the new path)
    auto spec = RefinementSpec::makePeriodic("Sphere", {0.15}, {0.5}, {0.5}, {0.5});
    auto refEvent = new Event(Event::makeRefinement(std::move(spec), 1u));
    REQUIRE(refEvent->getRefinementSpec() != nullptr);

    std::list<Event*> refEvents = {refEvent};

    Geometries geometries;
    geometries.geom_orig = &setup.originalGeometry;
    geometries.geom = &setup.workingGeometry;

    SolutionVectors solutionVectors;
    solutionVectors.q = &qCurrent;
    solutionVectors.v = &potential;

    // ACT - should not crash; may or may not refine depending on geometry
    REQUIRE_NOTHROW(handleMeshRefinement(refEvents, geometries, solutionVectors,
                                         *setup.simu, setup.simulationState, setup.alignment, setup.electrodes,
                                         0.5, setup.regGrid));

    // refEvents ownership was transferred to handleMeshRefinement, which deletes them
    // The important thing is no crash and no memory corruption
}

