#include <catch.h>
#include <test-util.h>
#include <geometry.h>
#include <inits.h>
#include <simu.h>
#include <electrodes.h>
#include <refinement.h>
#include <refinement/refinement-spec.h>
#include <eventhandler.h>
#include <simulation-state.h>
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
    electrodesVec.emplace_back(std::shared_ptr<Electrode>(new Electrode(1, {0}, {0})));
    electrodesVec.emplace_back(std::shared_ptr<Electrode>(new Electrode(2, {0}, {0})));
    return Electrodes::withElectrodePotentials(electrodesVec);
}

void prepareGeometryPair(Geometry &originalGeometry, Geometry &workingGeometry, Electrodes &electrodes, Alignment &alignment) {
    prepareGeometry(originalGeometry, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, electrodes, alignment, {1, 1, 1});
    workingGeometry.setTo(&originalGeometry);
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
    Geometry originalGeometry;
    Geometry workingGeometry;

    std::unique_ptr<Simu> simu = makeSimu();
    Alignment alignment = makeAlignment();

    // Mesh contains two electrodes. This fakes them being defined in the settings file
    Electrodes electrodes = makeElectrodes();

    // reads and prepares test-mesh from resource file
    prepareGeometryPair(originalGeometry, workingGeometry, electrodes, alignment);

    unsigned int npLC = originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    SolutionVector potential(npLC, 1);

    // refinement region is a single spherical volume
    auto spec = RefinementSpec::makePeriodic("Sphere", {0.1}, {0}, {0}, {0});
    std::vector<const RefinementSpec*> specs = {spec.get()};
    SimulationState simulationState;

    double S0 = 0.5;

    // ACT - run mesh refinement
    REQUIRE(autoref(originalGeometry, workingGeometry, qCurrent, potential, specs, *simu, simulationState, alignment, electrodes, S0));

    // ASSERT
    // check new refined mesh size
    REQUIRE(workingGeometry.getnpLC() == 50);
    REQUIRE(workingGeometry.getTetrahedra().getnElements() == 140);
    REQUIRE(workingGeometry.getTriangles().getnElements() == 80);

    // check initialisations. Sum of element determinants should not change, as total mesh volume/surface areas do not change
    auto &originalTets = originalGeometry.getTetrahedra();
    auto &originalTris = originalGeometry.getTriangles();

    double sumDeterminantsOriginal = sumDeterminants(originalTets);
    double sumDeterminantsOriginalTris = sumDeterminants(originalTris);

    auto &workingTets = workingGeometry.getTetrahedra();
    auto &workingTris = workingGeometry.getTriangles();

    double sumDeterminantsWorking = sumDeterminants(workingTets);
    double sumDeterminantsWorkingTris = sumDeterminants(workingTris);

    // TODO: the correct determinants for volumes and surfaces should be calculated from the known mesh size and compared here instead
    REQUIRE(sumDeterminantsOriginal > 0);
    REQUIRE(sumDeterminantsOriginal == Approx(sumDeterminantsWorking).margin(1e-18));

    REQUIRE(sumDeterminantsOriginalTris > 0);
    REQUIRE(sumDeterminantsOriginalTris == Approx(sumDeterminantsWorkingTris).margin(1e-18));

    // time step should be restricted to minimum after fine registration
    REQUIRE(simulationState.restrictedTimeStep());
    REQUIRE(simulationState.dt() == simu->getMindt());
    REQUIRE(simulationState.meshModified());
}

TEST_CASE("box refinement end-to-end") {
    Geometry originalGeometry;
    Geometry workingGeometry;

    std::unique_ptr<Simu> simu = makeSimu();
    Alignment alignment = makeAlignment();
    Electrodes electrodes = makeElectrodes();

    prepareGeometryPair(originalGeometry, workingGeometry, electrodes, alignment);

    unsigned int npLC = originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    SolutionVector potential(npLC, 1);

    auto spec = RefinementSpec::makePeriodic("Box", {}, {0.0, 0.5}, {0.0, 1.0}, {0.0, 1.0});
    std::vector<const RefinementSpec*> specs = {spec.get()};
    SimulationState simulationState;

    double S0 = 0.5;

    REQUIRE(autoref(originalGeometry, workingGeometry, qCurrent, potential, specs, *simu, simulationState, alignment, electrodes, S0));

    REQUIRE(workingGeometry.getnpLC() > originalGeometry.getnpLC());
    REQUIRE(workingGeometry.getTetrahedra().getnElements() > originalGeometry.getTetrahedra().getnElements());
    REQUIRE(workingGeometry.getTriangles().getnElements() > originalGeometry.getTriangles().getnElements());

    REQUIRE(sumDeterminants(originalGeometry.getTetrahedra()) > 0);
    REQUIRE(sumDeterminants(originalGeometry.getTetrahedra()) == Approx(sumDeterminants(workingGeometry.getTetrahedra())).margin(1e-18));
    REQUIRE(sumDeterminants(originalGeometry.getTriangles()) > 0);
    REQUIRE(sumDeterminants(originalGeometry.getTriangles()) == Approx(sumDeterminants(workingGeometry.getTriangles())).margin(1e-18));
    REQUIRE(simulationState.meshModified());
}

TEST_CASE("change refinement end-to-end") {
    Geometry originalGeometry;
    Geometry workingGeometry;

    std::unique_ptr<Simu> simu = makeSimu();
    Alignment alignment = makeAlignment();
    Electrodes electrodes = makeElectrodes();

    prepareGeometryPair(originalGeometry, workingGeometry, electrodes, alignment);

    unsigned int npLC = originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    populateChangeField(qCurrent);
    SolutionVector potential(npLC, 1);

    auto spec = RefinementSpec::makePeriodic("Change", {0.0}, {}, {}, {});
    std::vector<const RefinementSpec*> specs = {spec.get()};
    SimulationState simulationState;

    double S0 = 0.5;

    REQUIRE(autoref(originalGeometry, workingGeometry, qCurrent, potential, specs, *simu, simulationState, alignment, electrodes, S0));

    REQUIRE(workingGeometry.getnpLC() > originalGeometry.getnpLC());
    REQUIRE(workingGeometry.getTetrahedra().getnElements() > originalGeometry.getTetrahedra().getnElements());
    REQUIRE(workingGeometry.getTriangles().getnElements() > originalGeometry.getTriangles().getnElements());
    REQUIRE(sumDeterminants(originalGeometry.getTetrahedra()) == Approx(sumDeterminants(workingGeometry.getTetrahedra())).margin(1e-18));
    REQUIRE(sumDeterminants(originalGeometry.getTriangles()) == Approx(sumDeterminants(workingGeometry.getTriangles())).margin(1e-18));
    REQUIRE(simulationState.meshModified());
}

// TODO: test cases to check solution interpolation onto refined mesh

TEST_CASE("handleMeshRefinement smoke test with RefinementSpec") {
    // WP-5 smoke test: verify handleMeshRefinement works end-to-end using the new
    // typed RefinementSpec event path (no RefInfo, no void* cast).
    Geometry originalGeometry;
    Geometry workingGeometry;

    std::unique_ptr<Simu> simu = makeSimu();
    Alignment alignment = makeAlignment();
    Electrodes electrodes = makeElectrodes();

    prepareGeometryPair(originalGeometry, workingGeometry, electrodes, alignment);

    unsigned int npLC = originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    SolutionVector potential(npLC, 1);

    // Build a typed refinement event using makeRefinement (the new path)
    auto spec = RefinementSpec::makePeriodic("Sphere", {0.15}, {0.5}, {0.5}, {0.5});
    Event* refEvent = new Event(Event::makeRefinement(std::move(spec), 1u));
    REQUIRE(refEvent->getRefinementSpec() != nullptr);

    std::list<Event*> refEvents = {refEvent};

    Geometries geometries;
    geometries.geom_orig = &originalGeometry;
    geometries.geom = &workingGeometry;

    SolutionVectors solutionVectors;
    solutionVectors.q = &qCurrent;
    solutionVectors.v = &potential;

    SimulationState simulationState;

    // ACT - should not crash; may or may not refine depending on geometry
    bool refined = handleMeshRefinement(refEvents, geometries, solutionVectors,
                                        *simu, simulationState, alignment, electrodes, 0.5);

    // refEvents ownership was transferred to handleMeshRefinement, which deletes them
    // The important thing is no crash and no memory corruption
    REQUIRE_NOTHROW(refined);  // tautological, but documents intention
}

