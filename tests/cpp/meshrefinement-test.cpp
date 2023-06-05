#include <catch.h>
#include <test-util.h>
#include <geometry.h>
#include <inits.h>
#include <simu.h>
#include <electrodes.h>
#include <refinement.h>
#include <simulation-state.h>

TEST_CASE("mesh refinement") {
    // ARRANGE - create required two geometries that are identical
    Geometry originalGeometry;
    Geometry workingGeometry;

    std::unique_ptr<Simu> simu(SimuBuilder().build());
    Alignment alignment;
    alignment.addSurface(1, "strong", 1, {1, 0, 0}, 1, 1, {});
    alignment.addSurface(2, "strong", 1, {1, 0, 0}, 1, 1, {});

    Electrodes electrodes;
    electrodes.setnElectrodes(2); // Mesh contains two electrodes. This fakes them being defined in the settings file

    // reads and prepares test-mesh from resource file
    prepareGeometry(originalGeometry, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, *simu, alignment, electrodes);

    workingGeometry.setTo(&originalGeometry);

    unsigned int npLC = originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    SolutionVector qPrevious(npLC, 5);
    SolutionVector potential(npLC, 1);

    // refinement region is a single spherical volume
    std::unique_ptr<RefInfo> refInfo(RefInfo::make("Sphere", 0, 0, {0.1}, {0}, {0}, {0}));
    std::list<RefInfo> refInfos = {*refInfo};
    SimulationState simulationState;

    double S0 = 0.5;

    // ACT - run mesh refinement
    autoref(originalGeometry, workingGeometry, qCurrent, qPrevious, potential, refInfos, *simu, simulationState, alignment, electrodes, S0);

    // ASSERT
    // check new refined mesh size
    REQUIRE(workingGeometry.getnpLC() == 50);
    REQUIRE(workingGeometry.getTetrahedra().getnElements() == 140);
    REQUIRE(workingGeometry.getTriangles().getnElements() == 80);

    // check initialisations. Sum of element determinants should not change, as total mesh volume/surface areas do not change
    auto &originalTets = originalGeometry.getTetrahedra();
    double sumDeterminantsOriginal = std::accumulate(
            originalTets.getPtrToDeterminant(0),
            originalTets.getPtrToDeterminant(originalTets.getnElements()),
            0.);

    auto &workingTets = workingGeometry.getTetrahedra();
    double sumDeterminantsWorking = std::accumulate(
            workingTets.getPtrToDeterminant(0),
            workingTets.getPtrToDeterminant(workingTets.getnElements()),
            0.);

    REQUIRE(sumDeterminantsOriginal > 0);
    REQUIRE(sumDeterminantsOriginal == Approx(sumDeterminantsWorking).margin(1e-18));

    // time step should be restricted to minimum after fine registration
    REQUIRE(simulationState.restrictedTimeStep());
    REQUIRE(simulationState.dt() == simu->getMindt());
}

// TODO: test cases to check solution interpolation onto refined mesh