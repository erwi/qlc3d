#include <catch.h>
#include <test-util.h>
#include <geometry.h>
#include <inits.h>
#include <simu.h>
#include <electrodes.h>
#include <refinement.h>
#include <simulation-state.h>
#include <memory>

TEST_CASE("mesh refinement") {
    // ARRANGE - create required two geometries that are identical
    Geometry originalGeometry;
    Geometry workingGeometry;

    std::unique_ptr<Simu> simu(SimuBuilder().build());
    Alignment alignment;
    alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
    alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));

    std::vector<std::shared_ptr<Electrode>> electrodesVec;
    electrodesVec.emplace_back(std::shared_ptr<Electrode>(new Electrode(1, {0}, {0})));
    electrodesVec.emplace_back(std::shared_ptr<Electrode>(new Electrode(2, {0}, {0})));

    // Mesh contains two electrodes. This fakes them being defined in the settings file
    Electrodes electrodes = Electrodes::withElectrodePotentials(electrodesVec);

    // reads and prepares test-mesh from resource file
    prepareGeometry(originalGeometry, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, electrodes, alignment, {1, 1, 1});

    workingGeometry.setTo(&originalGeometry);

    unsigned int npLC = originalGeometry.getnpLC();
    SolutionVector qCurrent(npLC, 5);
    SolutionVector potential(npLC, 1);

    // refinement region is a single spherical volume
    std::unique_ptr<RefInfo> refInfo(RefInfo::make("Sphere", 0, 0, {0.1}, {0}, {0}, {0}));
    std::list<RefInfo> refInfos = {*refInfo};
    SimulationState simulationState;

    double S0 = 0.5;

    // ACT - run mesh refinement
    autoref(originalGeometry, workingGeometry, qCurrent, potential, refInfos, *simu, simulationState, alignment, electrodes, S0);

    // ASSERT
    // check new refined mesh size
    REQUIRE(workingGeometry.getnpLC() == 50);
    REQUIRE(workingGeometry.getTetrahedra().getnElements() == 140);
    REQUIRE(workingGeometry.getTriangles().getnElements() == 80);

    // check initialisations. Sum of element determinants should not change, as total mesh volume/surface areas do not change
    auto &originalTets = originalGeometry.getTetrahedra();
    auto &originalTris = originalGeometry.getTriangles();

    double sumDeterminantsOriginal = 0;
    for (int i = 0; i < originalTets.getnElements(); i++) {
      sumDeterminantsOriginal += originalTets.getDeterminant(i);
    }

    double sumDeterminantsOriginalTris = 0;
    for (int i = 0; i < originalTris.getnElements(); i++) {
      sumDeterminantsOriginalTris += originalTris.getDeterminant(i);
    }

    auto &workingTets = workingGeometry.getTetrahedra();
    auto &workingTris = workingGeometry.getTriangles();

    double sumDeterminantsWorking = 0;
    for (int i = 0; i < workingTets.getnElements(); i++) {
      sumDeterminantsWorking += workingTets.getDeterminant(i);
    }

    double sumDeterminantsWorkingTris = 0;
    for (int i = 0; i < workingTris.getnElements(); i++) {
      sumDeterminantsWorkingTris += workingTris.getDeterminant(i);
    }

    // TODO: the correct determinants for volumes and surfaces should be calculated from the known mesh size and compared here instead
    REQUIRE(sumDeterminantsOriginal > 0);
    REQUIRE(sumDeterminantsOriginal == Approx(sumDeterminantsWorking).margin(1e-18));

    REQUIRE(sumDeterminantsOriginalTris > 0);
    REQUIRE(sumDeterminantsOriginalTris == Approx(sumDeterminantsWorkingTris).margin(1e-18));

    // time step should be restricted to minimum after fine registration
    REQUIRE(simulationState.restrictedTimeStep());
    REQUIRE(simulationState.dt() == simu->getMindt());
}

// TODO: test cases to check solution interpolation onto refined mesh