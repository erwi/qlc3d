#include <catch.h>
#include <io/lcview-result-output.h>
#include <util/stringutil.h>
#include <test-util.h>
#include <simulation-state.h>
#include <geometry.h>
#include <solutionvector.h>

namespace fs = std::filesystem;

std::shared_ptr<Geometry> createSingleTetGeometry() {
  std::shared_ptr<Geometry> geom = std::make_shared<Geometry>();
  int np = 4;
  double coords[] = {0, 0, 0,
                     1, 0, 0,
                     0,  1, 0,
                     0, 0, 1};
  geom->setCoordinates(coords, np);
  Mesh *tetrahedra = new Mesh(1, np); // NOTE: geom takes ownership of tetrahedra and triangles and frees them in its destructor
  Mesh *triangles = new Mesh(1, 3);

  tetrahedra->setnNodes(np);
  idx tetNodes[] = {0, 1, 2, 3};
  tetrahedra->setAllNodes(tetNodes);

  geom->t = tetrahedra;


  triangles->setnNodes(3);
  idx triNodes[] = {0, 1, 2};
  triangles->setAllNodes(triNodes);
  geom->e = triangles;
  return geom;
}

TEST_CASE("Write binary LCView result file") {
  const std::string meshName = "path/to/mesh.msh";
  const double S0 = 0.5;
  const auto geom = createSingleTetGeometry();

  const SolutionVector potential(geom->getnpLC(), 1);
  const SolutionVector qTensor(geom->getnpLC(), 5);

  SimulationState simulationState;
  simulationState.state(RunningState::RUNNING);
  simulationState.currentIteration(0);

  SECTION("Writes result and mesh file to same directory") {
    TestUtil::TemporaryDirectory resDir;
    LcViewBinaryResultFormatWriter writer(resDir.path(), meshName, S0);

    REQUIRE_FALSE(writer.isDirectorRequired()); // LCView format does not require director file

    writer.setPotential(potential);
    writer.setQTensor(qTensor);
    writer.writeResult(*geom, simulationState);

    // Both files should be written to the same directory
    REQUIRE(fs::exists(resDir.path() / "mesh0.msh"));
    REQUIRE(fs::exists(resDir.path() / "result00000.dat"));
  }

  SECTION("New mesh file should be created when mesh number changes") {
    TestUtil::TemporaryDirectory resDir;
    simulationState.currentIteration(0);
    LcViewBinaryResultFormatWriter writer(resDir.path(), meshName, S0);

    writer.setPotential(potential);
    writer.setQTensor(qTensor);
    writer.writeResult(*geom, simulationState);

    REQUIRE(fs::exists(resDir.path() / "mesh0.msh"));
    REQUIRE(fs::exists(resDir.path() / "result00000.dat"));
    REQUIRE(resDir.listFiles().size() == 2); // 1 mesh file, 1 result file

    // WHEN:
    // Next iteration, including mesh refinement
    simulationState.incrementMeshNumber();
    simulationState.currentIteration(1);
    writer.writeResult(*geom, simulationState);

    // THEN:
    REQUIRE(fs::exists(resDir.path() / "mesh1.msh"));
    REQUIRE(fs::exists(resDir.path() / "result00001.dat"));
    REQUIRE(resDir.listFiles().size() == 4); // 2 mesh files, 2 result files

    // WHEN:
    // Simulation is completed
    simulationState.state(RunningState::COMPLETED);
    writer.writeResult(*geom, simulationState);

    // THEN:
    REQUIRE(fs::exists(resDir.path() / "result-final.dat"));
  }
}

TEST_CASE("Write text LCViewTxt result file") {
  const std::string meshName = "path/to/mesh.msh";
  const double S0 = 0.5;
  const auto geom = createSingleTetGeometry();

  const SolutionVector potential(geom->getnpLC(), 1);
  const SolutionVector qTensor(geom->getnpLC(), 5);
  const double director[] = {1, 0, 0,
                             1, 0, 0,
                             1, 0, 0,
                             1, 0, 0};
  SimulationState simulationState;
  simulationState.state(RunningState::RUNNING);
  simulationState.currentIteration(0);

  SECTION("Writes result and mesh file to same directory") {
    TestUtil::TemporaryDirectory resDir;
    LcViewTxtResultFormatWriter writer(resDir.path(), meshName, S0);

    REQUIRE(writer.isDirectorRequired()); // LCViewTxt format does require director file

    writer.setPotential(potential);
    writer.setQTensor(qTensor);
    writer.setDirector(director);
    writer.writeResult(*geom, simulationState);

    // Both files should be written to the same directory
    REQUIRE(fs::exists(resDir.path() / "mesh0.msh"));
    REQUIRE(fs::exists(resDir.path() / "result-t-00000.dat"));

    // WHNEN:
    // Simulation is completed
    simulationState.state(RunningState::COMPLETED);
    writer.writeResult(*geom, simulationState);
    REQUIRE(fs::exists(resDir.path() / "result-t-final.dat"));
  }

  SECTION("New mesh file should be created when mesh number changes") {
    TestUtil::TemporaryDirectory resDir;
    simulationState.currentIteration(0);
    LcViewTxtResultFormatWriter writer(resDir.path(), meshName, S0);

    writer.setPotential(potential);
    writer.setQTensor(qTensor);
    writer.setDirector(director);
    writer.writeResult(*geom, simulationState);

    REQUIRE(fs::exists(resDir.path() / "mesh0.msh"));
    REQUIRE(fs::exists(resDir.path() / "result-t-00000.dat"));
    REQUIRE(resDir.listFiles().size() == 2); // 1 mesh file, 1 result file

    // WHEN:
    // Next iteration, including mesh refinement
    simulationState.incrementMeshNumber();
    simulationState.currentIteration(1);
    writer.writeResult(*geom, simulationState);

    // THEN:
    REQUIRE(fs::exists(resDir.path() / "mesh1.msh"));
    REQUIRE(fs::exists(resDir.path() / "result-t-00001.dat"));
    REQUIRE(resDir.listFiles().size() == 4); // 2 mesh files, 2 result files
  }
}