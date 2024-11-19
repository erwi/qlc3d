#include <catch.h>
#include <memory>
#include <io/lcview-result-output.h>
#include <resultio.h>
#include <util/stringutil.h>
#include <test-util.h>
#include <simulation-state.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <solutionvector.h>
#include <lc-representation.h>

namespace fs = std::filesystem;

const unsigned int LC_MATERIAL = 4;
const unsigned int PERIODIC_BC = 3;

std::shared_ptr<Geometry> createSingleTetGeometry() {
  std::shared_ptr<Geometry> geom = std::make_shared<Geometry>();
  int np = 4;
  std::vector<Vec3> points = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
  std::shared_ptr<Coordinates> coordinates = std::make_shared<Coordinates>(std::move(points));
  geom->setCoordinates(coordinates);

  std::shared_ptr<Mesh> tetrahedra = Mesh::tetMesh();
  tetrahedra->setElementData(1, {0, 1, 2, 3}, {LC_MATERIAL});
  geom->setTetrahedra(tetrahedra);

  std::shared_ptr<Mesh> triangles = Mesh::triangleMesh();
  triangles->setElementData(1, {0, 1, 2}, {PERIODIC_BC});
  geom->setTriangles(triangles);

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

  std::vector<qlc3d::Director> director = {{1., 0., 0., 0.5},
                                           {1., 0., 0., 0.5},
                                           {1., 0., 0., 0.5},
                                           {1., 0., 0., 0.5}};

  SimulationState simulationState;
  simulationState.state(RunningState::RUNNING);
  simulationState.currentIteration(0);

  SECTION("Writes result and mesh file to same directory") {
    TestUtil::TemporaryDirectory resDir;
    LcViewTxtResultFormatWriter writer(resDir.path(), meshName, S0);

    REQUIRE(writer.isDirectorRequired()); // LCViewTxt format does require director file

    writer.setPotential(potential);
    writer.setQTensor(qTensor);
    writer.setDirector(&director);
    writer.writeResult(*geom, simulationState);

    // Both files should be written to the same directory
    REQUIRE(fs::exists(resDir.path() / "mesh0.msh"));
    REQUIRE(fs::exists(resDir.path() / "result-t-00000.dat"));

    // WHEN:
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
    writer.setDirector(&director);
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

void shouldEqual(const SolutionVector &q1, const SolutionVector &q2) {
  REQUIRE(q1.getnDoF() == q2.getnDoF());
  REQUIRE(q1.getnDimensions() == q2.getnDimensions());
  // check that the read values equal the written values
  for (int i = 0; i < q1.getnDoF(); i++) {
    for (int j = 0; j < q1.getnDimensions(); j++) {
      REQUIRE(q1.getValue(i, j) == Approx(q2.getValue(i, j)).margin(1e-6)); // saved using float precision
    }
  }
}

TEST_CASE("Write and read back Q-tensor as LCView format") {
  const std::string meshName = "path/to/mesh.msh";
  const double S0 = 0.5;
  const auto geom = createSingleTetGeometry();

  const SolutionVector potential(geom->getnpLC(), 1);
  SolutionVector qTensor(geom->getnpLC(), 5);

  // create q-tensor with arbitrary but valid values
  for (int i = 0; i < geom->getnpLC(); i++) {
    auto dir = qlc3d::Director::fromDegreeAngles(i * 10, i * 13, 0.5);
    qTensor.setValue(i, dir);
  }

  SimulationState simulationState;
  simulationState.state(RunningState::RUNNING);
  simulationState.currentIteration(0);
  TestUtil::TemporaryDirectory resDir;

  SECTION("Binary LCView file format") {
    LcViewBinaryResultFormatWriter writer(resDir.path(), meshName, S0);

    writer.setPotential(potential);
    writer.setQTensor(qTensor);
    writer.writeResult(*geom, simulationState);

    // WHEN:
    // Read back the Q-tensor from the binary file
    SolutionVector qTensorRead(geom->getnpLC(), 5);

    const std::string resultFile = (resDir.path() / "result00000.dat").string();
    REQUIRE(fs::exists(resultFile));
    ResultIO::ReadResult(resultFile, qTensorRead);

    // THEN:
    REQUIRE(qTensorRead.getnDoF() == geom->getnpLC());

    // check that the read values equal the written values
    shouldEqual(qTensor, qTensorRead);
  }

  SECTION("Text LCView file format") {
    LcViewTxtResultFormatWriter writer(resDir.path(), meshName, S0);

    // Text LCView requires director, not q-tensor
    std::vector<qlc3d::Director> director;
    for (int i = 0; i < geom->getnpLC(); i++) {
      director.push_back(qTensor.getDirector(i));
    }

    writer.setPotential(potential);
    writer.setDirector(&director);
    writer.writeResult(*geom, simulationState);

    const std::string resultFile = (resDir.path() / "result-t-00000.dat").string();
    REQUIRE(fs::exists(resultFile));
    // WHEN:
    // Read back the Q-tensor from the text file
    SolutionVector qTensorRead(geom->getnpLC(), 5);
    ResultIO::ReadResult(resultFile, qTensorRead);

    // THEN:
    // Check that the read values equal the written values
    shouldEqual(qTensor, qTensorRead);
  }
}

