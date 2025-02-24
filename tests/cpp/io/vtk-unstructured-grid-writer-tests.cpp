#include <catch.h>
#include <io/vtkiofun.h>
#include <mesh/mesh.h>
#include <test-util.h>
#include <solutionvector.h>
#include <lc-representation.h>
#include <geom/coordinates.h>
#include "util/stringutil.h"

TEST_CASE("write VTK unstructured ascii grid - linear 4 node tetrahedra") {
    using namespace vtkIOFun;
    using namespace std;

    // ARRANGE:
    // Set-up create dummy data consisting of a single 4-node tetrahedral mesh
    // and define some potential values and director vectors for each of the
    // mesh node values

    size_t numPoints = 4;
    size_t numLcPoints = 4;

    // 4 mesh node coordinates
    Coordinates coordinates(
          {{0, 0, 0},
           {1, 0, 0},
           {0, 1, 0},
           {0, 0, 1}});

    // 4 directors
    SolutionVector q(4, 5);

    q.setValue(0, qlc3d::Director(1, 0, 0, 0.1));
    q.setValue(1, qlc3d::Director(0, 1, 0, 0.2));
    q.setValue(2, qlc3d::Director(0, 0, 1, 0.3));

    qlc3d::Director d = qlc3d::Director::fromDegreeAngles(45, 45, 0.4);
    q.setValue(3, d);

    SolutionVector potentials(4, 1);
    potentials[0] = 0;
    potentials[1] = 0.1;
    potentials[2] = 0.2;
    potentials[3] = 0.3;

    // create mesh consisting of a single linear tetrahedron
    Mesh tetrahedra(3, ElementType::LINEAR_TETRAHEDRON);
    tetrahedra.setElementData(ElementType::LINEAR_TETRAHEDRON, {0, 1, 2, 3}, {4});

    // ACT: write the result to a temporary file
    auto tempFile = TestUtil::TemporaryFile::empty();
    UnstructuredGridWriter writer;
    writer.write(tempFile.name(),
                 numLcPoints,
                 coordinates,
                 tetrahedra,
                 potentials,
                 q);

    // ASSERT: check the file contents
    // Read the file contents and
    vector<string> lines = tempFile.readContentsAsText();

    REQUIRE(lines[0] == "# vtk DataFile Version 2.0");

    // check point coordinates section
    int pointsStartLine = 5;
    REQUIRE(lines[pointsStartLine] == "POINTS 4 float");
    REQUIRE(lines[pointsStartLine + 1] == "0 0 0");
    REQUIRE(lines[pointsStartLine + 2] == "1 0 0");
    REQUIRE(lines[pointsStartLine + 3] == "0 1 0");
    REQUIRE(lines[pointsStartLine + 4] == "0 0 1");

    // check tetrahedral cells section
    int tetsStartLine = 11;
    REQUIRE(lines[tetsStartLine] == "CELLS 1 5");
    REQUIRE(lines[tetsStartLine + 1] == "4 0 1 2 3");
    REQUIRE(lines[tetsStartLine + 3] == "CELL_TYPES 1");
    REQUIRE(lines[tetsStartLine + 4] == "10"); // 10 is linear tet, 24 is quadratic tet

    // check point data section
    int potentialsStartLine = 18;
    REQUIRE(lines[potentialsStartLine] == "SCALARS potential float 1");
    REQUIRE(lines[potentialsStartLine + 1] == "LOOKUP_TABLE default");
    REQUIRE(lines[potentialsStartLine + 2] == "0");
    REQUIRE(lines[potentialsStartLine + 3] == "0.1");
    REQUIRE(lines[potentialsStartLine + 4] == "0.2");
    REQUIRE(lines[potentialsStartLine + 5] == "0.3");

    // check director vector data section
    int directorStartLine = 25;
    REQUIRE(lines[directorStartLine] == "VECTORS director float");
    REQUIRE(lines[directorStartLine + 1] == "1 0 0");
    REQUIRE(lines[directorStartLine + 2] == "0 1 0");
    REQUIRE(lines[directorStartLine + 3] == "0 0 1");
    REQUIRE(lines[directorStartLine + 4] == "0.5 0.5 0.707107"); // 45 degrees tilt, twist

    // check order parameter data section
    int orderParamStartLine = 31;
    REQUIRE(lines[orderParamStartLine] == "SCALARS S float 1");
    REQUIRE(lines[orderParamStartLine + 1] == "LOOKUP_TABLE default");
    REQUIRE(lines[orderParamStartLine + 2] == "0.1");
    REQUIRE(lines[orderParamStartLine + 3] == "0.2");
    REQUIRE(lines[orderParamStartLine + 4] == "0.3");
    REQUIRE(lines[orderParamStartLine + 5] == "0.4");
}

TEST_CASE("write VTK unstructured ascii grid - quadratic 10 node tetrahedra") {
  using namespace vtkIOFun;
  using namespace std;
  unsigned int numLcPoints = 10;

  Coordinates coordinates(
          {{0, 0, 0},
           {1, 0, 0},
           {0, 1, 0},
           {0, 0, 1},
           {0.5, 0, 0},
           {0.5, 0.5, 0},
           {0, 0.5, 0},
           {0, 0, 0.5},
           {0.5, 0, 0.5},
           {0, 0.5, 0.5}});

  SolutionVector q(10, 5);
  SolutionVector potentials(10, 1);
  // create q-tensors whose tilt angle is proportional to the z-coordinate and potentials whose value equals the z-coordinate
  for (int i = 0; i < coordinates.size(); i++) {
    double z = coordinates.getPoint(i).z();
    q.setValue(i, qlc3d::Director::fromDegreeAngles(90 * z, 0, 0.5));
    potentials[i] = z;
  }

  Mesh tetrahedra(3, ElementType::QUADRATIC_TETRAHEDRON);
  tetrahedra.setElementData(ElementType::QUADRATIC_TETRAHEDRON, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {10});



  auto tempFile = TestUtil::TemporaryFile::empty();

  // WHEN
  UnstructuredGridWriter writer;
  writer.write(tempFile.name(),
               numLcPoints,
               coordinates,
               tetrahedra,
               potentials,
               q);

  // THEN
  vector<string> lines = tempFile.readContentsAsText();

  REQUIRE(lines[0] == "# vtk DataFile Version 2.0");

  int pointsStartLine = 5;
  REQUIRE(lines[pointsStartLine] == "POINTS 10 float");
  REQUIRE(lines[pointsStartLine + 1] == "0 0 0");
  REQUIRE(lines[pointsStartLine + 2] == "1 0 0");
  REQUIRE(lines[pointsStartLine + 3] == "0 1 0");
  REQUIRE(lines[pointsStartLine + 4] == "0 0 1");
  REQUIRE(lines[pointsStartLine + 5] == "0.5 0 0");
  REQUIRE(lines[pointsStartLine + 6] == "0.5 0.5 0");
  REQUIRE(lines[pointsStartLine + 7] == "0 0.5 0");
  REQUIRE(lines[pointsStartLine + 8] == "0 0 0.5");
  REQUIRE(lines[pointsStartLine + 9] == "0.5 0 0.5");
  REQUIRE(lines[pointsStartLine + 10] == "0 0.5 0.5");

  int tetsStartLine = 17;
  REQUIRE(lines[tetsStartLine] == "CELLS 1 11");
  REQUIRE(lines[tetsStartLine + 1] == "10 0 1 2 3 4 5 6 7 8 9");
  REQUIRE(lines[tetsStartLine + 3] == "CELL_TYPES 1"); // a singe cell type follows
  REQUIRE(lines[tetsStartLine + 4] == "24"); // 24 is quadratic tet type

  // Check potential data
  int potentialsStartLine = 24;
  REQUIRE(lines[potentialsStartLine] == "SCALARS potential float 1");
  REQUIRE(lines[potentialsStartLine + 1] == "LOOKUP_TABLE default");
  for (int i = 2; i <= 11; i++) {
    double expected = potentials.getValue(i - 2);
    double actual = std::stod(lines[potentialsStartLine + i]);
    REQUIRE(actual == Approx(expected).margin(1e-12));
  }

  // Check director data
  int directorStartLine = 37;
  REQUIRE(lines[directorStartLine] == "VECTORS director float");
  for (int i = 1; i <= 10; i++) {
    auto expectedDirector = q.getDirector(i - 1);
    auto l = lines[directorStartLine + i];
    auto splitLine = StringUtil::split(l, " ");

    REQUIRE(expectedDirector.nx() == Approx(std::stod(splitLine[0])).margin(1e-12));
    REQUIRE(expectedDirector.ny() == Approx(std::stod(splitLine[1])).margin(1e-12));
    REQUIRE(expectedDirector.nz() == Approx(std::stod(splitLine[2])).margin(1e-12));
  }

  int orderParamStartLine = 49;
  REQUIRE(lines[orderParamStartLine] == "SCALARS S float 1");
  REQUIRE(lines[orderParamStartLine + 1] == "LOOKUP_TABLE default");
  REQUIRE(lines[orderParamStartLine + 2] == "0.5");
  REQUIRE(lines[orderParamStartLine + 3] == "0.5");
  REQUIRE(lines[orderParamStartLine + 4] == "0.5");
  REQUIRE(lines[orderParamStartLine + 5] == "0.5");
  REQUIRE(lines[orderParamStartLine + 6] == "0.5");
  REQUIRE(lines[orderParamStartLine + 7] == "0.5");
  REQUIRE(lines[orderParamStartLine + 8] == "0.5");
  REQUIRE(lines[orderParamStartLine + 9] == "0.5");
  REQUIRE(lines[orderParamStartLine + 10] == "0.5");
}