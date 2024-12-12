#include <catch.h>
#include <io/vtkiofun.h>
#include <mesh.h>
#include <test-util.h>
#include <solutionvector.h>
#include <lc-representation.h>
#include <geom/coordinates.h>

TEST_CASE("write VTK unstructured ascii grid") {
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
    REQUIRE(lines[tetsStartLine + 4] == "10");

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