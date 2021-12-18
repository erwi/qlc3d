#include <catch.h>
#include <io/vtkiofun.h>
#include <mesh.h>
#include <test-util.h>

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
    double points[] = {
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
    };

    // 4 directors
    double director[] = {
            1, 0, -1, 0, // 4x nx
            0, 1, 0, -1, // 4x ny
            0, 0, 0, 0   // 4x nz
    };

    // 4 order parameter values
    double S[] = {
            0.1, 0.2, 0.3, 0.4
    };

    double potentials[] = {0, 0.1, 0.2, 0.3};

    // create mesh consisting of a single tetrahedron
    Mesh tetrahedra(1, numPoints);
    idx tetNodes[] = {0, 1, 2, 3};
    tetrahedra.setAllNodes(&tetNodes[0]);
    tetrahedra.setnNodes(numPoints);

    // ACT: write the result to a temporary file
    auto tempFile = TestUtil::TemporaryFile::empty();
    UnstructuredGridWriter writer;
    writer.write(tempFile.name(),
                 numPoints,
                 numLcPoints,
                 &points[0],
                 tetrahedra,
                 &potentials[0],
                 &director[0],
                 &director[4],
                 &director[8],
                 &S[0]);

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
    REQUIRE(lines[directorStartLine + 3] == "-1 0 0");
    REQUIRE(lines[directorStartLine + 4] == "0 -1 0");

    // check order parameter data section
    int orderParamStartLine = 31;
    REQUIRE(lines[orderParamStartLine] == "SCALARS S float 1");
    REQUIRE(lines[orderParamStartLine + 1] == "LOOKUP_TABLE default");
    REQUIRE(lines[orderParamStartLine + 2] == "0.1");
    REQUIRE(lines[orderParamStartLine + 3] == "0.2");
    REQUIRE(lines[orderParamStartLine + 4] == "0.3");
    REQUIRE(lines[orderParamStartLine + 5] == "0.4");
}