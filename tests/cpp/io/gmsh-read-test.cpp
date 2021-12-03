#include <catch.h>
#include <gmsh-read.h>

TEST_CASE("throw exception when") {
    GmshFileReader reader;
    REQUIRE_THROWS_AS(reader.readGmsh("path/to/missing/file"), std::invalid_argument);
}

TEST_CASE("read gmsh mesh") {
    using namespace std;
    GmshFileReader reader;
    //string fileName = "/home/eero/projects/qlc3d-tests/gmsh.msh";
    string fileName = "/home/eero/projects/qlc3d-tests/untitled.msh";
    auto meshData = reader.readGmsh(fileName);

    auto meshFormat = meshData->getMeshFormat();
    REQUIRE(4.1 == meshFormat->_versionNumber);
    REQUIRE(0 == meshFormat->_fileType); // ASCII
    REQUIRE(8 == meshFormat->_dataSize); // sizeof <floating point number>

    auto materials = meshData->getPhysicalNames();
    REQUIRE(materials != nullptr);

    auto entities = meshData->getEntities();
    REQUIRE(entities != nullptr);
}
