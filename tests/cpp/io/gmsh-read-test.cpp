#include <catch.h>
#include <test-util.h>
#include <io/gmsh-read.h>
#include <iostream>

TEST_CASE("throw exception when") {
    GmshFileReader reader;
    REQUIRE_THROWS_AS(reader.readGmsh("path/to/missing/file"), std::runtime_error);
}

TEST_CASE("read Gmsh mesh raw data") {
    using namespace std;
    GmshFileReader reader;
    auto meshData = reader.readGmsh(TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH);

    auto meshFormat = meshData->getMeshFormat();
    REQUIRE(4.1 == meshFormat->_versionNumber);
    REQUIRE(0 == meshFormat->_fileType); // ASCII
    REQUIRE(8 == meshFormat->_dataSize); // sizeof <floating point number>

    auto materials = meshData->getPhysicalNames();
    REQUIRE(materials != nullptr);

    auto entities = meshData->getEntities();
    REQUIRE(entities != nullptr);

    auto elements = meshData->getElements();
    REQUIRE(elements != nullptr);
    REQUIRE(elements->_numTriangles == elements->_triangleIndices.size() / 3);
    REQUIRE(elements->_numTetrahedra == elements->_tetrahedraIndices.size() / 4);

    REQUIRE(elements->_triangleEntityTags.size() == elements->_numTriangles);
    REQUIRE(elements->_triangleEntityTags.size() == elements->_numTetrahedra);

    auto nodes = meshData->getNodes();
    SECTION("coordinate values are set") {
        REQUIRE(14 == nodes->_numNodes);
        REQUIRE(3 * 14 == nodes->_coordinates.size());
        size_t numNaNs = std::count(nodes->_coordinates.begin(), nodes->_coordinates.end(), std::numeric_limits<double>::quiet_NaN());
        REQUIRE(0 == numNaNs);
    }

    SECTION("tet mesh should have 0-based node indexing") {
        bool numZeros = std::count(elements->_tetrahedraIndices.begin(), elements->_tetrahedraIndices.end(), 0);
        bool numMaxIndex = std::count(elements->_tetrahedraIndices.begin(), elements->_tetrahedraIndices.end(), nodes->_numNodes - 1);
        REQUIRE(0 < numZeros);
        REQUIRE(0 < numMaxIndex);
    }

    GmshPhysicalNamesMapper mapper(materials, entities, elements);
    SECTION("map Gmsh physical names to qlc3d surface material numbers") {
        vector<unsigned int> surfaceMaterialNumbers = mapper.mapTriangleNamesToMaterialNumbers();
        REQUIRE(surfaceMaterialNumbers.size() == elements->_numTriangles);
    }

    SECTION("map Gmsh physical names qlc3d volume material numbers") {
        vector<unsigned int> volumeMaterialNumbers = mapper.maptTetrahedraNamesToMaterialNumbers();
        REQUIRE(volumeMaterialNumbers.size() == elements->_numTetrahedra);
    }
}
