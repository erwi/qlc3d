#include <catch.h>
#include <test-util.h>
#include <io/gmsh-read.h>
#include <iostream>
#include <algorithm>
#include <set>
#include "material_numbers.h"

TEST_CASE("throw exception when") {
    GmshFileReader reader;
    REQUIRE_THROWS_AS(reader.readGmsh("path/to/missing/file"), std::runtime_error);
}

TEST_CASE("read Gmsh mesh raw data - linear elements") {
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

    REQUIRE(meshData->getElements()->_elementOrder == 1); // Linear elements

    auto elements = meshData->getElements();
    REQUIRE(elements != nullptr);
    REQUIRE(elements->_numTriangles == elements->_triangleIndices.size() / 3);
    REQUIRE(elements->_numTetrahedra == elements->_tetrahedraIndices.size() / 4);

    REQUIRE(elements->_triangleEntityTags.size() == elements->_numTriangles);
    REQUIRE(elements->_triangleEntityTags.size() == elements->_numTetrahedra);

    auto nodes = meshData->getNodes();
    SECTION("coordinate values are set") {
        REQUIRE(14 == nodes->_numNodes);
        REQUIRE(14 == nodes->_coordinates.size());
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

TEST_CASE("read Gmsh mesh raw data - quadratic elements") {
  using namespace std;
  GmshFileReader reader;
  auto meshData = reader.readGmsh(TestUtil::RESOURCE_SMALL_CUBE_QUADRATIC_GMSH_MESH);

  auto meshFormat = meshData->getMeshFormat();
  REQUIRE(4.1 == meshFormat->_versionNumber);
  REQUIRE(0 == meshFormat->_fileType); // ASCII
  REQUIRE(8 == meshFormat->_dataSize); // sizeof <floating point number>

  auto materials = meshData->getPhysicalNames();
  REQUIRE(materials != nullptr);
  REQUIRE(materials->_numNames == 6);

  auto names = materials->_physicalNames;
  REQUIRE(names.size() == 6);
  REQUIRE(names[1] == "fixlc1");
  REQUIRE(names[2] == "electrode1");
  REQUIRE(names[3] == "fixlc2");
  REQUIRE(names[4] == "electrode2");
  REQUIRE(names[5] == "periodic");
  REQUIRE(names[7] == "domain1");

  REQUIRE(meshData->getElements()->_elementOrder == 2); // Quadratic elements

  auto entities = meshData->getEntities();
  REQUIRE(entities != nullptr);

  auto elements = meshData->getElements();
  REQUIRE(elements != nullptr);
  REQUIRE(elements->_numTriangles == 24);
  REQUIRE(elements->_numTetrahedra == 24);
  REQUIRE(elements->_triangleIndices.size() == 24 * 6);
  REQUIRE(elements->_tetrahedraIndices.size() == 24 * 10);

  REQUIRE(elements->_triangleEntityTags.size() == elements->_numTriangles);
  REQUIRE(elements->_triangleEntityTags.size() == elements->_numTetrahedra);

  auto nodes = meshData->getNodes();
  SECTION("coordinate values are set") {
    REQUIRE(63 == nodes->_numNodes);
    REQUIRE(63 == nodes->_coordinates.size());
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



TEST_CASE("Read Gmsh mesh with dielectric and LC volumes") {
  GmshFileReader reader;
  auto fileData = reader.readGmsh(TestUtil::RESOURCE_UNIT_CUBE_DIELECTRIC_NEUMAN_GMSH_MESH);

  GmshPhysicalNamesMapper mapper(fileData->getPhysicalNames(), fileData->getEntities(), fileData->getElements());

  auto surfaceMaterialNumbers = mapper.mapTriangleNamesToMaterialNumbers();
  auto volumeMaterialNumbers = mapper.maptTetrahedraNamesToMaterialNumbers();

  std::set<unsigned int> surfaceMaterialsSet(surfaceMaterialNumbers.begin(), surfaceMaterialNumbers.end());
  std::set<unsigned int> volumeMaterialsSet(volumeMaterialNumbers.begin(), volumeMaterialNumbers.end());

  REQUIRE(4 == surfaceMaterialsSet.size());
  REQUIRE(1 == surfaceMaterialsSet.count(MAT_NEUMANN)); // side boundaries
  REQUIRE(1 == surfaceMaterialsSet.count(MAT_ELECTRODE2)); // bottom boundary
  REQUIRE(1 == surfaceMaterialsSet.count(MAT_FIXLC2)); // middle boundary between LC and dielectric region
  REQUIRE(1 == surfaceMaterialsSet.count(MAT_ELECTRODE1 | MAT_FIXLC1)); // top boundary

  REQUIRE(2 == volumeMaterialsSet.size());
  REQUIRE(1 == volumeMaterialsSet.count(MAT_DOMAIN1));
  REQUIRE(1 == volumeMaterialsSet.count(MAT_DIELECTRIC1));
}