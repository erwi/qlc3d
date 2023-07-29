#include <algorithm>
#include <catch.h>
#include <io/meshreader.h>
#include <material_numbers.h>
#include <test-util.h>

using namespace std;
TEST_CASE("Detect input mesh file format") {
    SECTION("Detect GiD mesh format") {
        MeshFormat format = MeshFormatInspector::inspectFormat(TestUtil::RESOURCE_THIN_GID_MESH);
        REQUIRE(MeshFormat::GID_MESH == format);
    }

    SECTION("Detect ASCII Gmsh format") {
        MeshFormat format = MeshFormatInspector::inspectFormat(TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH);
        REQUIRE(MeshFormat::GMSH_ASCII == format);
    }
}

TEST_CASE("Read Gmsh mesh without knowing format") {
    auto meshData = MeshReader::readMesh(TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH);

    unsigned int numPoints = meshData.points.size();

    SECTION("mesh should have expected number of points, tets, and tris") {
        REQUIRE(14 == numPoints);
        REQUIRE(24 == meshData.tetMaterials.size()); // 24 tets
        REQUIRE(24 * 4 == meshData.tetNodes.size()); // 4 nodes per tet
        REQUIRE(24 == meshData.triMaterials.size()); // 24 tris
        REQUIRE(24 * 3 == meshData.triNodes.size()); // 3 nodes per tri
    }

    SECTION("all tets should be LC material") {
        bool areAllDomain1 = std::all_of(meshData.tetMaterials.begin(), meshData.tetMaterials.end(),
                                         [](unsigned int v){ return v == MAT_DOMAIN1;});
        REQUIRE(areAllDomain1);
    }

    SECTION("mesh should contain 16 periodic material triangles") {
        size_t numPeriodic = std::count_if(meshData.triMaterials.begin(), meshData.triMaterials.end(),
                      [](idx v){ return v == MAT_PERIODIC; });
        REQUIRE(16 == numPeriodic);
    }

    SECTION("mesh should contain 4 FixLC1 Electrode1 triangles") {
        idx material = MAT_FIXLC1 + MAT_ELECTRODE1;
        size_t numMaterialTriangles = std::count_if(meshData.triMaterials.begin(), meshData.triMaterials.end(),
                                                    [material](idx v) { return v == material; });
        REQUIRE(4 == numMaterialTriangles);
    }

    SECTION("mesh should contain 4 FixLC2 Electrode2 triangles") {
        idx material = MAT_FIXLC2 + MAT_ELECTRODE2;
        size_t numMaterialTriangles = std::count_if(meshData.triMaterials.begin(), meshData.triMaterials.end(),
                                                    [material](idx v) { return v == material; });
        REQUIRE(4 == numMaterialTriangles);
    }
}