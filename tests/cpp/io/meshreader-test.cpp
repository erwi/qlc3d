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
    double *p;
    idx numPoints, numTets, numTris;
    idx *tets = nullptr, *tris = nullptr, *tetMaterials = nullptr, *triMaterials = nullptr;

    MeshReader::readMesh(TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH,
                         &p, &numPoints, &tets, &numTets, &tris, &numTris, &tetMaterials, &triMaterials);
    SECTION("mesh should have expected number of points, tets, and tris") {
        REQUIRE(14 == numPoints);
        REQUIRE(24 == numTets);
        REQUIRE(24 == numTris);
    }

    SECTION("all tets should be LC material") {
        bool areAllDomain1 = std::all_of(tetMaterials, tetMaterials + numTets,
                                         [](unsigned int v){ return v == MAT_DOMAIN1;});
        REQUIRE(areAllDomain1);
    }

    SECTION("mesh should contain 16 periodic material triangles") {
        size_t numPeriodic = std::count_if(triMaterials, triMaterials + numTris,
                      [](idx v){ return v == MAT_PERIODIC; });
        REQUIRE(16 == numPeriodic);
    }

    SECTION("mesh should contain 4 FixLC1 Electrode1 triangles") {
        idx material = MAT_FIXLC1 + MAT_ELECTRODE1;
        size_t numMaterialTriangles = std::count_if(triMaterials, triMaterials + numTris,
                                                    [material](idx v) { return v == material; });
        REQUIRE(4 == numMaterialTriangles);
    }

    SECTION("mesh should contain 4 FixLC2 Electrode2 triangles") {
        idx material = MAT_FIXLC2 + MAT_ELECTRODE2;
        size_t numMaterialTriangles = std::count_if(triMaterials, triMaterials + numTris,
                                                    [material](idx v) { return v == material; });
        REQUIRE(4 == numMaterialTriangles);
    }

    // TODO: noo....
    free((void*) p);
    free((void*) tets);
    free((void*) tris);
    free((void*) tetMaterials);
    free((void*) triMaterials);
}