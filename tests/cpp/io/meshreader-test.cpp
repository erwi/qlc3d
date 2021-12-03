#include <iostream>
#include <catch.h>
#include <io/meshreader.h>
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
