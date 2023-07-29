
#include <catch.h>
#include <io/meshreader.h>
#include <test-util.h>
#include <geom/vec3.h>

TEST_CASE("Read GiD mesh file") {
    std::vector<Vec3> points;
    std::vector<idx> tetNodes;
    std::vector<idx> tetMaterials;
    std::vector<idx> triNodes;
    std::vector<idx> triMaterials;

    // TODO: this fails when running target "All CTest", the mesh is not found
    // Create a ResourceFile class that manages the resource file paths?
    ReadGiDMesh3D(TestUtil::RESOURCE_THIN_GID_MESH, points, tetNodes, tetMaterials, triNodes, triMaterials);

    SECTION("Correct number of vertices and elements should have been read") {
        REQUIRE(points.size() == 2092);
        REQUIRE(tetNodes.size() == 4 * 3465);
        REQUIRE(tetMaterials.size() == 3465);
        REQUIRE(triMaterials.size() == 4004);
        REQUIRE(triNodes.size() == 3 * 4004);
    }

    // verify some hand-picked values from the file
    SECTION("Point coordinate values should be as expected") {
        REQUIRE(points[0].x() == 0.002);
        REQUIRE(points[0].y() == 0);
        REQUIRE(points[0].z() == 1);
    }

    SECTION("Tetrahedral element indices should be 0-based") {
        idx minIndex = points.size();
        idx nt = tetMaterials.size();
        for (int i = 0; i < nt; i++) {
            minIndex = std::min(minIndex, tetNodes[i]);
        }
        REQUIRE(minIndex == 0);
    }

    SECTION("Triangle element indices should be 0-based") {
      idx minIndex = *std::min_element(triNodes.begin(), triNodes.end());
        REQUIRE(minIndex == 0);
    }

    // TODO: check element material numbers for a few elements
}