//
// Created by eero on 07/04/2021.
//
#include <catch.h>
#include <resultio.h>
#include <qlc3d.h>
#include <test-util.h>

TEST_CASE("Read GiD mesh file") {
    double *p;
    idx np;
    idx *t;
    idx nt;
    idx *e;
    idx ne;
    idx *matt;
    idx *mate;

    // TODO: this fails when running target "All CTest", the mesh is not found
    // Create a ResourceFile class that manages the resource file paths?
    ReadGiDMesh3D(TestUtil::RESOURCE_THIN_GID_MESH, &p, &np, &t, &nt, &e, &ne, &matt, &mate);

    SECTION("Correct number of vertices and elements should have been read") {
        REQUIRE(np == 2092);
        REQUIRE(nt == 3465);
        REQUIRE(ne == 4004);
    }

    // verify some hand-picked values from the file
    SECTION("Point coordinate values should be as expected") {
        REQUIRE(p[0] == 0.002);
        REQUIRE(p[1] == 0);
        REQUIRE(p[2] == 1);
    }

    SECTION("Tetrahedral element indices should be 0-based") {
        idx minIndex = np;
        for (int i = 0; i < nt * 4; i++) {
            minIndex = std::min(minIndex, t[i]);
        }
        REQUIRE(minIndex == 0);
    }

    SECTION("Triangle element indices should be 0-based") {
        idx minIndex = ne;
        for (int i = 0; i < ne * 3; i++) {
            minIndex = std::min(minIndex, e[i]);
        }
        REQUIRE(minIndex == 0);
    }

    // TODO: check element material numbers for a few elements
}