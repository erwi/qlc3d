//
// Created by eero on 07/04/2021.
//
#include <catch.h>
#include <resultoutput.h>
#include <qlc3d.h>
TEST_CASE("Read resource mesh file") {
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
    ReadGiDMesh3D("./tests/resources/thin.msh", &p, &np, &t, &nt, &e, &ne, &matt, &mate);

    REQUIRE(np == 2092);
    REQUIRE(nt == 3465);
    REQUIRE(ne == 4004);
}
