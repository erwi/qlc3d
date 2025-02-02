#include <catch.h>
#include <inits.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <electrodes.h>
#include <alignment.h>
#include <test-util.h>


TEST_CASE("Reorder quadratic element node order") {
  // GMSH quadratic tet element node ordering is different from the usual used in literature.
  // The difference is that the last two nodes are swapped. When reading the GMSH mesh and preparing
  // the geometry, we should swap the last two nodes if they are in the GMSH order.

  Geometry geom;
  Electrodes electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));

  prepareGeometry(geom,
                  TestUtil::RESOURCE_SMALL_CUBE_QUADRATIC_GMSH_MESH,
                  electrodes,
                  alignment);

  Mesh &tets = geom.getTetrahedra();
  const Coordinates &coords = geom.getCoordinates();
  REQUIRE(tets.getnElements() == 24);

  // check node order in each tetrahedra
  for (unsigned int i = 0; i < tets.getnElements(); i++) {
    unsigned int tetNodes[10] = {};
    Vec3 coordNodes[10] = {};

    tets.loadNodes(i, tetNodes);
    coords.loadCoordinates(tetNodes, tetNodes + 10, coordNodes);

    // assume that the first 4 nodes are the same as in linear tetrahedra, i.e. corner nodes
    // Check the mid-point nodes positions are as expected w.r.t. to the corner nodes.
    // Not GMSH node ordering.

    // Node 4 should be the mid-point of edge 0-1
    Vec3 expectedPos = (coordNodes[0] + coordNodes[1]) * 0.5;
    REQUIRE(coordNodes[4].equals(expectedPos, 1e-9));

    // Node 5 should be the mid-point of edge 1-2
    expectedPos = (coordNodes[1] + coordNodes[2]) * 0.5;
    REQUIRE(coordNodes[5].equals(expectedPos, 1e-9));

    // Node 6 should be the mid-point of edge 0-2
    expectedPos = (coordNodes[2] + coordNodes[0]) * 0.5;
    REQUIRE(coordNodes[6].equals(expectedPos, 1e-9));

    // Node 7 should be the mid-point of edge 0-3
    expectedPos = (coordNodes[0] + coordNodes[3]) * 0.5;
    REQUIRE(coordNodes[7].equals(expectedPos, 1e-9));

    // Node 8 should be the mid-point of edge 1-3
    expectedPos = (coordNodes[1] + coordNodes[3]) * 0.5;
    REQUIRE(coordNodes[8].equals(expectedPos, 1e-9));

    // Node 9 should be the mid-point of edge 2-3
    expectedPos = (coordNodes[2] + coordNodes[3]) * 0.5;
    REQUIRE(coordNodes[9].equals(expectedPos, 1e-9));
  }
}