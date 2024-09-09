#include <catch.h>
#include <qlc3d.h>
#include <lc-representation.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <inits.h>
#include <test-util.h>
#include <memory>

const double MARGIN = 1e-12;

TEST_CASE("Set initial LC orientation") {
  double S0 = 0.55;
  Boxes boxes;

  SECTION("Set director to uniform 90 degrees inside Normal box") {
    SolutionVector q(2, 5); // 2 nodes

    Coordinates coordinates({
                                    {0, 0, 0}, // inside box
                                    {1, 0, 0}  // outside box
                            });

    // box with 90-degree tilt
    boxes.addBox(1, "Normal", {}, {0, 0.5}, {0, 1}, {0, 1}, {90, 0}, {0, 0});

    setVolumeQ(q, S0, boxes, coordinates);

    qlc3d::Director d1 = q.getDirector(0);
    qlc3d::Director d2 = q.getDirector(1);

    // d1 should be {0, 0, 1}
    REQUIRE(d1.nx() == Approx(0).margin(MARGIN));
    REQUIRE(d1.ny() == Approx(0).margin(MARGIN));
    REQUIRE(d1.nz() == Approx(1).margin(MARGIN));
    REQUIRE(d1.S() == Approx(S0).margin(MARGIN));

    // d2 should be {the default 1, 0, 0}
    REQUIRE(d2.nx() == Approx(1).margin(MARGIN));
    REQUIRE(d2.ny() == Approx(0).margin(MARGIN));
    REQUIRE(d2.nz() == Approx(0).margin(MARGIN));
    REQUIRE(d2.S() == Approx(S0).margin(MARGIN));
  }

  SECTION("Set director to linearly increasing 0-90 degree tilt, twist inside Normal box") {
    SolutionVector q(3, 5); // 3 nodes: bottom, mid, top

    Coordinates coordinates({{0, 0, 0}, {0, 0, 0.5}, {0, 0, 1}});

    boxes.addBox(1, "Normal", {},
                 {0, 1},
                 {0, 1},
                 {0, 1},
                 {0, 90}, // tilt increasing from 0 to 90
                 {0, 90});  // twist increasing from 0 to 90
    setVolumeQ(q, S0, boxes, coordinates);

    qlc3d::Director bottom = q.getDirector(0);
    qlc3d::Director mid = q.getDirector(1);
    qlc3d::Director top = q.getDirector(2);

    // bottom should be {1, 0, 0}
    REQUIRE(bottom.S() == Approx(S0).margin(MARGIN));
    REQUIRE(bottom.nx() == Approx(1).margin(MARGIN));
    REQUIRE(bottom.ny() == Approx(0).margin(MARGIN));
    REQUIRE(bottom.nz() == Approx(0).margin(MARGIN));

    // mid should have 45-degree tilt,twist
    REQUIRE(mid.S() == Approx(S0).margin(MARGIN));
    REQUIRE(mid.tiltDegrees() == Approx(45).margin(MARGIN));
    REQUIRE(mid.twistDegrees() == Approx(45).margin(MARGIN));

    // top should have 90-degree tilt,twist
    REQUIRE(top.S() == Approx(S0).margin(MARGIN));
    REQUIRE(top.nz() == Approx(1).margin(MARGIN));
    REQUIRE(top.tiltDegrees() == Approx(90).margin(MARGIN));
  }
}

TEST_CASE("Initial LC surface orientations") {
  Geometry geom;
  auto simu = std::unique_ptr<Simu>(SimuBuilder().build());
  auto lc = std::unique_ptr<LC>(LCBuilder().build());
  Boxes boxes; // volume orientations - leave empty, don't care in this test

  std::vector<std::shared_ptr<Electrode>> el;
  el.push_back(std::make_shared<Electrode>(1, std::vector<double>(), std::vector<double>()));
  el.push_back(std::make_shared<Electrode>(2, std::vector<double>(), std::vector<double>()));
  Electrodes electrodes(el);

  Alignment alignment;
  std::vector<idx> surfaceNodesIndex;
  geom.getTriangles().listFixLCSurfaces(surfaceNodesIndex, 1);

    SECTION("Weak homeotropic anchoring") {
      // expected to fail until weak anchoring is re-enabled
    alignment.addSurface(Surface::ofPlanarDegenerate(1, -1e-3));
    alignment.addSurface(Surface::ofPlanarDegenerate(2, -1e-3));
    prepareGeometry(geom, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, electrodes, alignment, {1, 1, 1});
    SolutionVector q(geom.getnpLC(), 5);

    // ACT
    initialiseLcSolutionVector(q, *simu, *lc, boxes, alignment, geom);

    // ASSERT
    // director should be parallel to surface normal
    for (idx i: surfaceNodesIndex) {
      auto director = q.getDirector(i);
      Vec3 d = director.vector();
      Vec3 surfaceNormal = geom.getNodeNormal(i);

      double dot = surfaceNormal.dot(d);
      REQUIRE(std::abs(dot) == Approx(1).margin(MARGIN));
      REQUIRE(director.S() == Approx(lc->S0()).margin(MARGIN));
    }
  }

  SECTION("Strong homeotropic anchoring") {
    alignment.addSurface(Surface::ofHomeotropic(1));
    alignment.addSurface(Surface::ofHomeotropic(2));
    prepareGeometry(geom, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, electrodes, alignment, {1, 1, 1});
    SolutionVector q(geom.getnpLC(), 5);
    // ACT
    initialiseLcSolutionVector(q, *simu, *lc, boxes, alignment, geom);

    // ASSERT
    // director should be parallel to surface normal
    for (idx i: surfaceNodesIndex) {
      auto director = q.getDirector(i);
      Vec3 d = director.vector();
      Vec3 surfaceNormal = geom.getNodeNormal(i);

      double dot = surfaceNormal.dot(d);

      REQUIRE(std::abs(dot) == Approx(1).margin(MARGIN));
      REQUIRE(director.S() == Approx(lc->S0()).margin(MARGIN));
    }
  }

  SECTION("Weak with pre-tilt and pre-twist") {
    // expected to fail until weak anchoring is re-enabled
    double tiltDegrees = 5;
    double twistDegrees = 12;
    alignment.addSurface(1, "Weak", 1e-3, {tiltDegrees, twistDegrees, 0}, 1., 1., {});
    alignment.addSurface(2, "Weak", 1e-3, {tiltDegrees, twistDegrees, 0}, 1., 1., {});
    prepareGeometry(geom, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, electrodes, alignment, {1, 1, 1});
    SolutionVector q(geom.getnpLC(), 5);
    // ACT
    initialiseLcSolutionVector(q, *simu, *lc, boxes, alignment, geom);

    // ASSERT
    // director should be parallel to expected value
    auto expectedDirector = qlc3d::Director::fromDegreeAngles(tiltDegrees, twistDegrees, lc->S0()).vector();
    for (idx i: surfaceNodesIndex) {
      auto director = q.getDirector(i);

      double dot = director.vector().dot(expectedDirector);

      REQUIRE(std::abs(dot) == Approx(1).margin(MARGIN));
      REQUIRE(director.S() == Approx(lc->S0()).margin(MARGIN));
    }
  }

  SECTION("Dont set surface orientation when enforce is false") {
    // expected to fail until weak anchoring is re-enabled

    // ARRANGE
    // Set up surface with 90-degree tilt, 0 degree twist and enforce flag set to false
    // This means that the surface should not override the volume LC orientation and the LC director should be
    // equal to (1, 0, 0) everywhere
    double tiltDegrees = 90;
    double twistDegrees = 0;
    alignment.addSurface(1, "Weak", 1e-3, {tiltDegrees, twistDegrees, 0}, 1., 1., {}, false);
    alignment.addSurface(2, "Weak", 1e-3, {tiltDegrees, twistDegrees, 0}, 1., 1., {}, false);
    prepareGeometry(geom, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, electrodes, alignment, {1, 1, 1});
    SolutionVector q(geom.getnpLC(), 5);
    // ACT
    initialiseLcSolutionVector(q, *simu, *lc, boxes, alignment, geom);

    // ASSERT
    // director should be parallel to surface normal
    Vec3 expectedDirector = {1, 0, 0};
    for (idx i: surfaceNodesIndex) {
      auto director = q.getDirector(i);

      double dot = director.vector().dot(expectedDirector);

      REQUIRE(std::abs(dot) == Approx(1).margin(MARGIN));
      REQUIRE(director.S() == Approx(lc->S0()).margin(MARGIN));
    }
  }
}
