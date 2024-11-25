#include <catch.h>
#include <qlc3d.h>
#include <lc-representation.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <inits.h>
#include <test-util.h>
#include <memory>
#include <alignment.h>

const double MARGIN = 1e-12;

TEST_CASE("Set initial LC orientation") {
  double S0 = 0.55;
  InitialVolumeOrientation boxes;

  SECTION("Set director to uniform 90 degrees inside Normal box") {
    SolutionVector q(2, 5); // 2 nodes

    Coordinates coordinates({
                                    {0, 0, 0}, // inside box
                                    {1, 0, 0}  // outside box
                            });

    // box with 90-degree tilt
    boxes.addBox(1, "Normal", {}, {0, 0.5}, {0, 1}, {0, 1}, {90, 0}, {0, 0});

    // ACT
    boxes.setVolumeQ(q, S0, coordinates);

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

    // ACT
    boxes.setVolumeQ(q, S0, coordinates);

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

  SECTION("Set director in Hedgehog box") {
    SolutionVector q(6, 5); // 6 nodes on each side of box centre

    Coordinates coordinates({
      {-1, 0, 0},
      {1, 0, 0},
      {0, -1, 0},
      {0, 1, 0},
      {0, 0, -1},
      {0, 0, 1}
      });

    boxes.addBox(1, "Hedgehog",{}, {-2, 2}, {-2, 2}, {-2, 2}, {0, 0}, {0, 0});

    // ACT
    boxes.setVolumeQ(q, S0, coordinates);

    // ASSERT: director should be parallel to vector from box centre to node
    auto dirLeft = q.getDirector(0);
    auto dirRight = q.getDirector(1);
    auto dirFront = q.getDirector(2);
    auto dirBack = q.getDirector(3);
    auto dirBottom = q.getDirector(4);
    auto dirTop = q.getDirector(5);

    REQUIRE((dirLeft.vector().equals({-1, 0, 0}, MARGIN) || dirLeft.vector().equals({1, 0, 0}, MARGIN)));
    REQUIRE((dirRight.vector().equals({-1, 0, 0}, MARGIN) || dirRight.vector().equals({1, 0, 0}, MARGIN)));

    REQUIRE((dirFront.vector().equals({0, -1, 0}, MARGIN) || dirFront.vector().equals({0, 1, 0}, MARGIN)));
    REQUIRE((dirBack.vector().equals({0, -1, 0}, MARGIN) || dirBack.vector().equals({0, 1, 0}, MARGIN)));

    REQUIRE((dirBottom.vector().equals({0, 0, -1}, MARGIN) || dirBottom.vector().equals({0, 0, 1}, MARGIN)));
    REQUIRE((dirTop.vector().equals({0, 0, -1}, MARGIN) || dirTop.vector().equals({0, 0, 1}, MARGIN)));
  }
}

TEST_CASE("Box tilt and twist expressions") {
  BoxBuilder bb(1);
  bb.setX({0, 1})
          .setY({0, 1})
          .setZ({0, 1})
          .setTiltExpression("45 * Z")  // Tilt increases linearly with Z
          .setTwistExpression("90 * X"); // Twist increases linearly with X

  auto box = bb.build();

  // Check director at different points
  auto d1 = box->getDirectorAt(Vec3(0, 0, 0));
  REQUIRE(d1.tiltDegrees() == Approx(0).margin(1e-12));
  REQUIRE(d1.twistDegrees() == Approx(0).margin(1e-12));

  auto d2 = box->getDirectorAt(Vec3(0, 0, 1));
  REQUIRE(d2.tiltDegrees() == Approx(45).margin(1e-12));
  REQUIRE(d2.twistDegrees() == Approx(0).margin(1e-12));

  auto d3 = box->getDirectorAt(Vec3(1, 0, 0));
  REQUIRE(d3.tiltDegrees() == Approx(0).margin(1e-12));
  REQUIRE(d3.twistDegrees() == Approx(90).margin(1e-12));

  auto d4 = box->getDirectorAt(Vec3(0.5, 0, 0.5));
  REQUIRE(d4.tiltDegrees() == Approx(22.5).margin(1e-12));
  REQUIRE(d4.twistDegrees() == Approx(45).margin(1e-12));
}

TEST_CASE("Initial LC surface orientations") {
  Geometry geom;
  auto simu = std::unique_ptr<Simu>(SimuBuilder().build());
  auto lc = std::unique_ptr<LC>(LCBuilder().build());
  InitialVolumeOrientation boxes; // volume orientations - leave empty, don't care in this test

  std::vector<std::shared_ptr<Electrode>> el;
  el.push_back(std::make_shared<Electrode>(1, std::vector<double>(), std::vector<double>()));
  el.push_back(std::make_shared<Electrode>(2, std::vector<double>(), std::vector<double>()));
  Electrodes electrodes(el);

  Alignment alignment;
  std::vector<idx> surfaceNodesIndex;
  geom.getTriangles().listFixLCSurfaces(surfaceNodesIndex, 1);

  SECTION("Weak homeotropic anchoring") {
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
    alignment.addSurface(Surface::ofStrongHomeotropic(1));
    alignment.addSurface(Surface::ofStrongHomeotropic(2));
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

  SECTION("Strong anchoring but overrideVolumes is false") {
    alignment.addSurface(Surface::ofStrongAnchoring(1, 90, 0, false));
    alignment.addSurface(Surface::ofStrongAnchoring(2, 90, 0, false));
    prepareGeometry(geom, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, electrodes, alignment, {1, 1, 1});
    SolutionVector q(geom.getnpLC(), 5);
    // ACT
    initialiseLcSolutionVector(q, *simu, *lc, boxes, alignment, geom);

    // ASSERT
    // director should be [1, 0, 0] even though anchoring easy direction is defined as [0, 0, 1]
    for (idx i: surfaceNodesIndex) {
      auto director = q.getDirector(i);
      Vec3 d = director.vector();
      double dot = Vec3(1, 0, 0).dot(d);

      REQUIRE(std::abs(dot) == Approx(1).margin(MARGIN));
      REQUIRE(director.S() == Approx(lc->S0()).margin(MARGIN));
    }
  }

  SECTION("Weak with pre-tilt and pre-twist") {
    // expected to fail until weak anchoring is re-enabled
    double tiltDegrees = 5;
    double twistDegrees = 12;
    alignment.addSurface(Surface::ofWeakAnchoring(1, tiltDegrees, twistDegrees, 1e-3, 1, 1));
    alignment.addSurface(Surface::ofWeakAnchoring(2, tiltDegrees, twistDegrees, 1e-3, 1, 1));
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
    alignment.addSurface(Surface::ofWeakAnchoring(1, tiltDegrees, twistDegrees, 1e-3, 1, 1));
    alignment.addSurface(Surface::ofWeakAnchoring(2, tiltDegrees, twistDegrees, 1e-3, 1, 1));
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

TEST_CASE("Surface strong anchoring with analytic expressions for tilt and twist") {
  // Create a Surface with strong anchoring and analytic expressions for tilt and twist
  Surface s = Surface::ofStrongAnchoring(1, "45 * Z", "90 * X");

  // Check director at different points
  auto d1 = s.getEasyDirectionAt(Vec3(0, 0, 0));
  REQUIRE(d1.equals({1, 0, 0}, 1e-15));

  auto d2 = s.getEasyDirectionAt(Vec3(0, 0, 2));
  REQUIRE(d2.equals({0, 0, 1}, 1e-15));

  auto d3 = s.getEasyDirectionAt(Vec3(1, 0, 0));
  REQUIRE(d3.equals({0, 1, 0}, 1e-15));
}

TEST_CASE("Surface easy direction calculation from angles should match director definition from angles") {
  double tiltDegrees = 45;
  double twistDegrees = 45;

  auto n = qlc3d::Director::fromDegreeAngles(tiltDegrees, twistDegrees, 1);

  Surface s = Surface::ofStrongAnchoring(1, tiltDegrees, twistDegrees);

  auto easy = s.getEasyVector();

  REQUIRE(easy.equals(n.vector(), 1e-15));
  REQUIRE(easy.equals(s.getV1().cross(s.getV2()), 1e-15));
}

TEST_CASE("Surface easy direction vectors for 0 tilt 0 twist") {
  double tiltDegrees = 0;
  double twistDegrees = 0;

  Surface s = Surface::ofStrongAnchoring(1, tiltDegrees, twistDegrees);

  auto easy = s.getEasyVector();

  // we require easy axis is v1 x v2
  REQUIRE(easy.equals({1, 0, 0}, 1e-15));
  REQUIRE(s.getV1().equals({0, 1, 0}, 1e-15));
  REQUIRE(s.getV2().equals({0, 0, 1}, 1e-15));
}