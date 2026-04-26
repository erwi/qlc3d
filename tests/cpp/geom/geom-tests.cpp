#include <catch.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <geometry.h>
#include <globals.h>
#include <material_numbers.h>
#include <solutionvector.h>
#include <mesh/mesh.h>

#include <vector>
#include <memory>

// ============================================================================
// Helper: a single straight-edged TET10 element.
//
// Corner nodes:
//   n0 = (0, 0, 0),  n1 = (1, 0, 0),  n2 = (0, 1, 0),  n3 = (0, 0, 1)
//
// Mid-edge nodes (at exact edge midpoints):
//   n4 = mid(n0,n1) = (0.5, 0,   0  )
//   n5 = mid(n1,n2) = (0.5, 0.5, 0  )
//   n6 = mid(n0,n2) = (0,   0.5, 0  )
//   n7 = mid(n0,n3) = (0,   0,   0.5)
//   n8 = mid(n1,n3) = (0.5, 0,   0.5)
//   n9 = mid(n2,n3) = (0,   0.5, 0.5)
// ============================================================================

struct Tet10Fixture {
  std::shared_ptr<Coordinates> coords;
  std::shared_ptr<Mesh> mesh;

  Tet10Fixture() {
    std::vector<Vec3> points = {
        Vec3(0.0, 0.0, 0.0),   // n0
        Vec3(1.0, 0.0, 0.0),   // n1
        Vec3(0.0, 1.0, 0.0),   // n2
        Vec3(0.0, 0.0, 1.0),   // n3
        Vec3(0.5, 0.0, 0.0),   // n4 = mid(n0,n1)
        Vec3(0.5, 0.5, 0.0),   // n5 = mid(n1,n2)
        Vec3(0.0, 0.5, 0.0),   // n6 = mid(n0,n2)
        Vec3(0.0, 0.0, 0.5),   // n7 = mid(n0,n3)
        Vec3(0.5, 0.0, 0.5),   // n8 = mid(n1,n3)
        Vec3(0.0, 0.5, 0.5),   // n9 = mid(n2,n3)
    };
    coords = std::make_shared<Coordinates>(std::move(points));
    mesh = std::make_shared<Mesh>(3, ElementType::QUADRATIC_TETRAHEDRON);
    mesh->setElementData(ElementType::QUADRATIC_TETRAHEDRON,
                         {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                         {MAT_DOMAIN1});
    mesh->calculateDeterminants3D(*coords);
    // Scale determinants the same way setMeshData does, so that calcLocCoords
    // (which multiplies the stored value back by CUBIC_METER_TO_CUBIC_MICROMETER)
    // recovers the correct volume in coordinate units.
    mesh->ScaleDeterminants(qlc3d::units::CUBIC_MICROMETER_TO_CUBIC_METER);
  }
};


TEST_CASE("Create regular grid") {
// ARRANGE:
// Set up geometry with a single tetrahedron whose centre point is (0.5, 0.5, 0.5);
  auto coord = std::shared_ptr<Coordinates>(new Coordinates({
    {0, 0, 0},
    {1, 0, 0},
    {0.5, 1, 0},
    {0.5, 0.5, 1}}));

  Geometry geom;
  geom.setMeshData(1, coord,
                   {0, 1, 2, 3}, {MAT_DOMAIN1},
                   {0, 1, 2}, {MAT_ELECTRODE1});


  // ACT:
  // Create a regular grid with one point at the centre of the tetrahedron
  geom.makeRegularGrid(1, 1, 1);

  // ASSERT:
  // Check that the regular grid weight are correct by interpolating the x,y,z coordinate values
  // at each mesh node. The interpolated x, y, z values should be at (0.5, 0.5, 0.5).

  auto grid = *geom.getRegularGrid();
  SolutionVector xValues(coord->size(), 1);
  SolutionVector yValues(coord->size(), 1);
  SolutionVector zValues(coord->size(), 1);

  for (idx i = 0; i < coord->size(); ++i) {
    auto p = coord->getPoint(i);
    xValues.setValue(i, 0, p.x());
    yValues.setValue(i, 0, p.y());
    zValues.setValue(i, 0, p.z());
  }

  auto regularX = grid.interpolateToRegular(xValues);
  auto regularY = grid.interpolateToRegular(yValues);
  auto regularZ = grid.interpolateToRegular(zValues);

  REQUIRE(regularX[0] == Approx(0.5).epsilon(1e-12));
  REQUIRE(regularY[0] == Approx(0.5).epsilon(1e-12));
  REQUIRE(regularZ[0] == Approx(0.5).epsilon(1e-12));
}

// ============================================================================
// Audit spatial-lookup primitives for TET10
// ============================================================================

TEST_CASE("TET10 containsCoordinate returns true for interior point") {
  Tet10Fixture f;
  // (0.1, 0.1, 0.1) is well inside the reference tet: x+y+z = 0.3 < 1
  REQUIRE(f.mesh->containsCoordinate(0, *f.coords, Vec3(0.1, 0.1, 0.1)) == true);
}

TEST_CASE("TET10 containsCoordinate returns false for exterior point") {
  Tet10Fixture f;
  // (0.8, 0.8, 0.8): x+y+z = 2.4 >> 1, clearly outside
  REQUIRE(f.mesh->containsCoordinate(0, *f.coords, Vec3(0.8, 0.8, 0.8)) == false);
}

TEST_CASE("TET10 containsCoordinate returns true for a corner node") {
  Tet10Fixture f;
  // n0 is a corner of the element
  REQUIRE(f.mesh->containsCoordinate(0, *f.coords, Vec3(0.0, 0.0, 0.0)) == true);
}

TEST_CASE("TET10 containsCoordinate returns true for a mid-edge node position") {
  Tet10Fixture f;
  // n4 = (0.5, 0, 0) lies on the edge between n0 and n1 (boundary of the tet)
  REQUIRE(f.mesh->containsCoordinate(0, *f.coords, Vec3(0.5, 0.0, 0.0)) == true);
}

TEST_CASE("TET10 calcLocCoords sum to 1 for an interior point") {
  Tet10Fixture f;
  double loc[4] = {};
  Vec3 p(0.1, 0.2, 0.3);
  f.mesh->calcLocCoords(0, *f.coords, p, loc);
  double sum = loc[0] + loc[1] + loc[2] + loc[3];
  REQUIRE(sum == Approx(1.0).epsilon(1e-12));
}

TEST_CASE("TET10 calcLocCoords reproduces point from barycentric weights") {
  Tet10Fixture f;
  double loc[4] = {};
  Vec3 p(0.1, 0.2, 0.3);
  f.mesh->calcLocCoords(0, *f.coords, p, loc);

  // Corner coordinates
  Vec3 corners[4] = {
      f.coords->getPoint(0),  // n0
      f.coords->getPoint(1),  // n1
      f.coords->getPoint(2),  // n2
      f.coords->getPoint(3),  // n3
  };

  // Reconstruct: p_interp = sum_i( loc[i] * corner[i] )
  Vec3 pInterp(0.0, 0.0, 0.0);
  for (int i = 0; i < 4; ++i) {
    pInterp += corners[i] * loc[i];
  }

  REQUIRE(pInterp.x() == Approx(p.x()).epsilon(1e-12));
  REQUIRE(pInterp.y() == Approx(p.y()).epsilon(1e-12));
  REQUIRE(pInterp.z() == Approx(p.z()).epsilon(1e-12));
}

TEST_CASE("TET10 elementCentroid equals TET4 corner average for straight-edged element") {
  // TET4 mesh with only the 4 corner nodes
  std::vector<Vec3> cornerPts = {
      Vec3(0.0, 0.0, 0.0),
      Vec3(1.0, 0.0, 0.0),
      Vec3(0.0, 1.0, 0.0),
      Vec3(0.0, 0.0, 1.0),
  };
  auto coordsLinear = std::make_shared<Coordinates>(std::vector<Vec3>(cornerPts));
  Mesh tet4(3, ElementType::LINEAR_TETRAHEDRON);
  tet4.setElementData(ElementType::LINEAR_TETRAHEDRON, {0, 1, 2, 3}, {MAT_DOMAIN1});

  // TET10 mesh with 10 nodes (corners + midpoints)
  Tet10Fixture f;

  Vec3 c4 = tet4.elementCentroid(0, *coordsLinear);
  Vec3 c10 = f.mesh->elementCentroid(0, *f.coords);

  REQUIRE(c10.x() == Approx(c4.x()).epsilon(1e-12));
  REQUIRE(c10.y() == Approx(c4.y()).epsilon(1e-12));
  REQUIRE(c10.z() == Approx(c4.z()).epsilon(1e-12));
}

