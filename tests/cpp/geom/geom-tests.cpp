#include <catch.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <geometry.h>
#include <solutionvector.h>

#include <vector>
#include <memory>


TEST_CASE("Create regular grid") {
// ARRANGE:
// Set up geometry with a single tetrahedron whose centre point is (0.5, 0.5, 0.5);
  auto coord = std::shared_ptr<Coordinates>(new Coordinates({
    {0, 0, 0},
    {1, 0, 0},
    {0.5, 1, 0},
    {0.5, 0.5, 1}}));

  Geometry geom;
  geom.setMeshData(coord,
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