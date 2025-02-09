#include <catch.h>
#include <test-util.h>
#include <sparsematrix.h>
#include <material_numbers.h>
#include <geometry.h>
#include <solutionvector.h>
#include <inits.h>
#include <spamtrix_ircmatrix.hpp>

TEST_CASE("Create SpamTrix - linear and quadratic meshes") {
  // GIVEN
  std::string meshName = GENERATE(TestUtil::RESOURCE_THIN_GID_MESH, TestUtil::RESOURCE_THIN_QUADRATIC_GMSH_MESH);

  Geometry geom;
  prepareGeometryWithDefaultBoundaries(geom, meshName);

  SolutionVector q(geom.getnpLC(), 5);
  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));

  setSurfacesQ(q, alignment, 0.5, geom);
  q.initialiseLcBoundaries(geom, alignment);

  auto &dofMap = q.getDofMap();

  // WHEN
  auto m = createQMatrix(geom, q, MAT_DOMAIN1);

  // THEN
  unsigned int expectedSize = 5 * dofMap.getnFreeNodes();
  REQUIRE(m->getNumRows() == expectedSize);
  REQUIRE(m->getNumCols() == expectedSize);
  REQUIRE(m->getnnz() > expectedSize); // at least diagonal non-zeros, but exact number depends on mesh connectivity
}
