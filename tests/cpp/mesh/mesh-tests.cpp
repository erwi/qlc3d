#include <catch.h>
#include <mesh/mesh.h>

TEST_CASE("Mesh dimension follows the element type") {
  auto tri = Mesh::triangleMesh();
  REQUIRE(tri->getElementType() == ElementType::LINEAR_TRIANGLE);
  REQUIRE(tri->getDimension() == 2);
  REQUIRE(tri->getnNodes() == 3);

  auto tet = Mesh::tetMesh();
  REQUIRE(tet->getElementType() == ElementType::LINEAR_TETRAHEDRON);
  REQUIRE(tet->getDimension() == 3);
  REQUIRE(tet->getnNodes() == 4);
}

TEST_CASE("Mesh dimension remains consistent after mesh mutations") {
  Mesh triMesh(2, ElementType::LINEAR_TRIANGLE);
  triMesh.setElementData(ElementType::QUADRATIC_TRIANGLE,
                         std::vector<unsigned int>{0, 1, 2, 3, 4, 5},
                         std::vector<unsigned int>{1});
  REQUIRE(triMesh.getElementType() == ElementType::QUADRATIC_TRIANGLE);
  REQUIRE(triMesh.getDimension() == 2);
  REQUIRE(triMesh.getnNodes() == 6);

  triMesh.ClearMesh();
  REQUIRE(triMesh.getElementType() == ElementType::QUADRATIC_TRIANGLE);
  REQUIRE(triMesh.getDimension() == 2);
  REQUIRE(triMesh.getnElements() == 0);

  Mesh tetMesh(3, ElementType::LINEAR_TETRAHEDRON);
  tetMesh.setElementData(ElementType::QUADRATIC_TETRAHEDRON,
                         std::vector<unsigned int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                         std::vector<unsigned int>{1});
  REQUIRE(tetMesh.getElementType() == ElementType::QUADRATIC_TETRAHEDRON);
  REQUIRE(tetMesh.getDimension() == 3);
  REQUIRE(tetMesh.getnNodes() == 10);

  tetMesh.ClearMesh();
  REQUIRE(tetMesh.getElementType() == ElementType::QUADRATIC_TETRAHEDRON);
  REQUIRE(tetMesh.getDimension() == 3);
  REQUIRE(tetMesh.getnElements() == 0);
}

