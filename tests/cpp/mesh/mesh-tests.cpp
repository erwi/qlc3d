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

TEST_CASE("Mesh element counts are derived from node storage after mutations") {
  Mesh triMesh(2, ElementType::LINEAR_TRIANGLE);
  triMesh.setElementData(ElementType::LINEAR_TRIANGLE,
                         std::vector<unsigned int>{0, 1, 2, 3, 4, 5},
                         std::vector<unsigned int>{11, 12});
  REQUIRE(triMesh.getnElements() == 2);

  triMesh.appendElements(std::vector<unsigned int>{6, 7, 8}, std::vector<unsigned int>{13});
  REQUIRE(triMesh.getnElements() == 3);
  REQUIRE(triMesh.getNode(2, 2) == 8);

  Mesh triCopy(2, ElementType::LINEAR_TRIANGLE);
  triCopy.CopyMesh(&triMesh);
  REQUIRE(triCopy.getElementType() == ElementType::LINEAR_TRIANGLE);
  REQUIRE(triCopy.getnElements() == 3);
  REQUIRE(triCopy.getNode(2, 2) == 8);

  triMesh.ClearMesh();
  REQUIRE(triMesh.getnElements() == 0);
  REQUIRE(triCopy.getnElements() == 3);

  Mesh tetMesh(3, ElementType::QUADRATIC_TETRAHEDRON);
  tetMesh.setElementData(ElementType::QUADRATIC_TETRAHEDRON,
                         std::vector<unsigned int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                                   10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
                         std::vector<unsigned int>{21, 22});
  REQUIRE(tetMesh.getnElements() == 2);

  tetMesh.appendElements(std::vector<unsigned int>{20, 21, 22, 23, 24, 25, 26, 27, 28, 29},
                         std::vector<unsigned int>{23});
  REQUIRE(tetMesh.getnElements() == 3);
  REQUIRE(tetMesh.getNode(2, 9) == 29);

  Mesh tetCopy(3, ElementType::QUADRATIC_TETRAHEDRON);
  tetCopy.CopyMesh(&tetMesh);
  REQUIRE(tetCopy.getElementType() == ElementType::QUADRATIC_TETRAHEDRON);
  REQUIRE(tetCopy.getnElements() == 3);
  REQUIRE(tetCopy.getNode(2, 9) == 29);
}

