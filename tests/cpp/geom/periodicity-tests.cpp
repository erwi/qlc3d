#include <catch.h>
#include <io/gmsh-read.h>
#include <test-util.h>
#include <io/meshreader.h>
#include <mesh.h>
#include <geom/coordinates.h>
#include <geom/periodicity.h>
#include "material_numbers.h"
#include "util/containerutil.h"
#include "geom/aabox.h"

struct TestMeshData {
  Mesh triMesh;
  Coordinates coordinates;
};

TestMeshData createTestMeshData(const std::string &meshFile) {
  auto meshData = MeshReader::readMesh(meshFile);

  Mesh triMesh(2, ElementType::LINEAR_TRIANGLE);
  triMesh.setElementData(ElementType::LINEAR_TRIANGLE, std::move(meshData.triNodes), std::move(meshData.triMaterials));

  Mesh tetMesh(3, ElementType::LINEAR_TETRAHEDRON);
  tetMesh.setElementData(ElementType::LINEAR_TETRAHEDRON, std::move(meshData.tetNodes), std::move(meshData.tetMaterials));

  Coordinates coordinates(std::move(meshData.points));

  triMesh.setConnectedVolume(&tetMesh);
  triMesh.calculateSurfaceNormals(coordinates, &tetMesh);

  return {triMesh, coordinates};
}

TEST_CASE("No periodic boundaries in mesh") {
  // ARRANGE
  auto meshData = createTestMeshData(TestUtil::RESOURCE_UNIT_CUBE_DIELECTRIC_NEUMAN_GMSH_MESH);
  auto &triMesh = meshData.triMesh;
  auto &coordinates = meshData.coordinates;
  auto boundingBox = coordinates.findBoundingBox();

  // ACT
  auto peri = PeriodicNodesMapping(triMesh, coordinates);

  // ASSERT
  auto pType = peri.getPeriodicityType();
  REQUIRE_FALSE(pType.isAnyPeriodic());
  REQUIRE_FALSE(pType.isLeftRightPeriodic());
  REQUIRE_FALSE(pType.isFrontBackPeriodic());
  REQUIRE_FALSE(pType.isTopBottomPeriodic());

  auto &mapping = peri.getMapping();
  REQUIRE(mapping.size() == coordinates.size());

  // each node should be mapped to itself
  for (unsigned int i = 0; i < coordinates.size(); i++) {
    REQUIRE(mapping[i] == i);
  }
}

TEST_CASE("Periodic front-back mesh") {
  // ARRANGE
  auto meshData = createTestMeshData(TestUtil::RESOURCE_PSEUDO_2D_NEUMANN_GMSH_MESH);
  auto &triMesh = meshData.triMesh;
  auto &coordinates = meshData.coordinates;
  auto boundingBox = coordinates.findBoundingBox();

  // ACT
  auto peri = PeriodicNodesMapping(triMesh, coordinates);

  // ASSERT
  // Check that periodicity type is correctly detected
  auto pType = peri.getPeriodicityType();
  REQUIRE(pType.isAnyPeriodic());
  REQUIRE(pType.isFrontBackPeriodic());
  REQUIRE_FALSE(pType.isLeftRightPeriodic());
  REQUIRE_FALSE(pType.isTopBottomPeriodic());

  auto &mapping = peri.getMapping();
  // Check that nodes on back surface are correctly mapped to front surface, but all other surface nodes
  // should be independent.

  for (unsigned int i = 0; i < triMesh.getnElements(); i++) {

    unsigned int nodes[3] = {0, 0, 0};
    triMesh.loadNodes(i, nodes);

    for (auto n1 : nodes) {
      auto n2 = mapping[n1];
      auto p1 = coordinates.getPoint(n1);
      auto p2 = coordinates.getPoint(n2);

      if (boundingBox.frontFaceContains(p1)) {
        REQUIRE(n1 == n2);
      }
      else if (boundingBox.backFaceContains(p1)) {
        REQUIRE(p2.x() == Approx(p1.x()).margin(1e-6));
        REQUIRE(p2.y() == Approx(boundingBox.getYMin()).margin(1e-6));
        REQUIRE(p2.z() == Approx(p1.z()).margin(1e-6));
      }
      else {
        REQUIRE(triMesh.getMaterialNumber(i) != MAT_PERIODIC);
        REQUIRE(n1 == n2);
      }
    }
  }
}

TEST_CASE("Periodic left-right-front-back mesh") {
  // ARRANGE
  auto meshData = createTestMeshData(TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH);
  auto &triMesh = meshData.triMesh;
  auto &coordinates = meshData.coordinates;
  auto boundingBox = coordinates.findBoundingBox();

  // ACT
  PeriodicNodesMapping pnm(triMesh, coordinates);

  // ASSERT
  auto pType = pnm.getPeriodicityType();
  REQUIRE(pType.isAnyPeriodic());
  REQUIRE(pType.isLeftRightPeriodic());
  REQUIRE(pType.isFrontBackPeriodic());
  REQUIRE_FALSE(pType.isTopBottomPeriodic());

  auto nodeIsOnFrontLeftEdge = [&boundingBox](Vec3 &p) { return boundingBox.frontFaceContains(p) && boundingBox.leftFaceContains(p); };
  auto nodeIsOnFrontRightEdge = [&boundingBox](Vec3 &p) { return boundingBox.frontFaceContains(p) && boundingBox.rightFaceContains(p); };
  auto nodeIsOnBackLeftEdge = [&boundingBox](Vec3 &p) { return boundingBox.backFaceContains(p) && boundingBox.leftFaceContains(p); };
  auto nodeIsOnBackRightEdge = [&boundingBox](Vec3 &p) { return boundingBox.backFaceContains(p) && boundingBox.rightFaceContains(p); };

  auto &mapping = pnm.getMapping();
  for (unsigned int i = 0; i < triMesh.getnElements(); i++) {
    unsigned int nodes[3] = {0, 0, 0};
    triMesh.loadNodes(i, nodes);

    // Check that nodes map to front or left surfaces on all vertical surfaces, but are independent on horizontal surfaces.
    // Note the order of checking is important here as edges nodes exist in multiple faces
    for (auto n1 : nodes) {
      auto n2 = mapping[n1];

      auto p1 = coordinates.getPoint(n1);
      auto p2 = coordinates.getPoint(n2);

      if (nodeIsOnFrontLeftEdge(p1)) {
        // xmin ymin vertical corner is independent
        REQUIRE(n1 == n2);
      }
      else if (nodeIsOnBackLeftEdge(p1) || nodeIsOnBackRightEdge(p1) || nodeIsOnFrontRightEdge(p1)) {
        // other vertical edges map to xmin-ymin corner
        REQUIRE(p2.x() == Approx(boundingBox.getXMin()).margin(1e-6));
        REQUIRE(p2.y() == Approx(boundingBox.getYMin()).margin(1e-6));
        REQUIRE(p1.z() == Approx(p2.z()).margin(1e-6));
      }
      else if (boundingBox.leftFaceContains(p1) || boundingBox.frontFaceContains(p1)) {
        // independent front and left surfaces
        REQUIRE(n1 == n2);
      }
      else if (boundingBox.rightFaceContains(p1)) {
        // right face maps to left surface
        REQUIRE(p2.x() == Approx(boundingBox.getXMin()).margin(1e-6));
        REQUIRE(p2.y() == Approx(p1.y()).margin(1e-6));
        REQUIRE(p2.z() == Approx(p1.z()).margin(1e-6));
      }
      else if (boundingBox.backFaceContains(p1)) {
        // back face maps to front surface
        REQUIRE(p2.x() == Approx(p1.x()).margin(1e-6));
        REQUIRE(p2.y() == Approx(boundingBox.getYMin()).margin(1e-6));
        REQUIRE(p2.z() == Approx(p1.z()).margin(1e-6));
      }
      else if (boundingBox.topFaceContains(p1) || boundingBox.bottomFaceContains(p1)) {
        // both top and bottom surfaces are independent
        REQUIRE(n1 == n2);
      }
    }
  }
}

TEST_CASE("Periodic left-right-top-bottom-front-back mesh") {
  auto meshData = createTestMeshData(TestUtil::RESOURCE_UNIT_CUBE_FULLY_PERIODIC_LINEAR_GMSH_MESH);
  auto &triMesh = meshData.triMesh;
  auto &coordinates = meshData.coordinates;
  auto boundingBox = coordinates.findBoundingBox();

  // ACT
  PeriodicNodesMapping pnm(triMesh, coordinates);

  // ASSERT
  auto pType = pnm.getPeriodicityType();
  REQUIRE(pType.isAnyPeriodic());
  REQUIRE(pType.isLeftRightPeriodic());
  REQUIRE(pType.isFrontBackPeriodic());
  REQUIRE(pType.isTopBottomPeriodic());

  // Check that all faces map to front, left and bottom surfaces.
  // All corners should map to bottom-front-left corner.
  // All edges should map to edges those joining at bottom-front-left corner.

  auto nodeIsOnFrontLeftBottomCorner = [&boundingBox](Vec3 &p) {
    return boundingBox.frontFaceContains(p) && boundingBox.leftFaceContains(p) && boundingBox.bottomFaceContains(p);
  };

  auto nodeIsAnyCorner = [&boundingBox](Vec3 &p) {
    int onFaceCount =
            (boundingBox.frontFaceContains(p) ? 1 : 0) +
            (boundingBox.backFaceContains(p) ? 1 : 0) +
            (boundingBox.leftFaceContains(p) ? 1 : 0) +
            (boundingBox.rightFaceContains(p) ? 1 : 0) +
            (boundingBox.topFaceContains(p) ? 1 : 0) +
            (boundingBox.bottomFaceContains(p) ? 1 : 0);
    assert(onFaceCount <= 3);
    return onFaceCount == 3;
  };

  auto nodeIsFrontLeftEdge = [&boundingBox](Vec3 &p) { return boundingBox.frontFaceContains(p) && boundingBox.leftFaceContains(p); };
  auto nodeIsFrontBottomEdge = [&boundingBox](Vec3 &p) { return boundingBox.frontFaceContains(p) && boundingBox.bottomFaceContains(p); };
  auto nodeIsLeftBottomEdge = [&boundingBox](Vec3 &p) { return boundingBox.leftFaceContains(p) && boundingBox.bottomFaceContains(p); };
  auto nodeIsAnyEdge = [&boundingBox](Vec3 &p) {
    int onFaceCount =
            (boundingBox.frontFaceContains(p) ? 1 : 0) +
            (boundingBox.backFaceContains(p) ? 1 : 0) +
            (boundingBox.leftFaceContains(p) ? 1 : 0) +
            (boundingBox.rightFaceContains(p) ? 1 : 0) +
            (boundingBox.topFaceContains(p) ? 1 : 0) +
            (boundingBox.bottomFaceContains(p) ? 1 : 0);
    assert(onFaceCount <= 2);
    return onFaceCount == 2;
  };

  auto &mapping = pnm.getMapping();

  for (unsigned int indTri = 0; indTri < triMesh.getnElements(); indTri++) {
    unsigned int nodes[3] = {0, 0, 0};
    triMesh.loadNodes(indTri, nodes);

    for (auto n1 : nodes) {
      auto n2 = mapping[n1];
      auto p1 = coordinates.getPoint(n1);
      auto p2 = coordinates.getPoint(n2);

      if (nodeIsAnyCorner(p1)) {
        if (nodeIsOnFrontLeftBottomCorner(p1)) {
          // front-left-bottom corner is independent
          REQUIRE(n1 == n2);
        } else {
          // other corners map to bottom-front-left corner
          REQUIRE(n1 != n2);
          REQUIRE(nodeIsOnFrontLeftBottomCorner(p2));
        }
      }
      else if (nodeIsAnyEdge(p1)) {
        if (nodeIsFrontLeftEdge(p1) || nodeIsFrontBottomEdge(p1) || nodeIsLeftBottomEdge(p1)) {
          // three independent edges
          REQUIRE(n1 == n2);
        } else {
          REQUIRE(n1 != n2);
          if (!boundingBox.frontFaceContains(p1) && !boundingBox.backFaceContains(p1)) {
            // all edges along y-axis map to bottom left edge
            REQUIRE(p1.y() == Approx(p2.y()).margin(1e-6));
            REQUIRE(nodeIsLeftBottomEdge(p2));
          }
          else if(!boundingBox.leftFaceContains(p1) && !boundingBox.rightFaceContains(p1)) {
            // all edges along x-axis map to bottom front edge
            REQUIRE(p1.x() == Approx(p2.x()).margin(1e-6));
            REQUIRE(nodeIsFrontBottomEdge(p2));
          }
          else if (!boundingBox.bottomFaceContains(p1) && !boundingBox.topFaceContains(p1)) {
            // all edges along z-axis map to front left edge
            REQUIRE(p1.z() == Approx(p2.z()).margin(1e-6));
            REQUIRE(nodeIsFrontLeftEdge(p2));
          } else {
            FAIL("Could not determine which edge node is on");
          }
        }
      }
      else if (boundingBox.rightFaceContains(p1)) {
        // right face maps to left surface
        REQUIRE(p2.x() == Approx(boundingBox.getXMin()).margin(1e-6));
        REQUIRE(p2.y() == Approx(p1.y()).margin(1e-6));
        REQUIRE(p2.z() == Approx(p1.z()).margin(1e-6));
      }
      else if (boundingBox.topFaceContains(p1)) {
        // top face maps to bottom surface
        REQUIRE(p2.x() == Approx(p1.x()).margin(1e-6));
        REQUIRE(p2.y() == Approx(p1.y()).margin(1e-6));
        REQUIRE(p2.z() == Approx(boundingBox.getZMin()).margin(1e-6));
      }
      else if (boundingBox.backFaceContains(p1)) {
        // back face maps to front surface
        REQUIRE(p2.x() == Approx(p1.x()).margin(1e-6));
        REQUIRE(p2.y() == Approx(boundingBox.getYMin()).margin(1e-6));
        REQUIRE(p2.z() == Approx(p1.z()).margin(1e-6));
      }
      else if (boundingBox.leftFaceContains(p1) || boundingBox.frontFaceContains(p1) || boundingBox.bottomFaceContains(p1)) {
        REQUIRE(n1 == n2);
      }
      else {
        FAIL("Node is not on any face");
      }
    }
  }
}
