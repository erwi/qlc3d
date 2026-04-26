#include <catch.h>
#include <geom/tet-mesh-search.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <geom/aabox.h>
#include <mesh/mesh.h>
#include <material_numbers.h>
#include <globals.h>

#include <memory>
#include <vector>

// ============================================================================
// Helper: build a single-element LINEAR_TETRAHEDRON mesh that encloses a unit-ish tet.
// Corner nodes: (0,0,0), (1,0,0), (0.5,1,0), (0.5,0.5,1)
// Centroid = (0.5, 0.5, 0.5)
// ============================================================================
struct SingleTetFixture {
  std::shared_ptr<Coordinates> coords;
  std::shared_ptr<Mesh> mesh;
  AABox box;

  SingleTetFixture() {
    std::vector<Vec3> pts = {
        {0, 0, 0}, {1, 0, 0}, {0.5, 1, 0}, {0.5, 0.5, 1}
    };
    coords = std::make_shared<Coordinates>(std::move(pts));
    mesh = Mesh::tetMesh();
    mesh->setElementData(ElementType::LINEAR_TETRAHEDRON, {0, 1, 2, 3}, {MAT_DOMAIN1});
    mesh->calculateDeterminants3D(*coords);
    mesh->ScaleDeterminants(qlc3d::units::CUBIC_MICROMETER_TO_CUBIC_METER);
    box = coords->findBoundingBox();
  }
};

// ============================================================================
// Helper: build a two-tet mesh sharing a face.
// Tet 0: nodes (0,0,0) (1,0,0) (0.5,1,0) (0.5,0.5,1)
// Tet 1: nodes (1,0,0) (2,0,0) (1.5,1,0) (1.5,0.5,1)
// ============================================================================
struct TwoTetFixture {
  std::shared_ptr<Coordinates> coords;
  std::shared_ptr<Mesh> mesh;
  AABox box;

  TwoTetFixture() {
    std::vector<Vec3> pts = {
        {0, 0, 0}, {1, 0, 0}, {0.5, 1, 0}, {0.5, 0.5, 1},
        {2, 0, 0}, {1.5, 1, 0}, {1.5, 0.5, 1}
    };
    coords = std::make_shared<Coordinates>(std::move(pts));
    mesh = Mesh::tetMesh();
    // Tet 0: 0,1,2,3;  Tet 1: 1,4,5,6
    mesh->setElementData(ElementType::LINEAR_TETRAHEDRON,
                         {0, 1, 2, 3,
                          1, 4, 5, 6},
                         {MAT_DOMAIN1, MAT_DOMAIN1});
    mesh->calculateDeterminants3D(*coords);
    mesh->ScaleDeterminants(qlc3d::units::CUBIC_MICROMETER_TO_CUBIC_METER);
    box = coords->findBoundingBox();
  }
};

// ============================================================================
// Helper: one LC tet + one dielectric tet
// ============================================================================
struct MixedMatFixture {
  std::shared_ptr<Coordinates> coords;
  std::shared_ptr<Mesh> mesh;
  AABox box;
  Vec3 lcCentroid;
  Vec3 dielCentroid;

  MixedMatFixture() {
    // LC tet at x in [0,1], dielectric tet at x in [1,2]
    std::vector<Vec3> pts = {
        {0, 0, 0}, {1, 0, 0}, {0.5, 1, 0}, {0.5, 0.5, 1},
        {2, 0, 0}, {1.5, 1, 0}, {1.5, 0.5, 1}
    };
    coords = std::make_shared<Coordinates>(std::move(pts));
    mesh = Mesh::tetMesh();
    mesh->setElementData(ElementType::LINEAR_TETRAHEDRON,
                         {0, 1, 2, 3,
                          1, 4, 5, 6},
                         {MAT_DOMAIN1, MAT_DIELECTRIC1});
    mesh->calculateDeterminants3D(*coords);
    mesh->ScaleDeterminants(qlc3d::units::CUBIC_MICROMETER_TO_CUBIC_METER);
    box = coords->findBoundingBox();

    // Approximate centroids
    lcCentroid   = {(0+1+0.5+0.5)/4.0, (0+0+1+0.5)/4.0, (0+0+0+1)/4.0};
    dielCentroid = {(1+2+1.5+1.5)/4.0, (0+0+1+0.5)/4.0, (0+0+0+1)/4.0};
  }
};

// ============================================================================
// Test 1 — bruteForceSearch finds the only tet (LINEAR_TETRAHEDRON)
// ============================================================================
TEST_CASE("TetMeshSearch: bruteForceSearch finds single tet for centroid") {
  SingleTetFixture f;
  TetMeshSearch search(*f.mesh, *f.coords, f.box);

  Vec3 centroid(0.5, 0.5, 0.5);
  unsigned int ind = 99;
  bool found = search.bruteForceSearch(ind, centroid, false, false);

  REQUIRE(found == true);
  REQUIRE(ind == 0u);
}

// ============================================================================
// Test 2 — bruteForceSearch returns false for a point outside
// ============================================================================
TEST_CASE("TetMeshSearch: bruteForceSearch returns false for exterior point") {
  SingleTetFixture f;
  TetMeshSearch search(*f.mesh, *f.coords, f.box);

  Vec3 outside(5, 5, 5);
  unsigned int ind = 99;
  bool found = search.bruteForceSearch(ind, outside, /*terminateOnError=*/false, false);

  REQUIRE(found == false);
  REQUIRE(ind == 99u); // unchanged
}

// ============================================================================
// Test 3 — genIndToTetsByCoords with two query points in different tets
// ============================================================================
TEST_CASE("TetMeshSearch: genIndToTetsByCoords locates points in correct tets") {
  TwoTetFixture f;
  TetMeshSearch search(*f.mesh, *f.coords, f.box);

  // Pick points clearly inside each tet (slightly off centroid toward interior)
  std::vector<Vec3> pts = {
      {0.5, 0.4, 0.4},   // inside tet 0
      {1.5, 0.4, 0.4}    // inside tet 1
  };
  Coordinates targets{std::vector<Vec3>(pts)};

  std::vector<unsigned int> result;
  REQUIRE_NOTHROW(search.genIndToTetsByCoords(result, targets, false, false));
  REQUIRE(result.size() == 2u);
  REQUIRE(result[0] == 0u);
  REQUIRE(result[1] == 1u);
}

// ============================================================================
// Test 4 — genIndToTetsByCoords marks exterior points as NOT_AN_INDEX
// ============================================================================
TEST_CASE("TetMeshSearch: genIndToTetsByCoords marks out-of-mesh points as NOT_AN_INDEX") {
  SingleTetFixture f;
  TetMeshSearch search(*f.mesh, *f.coords, f.box);

  std::vector<Vec3> pts = {{100, 100, 100}};
  Coordinates targets{std::vector<Vec3>(pts)};

  std::vector<unsigned int> result;
  search.genIndToTetsByCoords(result, targets, /*terminateOnError=*/false, false);

  REQUIRE(result.size() == 1u);
  REQUIRE(result[0] == TetMeshSearch::NOT_AN_INDEX);
}

// ============================================================================
// Test 5 — requireLCElement excludes dielectric elements
// ============================================================================
TEST_CASE("TetMeshSearch: requireLCElement excludes dielectric elements") {
  MixedMatFixture f;
  TetMeshSearch search(*f.mesh, *f.coords, f.box);

  // Query the centroid of the dielectric tet with requireLCElement=true
  std::vector<Vec3> pts = {f.dielCentroid};
  Coordinates targets{std::vector<Vec3>(pts)};

  // With LC requirement: should not find (or return NOT_AN_INDEX)
  std::vector<unsigned int> resultLC;
  search.genIndToTetsByCoords(resultLC, targets, false, true);
  REQUIRE(resultLC[0] == TetMeshSearch::NOT_AN_INDEX);

  // Without LC requirement: should find the dielectric tet
  std::vector<unsigned int> resultAny;
  search.genIndToTetsByCoords(resultAny, targets, false, false);
  REQUIRE(resultAny[0] != TetMeshSearch::NOT_AN_INDEX);
  REQUIRE(f.mesh->getMaterialNumber(resultAny[0]) == MAT_DIELECTRIC1);
}




