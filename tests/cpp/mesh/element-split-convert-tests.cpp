#include <test-util.h>
#include <io/gmsh-read.h>

#include <geometry.h>
#include <material_numbers.h>

#include "catch.h"
#include "mesh/element-split-convert.h"

namespace {

  // Return labeled coordinates for the quadratic tetrahedron used in the tests.
  // Labels follow the convention used in the tests: A, B, C, D, ab, bc, ac, ad, bd, cd
  static std::vector<Vec3> makeQuadraticTestCoordinates() {
    return std::vector<Vec3> {
      Vec3(0, 0, 0), // A
      Vec3(1, 0, 0), // B
      Vec3(0, 1, 0), // C
      Vec3(0, 0, 1), // D
      Vec3(0.5, 0, 0), // ab
      Vec3(0.5, 0.5, 0), // bc
      Vec3(0, 0.5, 0), // ac
      Vec3(0, 0, 0.5), // ad
      Vec3(0.5, 0, 0.5), // bd
      Vec3(0, 0.5, 0.5) // cd
    };
  }


RawMeshData makeQuadraticTestRawMeshData() {
  std::vector<Vec3> points = makeQuadraticTestCoordinates();

  return RawMeshData(2, std::move(points),
                     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {MAT_DOMAIN1},
                     {0, 1, 2, 4, 5, 6,
                      0, 1, 3, 4, 8, 7,
                      0, 2, 3, 6, 9, 7,
                      1, 2, 3, 5, 9, 8},
                     {MAT_FIXLC1, MAT_FIXLC1, MAT_FIXLC1, MAT_FIXLC1});
}

Geometry makeQuadraticTestGeometry() {
  return Geometry::fromRawMeshData(makeQuadraticTestRawMeshData());
}

} // namespace



TEST_CASE("Split and recombine tets") {
  // Tet node numbering as in fig 6.1 (b) in Eero's thesis
  unsigned int A = 0;
  unsigned int B = 1;
  unsigned int C = 2;
  unsigned int D = 3;

  unsigned int ab = 4;
  unsigned int bc = 5;
  unsigned int ac = 6;
  unsigned int ad = 7;
  unsigned int bd = 8;
  unsigned int cd = 9;

  SECTION("Split quadratic tetrahedron to 8 linear tets") {
    std::vector<unsigned int> quadraticTetrahedron = {A, B, C, D, ab, bc, ac, ad, bd, cd};

    // ACT
    std::vector<std::vector<unsigned int>> linearTets = splitQuadraticTetrahedronToLinear(quadraticTetrahedron);
    REQUIRE(linearTets.size() == 8);
    // corner tets
    REQUIRE(linearTets[0] == std::vector<unsigned int>({A, ab, ac, ad}));
    REQUIRE(linearTets[1] == std::vector<unsigned int>({B, bc, ab, bd}));
    REQUIRE(linearTets[2] == std::vector<unsigned int>({C, cd, ac, bc}));
    REQUIRE(linearTets[3] == std::vector<unsigned int>({D, cd, bd, ad}));

    // face tets
    REQUIRE(linearTets[4] == std::vector<unsigned int>({bd, ab, ac, ad})); // A-face
    REQUIRE(linearTets[5] == std::vector<unsigned int>({ac, bc, ab, bd})); // B-face
    REQUIRE(linearTets[6] == std::vector<unsigned int>({bd, cd, ac, bc})); // C-face
    REQUIRE(linearTets[7] == std::vector<unsigned int>({ac, cd, bd, ad})); // D-face
  }

  SECTION("Recombine 8 linear tets to quadratic tetrahedron") {
    std::vector<std::vector<unsigned int>> linearTets = {
      {A, ab, ac, ad},
      {B, bc, ab, bd},
      {C, cd, ac, bc},
      {D, cd, bd, ad},
      {bd, ab, ac, ad},
      {ac, bc, ab, bd},
      {bd, cd, ac, bc},
      {ac, cd, bd, ad}
    };
    auto coords = makeQuadraticTestCoordinates();
    std::vector<unsigned int> recombinedTet = recombineLinearTetsToQuadratic(linearTets, coords);
    REQUIRE(recombinedTet.size() == 10);
    REQUIRE(recombinedTet[0] == A);
    REQUIRE(recombinedTet[1] == B);
    REQUIRE(recombinedTet[2] == C);
    REQUIRE(recombinedTet[3] == D);
    REQUIRE(recombinedTet[4] == ab);
    REQUIRE(recombinedTet[5] == bc);
    REQUIRE(recombinedTet[6] == ac);
    REQUIRE(recombinedTet[7] == ad);
    REQUIRE(recombinedTet[8] == bd);
    REQUIRE(recombinedTet[9] == cd);

    // the determinant of the tet defined by the corners should be positive
    double jDet = det3D(coords[recombinedTet[0]], coords[recombinedTet[1]], coords[recombinedTet[2]], coords[recombinedTet[3]]);
    REQUIRE(jDet == Approx(1.0).epsilon(1e-9));
  }

  SECTION("Recombine 8 linear tets to quadratic tetrahedron when node order is swapped") {
    std::vector<std::vector<unsigned int>> linearTets = {
      {A, ab, ac, ad},
      {B, bc, ab, bd},
      {D, cd, bd, ad}, // <--- C and D order is swapped = negative jacobian determinant
      {C, cd, ac, bc}, // <---
      {bd, ab, ac, ad},
      {ac, bc, ab, bd},
      {bd, cd, ac, bc},
      {ac, cd, bd, ad}
    };
    auto coords = makeQuadraticTestCoordinates();
    std::vector<unsigned int> recombinedTet = recombineLinearTetsToQuadratic(linearTets, coords);


    // the determinant of the tet defined by the corners should be positive
    double jDet = det3D(coords[recombinedTet[0]], coords[recombinedTet[1]], coords[recombinedTet[2]], coords[recombinedTet[3]]);
    REQUIRE(jDet == Approx(1.0).epsilon(1e-9));
  }
}

TEST_CASE("Split and recombine triangles") {
  unsigned int A = 0;
  unsigned int B = 1;
  unsigned int C = 2;
  unsigned int ab = 3;
  unsigned int bc = 4;
  unsigned int ac = 5;

  SECTION("Split quadratic triangle to 4 linear triangles") {
    std::vector<unsigned int> quadraticTriangle = {A, B, C, ab, bc, ac};

    std::vector<std::vector<unsigned int>> linearTriangles = splitQuadraticTriangleToLinear(quadraticTriangle);
    REQUIRE(linearTriangles.size() == 4);
    REQUIRE(linearTriangles[0] == std::vector<unsigned int>({A, ab, ac}));
    REQUIRE(linearTriangles[1] == std::vector<unsigned int>({B, bc, ab}));
    REQUIRE(linearTriangles[2] == std::vector<unsigned int>({C, ac, bc}));
    REQUIRE(linearTriangles[3] == std::vector<unsigned int>({ab, bc, ac}));
  }

  SECTION("Recombine 4 linear triangles to 1 quadratic triangle") {
    std::vector<std::vector<unsigned int>> linearTris = {
      {A, ab, ac},
      {B, bc, ab},
      {C, ac, bc},
      {ab, bc, ac}
    };

    std::vector<unsigned int> recombinedTri = recombineLinearTrianglesToQuadratic(linearTris);
    REQUIRE(recombinedTri.size() == 6);
    REQUIRE(recombinedTri[0] == A);
    REQUIRE(recombinedTri[1] == B);
    REQUIRE(recombinedTri[2] == C);
    REQUIRE(recombinedTri[3] == ab);
    REQUIRE(recombinedTri[4] == bc);
    REQUIRE(recombinedTri[5] == ac);
  }
}

TEST_CASE("Convert single linear tet mesh to quadratic tet mesh") {
  // create raw mesh data for a single linear tetrahedron and triangle
  std::vector<Vec3> points = {
    Vec3(0, 0, 0),
    Vec3(1, 0, 0),
    Vec3(0, 1, 0),
    Vec3(0, 0, 1)
  };

  RawMeshData meshData(1, points, {0, 1, 2, 3}, {99}, {0, 1, 2}, {100});

  // ACT
  convertLinearMeshDataToQuadratic(meshData);

  // ASSERT
  REQUIRE(meshData.getElementOrder() == 2);
  REQUIRE(meshData.points.size() == 10);
  REQUIRE(meshData.tetNodes.size() == 10);
  REQUIRE(meshData.tetMaterials.size() == 1);
  REQUIRE(meshData.triNodes.size() == 6);
  REQUIRE(meshData.triMaterials.size() == 1);
}

TEST_CASE("Convert mesh with two linear tets sharing a face into quadratic mesh") {
  // create raw mesh data for two linear tets sharing a face
  std::vector<Vec3> points = {
    Vec3(0, 0, 0), // tet 1 and 2
    Vec3(1, 0, 0), // tet 1 and 2
    Vec3(0, 1, 0), // tet 1 and 2
    Vec3(0, 0, 1), // tet 1
    Vec3(0, 0, -1) // tet 2
  };

  RawMeshData meshData(1, points, {0, 1, 2, 3, 0, 1, 2, 4}, {99, 99}, {0, 1, 2}, {100});

  // ACT
  convertLinearMeshDataToQuadratic(meshData);

  // ASSERT
  REQUIRE(meshData.getElementOrder() == 2);
  REQUIRE(meshData.points.size() == 14);

  auto expectedTetNodes = std::vector<unsigned int>(
    {
      0, 1, 2, 3, 5, 6, 7, 8, 9, 10,
      0, 1, 2, 4, 5, 6, 7, 11, 12, 13
    });

  REQUIRE(meshData.tetNodes == expectedTetNodes); // two 2nd order tets that share a face

  auto expectedTriNodes = std::vector<unsigned int>({0, 1, 2, 5, 6, 7});
  REQUIRE(meshData.triNodes == expectedTriNodes); // single 2nd order triangle

  REQUIRE(meshData.tetMaterials.size() == 2); // unchanged
  REQUIRE(meshData.triMaterials.size() == 1); // unchanged
}


TEST_CASE("Convert linear Gmsh mesh to quadratic mesh") {
  auto meshData = MeshReader::readMesh(TestUtil::RESOURCE_UNIT_CUBE_NEUMANN_GMSH_MESH);

  auto numTetsOld = meshData.tetMaterials.size();
  auto numTrisOld = meshData.triMaterials.size();

  // ACT
  convertLinearMeshDataToQuadratic(meshData);

  // ASSERT
  REQUIRE(meshData.getElementOrder() == 2);
  REQUIRE(meshData.tetMaterials.size() == numTetsOld);
  REQUIRE(meshData.triMaterials.size() == numTrisOld);

  REQUIRE(meshData.tetNodes.size() == numTetsOld * 10);

}

TEST_CASE("Split quadratic geometry to linear mesh data") {
  Geometry geom = makeQuadraticTestGeometry();
  RawMeshData linearMeshData(1, {}, {}, {}, {}, {});

  splitQuadraticGeometryToLinear(geom, linearMeshData);

  REQUIRE(linearMeshData.getElementOrder() == 1);
  REQUIRE(linearMeshData.points.size() == geom.getnp());
  REQUIRE(linearMeshData.tetNodes.size() == 8 * 4);
  REQUIRE(linearMeshData.tetMaterials.size() == 8);
  REQUIRE(linearMeshData.triNodes.size() == 4 * 4 * 3);
  REQUIRE(linearMeshData.triMaterials.size() == 4 * 4);
}

TEST_CASE("Geometry::fromRawMeshData builds a quadratic geometry") {
  Geometry geom = Geometry::fromRawMeshData(makeQuadraticTestRawMeshData());

  REQUIRE(geom.getTetrahedra().getElementType() == ElementType::QUADRATIC_TETRAHEDRON);
  REQUIRE(geom.getTriangles().getElementType() == ElementType::QUADRATIC_TRIANGLE);
  REQUIRE(geom.getTetrahedra().getnElements() == 1);
  REQUIRE(geom.getTriangles().getnElements() == 4);
}

TEST_CASE("Split quadratic geometry and re-promote it to quadratic mesh data") {
  Geometry geom = makeQuadraticTestGeometry();
  RawMeshData linearMeshData(1, {}, {}, {}, {}, {});

  splitQuadraticGeometryToLinear(geom, linearMeshData);
  convertLinearMeshDataToQuadratic(linearMeshData);

  REQUIRE(linearMeshData.getElementOrder() == 2);
  REQUIRE(linearMeshData.points.size() > geom.getnp());
  REQUIRE(linearMeshData.tetNodes.size() == 8 * 10);
  REQUIRE(linearMeshData.tetMaterials.size() == 8);
  REQUIRE(linearMeshData.triNodes.size() == 16 * 6);
  REQUIRE(linearMeshData.triMaterials.size() == 16);
}

TEST_CASE("split and recombine quadratic gmsh mesh") {
  auto rawMeshData = MeshReader::readMesh(TestUtil::RESOURCE_SMALL_CUBE_QUADRATIC_GMSH_MESH);

  auto geom = Geometry::fromRawMeshData(rawMeshData);
  RawMeshData linearMeshData(1, {}, {}, {}, {}, {});

  splitQuadraticGeometryToLinear(geom, linearMeshData);

  REQUIRE(recombineLinearisedMeshToQuadratic(linearMeshData) == true);

  // ASSERT: linearMeshData should match the original quadratic mesh data
  // ARRANGE/ACT above: rawMeshData -> Geometry -> split -> recombine

  // element order
  REQUIRE(linearMeshData.getElementOrder() == rawMeshData.getElementOrder());

  // points: same count and coordinates (compare via distance to avoid relying on Vec3 internals)
  REQUIRE(linearMeshData.points.size() == rawMeshData.points.size());
  for (size_t i = 0; i < rawMeshData.points.size(); ++i) {
    REQUIRE(linearMeshData.points[i].distance(rawMeshData.points[i]) == Approx(0.0).margin(1e-9));
  }

  // tetrahedral elements and materials
  REQUIRE(linearMeshData.tetNodes.size() == rawMeshData.tetNodes.size());
  if (linearMeshData.tetNodes != rawMeshData.tetNodes) {
    // Log detailed mismatch information to help debug swapped node ordering.
    Log::info("Tet nodes differ: linearMeshData.tetNodes.size()={}, rawMeshData.tetNodes.size()={}", linearMeshData.tetNodes.size(), rawMeshData.tetNodes.size());
    const size_t nodesPerTet = 10;
    const size_t nTets = rawMeshData.tetNodes.size() / nodesPerTet;
    for (size_t t = 0; t < nTets; ++t) {
      size_t base = t * nodesPerTet;
      bool equal = true;
      for (size_t k = 0; k < nodesPerTet; ++k) {
        if (linearMeshData.tetNodes[base + k] != rawMeshData.tetNodes[base + k]) { equal = false; break; }
      }
      if (!equal) {
        std::string lstr;
        std::string rstr;
        for (size_t k = 0; k < nodesPerTet; ++k) {
          lstr += std::to_string(linearMeshData.tetNodes[base + k]);
          rstr += std::to_string(rawMeshData.tetNodes[base + k]);
          if (k + 1 < nodesPerTet) { lstr += ", "; rstr += ", "; }
        }
        Log::info("Tet {}: linear = [{}]", t, lstr);
        Log::info("Tet {}: raw    = [{}]", t, rstr);
      }
    }
  }
  REQUIRE(linearMeshData.tetNodes == rawMeshData.tetNodes);
  REQUIRE(linearMeshData.tetMaterials.size() == rawMeshData.tetMaterials.size());
  REQUIRE(linearMeshData.tetMaterials == rawMeshData.tetMaterials);

  // triangular elements and materials
  REQUIRE(linearMeshData.triNodes.size() == rawMeshData.triNodes.size());
  REQUIRE(linearMeshData.triNodes == rawMeshData.triNodes);
  REQUIRE(linearMeshData.triMaterials.size() == rawMeshData.triMaterials.size());
  REQUIRE(linearMeshData.triMaterials == rawMeshData.triMaterials);
}
