#include <test-util.h>
#include <io/gmsh-read.h>

#include <geometry.h>
#include <material_numbers.h>

#include "catch.h"
#include "mesh/element-split-convert.h"

namespace {

  // Return labeled coordinates for the quadratic tetrahedron used in the tests.
  // Coordinates are in Gmsh TET10 ordering: [0]=A, [1]=B, [2]=C, [3]=D,
  // [4]=AB, [5]=BC, [6]=AC, [7]=AD, [8]=CD, [9]=BD
  static std::vector<Vec3> makeQuadraticTestCoordinates() {
    auto a = Vec3(0, 0, 0);
    auto b = Vec3(1, 0, 0);
    auto c = Vec3(0, 1, 0);
    auto d = Vec3(0, 0, 1);

    return std::vector<Vec3> {
      a, b, c, d,           // corners: A=0, B=1, C=2, D=3
      Vec3::mean(a, b),     // [4] = AB
      Vec3::mean(b, c), // [5] = BC
      Vec3::mean(a, c),   // [6] = AC
      Vec3::mean(a, d),   // [7] = AD
      Vec3::mean(c, d), // [8] = CD
      Vec3::mean(b, d)  // [9] = BD
    };
  }


RawMeshData makeQuadraticTestRawMeshData() {
  std::vector<Vec3> points = makeQuadraticTestCoordinates();
  // Gmsh TET10: A,B,C,D, AB=4, BC=5, AC=6, AD=7, CD=8, BD=9
  // Gmsh TRI6 faces of the tet: corners + AB=edge(0,1), BC=edge(1,2), AC=edge(0,2)
  //   Face ABC: {0,1,2, AB=4, BC=5, AC=6}
  //   Face ABD: {0,1,3, AB=4, BD=9, AD=7}
  //   Face ACD: {0,2,3, AC=6, CD=8, AD=7}
  //   Face BCD: {1,2,3, BC=5, CD=8, BD=9}
  return RawMeshData(2, std::move(points),
                     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {MAT_DOMAIN1},
                     {0, 1, 2, 4, 5, 6,
                      0, 1, 3, 4, 9, 7,
                      0, 2, 3, 6, 8, 7,
                      1, 2, 3, 5, 8, 9},
                     {MAT_FIXLC1, MAT_FIXLC1, MAT_FIXLC1, MAT_FIXLC1});
}

Geometry makeQuadraticTestGeometry() {
  return Geometry::fromRawMeshData(makeQuadraticTestRawMeshData());
}

} // namespace


double determinant(const Coordinates &coords, const vector<unsigned int> &tet) {
  return det3D(
    coords.getPoint(tet[0]),
    coords.getPoint(tet[1]),
    coords.getPoint(tet[2]),
    coords.getPoint(tet[3]));
}


TEST_CASE("Split and recombine tets") {
  // Gmsh TET10 mid-edge node indices (assigned sequentially after corners 0-3):
  // [4]=AB, [5]=BC, [6]=AC, [7]=AD, [8]=CD, [9]=BD
  unsigned int A = 0;
  unsigned int B = 1;
  unsigned int C = 2;
  unsigned int D = 3;

  unsigned int ab = 4;  // [4] = AB mid-edge
  unsigned int bc = 5;  // [5] = BC mid-edge
  unsigned int ac = 6;  // [6] = AC mid-edge
  unsigned int ad = 7;  // [7] = AD mid-edge
  unsigned int cd = 8;  // [8] = CD mid-edge
  unsigned int bd = 9;  // [9] = BD mid-edge

  //const auto coords = makeQuadraticTestCoordinates();
  Coordinates coords = Coordinates(makeQuadraticTestCoordinates());
  SECTION("Split quadratic tetrahedron to 8 linear tets") {
    // Input is in Gmsh TET10 ordering: A,B,C,D, AB,BC,AC,AD,CD,BD
    std::vector<unsigned int> quadraticTetrahedron = {A, B, C, D, ab, bc, ac, ad, cd, bd};

    // ACT
    std::vector<std::vector<unsigned int>> linearTets = splitQuadraticTetrahedronToLinear(quadraticTetrahedron, coords);

    // ASSERT
    REQUIRE(linearTets.size() == 8);

    // corner tets
    REQUIRE(linearTets[0] == std::vector<unsigned int>({A, ab, ac, ad}));
    REQUIRE(linearTets[1] == std::vector<unsigned int>({B, bc, ab, bd}));
    REQUIRE(linearTets[2] == std::vector<unsigned int>({C, cd, ac, bc}));
    REQUIRE(linearTets[3] == std::vector<unsigned int>({D, cd, bd, ad}));

    // face tets
    REQUIRE(linearTets[4] == std::vector<unsigned int>({bd, ab, ad, ac})); // A-face
    REQUIRE(linearTets[5] == std::vector<unsigned int>({ac, bc, bd, ab})); // B-face
    REQUIRE(linearTets[6] == std::vector<unsigned int>({bd, cd, bc, ac})); // C-face
    REQUIRE(linearTets[7] == std::vector<unsigned int>({ac, cd, ad, bd})); // D-face

    // TODO: check that all the linear tets have a positive jacobian determinant. This is where negative tets may be introduced to the system
    for (int i = 0; i < 8; i++) {
      double jDet = determinant(coords, linearTets[i]);
      REQUIRE(jDet > 0);
    }
  }

  SECTION("Recombine 8 linear tets to quadratic tetrahedron") {
    // Linear tets derived from the split above (using Gmsh TET10 index assignments)
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
    // Gmsh TET10 mid-edge positions: [4]=AB, [5]=BC, [6]=AC, [7]=AD, [8]=CD, [9]=BD
    REQUIRE(recombinedTet[4] == ab);
    REQUIRE(recombinedTet[5] == bc);
    REQUIRE(recombinedTet[6] == ac);
    REQUIRE(recombinedTet[7] == ad);
    REQUIRE(recombinedTet[8] == cd);
    REQUIRE(recombinedTet[9] == bd);

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
  constexpr unsigned int A = 0, B = 1, C = 2, D = 3; // corner node indices
  // create raw mesh data for a single linear tetrahedron and triangle
  std::vector<Vec3> pIn = {
    Vec3(0, 0, 0),
    Vec3(1, 0, 0),
    Vec3(0, 1, 0),
    Vec3(0, 0, 1)
  };

  std::vector<unsigned int> linTet = {A, B, C, D};
  std::vector<unsigned int> linTri = {A, B, C};

  RawMeshData meshData(1, pIn, linTet, {99}, linTri, {100});

  // ACT
  convertLinearMeshDataToQuadratic(meshData);

  // ASSERT
  REQUIRE(meshData.getElementOrder() == 2);
  REQUIRE(meshData.points.size() == 10);
  REQUIRE(meshData.tetNodes.size() == 10);
  REQUIRE(meshData.tetMaterials.size() == 1);
  REQUIRE(meshData.triNodes.size() == 6);
  REQUIRE(meshData.triMaterials.size() == 1);

  // New node indices are assigned by createNewNodeNumbersMapping in order: AB=4, AC=5, AD=6, BC=7, BD=8, CD=9
  // createNewTetElementNodes then places them in Gmsh TET10 order: A,B,C,D, AB,BC,AC,AD,CD,BD
  // So tetNodes should be {A,B,C,D, 4=AB, 7=BC, 5=AC, 6=AD, 9=CD, 8=BD}
  REQUIRE(meshData.tetNodes == std::vector<unsigned int>{A, B, C, D, 4, 7, 5, 6, 9, 8});

  // Check Gmsh TRI6 node ordering: A,B,C, AB,BC,AC
  // For triangle ABC: AB=4, BC=7, AC=5
  REQUIRE(meshData.triNodes == std::vector<unsigned int>{A, B, C, 4, 7, 5});

  // Check that original points are unchanged
  constexpr double eps = 1e-9;
  REQUIRE(meshData.points[0].equals(pIn[0], eps));
  REQUIRE(meshData.points[1].equals(pIn[1], eps));
  REQUIRE(meshData.points[2].equals(pIn[2], eps));
  REQUIRE(meshData.points[3].equals(pIn[3], eps));

  // Node global indices: 4=AB, 5=AC, 6=AD, 7=BC, 8=BD, 9=CD (assigned by createNewNodeNumbersMapping)
  REQUIRE(meshData.points[4].equals(Vec3::mean(pIn[0], pIn[1]), eps)); // index 4 = AB midpoint
  REQUIRE(meshData.points[5].equals(Vec3::mean(pIn[0], pIn[2]), eps)); // index 5 = AC midpoint
  REQUIRE(meshData.points[6].equals(Vec3::mean(pIn[0], pIn[3]), eps)); // index 6 = AD midpoint
  REQUIRE(meshData.points[7].equals(Vec3::mean(pIn[1], pIn[2]), eps)); // index 7 = BC midpoint
  REQUIRE(meshData.points[8].equals(Vec3::mean(pIn[1], pIn[3]), eps)); // index 8 = BD midpoint
  REQUIRE(meshData.points[9].equals(Vec3::mean(pIn[2], pIn[3]), eps)); // index 9 = CD midpoint
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

  constexpr unsigned int A = 0, B = 1, C = 2, D = 3, E = 4;
  RawMeshData meshData(1, points,
    {A, B, C, D,
                A, B, C, E}, {99, 99},
{A, B, C}, {100});

  // ACT
  convertLinearMeshDataToQuadratic(meshData);

  // ASSERT
  REQUIRE(meshData.getElementOrder() == 2);
  REQUIRE(meshData.points.size() == 14);

  // createNewNodeNumbersMapping assigns new indices in traversal order:
  // Processing tet1(A,B,C,D): AB=5, AC=6, AD=7, BC=8, BD=9, CD=10
  // Processing tet2(A,B,C,E): AB=5(shared), AC=6(shared), AE=11, BC=8(shared), BE=12, CE=13
  // createNewTetElementNodes places them in Gmsh TET10 order: A,B,C,D, AB,BC,AC,AD,CD,BD
  auto expectedTetNodes = std::vector<unsigned int>(
    {
      A, B, C, D, 5, 8, 6, 7, 10, 9,   // Tet1: AB=5, BC=8, AC=6, AD=7, CD=10, BD=9
      A, B, C, E, 5, 8, 6, 11, 13, 12  // Tet2: AB=5, BC=8, AC=6, AE=11, CE=13, BE=12
    });

  REQUIRE(meshData.tetNodes == expectedTetNodes); // two 2nd order tets that share a face

  // Gmsh TRI6 for face ABC: A,B,C, AB=5, BC=8, AC=6
  auto expectedTriNodes = std::vector<unsigned int>({A, B, C, 5, 8, 6});
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
  RawMeshData workingMeshData(1, {}, {}, {}, {}, {});

  splitQuadraticGeometryToLinear(geom, workingMeshData);

  REQUIRE(recombineLinearisedMeshToQuadratic(workingMeshData) == true);

  // ASSERT: linearMeshData should match the original quadratic mesh data
  // ARRANGE/ACT above: rawMeshData -> Geometry -> split -> recombine


  auto geomRecovered = Geometry::fromRawMeshData(workingMeshData);

  // element order
  REQUIRE(workingMeshData.getElementOrder() == rawMeshData.getElementOrder());

  // points: same count and coordinates (compare via distance to avoid relying on Vec3 internals)
  REQUIRE(workingMeshData.points.size() == rawMeshData.points.size());
  for (size_t i = 0; i < rawMeshData.points.size(); ++i) {
    REQUIRE(workingMeshData.points[i].distance(rawMeshData.points[i]) == Approx(0.0).margin(1e-9));
  }

  // tetrahedral elements and materials
  REQUIRE(workingMeshData.tetNodes.size() == rawMeshData.tetNodes.size());
  if (workingMeshData.tetNodes != rawMeshData.tetNodes) {
    // Log detailed mismatch information to help debug swapped node ordering.
    Log::info("Tet nodes differ: linearMeshData.tetNodes.size()={}, rawMeshData.tetNodes.size()={}", workingMeshData.tetNodes.size(), rawMeshData.tetNodes.size());
    const size_t nodesPerTet = 10;
    const size_t nTets = rawMeshData.tetNodes.size() / nodesPerTet;
    for (size_t t = 0; t < nTets; ++t) {
      size_t base = t * nodesPerTet;
      bool equal = true;
      for (size_t k = 0; k < nodesPerTet; ++k) {
        if (workingMeshData.tetNodes[base + k] != rawMeshData.tetNodes[base + k]) { equal = false; break; }
      }
      if (!equal) {
        std::string lstr;
        std::string rstr;
        for (size_t k = 0; k < nodesPerTet; ++k) {
          lstr += std::to_string(workingMeshData.tetNodes[base + k]);
          rstr += std::to_string(rawMeshData.tetNodes[base + k]);
          if (k + 1 < nodesPerTet) { lstr += ", "; rstr += ", "; }
        }
        Log::info("Tet {}: linear = [{}]", t, lstr);
        Log::info("Tet {}: raw    = [{}]", t, rstr);
      }
    }
  }
  REQUIRE(workingMeshData.tetNodes == rawMeshData.tetNodes);
  REQUIRE(workingMeshData.tetMaterials.size() == rawMeshData.tetMaterials.size());
  REQUIRE(workingMeshData.tetMaterials == rawMeshData.tetMaterials);

  // triangular elements and materials
  REQUIRE(workingMeshData.triNodes.size() == rawMeshData.triNodes.size());
  REQUIRE(workingMeshData.triNodes == rawMeshData.triNodes);
  REQUIRE(workingMeshData.triMaterials.size() == rawMeshData.triMaterials.size());
  REQUIRE(workingMeshData.triMaterials == rawMeshData.triMaterials);
}
