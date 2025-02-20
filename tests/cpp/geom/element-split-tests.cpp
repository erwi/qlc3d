#include <catch.h>
#include <geom/element-split.h>

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

    std::vector<unsigned int> recombinedTet = recombineLinearTetsToQuadratic(linearTets);
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