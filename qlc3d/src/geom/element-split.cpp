#include <geom/element-split.h>
#include <util/exception.h>
#include <unordered_set>
#include <cassert>

std::vector<std::vector<unsigned int>> splitQuadraticTetrahedronToLinear(const std::vector<unsigned int> &quadraticTetrahedron) {
  if (quadraticTetrahedron.size() != 10) {
    RUNTIME_ERROR("Quadratic tetrahedron must have 10 nodes, got " + std::to_string(quadraticTetrahedron.size()));
  }

  std::vector<std::vector<unsigned int>> linearTets;
  linearTets.reserve(8);

  // node numbering w.r.t fig 6.1 (b) in Eero's thesis
  unsigned int A = quadraticTetrahedron[0];
  unsigned int B = quadraticTetrahedron[1];
  unsigned int C = quadraticTetrahedron[2];
  unsigned int D = quadraticTetrahedron[3];

  unsigned int ab = quadraticTetrahedron[4];
  unsigned int bc = quadraticTetrahedron[5];
  unsigned int ac = quadraticTetrahedron[6];
  unsigned int ad = quadraticTetrahedron[7];
  unsigned int bd = quadraticTetrahedron[8];
  unsigned int cd = quadraticTetrahedron[9];

  linearTets.push_back({A, ab, ac, ad}); // A-corner
  linearTets.push_back({B, bc, ab, bd}); // B-corner
  linearTets.push_back({C, cd, ac, bc}); // C-corner
  linearTets.push_back({D, cd, bd, ad}); // D-corner

  linearTets.push_back({bd, ab, ac, ad}); // A-face
  linearTets.push_back({ac, bc, ab, bd}); // B-face
  linearTets.push_back({bd, cd, ac, bc}); // C-face
  linearTets.push_back({ac, cd, bd, ad}); // D-face


  return linearTets;
}

std::vector<unsigned int> recombineLinearTetsToQuadratic(const std::vector<std::vector<unsigned int>> &linearTets) {
  if (linearTets.size() != 8) {
    RUNTIME_ERROR("Expected 8 linear tets, got " + std::to_string(linearTets.size()));
  }

  unsigned int A = linearTets[0][0];
  unsigned int B = linearTets[1][0];
  unsigned int C = linearTets[2][0];
  unsigned int D = linearTets[3][0];

  unsigned int ab = linearTets[0][1];
  unsigned int ac = linearTets[0][2];
  unsigned int ad = linearTets[0][3];
  unsigned int bc = linearTets[1][1];
  unsigned int bd = linearTets[1][3];
  unsigned int cd = linearTets[2][1];

#ifndef NDEBUG
  // check that all 10 nodes are unique
  std::unordered_set<unsigned int> nodes;
  nodes.insert(A);
  nodes.insert(B);
  nodes.insert(C);
  nodes.insert(D);
  nodes.insert(ab);
  nodes.insert(ac);
  nodes.insert(ad);
  nodes.insert(bc);
  nodes.insert(bd);
  nodes.insert(cd);
  if (nodes.size() != 10) {
    throw ElementSplitCombineException("Recombined tetrahedron has duplicate nodes");
  }
  // could add additional checks here to verify every node value is as expected. Should it run also in
  // release mode?
#endif

  return {A, B, C, D, ab, bc, ac, ad, bd, cd};
}

std::vector<std::vector<unsigned int>> splitQuadraticTriangleToLinear(const std::vector<unsigned int> &quadraticTriangle) {
  if (quadraticTriangle.size() != 6) {
    RUNTIME_ERROR("Quadratic triangle must have 6 nodes, got " + std::to_string(quadraticTriangle.size()));
  }

  std::vector<std::vector<unsigned int>> linearTriangles;
  linearTriangles.reserve(4);

  // node numbering w.r.t fig 6.1 (b) in Eero's thesis
  unsigned int A = quadraticTriangle[0];
  unsigned int B = quadraticTriangle[1];
  unsigned int C = quadraticTriangle[2];

  unsigned int ab = quadraticTriangle[3];
  unsigned int bc = quadraticTriangle[4];
  unsigned int ac = quadraticTriangle[5];

  linearTriangles.push_back({A, ab, ac});
  linearTriangles.push_back({B, bc, ab});
  linearTriangles.push_back({C, ac, bc});
  linearTriangles.push_back({ab, bc, ac});

  return linearTriangles;
}

std::vector<unsigned int> recombineLinearTrianglesToQuadratic(const std::vector<std::vector<unsigned int>> &linearTriangles) {
  if (linearTriangles.size() != 4) {
    RUNTIME_ERROR("Expected 4 linear triangles, got " + std::to_string(linearTriangles.size()));
  }

  unsigned int A = linearTriangles[0][0];
  unsigned int B = linearTriangles[1][0];
  unsigned int C = linearTriangles[2][0];

  unsigned int ab = linearTriangles[0][1];
  unsigned int bc = linearTriangles[1][1];
  unsigned int ac = linearTriangles[2][1];

  if (linearTriangles[0][2] != ac &&
    linearTriangles[1][2] != ab &&
    linearTriangles[2][2] != bc &&
    linearTriangles[3][0] != ab &&
    linearTriangles[3][1] != bc &&
    linearTriangles[3][2] != ac) {
    throw ElementSplitCombineException("Third node of first triangle must be ac");
  }

  return {A, B, C, ab, bc, ac};
}