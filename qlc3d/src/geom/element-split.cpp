#include "mesh/element-split.h"
#include <util/exception.h>
#include <util/logging.h>
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

bool recombineLinearisedMeshToQuadratic(RawMeshData &meshData) {

  const int elementOrder = meshData.getElementOrder();
  if (elementOrder != 1) {
    RUNTIME_ERROR("Can only recombine linearised mesh to quadratic if the mesh is linearised. Got mesh with elementOrder=" + std::to_string(elementOrder));
  }

  const size_t numTetsIn = meshData.tetNodes.size() / 4;
  const size_t numTrisIn = meshData.triNodes.size() / 3;
  const size_t numTetsOut = numTetsIn / 8;
  const size_t numTrisOut = numTrisIn / 4;

  Log::info("Attempting to recombine first order mesh with {} tetrahedra {} triangles to second order mesh with {} tetrahedra {} triangles.",
            numTetsIn, numTrisIn, numTetsOut, numTrisOut);

  if (numTetsIn % 8 != 0) {
    Log::info("Number of tetrahedra is not a multiple of 8, cannot recombine.");
    return false;
  }
  if (numTrisIn % 4 != 0) {
    Log::info("Number of triangles is not a multiple of 4, cannot recombine.");
    return false;
  }

  std::vector<unsigned int> tetMaterialsOut;
  tetMaterialsOut.reserve(numTetsOut);

  for (size_t i = 0; i < numTetsIn; i+= 8) {
    // materials for all 8 tets must match
    unsigned int tetMat = meshData.tetMaterials[i];
    for (size_t j = 0; j < 8; j++) {
      if (meshData.tetMaterials[i + j] != tetMat) {
        Log::info("Tetrahedron {} has different material than {}, cannot recombine.", i + j, i);
        return false;
      }
    }
    tetMaterialsOut.push_back(tetMat);
  }

  // Recombine tet elements
  auto &t = meshData.tetNodes;
  std::vector<unsigned int> tetNodesOut;
  tetNodesOut.reserve(numTetsOut * 10);

  try {
    for (size_t i = 0; i < numTetsOut; i++) {
      std::vector<std::vector<unsigned int>> linearTets;
      auto start = i * 8 * 4; // 8 tets, 4 nodes each
      linearTets.emplace_back(std::vector<unsigned int>({t[start + 0], t[start + 1], t[start + 2], t[start + 3]}));
      linearTets.emplace_back(std::vector<unsigned int>({t[start + 4], t[start + 5], t[start + 6], t[start + 7]}));
      linearTets.emplace_back(std::vector<unsigned int>({t[start + 8], t[start + 9], t[start + 10], t[start + 11]}));
      linearTets.emplace_back(std::vector<unsigned int>({t[start + 12], t[start + 13], t[start + 14], t[start + 15]}));
      linearTets.emplace_back(std::vector<unsigned int>({t[start + 16], t[start + 17], t[start + 18], t[start + 19]}));
      linearTets.emplace_back(std::vector<unsigned int>({t[start + 20], t[start + 21], t[start + 22], t[start + 23]}));
      linearTets.emplace_back(std::vector<unsigned int>({t[start + 24], t[start + 25], t[start + 26], t[start + 27]}));
      linearTets.emplace_back(std::vector<unsigned int>({t[start + 28], t[start + 29], t[start + 30], t[start + 31]}));

      auto recombinedNodes = recombineLinearTetsToQuadratic(linearTets);
      tetNodesOut.insert(tetNodesOut.end(), recombinedNodes.begin(), recombinedNodes.end());
    }
  } catch (ElementSplitCombineException &e) {
    Log::info("Can not combine quadratic tetrahedron from input");
    return false;
  }

  // Recombine triangle materials
  std::vector<unsigned int> triMaterialsOut;
  triMaterialsOut.reserve(numTrisOut);
  for (size_t i = 0; i < numTrisIn; i+= 4) {
    unsigned int triMat = meshData.triMaterials[i];
    for (size_t j = 0; j < 4; j++) { // materials for all 4 source tris must match
      if (meshData.triMaterials[i + j] != triMat) {
        Log::info("Triangle {} has different material than {}, cannot recombine.", i + j, i);
        return false;
      }
    }
    triMaterialsOut.push_back(triMat);
  }

  // Recombine triangle elements
  auto &e = meshData.triNodes;
  std::vector<unsigned int> triNodesOut;
  triNodesOut.reserve(numTrisOut * 6);
  try {
    for (size_t i = 0; i < numTrisOut; i++) {
      std::vector<std::vector<unsigned int>> linearTris;
      auto start = i * 4 * 3; // 4 tris, 3 nodes each
      linearTris.emplace_back(std::vector<unsigned int>({e[start + 0], e[start + 1], e[start + 2]}));
      linearTris.emplace_back(std::vector<unsigned int>({e[start + 3], e[start + 4], e[start + 5]}));
      linearTris.emplace_back(std::vector<unsigned int>({e[start + 6], e[start + 7], e[start + 8]}));
      linearTris.emplace_back(std::vector<unsigned int>({e[start + 9], e[start + 10], e[start + 11]}));

      auto recombinedNodes = recombineLinearTrianglesToQuadratic(linearTris);
      triNodesOut.insert(triNodesOut.end(), recombinedNodes.begin(), recombinedNodes.end());
    }
  } catch (ElementSplitCombineException &e) {
    Log::info("Can not combine quadratic triangle from input");
    return false;
  }

  // Update mesh data with quadratic element data
  meshData.setElementOrder(2);
  meshData.tetMaterials = std::move(tetMaterialsOut);
  meshData.tetNodes = std::move(tetNodesOut);
  meshData.triMaterials = std::move(triMaterialsOut);
  meshData.triNodes = std::move(triNodesOut);
  return true;
}
