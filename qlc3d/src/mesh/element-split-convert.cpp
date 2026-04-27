#include "mesh/element-split-convert.h"

#include <cassert>
#include <geometry.h>
#include <geom/coordinates.h>
#include <qlc3d.h>
#include <spamtrix_matrixmaker.hpp>
#include <array>

#include "util/exception.h"
#include "util/logging.h"
#include <unordered_set>

namespace {
constexpr double MIDPOINT_TOLERANCE_RATIO = 1e-3;

void validateAndSnapQuadraticTetrahedron(std::vector<Vec3> &points, const unsigned int *tetNodes, size_t tetIndex) {
  const auto validateEdge = [&](unsigned int cornerA, unsigned int cornerB, unsigned int midNode, const char *edgeName) {
    const Vec3 &a = points[cornerA];
    const Vec3 &b = points[cornerB];
    const Vec3 expected = Vec3::mean(a, b);
    Vec3 &actual = points[midNode];

    const double edgeLength = a.distance(b);
    if (edgeLength <= 0.0) {
      RUNTIME_ERROR(fmt::format("Quadratic tetrahedron {} has zero-length edge {} (nodes {} and {}).", tetIndex, edgeName, cornerA, cornerB));
    }

    const double deviation = actual.distance(expected);
    const double tolerance = edgeLength * MIDPOINT_TOLERANCE_RATIO;
    if (deviation > tolerance) {
      RUNTIME_ERROR(fmt::format("Quadratic tetrahedron {} mid-edge node {} on edge {} is too far from midpoint: deviation = {}, tolerance = {}, actual = {}, expected = {}.",
                                tetIndex, midNode, edgeName, deviation, tolerance, actual, expected));
    }

    actual = expected;
  };

  // Gmsh TET10 node ordering: [0-3]=corners A,B,C,D; [4]=AB, [5]=BC, [6]=AC, [7]=AD, [8]=CD, [9]=BD
  validateEdge(tetNodes[0], tetNodes[1], tetNodes[4], "AB");
  validateEdge(tetNodes[1], tetNodes[2], tetNodes[5], "BC");
  validateEdge(tetNodes[0], tetNodes[2], tetNodes[6], "AC");
  validateEdge(tetNodes[0], tetNodes[3], tetNodes[7], "AD");
  validateEdge(tetNodes[2], tetNodes[3], tetNodes[8], "CD");
  validateEdge(tetNodes[1], tetNodes[3], tetNodes[9], "BD");
}
} // namespace

void splitQuadraticGeometryToLinear(const Geometry &quadGeom, RawMeshData &outMeshData) {
  const Mesh &tetrahedra = quadGeom.getTetrahedra();
  const Mesh &triangles = quadGeom.getTriangles();

  if (tetrahedra.getElementType() != ElementType::QUADRATIC_TETRAHEDRON) {
    RUNTIME_ERROR(fmt::format("splitQuadraticGeometryToLinear expects quadratic tetrahedra, got {}.", tetrahedra.getElementType()));
  }
  if (triangles.getElementType() != ElementType::QUADRATIC_TRIANGLE) {
    RUNTIME_ERROR(fmt::format("splitQuadraticGeometryToLinear expects quadratic triangles, got {}.", triangles.getElementType()));
  }

  std::vector<Vec3> points;
  points.reserve(quadGeom.getnp());
  for (idx i = 0; i < quadGeom.getnp(); ++i) {
    points.push_back(quadGeom.getCoordinates().getPoint(i));
  }

  std::vector<idx> linearTetNodes;
  std::vector<idx> linearTetMaterials;
  linearTetNodes.reserve(tetrahedra.getnElements() * 8 * 4);
  linearTetMaterials.reserve(tetrahedra.getnElements() * 8);
  for (idx i = 0; i < tetrahedra.getnElements(); ++i) {
    std::vector<unsigned int> quadraticTet(tetrahedra.getnNodes());
    for (idx n = 0; n < tetrahedra.getnNodes(); ++n) {
      quadraticTet[n] = tetrahedra.getNode(i, n);
    }

    const auto linearTets = splitQuadraticTetrahedronToLinear(quadraticTet);
    for (const auto &linearTet : linearTets) {
      linearTetNodes.insert(linearTetNodes.end(), linearTet.begin(), linearTet.end());
      linearTetMaterials.push_back(tetrahedra.getMaterialNumber(i));
    }
  }

  std::vector<idx> linearTriNodes;
  std::vector<idx> linearTriMaterials;
  linearTriNodes.reserve(triangles.getnElements() * 4 * 3);
  linearTriMaterials.reserve(triangles.getnElements() * 4);
  for (idx i = 0; i < triangles.getnElements(); ++i) {
    std::vector<unsigned int> quadraticTri(triangles.getnNodes());
    for (idx n = 0; n < triangles.getnNodes(); ++n) {
      quadraticTri[n] = triangles.getNode(i, n);
    }

    const auto linearTris = splitQuadraticTriangleToLinear(quadraticTri);
    for (const auto &linearTri : linearTris) {
      linearTriNodes.insert(linearTriNodes.end(), linearTri.begin(), linearTri.end());
      linearTriMaterials.push_back(triangles.getMaterialNumber(i));
    }
  }

  outMeshData = RawMeshData(1, std::move(points), std::move(linearTetNodes), std::move(linearTetMaterials),
                            std::move(linearTriNodes), std::move(linearTriMaterials));
}

void validateAndSnapQuadraticTetrahedra(RawMeshData &meshData) {
  if (meshData.getElementOrder() != 2) {
    return;
  }

  const unsigned int nodesPerTetQuadratic = 10;
  const size_t numTetrahedra = meshData.tetNodes.size() / nodesPerTetQuadratic;
  for (size_t tetCounter = 0; tetCounter < numTetrahedra; tetCounter++) {
    validateAndSnapQuadraticTetrahedron(meshData.points, &meshData.tetNodes[tetCounter * nodesPerTetQuadratic], tetCounter);
  }
}

std::vector<std::vector<unsigned int>> splitQuadraticTetrahedronToLinear(const std::vector<unsigned int> &quadraticTetrahedron) {
  if (quadraticTetrahedron.size() != 10) {
    RUNTIME_ERROR("Quadratic tetrahedron must have 10 nodes, got " + std::to_string(quadraticTetrahedron.size()));
  }

  std::vector<std::vector<unsigned int>> linearTets;
  linearTets.reserve(8);

  // Gmsh TET10 ordering: [4]=AB, [5]=BC, [6]=AC, [7]=AD, [8]=CD, [9]=BD
  unsigned int A = quadraticTetrahedron[0];
  unsigned int B = quadraticTetrahedron[1];
  unsigned int C = quadraticTetrahedron[2];
  unsigned int D = quadraticTetrahedron[3];

  unsigned int ab = quadraticTetrahedron[4];  // AB mid-edge
  unsigned int bc = quadraticTetrahedron[5];  // BC mid-edge
  unsigned int ac = quadraticTetrahedron[6];  // AC mid-edge
  unsigned int ad = quadraticTetrahedron[7];  // AD mid-edge
  unsigned int cd = quadraticTetrahedron[8];  // CD mid-edge
  unsigned int bd = quadraticTetrahedron[9];  // BD mid-edge

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

std::vector<unsigned int> recombineLinearTetsToQuadratic(const std::vector<std::vector<unsigned int>> &linearTets,
  const std::vector<Vec3> &points) {
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

  // load coordinates for the 10 nodes and make sure the order is correct
  const double eps = 1e-9;
  auto pA = points[A];
  auto pB = points[B];
  auto pC = points[C];
  auto pD = points[D];


  double det = det3D(pA, pB, pC, pD);
  if (det < 0) {
    Log::info("swapping {} and {}", A, B);
    std::swap(C, D);
    pC = points[C];
    pD = points[D];
  }

  std::vector<std::pair<unsigned int, Vec3>> midNodes = {
    {ab, points[ab]},
    {bc, points[bc]},
    {ac, points[ac]},
    {ad, points[ad]},
    {bd, points[bd]},
    {cd, points[cd]}
  };

  // lambda which extracts the point index of the pair whose point is closest to expected point. This ensures
  auto findClosest = [](const std::vector<std::pair<unsigned int, Vec3>> &midNodes, const Vec3 &expected, double eps) {
     unsigned int closestNode = 0;
     double closestDistance = std::numeric_limits<double>::max();
     for (const auto &pair : midNodes) {
       double distance = pair.second.distance(expected);
       if (distance < closestDistance) {
         closestDistance = distance;
         closestNode = pair.first;
       }
     }
     if (closestDistance > eps) {
       throw ElementSplitCombineException("No mid-edge node between A and B is close enough to expected position");
     }
     return closestNode;
  };

  ab = findClosest(midNodes, Vec3::mean(pA, pB), eps);
  ac = findClosest(midNodes, Vec3::mean(pA, pC), eps);
  ad = findClosest(midNodes, Vec3::mean(pA, pD), eps);
  bc = findClosest(midNodes, Vec3::mean(pB, pC), eps);
  bd = findClosest(midNodes, Vec3::mean(pB, pD), eps);
  cd = findClosest(midNodes, Vec3::mean(pC, pD), eps);

  // sanity check, check that all 10 nodes are unique
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

  // Using Gmsh TET10 node ordering: A, B, C, D, AB, BC, AC, AD, CD, BD
  return {A, B, C, D, ab, bc, ac, ad, cd, bd};
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

  const unsigned int elementOrder = meshData.getElementOrder();
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

      auto recombinedNodes = recombineLinearTetsToQuadratic(linearTets, meshData.points);
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

  const unsigned int nodesPerTetQuadratic = 10;
  for (size_t tetCounter = 0; tetCounter < numTetsOut; tetCounter++) {
    validateAndSnapQuadraticTetrahedron(meshData.points, &meshData.tetNodes[tetCounter * nodesPerTetQuadratic], tetCounter);
  }
  return true;
}

// <editor-fold desc="convertLinearMeshDataToQuadratic">

/**
 * Calculate node indices for new nodes that are added when converting linear tetrahedra to quadratic.
 */
std::pair<SpaMtrix::IRCMatrix, unsigned int> createNewNodeNumbersMapping(std::vector<unsigned int> &tetNodes, unsigned int maxNodeNumber) {
  const unsigned int nodesPerTet = 4;
  const unsigned int numTetrahedra = tetNodes.size() / nodesPerTet;

  // Create sparse matrix for mapping tet corned node numbers to mid edge node numbers
  SpaMtrix::MatrixMaker tetNodeMaker(maxNodeNumber + 1, maxNodeNumber + 1);
  for (unsigned int tetCounter = 0; tetCounter < numTetrahedra; tetCounter++) {
    const unsigned int iBegin = tetCounter * nodesPerTet;
    const unsigned int A = tetNodes[iBegin + 0];
    const unsigned int B = tetNodes[iBegin + 1];
    const unsigned int C = tetNodes[iBegin + 2];
    const unsigned int D = tetNodes[iBegin + 3];

    tetNodeMaker.addNonZero(A, B);
    tetNodeMaker.addNonZero(B, A);

    tetNodeMaker.addNonZero(B, C);
    tetNodeMaker.addNonZero(C, B);

    tetNodeMaker.addNonZero(A, C);
    tetNodeMaker.addNonZero(C, A);

    tetNodeMaker.addNonZero(A, D);
    tetNodeMaker.addNonZero(D, A);

    tetNodeMaker.addNonZero(B, D);
    tetNodeMaker.addNonZero(D, B);

    tetNodeMaker.addNonZero(C, D);
    tetNodeMaker.addNonZero(D, C);
  }

  auto tetMap = tetNodeMaker.getIRCMatrix();

  unsigned int newNodeIndex = maxNodeNumber;
  for (unsigned int tetCounter = 0; tetCounter < numTetrahedra; tetCounter++) {
    const unsigned int iBegin = tetCounter * nodesPerTet;
    const unsigned int A = tetNodes[iBegin + 0];
    const unsigned int B = tetNodes[iBegin + 1];
    const unsigned int C = tetNodes[iBegin + 2];
    const unsigned int D = tetNodes[iBegin + 3];

    if (tetMap.isNonZero(A, B) && tetMap.getValue(A, B) == 0.) {
      newNodeIndex++;
      (*tetMap.getValuePtr(A, B)) = newNodeIndex;
      (*tetMap.getValuePtr(B, A)) = newNodeIndex;
    }

    if (tetMap.isNonZero(A, C) && tetMap.getValue(A, C) == 0.) {
      newNodeIndex++;
      (*tetMap.getValuePtr(A, C)) = newNodeIndex;
      (*tetMap.getValuePtr(C, A)) = newNodeIndex;
    }

    if (tetMap.isNonZero(A, D) && tetMap.getValue(A, D) == 0.) {
      newNodeIndex++;
      (*tetMap.getValuePtr(A, D)) = newNodeIndex;
      (*tetMap.getValuePtr(D, A)) = newNodeIndex;
    }

    if (tetMap.isNonZero(B, C) && tetMap.getValue(B, C) == 0.) {
      newNodeIndex++;
      (*tetMap.getValuePtr(B, C)) = newNodeIndex;
      (*tetMap.getValuePtr(C, B)) = newNodeIndex;
    }

    if (tetMap.isNonZero(B, D) && tetMap.getValue(B, D) == 0.) {
      newNodeIndex++;
      (*tetMap.getValuePtr(B, D)) = newNodeIndex;
      (*tetMap.getValuePtr(D, B)) = newNodeIndex;
    }

    if (tetMap.isNonZero(C, D) && tetMap.getValue(C, D) == 0.) {
      newNodeIndex++;
      (*tetMap.getValuePtr(C, D)) = newNodeIndex;
      (*tetMap.getValuePtr(D, C)) = newNodeIndex;
    }
  }

  return {tetMap, newNodeIndex};
}

/**
 * Calculate the locations of new nodes and appends them to the coordinates vector.
 */
void appendNewNodeCoordinates(std::vector<Vec3> &coordinates,
                              const SpaMtrix::IRCMatrix &newNodeNumbers,
                              const std::vector<unsigned int> &tetNodes,
                              unsigned int numTetrahedra,
                              unsigned int newNodeCount) {
  const unsigned int nodesPerTet = 4;
  assert(tetNodes.size() == numTetrahedra * nodesPerTet);

  coordinates.resize(newNodeCount, Vec3(0, 0, 0));

  auto addCoordinate = [&](unsigned int n1, unsigned int n2) {
    auto n12 = (unsigned int) newNodeNumbers.getValue(n1, n2);
    coordinates[n12] = Vec3::mean(coordinates[n1], coordinates[n2]);
  };


  for (unsigned int tetCounter = 0; tetCounter < numTetrahedra; tetCounter++) {
    const unsigned int A = tetNodes[tetCounter * nodesPerTet + 0];
    const unsigned int B = tetNodes[tetCounter * nodesPerTet + 1];
    const unsigned int C = tetNodes[tetCounter * nodesPerTet + 2];
    const unsigned int D = tetNodes[tetCounter * nodesPerTet + 3];

    // GMSH 10 node quadratic tet node ordering
    addCoordinate(A, B);
    addCoordinate(A, C);
    addCoordinate(A, D);
    addCoordinate(B, C);
    addCoordinate(B, D);
    addCoordinate(C, D);
  }
}

std::vector<unsigned int> createNewTetElementNodes(const std::vector<unsigned int> &tetNodes,
                                                   const SpaMtrix::IRCMatrix &newNodeNumbers,
                                                   unsigned int numTetrahedra) {
  const unsigned int nodesPerTet = 4;
  assert(tetNodes.size() == nodesPerTet * numTetrahedra);

  std::vector<unsigned int> newTets;
  newTets.reserve(numTetrahedra * 10);

  for (unsigned int tetCounter = 0; tetCounter < numTetrahedra; tetCounter ++) {
    const unsigned int A = tetNodes[tetCounter * nodesPerTet + 0];
    const unsigned int B = tetNodes[tetCounter * nodesPerTet + 1];
    const unsigned int C = tetNodes[tetCounter * nodesPerTet + 2];
    const unsigned int D = tetNodes[tetCounter * nodesPerTet + 3];

    // Gmsh TET10 node ordering: A, B, C, D, AB, BC, AC, AD, CD, BD
    newTets.push_back(A);
    newTets.push_back(B);
    newTets.push_back(C);
    newTets.push_back(D);
    newTets.push_back((unsigned int) newNodeNumbers.getValue(A, B));  // [4] AB
    newTets.push_back((unsigned int) newNodeNumbers.getValue(B, C));  // [5] BC
    newTets.push_back((unsigned int) newNodeNumbers.getValue(A, C));  // [6] AC
    newTets.push_back((unsigned int) newNodeNumbers.getValue(A, D));  // [7] AD
    newTets.push_back((unsigned int) newNodeNumbers.getValue(C, D));  // [8] CD
    newTets.push_back((unsigned int) newNodeNumbers.getValue(B, D));  // [9] BD
  }

  return newTets;
}

std::vector<unsigned int> createNewTriElementNodes(const std::vector<unsigned int> &triNodes,
                                                   const SpaMtrix::IRCMatrix &newNodeNumbers,
                                                   unsigned int numTriangles) {
  const unsigned int nodesPerTri = 3;
  assert(triNodes.size() == nodesPerTri * numTriangles);

  std::vector<unsigned int> newTris;
  newTris.reserve(numTriangles * 6);

  for (unsigned int triCounter = 0; triCounter < numTriangles; triCounter ++) {
    const unsigned int A = triNodes[triCounter * nodesPerTri + 0];
    const unsigned int B = triNodes[triCounter * nodesPerTri + 1];
    const unsigned int C = triNodes[triCounter * nodesPerTri + 2];

    // Use GMSH 6 node quadratic triangle node order
    newTris.push_back(A);
    newTris.push_back(B);
    newTris.push_back(C);
    newTris.push_back((unsigned int) newNodeNumbers.getValue(A, B));
    newTris.push_back((unsigned int) newNodeNumbers.getValue(B, C));
    newTris.push_back((unsigned int) newNodeNumbers.getValue(A, C));
  }

  return newTris;
}

void convertLinearMeshDataToQuadratic(RawMeshData &meshData) {
  assert(meshData.getElementOrder() == 1);

  Log::info("Converting mesh with linear elements to quadratic element mesh");
  const unsigned int nodesPerTet = 4;
  const unsigned int nodesPerTri = 3;
  const unsigned int numTetrahedra = meshData.tetNodes.size() / nodesPerTet;
  const unsigned int numTriangles = meshData.triNodes.size() / nodesPerTri;;

  // find largest current node index in tetrahedra. New nodes will be added after this index
  const unsigned int maxNodeIndex = *std::max_element(meshData.tetNodes.begin(), meshData.tetNodes.end());

  std::pair<SpaMtrix::IRCMatrix, unsigned int> nodeMapping = createNewNodeNumbersMapping(meshData.tetNodes, maxNodeIndex);

  SpaMtrix::IRCMatrix &newNodesMap = nodeMapping.first;
  unsigned int maxNodeNumber = nodeMapping.second;
  const unsigned int newNodeCount = maxNodeNumber + 1;
  Log::info("new node count: {}", newNodeCount);

  appendNewNodeCoordinates(meshData.points, newNodesMap, meshData.tetNodes, numTetrahedra, newNodeCount);
  assert(meshData.points.size() == newNodeCount);

  meshData.tetNodes = createNewTetElementNodes(meshData.tetNodes, newNodesMap, numTetrahedra);
  assert(meshData.tetNodes.size() == numTetrahedra * 10);

  meshData.triNodes = createNewTriElementNodes(meshData.triNodes, newNodesMap, numTriangles);
  assert(meshData.triNodes.size() == numTriangles * 6);

  meshData.setElementOrder(2);
  validateAndSnapQuadraticTetrahedra(meshData);
}

//</editor-fold>