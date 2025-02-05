#ifndef PROJECT_QLC3D_FE_UTIL_H
#define PROJECT_QLC3D_FE_UTIL_H
#include <algorithm>
#include <lc-representation.h>
#include <util/exception.h>
#include <unordered_set>
#include <geom/coordinates.h>
#include <geom/vec3.h>


inline bool containsIndex(const unsigned int* start, const unsigned int* end, const unsigned int value) {
  return std::find(start, end, value) != end;
}

inline void reorderBoundaryTetNodes(unsigned int tetNodes[4], const unsigned int triNodes[3]) {
  // find tetNode not in triNodes
  unsigned int tetNode =  tetNodes[0];
  for (int i = 0; i < 4; ++i) {
    tetNode = tetNodes[i];
    if (!containsIndex(triNodes, triNodes + 3, tetNode)) {
      tetNode = tetNodes[i];
      break;
    }
  }

  for (int i = 0; i < 3; i++) {
    tetNodes[i] = triNodes[i];
  }
  tetNodes[3] = tetNode;
}

inline void reorderQuadraticBoundaryTetNodes(std::vector<unsigned int> &tetNodes,
                                             const std::vector<unsigned int> &triNodes,
                                             const Coordinates &coords) {
  std::unordered_set<unsigned int> triNodesSet(triNodes.begin(), triNodes.end());
  std::unordered_set<unsigned int> tetNodesSet(tetNodes.begin(), tetNodes.end());
  // populate tetNodesOut with nodes in correct order
  std::vector<unsigned int> tetNodesOut;
  tetNodesOut.resize(tetNodes.size(), 0);

  // first find the corner nodes of the triangle, these will match the first 3 nodes of the tet
  tetNodesOut[0] = triNodes[0];
  tetNodesOut[1] = triNodes[1];
  tetNodesOut[2] = triNodes[2];

  // find which tetNode is not in triangle, this will be the 4th node in tetNodesOut
  //unsigned int tetNode = 0;
  for (unsigned int i = 0; i < tetNodes.size(); ++i) {
    if (triNodesSet.find(tetNodes[i]) == triNodesSet.end()) {
      tetNodesOut[3] = tetNodes[i];
      break;
    }
  }

  // 5, 6, 7 are the mid-edge nodes, same as for the triangle
  tetNodesOut[4] = triNodes[3];
  tetNodesOut[5] = triNodes[4];
  tetNodesOut[6] = triNodes[5];

  tetNodesSet.erase(tetNodesOut[0]);
  tetNodesSet.erase(tetNodesOut[1]);
  tetNodesSet.erase(tetNodesOut[2]);
  tetNodesSet.erase(tetNodesOut[3]);
  tetNodesSet.erase(tetNodesOut[4]);
  tetNodesSet.erase(tetNodesOut[5]);
  tetNodesSet.erase(tetNodesOut[6]);
  assert(tetNodesSet.size() == 3); // should contain the remaining mid-edge nodes, but we don't know the order

  // tetNodesOut[7] should be between tetNodesOut[0] and tetNodesOut[3]
  // tetNodesOut[8] should be between tetNodesOut[1] and tetNodesOut[3]
  // tetNodesOut[9] should be between tetNodesOut[2] and tetNodesOut[3]
  auto p3 = coords.getPoint(tetNodesOut[3]);
  auto p7 = 0.5 * (coords.getPoint(tetNodesOut[0]) + p3);
  auto p8 = 0.5 * (coords.getPoint(tetNodesOut[1]) + p3);
  auto p9 = 0.5 * (coords.getPoint(tetNodesOut[2]) + p3);

  for (auto i : tetNodesSet) {
    auto p = coords.getPoint(i);
    // check which mid-edge node it is closest to
    double distances[3] = {p.distance(p7), p.distance(p8), p.distance(p9)};
    auto minIndex = std::min_element(distances, distances + 3) - distances;

    assert(tetNodesOut[7 + minIndex] == 0);
    tetNodesOut[7 + minIndex] = i;
  }

  // finally copy back tetNodesOut to tetNodes. They should now be in correct order
  for (unsigned int i = 0; i < tetNodes.size(); ++i) {
    tetNodes[i] = tetNodesOut[i];
  }
}

inline void reorderBoundaryTetNodes(std::vector<unsigned int> &tetNodes, const std::vector<unsigned int> &triNodes) {
  const unsigned int numTetNodes = tetNodes.size();
  const unsigned int numTriNodes = triNodes.size();

  if (numTetNodes != 4) {
    throw NotYetImplementedException("Reordering only 4-node tetrahedra is supported");
  }

  // find tetNode not in triNodes
  unsigned int tetNode =  tetNodes[0];
  for (unsigned int i = 0; i < numTetNodes; ++i) {
    tetNode = tetNodes[i];
    if (!containsIndex(&triNodes[0], &triNodes[numTriNodes], tetNode)) {
      tetNode = tetNodes[i];
      break;
    }
  }

  for (unsigned int i = 0; i < numTriNodes; i++) {
    tetNodes[i] = triNodes[i];
  }
  tetNodes[numTriNodes] = tetNode; // last node is the one not in triNodes TODO, for higher order elements the first 6 nodes should be the tri nodes and the rest the internal tet nodes
}


/**
 * Redistributes the tensor components to separate destination arrays. It is assumed the lengths of the
 * destination arrays match the number of source tensors.
 */
inline void groupComponentsByNode(qlc3d::TTensor* start, qlc3d::TTensor* end,
                            double* q1, double* q2, double* q3, double* q4, double* q5) {
  unsigned int i = 0;
  for (qlc3d::TTensor* it = start; it != end; ++it, ++i) {
      q1[i] = it->t1();
      q2[i] = it->t2();
      q3[i] = it->t3();
      q4[i] = it->t4();
      q5[i] = it->t5();
  }
}

#endif //PROJECT_QLC3D_FE_UTIL_H
