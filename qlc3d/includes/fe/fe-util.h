#ifndef PROJECT_QLC3D_FE_UTIL_H
#define PROJECT_QLC3D_FE_UTIL_H
#include <algorithm>
#include <lc-representation.h>
#include "util/exception.h"

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

inline void reorderBoundaryTetNodes(std::vector<unsigned int> &tetNodes, const std::vector<unsigned int> &triNodes) {
  const unsigned int numTetNodes = tetNodes.size();
  const unsigned int numTriNodes = triNodes.size();

  if (numTetNodes > 4) {
    throw NotYetImplementedException("Higher order elements are not yet supported.");
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
