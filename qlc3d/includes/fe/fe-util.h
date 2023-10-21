#ifndef PROJECT_QLC3D_FE_UTIL_H
#define PROJECT_QLC3D_FE_UTIL_H
#include <algorithm>

inline bool containsIndex(const unsigned int* start, const unsigned int* end, const unsigned int index) {
  return std::find(start, end, index) != end;
}

inline void reorderTetNodes(unsigned int tetNodes[4], const unsigned int triNodes[3]) {
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

#endif //PROJECT_QLC3D_FE_UTIL_H
