#ifndef PROJECT_QLC3D_FIXEDNODES_H
#define PROJECT_QLC3D_FIXEDNODES_H
#include <vector>
#include <unordered_map>

class Mesh;

class FixedNodes {
  std::unordered_map<unsigned int, double> fixedValueByNodeIndex; // TODO: there should be no collisions, so could make this more efficient

public:

  void setFixedNodesPot(const Mesh &triangles, const std::unordered_map<unsigned int, double> &potentialByElectrode);
  void setFixedNodesAndValues(const std::unordered_map<unsigned int, double> &fixedNodesAndValues);
  void clear() { fixedValueByNodeIndex.clear(); }

  [[nodiscard]] bool isFixedNode(unsigned int i) const {
    //return whether i is a key in fixedValueByNodeIndex
    return fixedValueByNodeIndex.find(i) != fixedValueByNodeIndex.end();
  };

  [[nodiscard]] size_t getnFixedNodes() const {
    return fixedValueByNodeIndex.size();
  }

  [[nodiscard]] const std::unordered_map<unsigned int, double>& getFixedValueByNodeIndex() const {
    return fixedValueByNodeIndex;
  }
};
#endif //PROJECT_QLC3D_FIXEDNODES_H
