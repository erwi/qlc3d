#ifndef PROJECT_QLC3D_DOFMAP_H
#define PROJECT_QLC3D_DOFMAP_H
#include <vector>
#include <limits>
#include <unordered_set>

class Geometry;
class FixedNodes;

class DofMap {
  unsigned int nDof;
  unsigned int nDimensions;
  unsigned int nFreeNodes;
  std::vector<unsigned int> dofs;

public:

  static constexpr unsigned int NOT_DOF = std::numeric_limits<unsigned int>::max();

  DofMap(unsigned int nDof, unsigned int nDimensions);
  DofMap(const std::vector<unsigned int> map);
  DofMap(const DofMap &other);

  void calculateMapping(const Geometry &geom, const std::unordered_set<unsigned int> &fixedNodes);

  [[nodiscard]] unsigned int getDof(unsigned int index) const { return dofs[index]; }
  [[nodiscard]] bool isFixedNode(unsigned int index) const { return dofs[index] == NOT_DOF; }
  [[nodiscard]] bool isFreeNode(unsigned int index) const { return dofs[index] < NOT_DOF; }
  [[nodiscard]] unsigned int getnDof() const { return nDof; }
  [[nodiscard]] unsigned int getnDimensions() const { return nDimensions; }
  [[nodiscard]] unsigned int getnFreeNodes() const { return nFreeNodes; }
};

#endif //PROJECT_QLC3D_DOFMAP_H
