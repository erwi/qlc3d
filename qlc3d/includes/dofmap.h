#ifndef PROJECT_QLC3D_DOFMAP_H
#define PROJECT_QLC3D_DOFMAP_H
#include <vector>
#include <limits>

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

  void calculateMapping(const Geometry &geom, const FixedNodes &fixedNodules);

  [[nodiscard]] unsigned int getDof(unsigned int index) const { return dofs[index]; }
  [[nodiscard]] unsigned int getnFreeNodes() const { return nFreeNodes; }
  [[nodiscard]] unsigned int getnDof() const { return nDof; }
  [[nodiscard]] unsigned int getnDimensions() const { return nDimensions; }
};

#endif //PROJECT_QLC3D_DOFMAP_H
