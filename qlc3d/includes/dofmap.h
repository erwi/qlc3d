#ifndef PROJECT_QLC3D_DOFMAP_H
#define PROJECT_QLC3D_DOFMAP_H
#include <vector>
#include <limits>
#include <unordered_set>

class Geometry;
class FixedNodes;
class PeriodicNodesMapping;

class DofMap {
  unsigned int nDof;
  unsigned int nDimensions;
  unsigned int nFreeNodes;
  std::vector<unsigned int> dofs;
  static constexpr unsigned int NOT_DOF = std::numeric_limits<unsigned int>::max();

public:
  DofMap(unsigned int nDof, unsigned int nDimensions);

  void calculateMapping(const std::unordered_set<unsigned int> &fixedNodes,
                        const std::vector<unsigned int> &periodicNodesMapping);

  /**
   * Returns the mapped degree of freedom for the given index and dimension 0. If the node is fixed, i.e. not a degree of freedom, returns NOT_DOF.
   * @param index The index of the node.
   * @param dimension The dimension (0-based) for which to get the degree of freedom.
   * @return The mapped degree of freedom, or NOT_DOF if the node is fixed.
   */
  [[nodiscard]] unsigned int getDof(unsigned int index) const { return getDof(index, 0); }
  [[nodiscard]] unsigned int getDof(unsigned int index, unsigned int dimension) const { return dofs[index + dimension * nDof]; };

  /** Returns true if the given degree of freedom is fixed (i.e. not a degree of freedom). NOTE: operates on the mapped dof, i.e. the result of calling getDof(index) **/
  [[nodiscard]] bool isFixedDof(unsigned int dof) const { return dof == NOT_DOF; }
  /** Returns true if the given degree of freedom is free (i.e. a degree of freedom). NOTE: operates on the mapped dof, i.e. the result of calling getDof(index) **/
  [[nodiscard]] bool isFreeDof(unsigned int dof) const { return !isFixedDof(dof); }
  /** Returns true if the node at the given index is fixed (i.e. not a degree of freedom). NOTE: operates on the unmapped raw index, not the result of calling getDof(index) */
  [[nodiscard]] bool isFixedNode(unsigned int index) const { return dofs[index] == NOT_DOF; }
  /** Returns true if the node at the given index is free (i.e. a degree of freedom). NOTE: operates on the unmapped raw index, not the result of calling getDof(index) */
  [[nodiscard]] bool isFreeNode(unsigned int index) const { return !isFixedNode(index); }
  /** Number of degrees of freedom per dimension, including fixed and periodic nodes. */
  [[nodiscard]] unsigned int getnDof() const { return nDof; }
  [[nodiscard]] unsigned int getnDimensions() const { return nDimensions; }
  /** number of free nodes per dimension */
  [[nodiscard]] unsigned int getnFreeNodes() const { return nFreeNodes; }
};

#endif //PROJECT_QLC3D_DOFMAP_H
