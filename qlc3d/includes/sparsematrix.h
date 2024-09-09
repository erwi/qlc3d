#ifndef PROJECT_QLC3D_SPARSEMATRIX_H
#define PROJECT_QLC3D_SPARSEMATRIX_H
#include <memory>

namespace SpaMtrix {
  class IRCMatrix;
  class Vector;
}
class Geometry;
class SolutionVector;

// Creates sparse matrix for Q-tensor solution
std::unique_ptr<SpaMtrix::IRCMatrix> createQMatrix(const Geometry &geom,
                                  const SolutionVector &q,
                                  const int& materialNumber);
/**
 * Creates an empty mass-matrix for Q-tensor solution. Non-zeroes are expanded along the diagonal so that the matrix
 * contains 5 blocks of the same size, once per Q-tensor component.
 */
std::unique_ptr<SpaMtrix::IRCMatrix> createQMassMatrix(const Geometry &geom,
                                      const SolutionVector &q,
                                      const int& materialNumber);

#endif //PROJECT_QLC3D_SPARSEMATRIX_H
