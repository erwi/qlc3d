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

std::unique_ptr<SpaMtrix::IRCMatrix> createQMassMatrix(const Geometry &geom,
                                      const SolutionVector &q,
                                      const int& materialNumber);

#endif //PROJECT_QLC3D_SPARSEMATRIX_H
