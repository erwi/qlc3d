#ifndef PROJECT_QLC3D_LC_SOLVER_H
#define PROJECT_QLC3D_LC_SOLVER_H
#include <memory>
#include <fe/gaussian-quadrature.h>
class LC;
class Simu;
class SolverSettings;
class SolutionVector;
class Geometry;
class SimulationState;
class Alignment;

namespace SpaMtrix {
  class IRCMatrix;
  class Vector;
}

class LCSolver {
  const double rt2 = std::sqrt(2.0);
  const double rt3 = std::sqrt(3.0);
  const double rt6 = std::sqrt(6.0);
  const double A, B, C;
  const double L1, L2, L3, L6;
  const double deleps;
  const SolverSettings &solverSettings;

  bool isSymmetricMatrix;
  bool isThreeElasticConstants;

  std::unique_ptr<SpaMtrix::IRCMatrix> K;
  std::unique_ptr<SpaMtrix::Vector> L;
  std::unique_ptr<SpaMtrix::Vector> X;

  void initialiseMatrixSystem(const SolutionVector &q, const Geometry &geom);
  void assembleMatrixSystem(const SolutionVector &q, const SolutionVector &v, const Geometry &geom);
  void assembleVolumeTerms(const SolutionVector &q, const SolutionVector &v, const Geometry &geom);
  void solveMatrixSystem();
  double updateQ(SolutionVector &q);

  void assembleLocalVolumeMatrix(unsigned int indTet,
                                 double lK[20][20],
                                 double lL[20],
                                 unsigned int tetNodes[4],
                                 unsigned int tetDofs[4],
                                 GaussianQuadratureTet<11> shapes,
                                 const SolutionVector &q,
                                 const SolutionVector &v,
                                 const Geometry &geom);

  void addToGlobalMatrix(double lK[20][20],
                         double lL[20],
                         const SolutionVector &q,
                         const unsigned int tetNodes[4],
                         const unsigned int tetDofs[4]);

public:
  ~LCSolver();
  LCSolver(const LC &lc, const SolverSettings &solverSettings);

  double solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState);


};

#endif //PROJECT_QLC3D_LC_SOLVER_H
