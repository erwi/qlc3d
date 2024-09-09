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

enum class LCSolverType {
  STEADY_STATE,
  TIME_STEPPING
};

struct LCSolverResult {
  const LCSolverType solverType;
  const unsigned int iterations;
  const double dq;
  const bool converged;
};


class ILCSolver {
public:
  virtual ~ILCSolver() = default;
  virtual LCSolverResult solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState) = 0;
};

class ImplicitLCSolver {

protected:
  const LC &lc;
  const SolverSettings &solverSettings;
  const bool isSymmetricMatrix;
  const bool isThreeElasticConstants;

  /** Jacobian matrix */
  std::unique_ptr<SpaMtrix::IRCMatrix> K;

  /**
   * Solve the matrix system Kx = r. The result is written into the vector x.
   * Returns true/false depnding on whether converged to required tolerance
   */
  bool solveMatrixSystem(const SpaMtrix::IRCMatrix &Kmatrix, const SpaMtrix::Vector &r, SpaMtrix::Vector &x) const;
  double maxAbs(const SpaMtrix::Vector &v) const;
public:
  ImplicitLCSolver(const LC &lc, const SolverSettings &solverSettings);

};



class LCSolver : public ILCSolver, protected ImplicitLCSolver {
protected:
  const double rt2 = std::sqrt(2.0);
  const double rt3 = std::sqrt(3.0);
  const double rt6 = std::sqrt(6.0);
  const double A, B, C;
  const double L1, L2, L3, L6;
  const double deleps;

  //std::unique_ptr<SpaMtrix::IRCMatrix> K;
  std::unique_ptr<SpaMtrix::Vector> L;
  std::unique_ptr<SpaMtrix::Vector> X;

  void initialiseMatrixSystem(const SolutionVector &q, const Geometry &geom, double dt);
  void assembleMatrixSystem(const SolutionVector &q, const SolutionVector &v, const Geometry &geom);
  void assembleVolumeTerms(const SolutionVector &q, const SolutionVector &v, const Geometry &geom);
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
  virtual ~LCSolver();
  LCSolver(const LC &lc, const SolverSettings &solverSettings);

  LCSolverResult solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState) override;


};

#endif //PROJECT_QLC3D_LC_SOLVER_H
