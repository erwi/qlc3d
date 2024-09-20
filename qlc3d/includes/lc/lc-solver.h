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

struct LCSolverParams {
  double A;
  double B;
  double C;
  double L1;
  double L2;
  double L3;
  double L6;
  double deleps;
  double dt;
  double u1;
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
  /** Euler-Lagrange equations vector */
  std::unique_ptr<SpaMtrix::Vector> L;
  /** The unknown vector of solutions for the LC */
  std::unique_ptr<SpaMtrix::Vector> X;

  /**
   * Solve the matrix system Kx = r. The result is written into the vector x.
   * Returns true/false depnding on whether converged to required tolerance
   */
  bool solveMatrixSystem(const SpaMtrix::IRCMatrix &Kmatrix, const SpaMtrix::Vector &r, SpaMtrix::Vector &x) const;
  double maxAbs(const SpaMtrix::Vector &v) const;

  /**
   * This is a common local element matrix assembly used by both steady-state and time-stepping solvers. This assembles
   * the Jacobian matrix and the RHS vector (Euler-Lagrange equations) for a single tetrahedral element.
   *
   * In the time-stepping solver, the results of this method are further augmented by terms required by the
   * implicit time-stepping algorithm.
   */
  void assembleLocalVolumeMatrix(unsigned int indTet,
                                 double lK[20][20],
                                 double lL[20],
                                 unsigned int tetNodes[4],
                                 unsigned int tetDofs[4],
                                 GaussianQuadratureTet<11> shapes,
                                 const SolutionVector &q,
                                 const SolutionVector &v,
                                 const Geometry &geom,
                                 const LCSolverParams &params);

  void addToGlobalMatrix(double lK[20][20], double lL[20], const SolutionVector &q, const unsigned int tetNodes[4]);
  void assembleMatrixSystem(const SolutionVector &q, const SolutionVector &v, const Geometry &geom, const LCSolverParams &params);
public:
  ImplicitLCSolver(const LC &lc, const SolverSettings &solverSettings);
};

class SteadyStateLCSolver : public ILCSolver, protected ImplicitLCSolver {
protected:
  const double rt2 = std::sqrt(2.0);
  const double rt3 = std::sqrt(3.0);
  const double rt6 = std::sqrt(6.0);
  const double A, B, C;
  const double L1, L2, L3, L6;
  const double deleps;

  void initialiseMatrixSystem(const SolutionVector &q, const Geometry &geom, double dt);


public:
  virtual ~SteadyStateLCSolver();
  SteadyStateLCSolver(const LC &lc, const SolverSettings &solverSettings);

  LCSolverResult solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState) override;
};

/**
 * Time stepping LC solver implementing a non-linear implicit time stepping scheme where time derivatives are calculated
 * at mid-points between time steps.
 */
class TimeSteppingLCSolver : public ILCSolver, protected ImplicitLCSolver {
  const double maxError;
  /** Mass matrix, required by implicit time stepping */
  std::unique_ptr<SpaMtrix::IRCMatrix> M;
  /** Q-tensor at previous time step */
  std::unique_ptr<SpaMtrix::Vector> q1;
  std::unique_ptr<SpaMtrix::Vector> dqdt; // time derivative of Q-tensor, estimated at end of time step
  /** RHS vector at previous time step */
  std::unique_ptr<SpaMtrix::Vector> f_prev;
  bool isFirstRun = true;

  void initialiseMatrixSystem(const SolutionVector &q, const Geometry &geom);

public:
  TimeSteppingLCSolver(const LC &lc, const SolverSettings &solverSettings, double maxError);
  ~TimeSteppingLCSolver() = default;
  LCSolverResult solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState) override;
};

#endif //PROJECT_QLC3D_LC_SOLVER_H
