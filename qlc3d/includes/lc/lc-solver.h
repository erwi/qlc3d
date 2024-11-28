#ifndef PROJECT_QLC3D_LC_SOLVER_H
#define PROJECT_QLC3D_LC_SOLVER_H
#include <memory>
#include <fe/gaussian-quadrature.h>
#include <alignment.h>

class LC;
class Simu;
class SolverSettings;
class SolutionVector;
class Geometry;
class SimulationState;
class Alignment;
class DofMap;

namespace SpaMtrix {
  class IRCMatrix;
  class Vector;
}

enum class LCSolverType {
  STEADY_STATE,
  TIME_STEPPING
};

struct ElapsedTimes {
  const double assemblyTimeSeconds;
  const double solveTimeSeconds;
};

struct LCSolverResult {
  const LCSolverType solverType;
  const int iterations;
  const double dq;
  const bool converged;
  const bool maxIterationsReached;
  const ElapsedTimes elapsedTimes;
};

struct LCSolverParams {
  const double A;
  const double B;
  const double C;
  const double L1;
  const double L2;
  const double L3;
  const double L4;
  const double L6;
  const double deleps;
  const double dt;
  const double u1;
  const double S0;
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
  const Alignment &alignment;
  const bool isChiral;
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
  [[nodiscard]] double maxAbs(const SpaMtrix::Vector &v) const;

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

  void assembleLocalWeakAnchoringMatrix(unsigned int indTri, double lK[15][15], double lL[15],
                                        unsigned int triNodes[3], unsigned int triDofs[3],
                                        GaussianQuadratureTri<7> shapes, const SolutionVector &q,
                                        const Geometry &geom, const Surface &surface,
                                        double surfaceOrder);

  /**
   * Assemble the global matrix system from the local element matrices.
   * @param lK local element matrix. The size of the matrix is 20x20 and 15x15 for 1st order tets and tris.
   * @param lL local elemnt RHS vector. The size is 20 and 15 for 1st order tets and tris.
   * @param dofMap DofMap for Q-tensor.
   * @param elemNodes tet or triangle element nodes. Size is elemNodeCount.
   * @param elemNodeCount number of nodes in the current element
   */
  void addToGlobalMatrix(double* lK, double* lL, const DofMap &dofMap,
                         const unsigned int* elemNodes, int elemNodeCount);
  void assembleMatrixSystemVolumeTerms(const SolutionVector &q, const SolutionVector &v, const Geometry &geom, const LCSolverParams &params);
  void assembleMatrixSystemWeakAnchoring(const SolutionVector &q, const Geometry &geom, const LCSolverParams &params);
  void assembleMatrixSystem(const SolutionVector &q, const SolutionVector &v, const Geometry &geom, const LCSolverParams &params);
public:
  ImplicitLCSolver(const LC &lc, const SolverSettings &solverSettings, const Alignment &alignment);
};

class SteadyStateLCSolver : public ILCSolver, protected ImplicitLCSolver {
protected:
  const double rt2 = std::sqrt(2.0);
  const double rt3 = std::sqrt(3.0);
  const double rt6 = std::sqrt(6.0);
  const double A, B, C;
  const double L1, L2, L3, L4, L6;
  const double deleps;

  void initialiseMatrixSystem(const SolutionVector &q, const Geometry &geom, double dt);


public:
  virtual ~SteadyStateLCSolver();
  SteadyStateLCSolver(const LC &lc, const SolverSettings &solverSettings, const Alignment &alignment);

  LCSolverResult solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState) override;
};

/**
 * Time stepping LC solver implementing a non-linear implicit time stepping scheme where time derivatives are calculated
 * at mid-points between time steps.
 */
class TimeSteppingLCSolver : public ILCSolver, protected ImplicitLCSolver {
  const double maxError;
  const int maxNewtonIterations;
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
  TimeSteppingLCSolver(const LC &lc, const SolverSettings &solverSettings, double maxError, const Alignment &alignment, int maxNewtonIterations);
  ~TimeSteppingLCSolver() = default;
  LCSolverResult solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState) override;
};

#endif //PROJECT_QLC3D_LC_SOLVER_H
