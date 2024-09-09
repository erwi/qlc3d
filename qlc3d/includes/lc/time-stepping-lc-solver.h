#ifndef PROJECT_QLC3D_TIME_STEPPING_LC_SOLVER_H
#define PROJECT_QLC3D_TIME_STEPPING_LC_SOLVER_H
#include <lc/lc-solver.h>

struct Params {
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

class TimeSteppingLCSolver : public ILCSolver, protected ImplicitLCSolver {

  std::unique_ptr<SpaMtrix::IRCMatrix> M; // mass-matrix
  //std::unique_ptr<SpaMtrix::IRCMatrix> K;
  std::unique_ptr<SpaMtrix::Vector> L;
  std::unique_ptr<SpaMtrix::Vector> X;

  /** Q-tensor at previous time step */
  std::unique_ptr<SpaMtrix::Vector> q1;
  /** RHS vector at previous time step */
  std::unique_ptr<SpaMtrix::Vector> f_prev;
  bool isFirstRun = true;

  void initialiseMatrixSystem(const SolutionVector &q, const Geometry &geom);

  void addToGlobalMatrix(double lK[20][20],
                         double lL[20],
                         const SolutionVector &q,
                         const unsigned int tetNodes[4]);

  void assembleLocalVolumeMatrix_impl(unsigned int indTet,
                                 double lK[20][20],
                                 double lL[20],
                                 unsigned int tetNodes[4],
                                 unsigned int tetDofs[4],
                                 GaussianQuadratureTet<11> shapes,
                                 const SolutionVector &q,
                                 const SolutionVector &v,
                                 const Geometry &geom,
                                 const Params &params);

  void assembleVolumeTerms_impl(const SolutionVector &q, const SolutionVector &v, const Geometry &geom,
                                const Params &params);

  //void solveMatrixSystem(const SpaMtrix::IRCMatrix &K, const SpaMtrix::Vector &r, SpaMtrix::Vector &X);

public:
  TimeSteppingLCSolver(const LC &lc, const SolverSettings &solverSettings);
  ~TimeSteppingLCSolver() = default;
  LCSolverResult solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState) override;
};

#endif //PROJECT_QLC3D_TIME_STEPPING_LC_SOLVER_H
