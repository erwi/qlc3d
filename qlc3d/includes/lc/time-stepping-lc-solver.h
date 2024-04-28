#ifndef PROJECT_QLC3D_TIME_STEPPING_LC_SOLVER_H
#define PROJECT_QLC3D_TIME_STEPPING_LC_SOLVER_H
#include <lc/lc-solver.h>


class TimeSteppingLCSolver : public ILCSolver {
  const LC &lc;
  const SolverSettings &solverSettings;
  std::unique_ptr<SpaMtrix::IRCMatrix> K;
  std::unique_ptr<SpaMtrix::Vector> L;
  std::unique_ptr<SpaMtrix::Vector> X;

  void initialiseMatrixSystem(const SolutionVector &q, const Geometry &geom, double dt);

public:
  TimeSteppingLCSolver(const LC &lc, const SolverSettings &solverSettings);
  ~TimeSteppingLCSolver() = default;
  double solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState) override;
};

#endif //PROJECT_QLC3D_TIME_STEPPING_LC_SOLVER_H
