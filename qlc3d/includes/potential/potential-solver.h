#ifndef PROJECT_QLC3D_POTENTIAL_SOLVER_H
#define PROJECT_QLC3D_POTENTIAL_SOLVER_H
#include <memory>

class SolutionVector;
class Geometry;
class Electrodes;
class LC;
class SolverSettings;
namespace SpaMtrix {
  class IRCMatrix;
}

// todo: make private...
void calcpot3d(SpaMtrix::IRCMatrix &Kpot,
               SolutionVector &v,
               const SolutionVector &q,
               const LC &lc,
               const Geometry &geom,
               const SolverSettings &settings,
               const Electrodes &electrodes);

class PotentialSolver {
  std::shared_ptr<Electrodes> electrodes;
  std::shared_ptr<LC> lc;
  std::shared_ptr<SolverSettings> solverSettings;
  std::unique_ptr<SpaMtrix::IRCMatrix> K;

  void createPotentialMatrix(const Geometry &geom, const SolutionVector &sol);

public:
  PotentialSolver(std::shared_ptr<Electrodes> electrodes, std::shared_ptr<LC> lc, std::shared_ptr<SolverSettings> solverSettings) :
    electrodes(std::move(electrodes)), lc(std::move(lc)), solverSettings(std::move(solverSettings)), K(nullptr) {};

  /** Solves potential and stores it in vOut */
  void solvePotential(SolutionVector &vOut,
                      const SolutionVector &q,
                      const Geometry &geom);

  void onGeometryChanged();
};
#endif //PROJECT_QLC3D_POTENTIAL_SOLVER_H
