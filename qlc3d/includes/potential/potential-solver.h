#ifndef PROJECT_QLC3D_POTENTIAL_SOLVER_H
#define PROJECT_QLC3D_POTENTIAL_SOLVER_H
#include <memory>
#include <fe/gaussian-quadrature.h>
#include <spamtrix_densematrix.hpp>

class SolutionVector;
class Geometry;
class Mesh;
class Coordinates;
class Electrodes;
class LC;
class SolverSettings;

namespace SpaMtrix {
  class IRCMatrix;
  class Vector;
  class DenseMatrix;
}

class DofMap;

// todo: make private...
void calcpot3d(SpaMtrix::IRCMatrix &Kpot,
               SolutionVector &v,
               const SolutionVector &q,
               const LC &lc,
               const Geometry &geom,
               const SolverSettings &settings,
               const Electrodes &electrodes);

class PotentialSolver {

  const double rt2 = std::sqrt(2.0);
  const double rt3 = std::sqrt(3.0);
  const double rt6 = std::sqrt(6.0);

  double eper_lc, deleps;
  double efe, efe2;
  double S0;

  Electrodes &electrodes;
  std::shared_ptr<LC> lc;
  std::shared_ptr<SolverSettings> solverSettings;

  /** Sparse matrix for linear equations system */
  std::unique_ptr<SpaMtrix::IRCMatrix> K;
  /** RHS vector for linear equations system */
  std::unique_ptr<SpaMtrix::Vector> L;
  /** Solution vector for linear equations system */
  std::unique_ptr<SpaMtrix::Vector> V;

  void createPotentialMatrix(const Geometry &geom, const DofMap &dofMap);

  bool isPotentialSolutionRequired(const SolutionVector &v) const;
  void setUniformEField(SolutionVector &vOut, const Coordinates &coordinates) const;
  void initialiseMatrixSystem(const SolutionVector &vOut, const Geometry &geom);
  void assembleMatrixSystem(const SolutionVector &v, const SolutionVector &q, const Geometry &geom);
  void assembleVolume(const SolutionVector &v, const SolutionVector &q, const Geometry &geom);
  void assembleNeumann(const SolutionVector &v, const SolutionVector &q, const Geometry &geom);
  void addToGlobalMatrix(const SpaMtrix::DenseMatrix &lK, const std::vector<double> &lL,
                         const SolutionVector &v,
                         const std::vector<unsigned int> &tetNodes,
                         const std::vector<unsigned int> &tetDofs);

  void localKL(const Geometry &geom,
               SpaMtrix::DenseMatrix &lK,
               std::vector<double> &lL,
               unsigned int elementIndex,
               const SolutionVector &q,
               const LC &lc,
               const Electrodes &electrodes,
              TetShapeFunction &s);

  void localKLNeumann(
          const Coordinates &coordinates,
          SpaMtrix::DenseMatrix &lK,
          std::vector<double> &lL,
          const std::vector<unsigned int> &tetNodes,
          const SolutionVector &q,
          double triDet,
          double tetDet,
          const Vec3 &n,
          GaussianQuadratureTet<7> &shapes);

  void solveMatrixSystem(SolutionVector &v);

  [[nodiscard]] bool isFixedNode(unsigned int i);
  [[nodiscard]] bool isFreeNode(unsigned int i);

public:
  PotentialSolver(Electrodes &electrodes, std::shared_ptr<LC> lc, std::shared_ptr<SolverSettings> solverSettings);
  ~PotentialSolver(); // default implementation required to avoid complains from unique_ptr about incomplete type

  /** Solves potential and stores it in vOut */
  void solvePotential(SolutionVector &vOut,
                      const SolutionVector &q,
                      const Geometry &geom);

  void onGeometryChanged();



};
#endif //PROJECT_QLC3D_POTENTIAL_SOLVER_H
