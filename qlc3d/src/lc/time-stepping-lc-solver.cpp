#include <lc/time-stepping-lc-solver.h>
#include <lc/lc-energy-terms.h>
#include <lc.h>
#include <sparsematrix.h>
#include <material_numbers.h>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <simulation-state.h>
#include <solutionvector.h>
#include <memory>
#include <fe/gaussian-quadrature.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <spamtrix_blas.hpp>
#include <util/logging.h>
#include <solver-settings.h>

TimeSteppingLCSolver::TimeSteppingLCSolver(const LC &lc, const SolverSettings &solverSettings) : ImplicitLCSolver(lc, solverSettings) {
  Log::info("Creating time stepping solver for Q-tensor, isTreeElasticConstants={}, isSymmetricMatrix={}", isThreeElasticConstants, isSymmetricMatrix);
}

bool isFixedNode(idx i) {
  return i == NOT_AN_INDEX;
}

void addToGlobalMatrix(const double lK[4][4], idx tetDofs[4], SpaMtrix::IRCMatrix &M) {
  const unsigned int numFreeNodes = M.getNumRows() / 5;

  for (idx i = 0; i < 4; i++) {
    const idx iDof = tetDofs[i];
    if (isFixedNode(iDof)) { continue; } // this row should not even exist in the matrix

    for (idx j = 0; j < 4; j++) {
      const idx jDof = tetDofs[j];
      if (isFixedNode(jDof)) { continue; }

      for (int dof = 0; dof < 5; dof ++) {
        M.sparse_add(iDof + dof * numFreeNodes, jDof + dof * numFreeNodes, lK[i][j]);
      }
    }
  }
}

void assembleElementMassMatrix(double lK[4][4], idx tetNodes[4], double tetDeterminant, const Geometry &geom, GaussianQuadratureTet<11> &shapes) {
  memset(lK, 0, 4 * 4 * sizeof(double));
  Vec3 coordinates[4];
  geom.getCoordinates().loadCoordinates(tetNodes, tetNodes + 4, coordinates);
  coordinates[0] *= 1e-6;
  coordinates[1] *= 1e-6;
  coordinates[2] *= 1e-6;
  coordinates[3] *= 1e-6;

  shapes.initialiseElement(coordinates, tetDeterminant);

  for (; shapes.hasNextPoint(); shapes.nextPoint()) {
    const double mul = shapes.weight() * tetDeterminant;
    for (idx i = 0; i < 4; i++) {
      for (idx j = 0; j < 4; j++) {
        lK[i][j] += mul * shapes.N(i) * shapes.N(j);
      }
    }
  }
}

void assembleGlobalMassMatrix(SpaMtrix::IRCMatrix &M, const Geometry &geom, const SolutionVector &q) {
  GaussianQuadratureTet<11> shapes = gaussQuadratureTet4thOrder();
  const Mesh tetMesh = geom.getTetrahedra();
  const idx tetCount = tetMesh.getnElements();
  double lK[4][4];
  idx tetNodes[4];
  idx tetDofs[4];

  for (idx tetIndex = 0; tetIndex < tetCount; tetIndex++) {
    if (tetMesh.getMaterialNumber(tetIndex) != MAT_DOMAIN1) {
      continue;
    }

    tetMesh.loadNodes(tetIndex, tetNodes);
    q.loadEquNodes(tetNodes, tetNodes + 4, tetDofs);

    assembleElementMassMatrix(lK, tetNodes, tetMesh.getDeterminant(tetIndex), geom, shapes);

    addToGlobalMatrix(lK, tetDofs, M);
  }
}

void TimeSteppingLCSolver::addToGlobalMatrix(double lK[20][20], double lL[20], const SolutionVector &q, const unsigned int tetNodes[4]) {
  idx npLC = q.getnDoF();
  for (int i = 0; i < 20; i++) {
    idx ri = tetNodes[i % 4] + npLC * (i / 4);
    idx eqr = q.getEquNode(ri);

    if (eqr == NOT_AN_INDEX) { // this row is fixed, so doesn't even exist in the system
      continue;
    }

    (*L)[eqr] += lL[i];
    for (int j = 0 ; j < 20 ; j++) { // LOOP OVER COLUMNS
      idx rj = tetNodes[j % 4] + npLC * (j / 4);
      idx eqc = q.getEquNode(rj);

      // Note: we can simply ignore fixed columns as they correspond to fixed zeroes in the RHS vector L,
      // so that their contribution to other rows in L would be zero anyway. We only set the free rows/columns
      // to the global matrix
      if (eqc != NOT_AN_INDEX) { // it's free
        K->sparse_add(eqr, eqc, lK[i][j]);
      }
    }
  }
}

void TimeSteppingLCSolver::assembleLocalVolumeMatrix_impl(unsigned int indTet,
                                                          double (*lK)[20],
                                                          double *lL,
                                                          unsigned int *tetNodes,
                                                          unsigned int *tetDofs,
                                                          GaussianQuadratureTet<11> shapes,
                                                          const SolutionVector &q,
                                                          const SolutionVector &v,
                                                          const Geometry &geom,
                                                          const Params &params) {

  memset(lK, 0., 20 * 20 * sizeof(double));
  memset(lL, 0., 20 * sizeof(double));
  Vec3 tetCoords[4];
  geom.getCoordinates().loadCoordinates(tetNodes, tetNodes + 4, tetCoords);
  tetCoords[0] *= 1e-6;
  tetCoords[1] *= 1e-6;
  tetCoords[2] *= 1e-6;
  tetCoords[3] *= 1e-6;

  const double tetDeterminant = geom.getTetrahedra().getDeterminant(indTet);
  shapes.initialiseElement(tetCoords, tetDeterminant);

  qlc3d::TTensor qNodal[4];
  q.loadQtensorValues(tetNodes, tetNodes + 4, qNodal);

  double potential[4];
  v.loadValues(tetNodes, tetNodes + 4, potential);

  double q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z;
  shapes.sampleAllX(qNodal, q1x, q2x, q3x, q4x, q5x);
  shapes.sampleAllY(qNodal, q1y, q2y, q3y, q4y, q5y);
  shapes.sampleAllZ(qNodal, q1z, q2z, q3z, q4z, q5z);

  double Ml[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  const double timeMultiplier = 2. * params.u1 / params.dt; // effect of viscosity and time step size
  for (; shapes.hasNextPoint(); shapes.nextPoint()) {
    double q1, q2, q3, q4, q5;
    double Vx, Vy, Vz;
    shapes.sampleAll(qNodal, q1, q2, q3, q4, q5);
    Vx = shapes.sampleX(potential);
    Vy = shapes.sampleY(potential);
    Vz = shapes.sampleZ(potential);

    LcEnergyTerms::assembleThermotropic(lK, lL, shapes, tetDeterminant, q1, q2, q3, q4, q5, params.A, params.B, params.C);
    LcEnergyTerms::assembleElasticL1(lK, lL, shapes, tetDeterminant, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z, params.L1);
    LcEnergyTerms::assembleThreeElasticConstants(lK, lL, shapes, tetDeterminant, q1, q2, q3, q4, q5, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z, params.L2, params.L3, params.L6);

    if (Vx != 0.0 || Vy != 0.0 || Vz != 0.0) {
      LcEnergyTerms::assembleDielectric(lL, shapes, tetDeterminant, Vx, Vy, Vz, params.deleps);
    }

    LcEnergyTerms::assembleMassMatrix(Ml, shapes, tetDeterminant, timeMultiplier);
  }

  for (idx i = 0; i < 4; i++) {
    // add M to local Jacobian matrix
    for (idx j = 0; j < 4; j++) {
      const double m = Ml[i][j];
      lK[i][j] += m;
      lK[i + 4][j + 4] += m;
      lK[i + 8][j + 8] += m;
      lK[i + 12][j + 12] += m;
      lK[i + 16][j + 16] += m;
    }
  }
}

void TimeSteppingLCSolver::assembleVolumeTerms_impl(const SolutionVector &q, const SolutionVector &v, const Geometry &geom,
                                               const Params &params) {
  GaussianQuadratureTet<11> shapes = gaussQuadratureTet4thOrder();
  const unsigned int elementCount = geom.getTetrahedra().getnElements();
  const Mesh &tets = geom.getTetrahedra();
  double lK[4 * 5][4 * 5];
  double lL[4 * 5];
  idx tetNodes[4];
  idx tetDofs[4];

  for (unsigned int indTet = 0; indTet < elementCount; indTet++) {
    if (tets.getMaterialNumber(indTet) != MAT_DOMAIN1) {
      continue;
    }
    tets.loadNodes(indTet, tetNodes);
    q.loadEquNodes(tetNodes, tetNodes + 4, tetDofs);

    assembleLocalVolumeMatrix_impl(indTet, lK, lL, tetNodes, tetDofs, shapes, q, v, geom, params);

    addToGlobalMatrix(lK, lL, q, tetNodes);
  }
}

void TimeSteppingLCSolver::initialiseMatrixSystem(const SolutionVector &q, const Geometry &geom) {
  if (K != nullptr) {
    return;
  }
  const idx numFreeNodes = q.getnFreeNodes();
  const idx N = numFreeNodes * 5;

  M = createQMassMatrix(geom, q, MAT_DOMAIN1);
  assembleGlobalMassMatrix(*M, geom, q);

  K = createQMatrix(geom, q, MAT_DOMAIN1);
  L = std::make_unique<SpaMtrix::Vector>(N);
  X = std::make_unique<SpaMtrix::Vector>(N);

  q1 = std::make_unique<SpaMtrix::Vector>(N);
  f_prev = std::make_unique<SpaMtrix::Vector>(N);
  *f_prev = 0;
}

void calculateRhsResidualVector(SpaMtrix::Vector &r,
                                const SpaMtrix::Vector &q1,
                                const SpaMtrix::Vector &q2,
                                const SpaMtrix::Vector &f1,
                                const SpaMtrix::Vector &f2,
                                const SpaMtrix::IRCMatrix &M,
                                double timeMultiplier
                  ) {

  assert(r.getLength() == q1.getLength());
  assert(r.getLength() == q2.getLength());
  assert(r.getLength() == f1.getLength());
  assert(r.getLength() == f2.getLength());
  assert(r.getLength() == M.getNumRows());

  SpaMtrix::Vector temp(r.getLength());

  // timeMultiplier * M (q2 - q1)
  temp = q2;
  temp -= q1;
  temp *= timeMultiplier;
  SpaMtrix::multiply(M, temp, r);

  r += f1;
  r += f2;
}

LCSolverResult TimeSteppingLCSolver::solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom,
                                   SimulationState &simulationState) {

  if (K == nullptr) {
    initialiseMatrixSystem(q, geom);
  }

  Log::incrementIndent();

  Params params = {lc.A(), lc.B(), lc.C(), lc.L1(), lc.L2(), lc.L3(), lc.L6(), lc.deleps(), simulationState.dt(), lc.u1()};
  double maxDq = 0., firstDq = 0.;
  unsigned int iter = 0;

  //copyFreeNodes(q, *q1); // q1 is q from previous time step
  q.copyFreeDofsTo(*q1); // q1 is q from previous time step
  SpaMtrix::Vector r((*q1).getLength());
  SpaMtrix::Vector q2((*q1).getLength());

  SolutionVector q_prev = q;


  const double timeMultiplier = 2. * params.u1 / params.dt;
  bool solverConverged = true;
  string message = "Newton iterations: {}";
  do {
    //copyFreeNodes(q, q2); // make sure q2 is updated with the latest q values from previous iteration
    q.copyFreeDofsTo(q2); // make sure q2 is updated with the latest q values from previous iteration
    (*L) = 0;
    (*K) = 0;
    (*X) = 0;

    assembleVolumeTerms_impl(q, v, geom, params); // updates K and L using q2

    if (isFirstRun) {
      *f_prev = *L;
      isFirstRun = false;
    }
    calculateRhsResidualVector(r, *q1, q2, *f_prev, *L, *M, timeMultiplier);

    // X = K^-1 * r
    solverConverged = solverConverged && solveMatrixSystem(*K, r, *X);

    // q(m+1) = q(m) - dq(m)
    q.incrementFreeDofs(*X, -1.0); // decrement the result of the solver from the Q-tensor
    maxDq = maxAbs(*X);
    if (iter == 0) {
      firstDq = maxDq;
    }
    Log::info("Newton iteration={} maxDq={}", ++iter, maxDq);
  } while (maxDq > 1e-3);
  Log::decrementIndent();
  // save RHS vector for next time step
  *f_prev = *L;
  return {LCSolverType::TIME_STEPPING, iter, firstDq, solverConverged};
}

