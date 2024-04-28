#include <lc/lc-solver.h>
#include <lc/lc-energy-terms.h>
#include <qlc3d.h>
#include <solver-settings.h>
#include <simulation-state.h>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <sparsematrix.h>
#include <mesh.h>
#include <fe/gaussian-quadrature.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <spamtrix_diagpreconditioner.hpp>
#include <spamtrix_iterativesolvers.hpp>
#include <util/logging.h>


LCSolver::~LCSolver() = default;

LCSolver::LCSolver(const LC &lc, const SolverSettings &solverSettings) :
A{lc.A()}, B{lc.B()}, C{lc.C()},
L1{lc.L1()}, L2{lc.L2()}, L3{lc.L3()}, L6{lc.L6()},
deleps{lc.deleps()},
solverSettings{solverSettings}
{
  Log::info("Creating steady state solver for Q-tensor");
  // The matrix is symmetric if all elastic constants are equal and there is no chirality
  isThreeElasticConstants = lc.K11() != lc.K22() || lc.K11() != lc.K33();
  isSymmetricMatrix = lc.p0() == 0.0;
}



double LCSolver::solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState) {
  initialiseMatrixSystem(q, geom, 0);
  assembleMatrixSystem(q, v, geom);
  solveMatrixSystem();
  double dq = updateQ(q);
  return dq;
}

void LCSolver::solveMatrixSystem() {
  // Solve linear system
  SpaMtrix::DiagPreconditioner M(*K);
  idx maxiter 	= solverSettings.getQ_GMRES_Maxiter();
  idx restart 	= solverSettings.getQ_GMRES_Restart();
  real toler    = solverSettings.getQ_GMRES_Toler();
  SpaMtrix::IterativeSolvers solver(maxiter, restart, toler);

  if (isSymmetricMatrix) {
    if (!solver.pcg(*K, *X, *L, M)) {
      Log::warn("PCG did not converge in {} iterations. Tolerance achieved is {}.", solver.maxIter, solver.toler);
    }
  } else {
    if (!solver.gmres(*K, *X, *L, M)) {
      Log::warn("GMRES did not converge in {} iterations. Tolerance achieved is {}.", solver.maxIter, solver.toler);
    }
  }
}

double LCSolver::updateQ(SolutionVector &q) {
  const idx npLC = q.getnDoF();
  double maxdqOut = 0.0;
  for (unsigned int i = 0; i < 5; i++) {
    for (idx j = 0; j < q.getnDoF(); j++) {
      const idx n = j + i * npLC;
      const idx effDoF = q.getEquNode(n);

      if (effDoF != NOT_AN_INDEX) {
        const double dqj = (*X)[ effDoF ]; // TODO: should we multiply by "damping" here, it looks like it's never used !!
        q[n] -= dqj;
        // KEEP TRACK OF LARGEST CHANGE IN Q-TENSOR
        maxdqOut = fabs(dqj) > maxdqOut ? fabs(dqj) : maxdqOut;
      }
    }
  }

  return maxdqOut;
}

void LCSolver::initialiseMatrixSystem(const SolutionVector &q, const Geometry &geom, double dt) {
  if (K != nullptr) {
    return;
  }

  const idx N = q.getnFreeNodes() * 5;

  K = createQMatrix(geom, q, MAT_DOMAIN1);
  L = std::make_unique<SpaMtrix::Vector>(N);
  X = std::make_unique<SpaMtrix::Vector>(N);
}

void LCSolver::assembleMatrixSystem(const SolutionVector &q, const SolutionVector &v, const Geometry &geom) {
  *(K) = 0.;
  *(L) = 0.;
  *(X) = 0.;

  assembleVolumeTerms(q, v, geom);
}

void LCSolver::assembleVolumeTerms(const SolutionVector &q, const SolutionVector &v, const Geometry &geom) {
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

    assembleLocalVolumeMatrix(indTet, lK, lL, tetNodes, tetDofs, shapes, q, v, geom);

    addToGlobalMatrix(lK, lL, q, tetNodes, tetDofs);
  }
}

void LCSolver::assembleLocalVolumeMatrix(const unsigned int indTet, double lK[20][20], double lL[20], unsigned int tetNodes[4], unsigned int *tetDofs,
                                         GaussianQuadratureTet<11> shapes, const SolutionVector &q,
                                         const SolutionVector &v, const Geometry &geom) {
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

  for (; shapes.hasNextPoint(); shapes.nextPoint()) {
    double q1, q2, q3, q4, q5;
    double Vx, Vy, Vz;
    shapes.sampleAll(qNodal, q1, q2, q3, q4, q5);
    Vx = shapes.sampleX(potential);
    Vy = shapes.sampleY(potential);
    Vz = shapes.sampleZ(potential);

    LcEnergyTerms::assembleThermotropic(lK, lL, shapes, tetDeterminant, q1, q2, q3, q4, q5, A, B, C);
    LcEnergyTerms::assembleElasticL1(lK, lL, shapes, tetDeterminant, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z, L1);
    LcEnergyTerms::assembleThreeElasticConstants(lK, lL, shapes, tetDeterminant, q1, q2, q3, q4, q5, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z, L2, L3, L6);

    if (Vx != 0.0 || Vy != 0.0 || Vz != 0.0) {
      LcEnergyTerms::assembleDielectric(lL, shapes, tetDeterminant, Vx, Vy, Vz, deleps);
    }
  }
}

void LCSolver::addToGlobalMatrix(double lK[20][20], double lL[20], const SolutionVector &q, const unsigned int tetNodes[4],
                                 const unsigned int tetDofs[4]) {
  idx npLC = q.getnDoF();
  for (int i = 0; i < 20; i++) {
    idx ri = tetNodes[i % 4] + npLC * (i / 4);
    idx eqr = q.getEquNode(ri);

    if (eqr == NOT_AN_INDEX) { // this row is fixed, so doesn't even exist in the system
      continue;
    }

    (*L)[eqr] += lL[i] * 2e16;
    for (int j = 0 ; j < 20 ; j++) { // LOOP OVER COLUMNS
      idx rj = tetNodes[j % 4] + npLC * (j / 4);
      idx eqc = q.getEquNode(rj);

      // Note: we can simply ignore fixed columns as they correspond to fixed zeroes in the RHS vector L,
      // so that their contribution to other rows in L would be zero anyway. We only set the free rows/columns
      // to the global matrix
      if (eqc != NOT_AN_INDEX) { // it's free
        K->sparse_add(eqr, eqc, lK[i][j] * 2e16);
      }
    }
  }
}
