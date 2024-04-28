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
#include "spamtrix_diagpreconditioner.hpp"
#include "spamtrix_iterativesolvers.hpp"
#include "util/logging.h"

struct Params {
  double A;
  double B;
  double C;
  double L1;
  double L2;
  double L3;
  double L6;
  double deleps;
};

TimeSteppingLCSolver::TimeSteppingLCSolver(const LC &lc, const SolverSettings &solverSettings) :
  lc(lc),
  solverSettings(solverSettings) {
  Log::info("Creating time stepping solver for Q-tensor");
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

void assembleLocalMassMatrix(double lK[4][4], idx tetNodes[4], double tetDeterminant, const Geometry &geom, GaussianQuadratureTet<11> &shapes) {
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

void assembleMassMatrix(SpaMtrix::IRCMatrix &M, const Geometry &geom, const SolutionVector &q) {
  GaussianQuadratureTet<11> shapes = gaussQuadratureTet4thOrder();
  const Mesh tetMesh = geom.getTetrahedra();
  const idx tetCount = tetMesh.getnElements();
  //const idx numFreeNodes = q.getnFreeNodes();
  double lK[4][4];
  idx tetNodes[4];
  idx tetDofs[4];

  for (idx tetIndex = 0; tetIndex < tetCount; tetIndex++) {
    if (tetMesh.getMaterialNumber(tetIndex) != MAT_DOMAIN1) {
      continue;
    }

    tetMesh.loadNodes(tetIndex, tetNodes);
    q.loadEquNodes(tetNodes, tetNodes + 4, tetDofs);

    assembleLocalMassMatrix(lK, tetNodes, tetMesh.getDeterminant(tetIndex), geom, shapes);

    addToGlobalMatrix(lK, tetDofs, M);
  }
}

void assembleLocalRhsVector(double lL[20], const idx tetNodes[4], double tetDeterminant,
                             const Geometry &geom, const SolutionVector &q, const SolutionVector &v,
                             GaussianQuadratureTet<11> &shapes, const Params &params) {
  memset(lL, 0, 20 * sizeof(double));
  Vec3 tetCoords[4];
  geom.getCoordinates().loadCoordinates(tetNodes, tetNodes + 4, tetCoords);
  tetCoords[0] *= 1e-6;
  tetCoords[1] *= 1e-6;
  tetCoords[2] *= 1e-6;
  tetCoords[3] *= 1e-6;

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

    LcEnergyTerms::assembleThermotropic(nullptr, lL, shapes, tetDeterminant, q1, q2, q3, q4, q5, params.A, params.B, params.C);
    LcEnergyTerms::assembleElasticL1(nullptr, lL, shapes, tetDeterminant, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z, params.L1);
    //LcEnergyTerms::assembleThreeElasticConstants(nullptr, lL, shapes, tetDeterminant, q1, q2, q3, q4, q5, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z, params.L2, params.L3, params.L6);

    if (Vx != 0.0 || Vy != 0.0 || Vz != 0.0) {
      LcEnergyTerms::assembleDielectric(lL, shapes, tetDeterminant, Vx, Vy, Vz, params.deleps);
    }
  }


}


void assembleRhsVector(SpaMtrix::Vector &L, const SolutionVector &q, const SolutionVector &v, const Geometry &geom, const Params &params) {
  GaussianQuadratureTet<11> shapes = gaussQuadratureTet4thOrder();
  const Mesh tetMesh = geom.getTetrahedra();
  const idx tetCount = tetMesh.getnElements();
  double lL[20];
  idx tetNodes[4];
  const idx npLC = q.getnDoF();

  for (idx tetIndex = 0; tetIndex < tetCount; tetIndex++) {
    if (tetMesh.getMaterialNumber(tetIndex) != MAT_DOMAIN1) {
      continue;
    }

    tetMesh.loadNodes(tetIndex, tetNodes);

    assembleLocalRhsVector(lL, tetNodes, tetMesh.getDeterminant(tetIndex), geom, q, v, shapes, params);

    for (idx i = 0; i < 20; i++) {
      idx ri = tetNodes[i % 4] + npLC * (i / 4);
      idx eqr = q.getEquNode(ri);

      if (isFixedNode(eqr)) { continue; } // this row should not even exist in the matrix

      L[eqr] += lL[i];
    }
  }
}

void TimeSteppingLCSolver::initialiseMatrixSystem(const SolutionVector &q, const Geometry &geom, double dt) {
  if (K != nullptr) {
    return;
  }

  K = createQMassMatrix(geom, q, MAT_DOMAIN1);

  L = std::make_unique<SpaMtrix::Vector>(K->getNumRows());
  X = std::make_unique<SpaMtrix::Vector>(K->getNumRows());
}

void solveMatrixSystema(const SpaMtrix::IRCMatrix &M, const SpaMtrix::Vector &L, SpaMtrix::Vector &X) {
  SpaMtrix::DiagPreconditioner I(M);
  idx maxiter 	= M.getNumRows(); //solverSettings.getQ_GMRES_Maxiter();
  idx restart 	= 0; //solverSettings.getQ_GMRES_Restart();
  real toler    = 1e-9; //solverSettings.getQ_GMRES_Toler();
  SpaMtrix::IterativeSolvers solver(maxiter, restart, toler);

  if (!solver.pcg(M, X, L, I)) {
    Log::warn("PCG did not converge in {} iterations. Tolerance achieved is {}.", solver.maxIter, solver.toler);
  }
}

double updateSolutionVector(SolutionVector &q, const SpaMtrix::Vector &X, double dt, double viscosity) {
  const idx npLC = q.getnDoF();
  double maxdqOut = 0.0;
  for (unsigned int i = 0; i < 5; i++) {
    for (idx j = 0; j < q.getnDoF(); j++) {
      const idx n = j + i * npLC;
      const idx effDoF = q.getEquNode(n);

      if (effDoF != NOT_AN_INDEX) {
        const double dqj = X[ effDoF ]; // TODO: should we multiply by "damping" here, it looks like it's never used !!
        q[n] -= (dt / viscosity) * dqj;
        // KEEP TRACK OF LARGEST CHANGE IN Q-TENSOR
        maxdqOut = fabs(dqj) > maxdqOut ? fabs(dqj) : maxdqOut;
      }
    }
  }

  return maxdqOut;
}


double TimeSteppingLCSolver::solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom,
                                   SimulationState &simulationState) {

  if (K == nullptr) {
    initialiseMatrixSystem(q, geom, simulationState.dt());
    assembleMassMatrix(*K, geom, q);
  }

  L->setAllValuesTo(0);
  X->setAllValuesTo(0);

  Params params = {lc.A(), lc.B(), lc.C(), lc.L1(), lc.L2(), lc.L3(), lc.L6(), lc.deleps()};
  assembleRhsVector(*L, q, v, geom, params);

  solveMatrixSystema(*K, *L, *X);

  double dq = updateSolutionVector(q, *X, simulationState.dt(), 0.1);

  return dq;
}

