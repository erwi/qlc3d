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
#include <spamtrix_blas.hpp>
#include <util/logging.h>
#include <util/exception.h>

// <editor-fold ImplicitLCSolver>
ImplicitLCSolver::ImplicitLCSolver(const LC &lc, const SolverSettings &solverSettings, const Alignment &alignment) : lc(lc),
  solverSettings(solverSettings),
  alignment(alignment),
  isSymmetricMatrix(lc.p0() == 0.0),
  isThreeElasticConstants(lc.K11() != lc.K22() || lc.K11() != lc.K33()) {
  Log::info("Creating steady state solver for Q-tensor");
}

bool ImplicitLCSolver::solveMatrixSystem(const SpaMtrix::IRCMatrix &Kmatrix, const SpaMtrix::Vector &r, SpaMtrix::Vector &x) const {
  // check that matrix and vector sizes match
  if (Kmatrix.getNumRows() != r.getLength() || Kmatrix.getNumCols() != r.getLength() || Kmatrix.getNumCols() != x.getLength()) {
    RUNTIME_ERROR("Matrix and vector sizes do not match");
  }

  SpaMtrix::DiagPreconditioner I(Kmatrix);
  idx maxiter 	=  solverSettings.getQ_GMRES_Maxiter();
  idx restart 	=  solverSettings.getQ_GMRES_Restart();
  real toler    =  solverSettings.getQ_GMRES_Toler();
  SpaMtrix::IterativeSolvers solver(maxiter, restart, toler);

  bool converged;
  std::string method = "Unknown";
  if (isSymmetricMatrix) {
    method = "PCG";
    converged = solver.pcg(Kmatrix, x, r, I);
  } else {
    method = "GMRES";
    converged = solver.gmres(Kmatrix, x, r, I);
  }

  if (!converged) {
    Log::warn("{} Solver did not converge in {} iterations to required tolerance {}. Tolerance achieved is {}.", method, solver.maxIter, toler, solver.toler);
  }
  return converged;
}

void ImplicitLCSolver::assembleLocalVolumeMatrix(const unsigned int indTet, double lK[20][20], double lL[20], unsigned int tetNodes[4], unsigned int *tetDofs,
                                                 GaussianQuadratureTet<11> shapes, const SolutionVector &q,
                                                 const SolutionVector &v, const Geometry &geom, const LCSolverParams &params) {
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

    LcEnergyTerms::assembleThermotropic(lK, lL, shapes, tetDeterminant, q1, q2, q3, q4, q5, params.A, params.B, params.C);

    LcEnergyTerms::assembleElasticL1(lK, lL, shapes, tetDeterminant, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z, params.L1);

    if (isThreeElasticConstants) {
      LcEnergyTerms::assembleThreeElasticConstants(lK, lL, shapes, tetDeterminant, q1, q2, q3, q4, q5, q1x, q1y, q1z,
                                                   q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z,
                                                   params.L2, params.L3, params.L6);
    }

    if (Vx != 0.0 || Vy != 0.0 || Vz != 0.0) {
      LcEnergyTerms::assembleDielectric(lL, shapes, tetDeterminant, Vx, Vy, Vz, params.deleps);
    }
  }
}

void ImplicitLCSolver::assembleLocalWeakAnchoringMatrix(unsigned int indTri, double lK[15][15], double lL[15],
                                                        unsigned int triNodes[3], unsigned int triDofs[3],
                                                        GaussianQuadratureTri<7> shapes, const SolutionVector &q,
                                                        const Geometry &geom, const Surface &surface,
                                                        double surfaceOrder) {
  memset(lK, 0., 15 * 15 * sizeof(double));
  memset(lL, 0., 15 * sizeof(double));

  Vec3 triCoords[3];
  geom.getCoordinates().loadCoordinates(triNodes, triNodes + 3, triCoords);
  triCoords[0] *= 1e-6;
  triCoords[1] *= 1e-6;
  triCoords[2] *= 1e-6;

  const double triDeterminant = geom.getTriangles().getDeterminant(indTri);
  shapes.initialiseElement(triCoords, triDeterminant);

  Vec3 vec1[3];
  Vec3 vec2[3];
  if (surface.usesSurfaceNormal()) {
    assert(surface.getK1() == 0.);
    vec2[0] = geom.getNodeNormal(triNodes[0]);
    vec2[1] = geom.getNodeNormal(triNodes[1]);
    vec2[2] = geom.getNodeNormal(triNodes[2]);
  } else {
    // use the principal axes of the element anchoring at each node
    vec1[0] = surface.getV1();
    vec1[1] = surface.getV1();
    vec1[2] = surface.getV1();
    vec2[0] = surface.getV2();
    vec2[1] = surface.getV2();
    vec2[2] = surface.getV2();
  }

  const double W = surface.getAnchoringType() == WeakHomeotropic ? -surface.getStrength() : surface.getStrength();
  const double K1 = surface.getK1();
  const double K2 = surface.getK2();
  qlc3d::TTensor qNodal[3];
  q.loadQtensorValues(triNodes, triNodes + 3, qNodal);

  for (; shapes.hasNextPoint(); shapes.nextPoint()) {
    double q1, q2, q3, q4, q5;
    shapes.sampleAll(qNodal, q1, q2, q3, q4, q5);

    Vec3 v1, v2;
    if (!surface.usesSurfaceNormal()) {
      shapes.sample(vec1, v1);
    }
    shapes.sample(vec2, v2);

    LcEnergyTerms::addWeakAnchoring(lK, lL, shapes, triDeterminant, q1, q2, q3, q4, q5, v1, v2, W, K1, K2, surfaceOrder);
  }
}

void ImplicitLCSolver::addToGlobalMatrix(double* lK, double* lL, const SolutionVector &q,
                                         const unsigned int* elemNodes, int elemNodeCount) {
  //const int elemNodeCount = 4;
  idx npLC = q.getnDoF();
  for (int i = 0; i < 5 * elemNodeCount; i++) {
    idx ri = elemNodes[i % elemNodeCount] + npLC * (i / elemNodeCount);
    idx eqr = q.getEquNode(ri);

    if (eqr == NOT_AN_INDEX) { // this row is fixed, so doesn't even exist in the system
      continue;
    }

    (*L)[eqr] += lL[i];
    for (int j = 0 ; j < 5 * elemNodeCount ; j++) { // LOOP OVER COLUMNS
      idx rj = elemNodes[j % elemNodeCount] + npLC * (j / elemNodeCount);
      idx eqc = q.getEquNode(rj);

      // Note: we can simply ignore fixed columns as they correspond to fixed zeroes in the RHS vector L,
      // so that their contribution to other rows in L would be zero anyway. We only set the free rows/columns
      // to the global matrix
      if (eqc != NOT_AN_INDEX) { // it's free
        // [i][j] = j + 5 * elemNodeCount * i
        K->sparse_add(eqr, eqc, lK[j + 5 * elemNodeCount*i]);
      }
    }
  }
}

void ImplicitLCSolver::assembleMatrixSystemVolumeTerms(const SolutionVector &q, const SolutionVector &v, const Geometry &geom, const LCSolverParams &params) {
  const int elementNodeCount = 4;
  GaussianQuadratureTet<11> shapes = gaussQuadratureTet4thOrder();
  const unsigned int elementCount = geom.getTetrahedra().getnElements();
  const Mesh &tets = geom.getTetrahedra();
  double lK[elementNodeCount * 5][elementNodeCount * 5];
  double lL[elementNodeCount * 5];
  idx tetNodes[elementNodeCount];
  idx tetDofs[elementNodeCount];

  for (unsigned int indTet = 0; indTet < elementCount; indTet++) {
    if (tets.getMaterialNumber(indTet) != MAT_DOMAIN1) {
      continue;
    }
    tets.loadNodes(indTet, tetNodes);
    q.loadEquNodes(tetNodes, tetNodes + elementNodeCount, tetDofs);

    assembleLocalVolumeMatrix(indTet, lK, lL, tetNodes, tetDofs, shapes, q, v, geom, params);

    addToGlobalMatrix(&lK[0][0], lL, q, tetNodes, elementNodeCount);
  }
}

void ImplicitLCSolver::assembleMatrixSystemWeakAnchoring(const SolutionVector &q, const Geometry &geom,
                                                         const LCSolverParams &params) {
  const int elementNodeCount = 3;
  auto shapes = gaussianQuadratureTri4thOrder();
  auto &tris = geom.getTriangles();

  double lK[elementNodeCount * 5][elementNodeCount * 5];
  double lL[elementNodeCount * 5];
  idx triNodes[elementNodeCount];
  idx triDofs[elementNodeCount];

  std::unordered_map<unsigned int, Surface> weakSurfaces = alignment.getWeakSurfacesByFixLcNumber();

  for (unsigned int indTri = 0; indTri < tris.getnElements(); indTri++) {
    unsigned int fixLcNumber = tris.getFixLCNumber(indTri);

    // check if a weak surface exists by the FixLC number
    auto weakSurface = weakSurfaces.find(fixLcNumber);
    if (weakSurface == weakSurfaces.end()) {
      continue;
    }
    tris.loadNodes(indTri, triNodes);
    q.loadEquNodes(triNodes, triNodes + elementNodeCount, triDofs);

    // assemble local matrix for weak anchoring
    assembleLocalWeakAnchoringMatrix(indTri, lK, lL, triNodes, triDofs, shapes, q, geom, weakSurface->second, params.S0);
    // add to global matrix

    addToGlobalMatrix(&lK[0][0], lL, q, triNodes, elementNodeCount);
  }
}

void ImplicitLCSolver::assembleMatrixSystem(const SolutionVector &q, const SolutionVector &v, const Geometry &geom, const LCSolverParams &params) {
  *(K) = 0.;
  *(L) = 0.;
  *(X) = 0.;

  assembleMatrixSystemVolumeTerms(q,v, geom, params);

  if (alignment.weakSurfacesExist()) {
    assembleMatrixSystemWeakAnchoring(q, geom, params);
  }
}


double ImplicitLCSolver::maxAbs(const SpaMtrix::Vector &v) const {
  double max = 0.0;
  for (idx i = 0; i < v.getLength(); i++) {
    max = std::max(max, fabs(v[i]));
  }
  return max;
}

// </editor-fold>

// <editor-fold SteadyStateLCSolver>
SteadyStateLCSolver::~SteadyStateLCSolver() = default;

SteadyStateLCSolver::SteadyStateLCSolver(const LC &lc, const SolverSettings &solverSettings, const Alignment &alignment) :
        ImplicitLCSolver(lc, solverSettings, alignment),
        A{lc.A()}, B{lc.B()}, C{lc.C()},
        L1{lc.L1()}, L2{lc.L2()}, L3{lc.L3()}, L6{lc.L6()},
        deleps{lc.deleps()}
        {
  Log::info("Creating steady state solver for Q-tensor with isThreeElasticConstants={}, isSymmetricMatrix={}", isThreeElasticConstants, isSymmetricMatrix);
}

LCSolverResult SteadyStateLCSolver::solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState) {
  initialiseMatrixSystem(q, geom, 0);

  LCSolverParams params = {A, B, C, L1, L2, L3, L6, deleps, simulationState.dt(), 0, lc.S0()}; // todo: this is replication

  assembleMatrixSystem(q, v, geom, params);

  // X = K^-1 * r
  bool solverConverged = solveMatrixSystem(*K, *L, *X);

  // q(m+1) = q(m) - dq(m)
  q.incrementFreeDofs(*X, -1.0); // decrement the result of the solver from the Q-tensor

  return {LCSolverType::STEADY_STATE, 1, maxAbs(*X), solverConverged};
}

void SteadyStateLCSolver::initialiseMatrixSystem(const SolutionVector &q, const Geometry &geom, double dt) {
  if (K != nullptr) {
    return;
  }

  const idx N = q.getnFreeNodes() * 5;

  K = createQMatrix(geom, q, MAT_DOMAIN1);
  L = std::make_unique<SpaMtrix::Vector>(N);
  X = std::make_unique<SpaMtrix::Vector>(N);
}
// </editor-fold>

// <editor-fold TimeSteppingLCSolver>

TimeSteppingLCSolver::TimeSteppingLCSolver(const LC &lc, const SolverSettings &solverSettings, double maxError, const Alignment &alignment) :
ImplicitLCSolver(lc, solverSettings, alignment), maxError(maxError) {
  Log::info("Creating time stepping solver for Q-tensor, isTreeElasticConstants={}, isSymmetricMatrix={}", isThreeElasticConstants, isSymmetricMatrix);
}

bool isFixedNode(idx i) {
  return i == NOT_AN_INDEX;
}

void addToGlobalMassMatrix(const double lK[4][4], idx tetDofs[4], SpaMtrix::IRCMatrix &M) {
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

    addToGlobalMassMatrix(lK, tetDofs, M);
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
  dqdt = std::make_unique<SpaMtrix::Vector>(N);
  f_prev = std::make_unique<SpaMtrix::Vector>(N);
  *f_prev = 0;
}

/**
 * Modifies the matrix system with extra terms required by implicit time stepping
 * @param r the RHS modified vector is written  here
 * @param K the Jacobian matrix, is modified in place
 * @param q1 LC solution from previous time step
 * @param q2 most recent LC solution
 * @param f1 Euler-Lagrange function from previous time step (calculated from q1)
 * @param f2 most recent Euler-Lagrange function (calculated from q2)
 * @param M pre-computed mass matrix
 * @param timeMultiplier scalar multiplier, includes time step size (= 2 * u1 / dt)
 */
void modifySystemForTimeStepping(SpaMtrix::Vector &r,
                                 SpaMtrix::IRCMatrix &K,
                                 const SpaMtrix::Vector &q1,
                                 const SpaMtrix::Vector &q2,
                                 const SpaMtrix::Vector &f1,
                                 const SpaMtrix::Vector &f2,
                                 const SpaMtrix::IRCMatrix &M,
                                 double timeMultiplier) {

  assert(r.getLength() == q1.getLength());
  assert(r.getLength() == q2.getLength());
  assert(r.getLength() == f1.getLength());
  assert(r.getLength() == f2.getLength());
  assert(r.getLength() == M.getNumRows());

  // Modify RHS vector for time stepping
  SpaMtrix::Vector temp(r.getLength());
  // timeMultiplier * M (q2 - q1)
  temp = q2;
  temp -= q1;
  temp *= timeMultiplier;

  SpaMtrix::multiply(M, temp, r);

  r += f1;
  r += f2;

  // Modify Jacobian matrix for time stepping
  K.add(M, timeMultiplier);
}

LCSolverResult TimeSteppingLCSolver::solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom,
                                           SimulationState &simulationState) {

  if (K == nullptr) {
    initialiseMatrixSystem(q, geom);
  }

  // todo: params is replicated, move to super class
  LCSolverParams params = {lc.A(), lc.B(), lc.C(), lc.L1(), lc.L2(), lc.L3(), lc.L6(), lc.deleps(), simulationState.dt(), lc.u1(), lc.S0()};
  double maxDq = 0., firstDq = 0.;
  unsigned int iter = 0;

  q.copyFreeDofsTo(*q1); // q1 is q from previous time step
  SpaMtrix::Vector r((*q1).getLength());
  SpaMtrix::Vector q2((*q1).getLength());

  // prediction based on dqdt from previous time step: q = q + dt * dqdt
  *dqdt *= params.dt;
  q.incrementFreeDofs(*dqdt);

  const double timeMultiplier = 2. * params.u1 / params.dt;
  bool solverConverged = true;
  string message = "Newton iterations: {}";
  do {
    q.copyFreeDofsTo(q2); // make sure q2 is updated with the latest q values from previous iteration
    assembleMatrixSystem(q, v, geom, params); // updates K and L using q2

    if (isFirstRun) {
      *f_prev = *L;
      isFirstRun = false;
    }
    modifySystemForTimeStepping(r, *K, *q1, q2, *f_prev, *L, *M, timeMultiplier);

    // X = K^-1 * r
    solverConverged = solverConverged && solveMatrixSystem(*K, r, *X);

    // q(m+1) = q(m) - dq(m)
    q.incrementFreeDofs(*X, -1.0); // decrement the result of the solver from the Q-tensor
    maxDq = maxAbs(*X);
    if (iter == 0) {
      firstDq = maxDq;
      Log::enableInfoNewline(false);
      Log::info("Newton iterations. dQ={:.4e}", maxDq);
    } else {
      Log::append_info(", {:.4e}", maxDq);
    }
    iter++;
  } while (maxDq > maxError);
  Log::enableInfoNewline(true);

  // save RHS vector for next time step
  *f_prev = *L;

  // calculate dq/dt and save it for next time step
  q.copyFreeDofsTo(*dqdt);
  *dqdt -= *q1;
  *dqdt *= 1.0 / params.dt;
  double maxDqdt = maxAbs(*dqdt);
  Log::append_info(". dQ/dt={:.4e}", maxDqdt);

  return {LCSolverType::TIME_STEPPING, iter, firstDq, solverConverged};
}
// </editor-fold>