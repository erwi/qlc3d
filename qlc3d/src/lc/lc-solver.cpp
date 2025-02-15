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
#include <dofmap.h>
#include <util/logging.h>
#include <util/exception.h>
#include <util/stopwatch.h>
#include <spamtrix_luincpreconditioner.hpp>

// <editor-fold ImplicitLCSolver>
ImplicitLCSolver::ImplicitLCSolver(const LC &lc, const SolverSettings &solverSettings, const Alignment &alignment) : lc(lc),
  solverSettings(solverSettings),
  alignment(alignment),
  isChiral(lc.L4() != 0.0),
  isSymmetricMatrix(lc.p0() == 0.0),
  isThreeElasticConstants(lc.K11() != lc.K22() || lc.K11() != lc.K33()) {
  Log::info("Creating steady state solver for Q-tensor");
}

bool ImplicitLCSolver::solveMatrixSystem(const SpaMtrix::IRCMatrix &Kmatrix, const SpaMtrix::Vector &r, SpaMtrix::Vector &x) const {
  // check that matrix and vector sizes match
  if (Kmatrix.getNumRows() != r.getLength() || Kmatrix.getNumCols() != r.getLength() || Kmatrix.getNumCols() != x.getLength()) {
    RUNTIME_ERROR("Matrix and vector sizes do not match");
  }

  idx maxiter 	=  solverSettings.getQ_GMRES_Maxiter();
  idx restart 	=  solverSettings.getQ_GMRES_Restart();
  real toler    =  solverSettings.getQ_GMRES_Toler();
  SpaMtrix::IterativeSolvers solver(maxiter, restart, toler);

  bool converged;
  std::string method = "Unknown";
  if (isSymmetricMatrix) {
    method = "PCG";
    SpaMtrix::DiagPreconditioner I(Kmatrix);
    converged = solver.pcg(Kmatrix, x, r, I);
  } else {
    method = "GMRES";
    SpaMtrix::LUIncPreconditioner LU(Kmatrix);
    converged = solver.gmres(Kmatrix, x, r, LU);
  }

  if (!converged) {
    Log::warn("{} Solver did not converge in {} iterations to required tolerance {}. Tolerance achieved is {}.", method, solver.maxIter, toler, solver.toler);
  }
  return converged;
}

void ImplicitLCSolver::assembleLocalVolumeMatrix(const unsigned int indTet,
                                                 SpaMtrix::DenseMatrix &lK,
                                                 std::vector<double> &lL,
                                                 std::vector<unsigned int> &tetNodes,
                                                 std::vector<unsigned int> &tetDofs,
                                                 TetShapeFunction &shapes,
                                                 const SolutionVector &q,
                                                 const SolutionVector &v,
                                                 const Geometry &geom,
                                                 const LCSolverParams &params) {

  const unsigned int npe = shapes.getNumPointsPerElement();
  lK.setAllValuesTo(0.0);
  memset(&lL[0], 0., lL.size() * sizeof(double));

  std::vector<Vec3> tetCoords(npe, Vec3());
  geom.getCoordinates().loadCoordinates(tetNodes.data(), tetNodes.data() + npe, tetCoords.data());
  for (auto &coord : tetCoords) {
    coord *= 1e-6;
  }

  const double tetDeterminant = geom.getTetrahedra().getDeterminant(indTet);

#ifndef NDEBUG
  if (npe == 10) {
    // sanity checks
    // calculate the volume enclosed by the tetrahedron defined by the first 4 nodes in tetCoords
    auto vec1 = tetCoords[3] - tetCoords[0];
    auto vec2 = tetCoords[2] - tetCoords[0];
    auto vec3 = tetCoords[1] - tetCoords[0];
    const double volume = std::abs(vec1.dot(vec2.cross(vec3)));

    double diff = std::abs(tetDeterminant - volume);
    double eps = 1e-12;
    assert(diff < eps);

    auto p5 = 0.5 * (tetCoords[0] + tetCoords[1]);
    auto p6 = 0.5 * (tetCoords[1] + tetCoords[2]);
    auto p7 = 0.5 * (tetCoords[0] + tetCoords[2]);
    auto p8 = 0.5 * (tetCoords[0] + tetCoords[3]);
    auto p9 = 0.5 * (tetCoords[1] + tetCoords[3]);
    auto p10 = 0.5 * (tetCoords[2] + tetCoords[3]);

    assert(p5.equals(tetCoords[4], eps));
    assert(p6.equals(tetCoords[5], eps));
    assert(p7.equals(tetCoords[6], eps));
    assert(p8.equals(tetCoords[7], eps));
    assert(p9.equals(tetCoords[8], eps));
    assert(p10.equals(tetCoords[9], eps));
  }
#endif

  std::vector<qlc3d::TTensor> qNodal(npe, qlc3d::TTensor());
  q.loadQtensorValues(tetNodes.data(), tetNodes.data() + npe, qNodal.data());

  std::vector<double> potential(npe, 0.);
  v.loadValues(tetNodes.data(), tetNodes.data() + npe, potential.data());

  for (; shapes.hasNextPoint(); shapes.nextPoint()) {
    shapes.initialiseElement(tetCoords.data(), tetDeterminant);
    double q1, q2, q3, q4, q5;
    double Vx, Vy, Vz;
    double q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z;
    shapes.sampleQ(qNodal, q1, q2, q3, q4, q5);
    shapes.sampleQX(qNodal, q1x, q2x, q3x, q4x, q5x);
    shapes.sampleQY(qNodal, q1y, q2y, q3y, q4y, q5y);
    shapes.sampleQZ(qNodal, q1z, q2z, q3z, q4z, q5z);
    Vx = shapes.sampleX(potential.data());
    Vy = shapes.sampleY(potential.data());
    Vz = shapes.sampleZ(potential.data());

    LcEnergyTerms::assembleThermotropic(&lK, lL, shapes, tetDeterminant, q1, q2, q3, q4, q5, params.A, params.B, params.C);

    LcEnergyTerms::assembleElasticL1(&lK, lL, shapes, tetDeterminant, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z, params.L1);

    if (isThreeElasticConstants) {
      LcEnergyTerms::assembleThreeElasticConstants(&lK, lL, shapes, tetDeterminant, q1, q2, q3, q4, q5, q1x, q1y, q1z,
                                                   q2x, q2y, q2z, q3x, q3y, q3z, q4x, q4y, q4z, q5x, q5y, q5z,
                                                   params.L2, params.L3, params.L6);
    }

    if (isChiral) {
      LcEnergyTerms::assembleChiralTerms(&lK, lL, shapes, tetDeterminant,
                                         q1, q2, q3, q4, q5,
                                         q1x, q1y, q1z,
                                          q2x, q2y, q2z,
                                          q3x, q3y, q3z,
                                          q4x, q4y, q4z,
                                          q5x, q5y, q5z,
                                         params.L4);
    }

    if (Vx != 0.0 || Vy != 0.0 || Vz != 0.0) {
      LcEnergyTerms::assembleDielectric(lL, shapes, tetDeterminant, Vx, Vy, Vz, params.deleps);
    }
  }
}

void ImplicitLCSolver::assembleLocalWeakAnchoringMatrix(unsigned int indTri,
                                                        SpaMtrix::DenseMatrix &lK,
                                                        std::vector<double> &lL,
                                                        const std::vector<unsigned int> &triNodes,
                                                        const std::vector<unsigned int> &triDofs,
                                                        TriShapeFunction &shapes,
                                                        const SolutionVector &q,
                                                        const Geometry &geom,
                                                        const Surface &surface,
                                                        double surfaceOrder) {
  lK.setAllValuesTo(0);
  memset(lL.data(), 0., lL.size() * sizeof(double));

  const unsigned int npe = lL.size() / 5; // nodes per element

  std::vector<Vec3> triCoords(npe, Vec3());
  geom.getCoordinates().loadCoordinates(triNodes.data(), triNodes.data() + npe, triCoords.data());
  for(auto &coord : triCoords) {
    coord *= 1e-6;
  }

  const double triDeterminant = geom.getTriangles().getDeterminant(indTri);

  //shapes.initialiseElement(triCoords, triDeterminant);

  // anchoring principal axes vectors for each node
  std::vector<Vec3> vec1(npe, Vec3());
  std::vector<Vec3> vec2(npe, Vec3());
  if (surface.usesSurfaceNormal()) {
    assert(surface.getK1() == 0.);
    for (unsigned int i = 0; i < npe; i++) {
      vec2[i] = geom.getNodeNormal(triNodes[i]);
    }
  } else {
    // use the principal axes of the element anchoring at each node
    for(unsigned  int i = 0; i < npe; i++) {
      vec1[i] = surface.getV1();
      vec2[i] = surface.getV2();
    }
  }

  const double W = surface.getAnchoringType() == WeakHomeotropic ? -surface.getStrength() : surface.getStrength();
  const double K1 = surface.getK1();
  const double K2 = surface.getK2();
  std::vector<qlc3d::TTensor> qNodal(npe, qlc3d::TTensor());
  q.loadQtensorValues(triNodes.data(), triNodes.data() + npe, qNodal.data());

  for (; shapes.hasNextPoint(); shapes.nextPoint()) {
    double q1, q2, q3, q4, q5;
    shapes.sampleQ(qNodal.data(), q1, q2, q3, q4, q5);

    Vec3 v1, v2;
    if (!surface.usesSurfaceNormal()) {
      shapes.sample(vec1, v1);
    }
    shapes.sample(vec2, v2);

    LcEnergyTerms::addWeakAnchoring(&lK, lL, shapes, triDeterminant, q1, q2, q3, q4, q5, v1, v2, W, K1, K2, surfaceOrder);
  }
}

void ImplicitLCSolver::addToGlobalMatrix(const SpaMtrix::DenseMatrix &lK,
                                         const std::vector<double> &lL,
                                         const DofMap &dofMap,
                                         const std::vector<unsigned int> &elemNodes) {
  const idx npLC = dofMap.getnDof();
  const idx elemNodeCount = elemNodes.size();
  const idx len = lL.size();
  for (unsigned int i = 0; i < len; i++) {
    idx ri = elemNodes[i % elemNodeCount] + npLC * (i / elemNodeCount);

    if (dofMap.isFixedNode(ri)) {
      continue;
    }
    idx eqr = dofMap.getDof(ri);

    #pragma omp atomic
    (*L)[eqr] += lL[i];
    for (unsigned int j = 0 ; j < len ; j++) { // LOOP OVER COLUMNS
      idx rj = elemNodes[j % elemNodeCount] + npLC * (j / elemNodeCount);

      // Note: we can simply ignore fixed columns as they correspond to fixed zeroes in the RHS vector L,
      // so that their contribution to other rows in L would be zero anyway. We only set the free rows/columns
      // to the global matrix
      if (dofMap.isFixedNode(rj)) {
        continue;
      }
      idx eqc = dofMap.getDof(rj);

      double* value = K->getValuePtr(eqr, eqc);

      if (value == nullptr) {
        RUNTIME_ERROR(fmt::format("Value at row {} and column {} is not found in the matrix for LC solution", eqr, eqc));
      }

      #pragma omp atomic
      *value += lK(i, j);
    }
  }
}

void ImplicitLCSolver::assembleMatrixSystemVolumeTerms(const SolutionVector &q, const SolutionVector &v, const Geometry &geom, const LCSolverParams &params) {
  const unsigned int elementCount = geom.getTetrahedra().getnElements();
  const Mesh &tets = geom.getTetrahedra();

  const ElementType elementType = tets.getElementType();
  const unsigned int elementNodeCount = tets.getnNodes();

  SpaMtrix::DenseMatrix lK(elementNodeCount * 5, elementNodeCount * 5);
  std::vector<double> lL(elementNodeCount * 5, 0.0);
  std::vector<unsigned int> tetNodes(elementNodeCount, 0);
  std::vector<unsigned int> tetDofs(elementNodeCount, 0);

  TetShapeFunction shapes(getElementOrder(elementType));
  //shapes.setIntegrationPoints(Keast4); // TODO: need higher order integration points
  shapes.setIntegrationPoints(Keast8);

  #pragma omp parallel for firstprivate(shapes, lK, lL, tetNodes, tetDofs) schedule(guided)
  for (unsigned int indTet = 0; indTet < elementCount; indTet++) {
    if (tets.getMaterialNumber(indTet) != MAT_DOMAIN1) {
      continue;
    }
    tets.loadNodes(indTet, &tetNodes[0]);
    q.loadEquNodes(&tetNodes[0], &tetNodes[elementNodeCount], &tetDofs[0]);

    assembleLocalVolumeMatrix(indTet, lK, lL, tetNodes, tetDofs, shapes, q, v, geom, params);
    addToGlobalMatrix(lK, lL, q.getDofMap(), tetNodes);
  }
}

void ImplicitLCSolver::assembleMatrixSystemWeakAnchoring(const SolutionVector &q, const Geometry &geom,
                                                         const LCSolverParams &params) {
  auto &tris = geom.getTriangles();
  const auto elementType = tris.getElementType();
  const unsigned int elementCount = tris.getnElements();
  const unsigned int nodesPerElement = tris.getnNodes();

  TriShapeFunction shapes(getElementOrder(elementType));
  shapes.setIntegrationPoints(Tri4thOrder);

  SpaMtrix::DenseMatrix lK(nodesPerElement * 5, nodesPerElement * 5);
  std::vector<double> lL(nodesPerElement * 5, 0.0);
  std::vector<idx> triNodes(nodesPerElement, 0);
  std::vector<idx> triDofs(nodesPerElement, 0);

  std::unordered_map<unsigned int, Surface> weakSurfaces = alignment.getWeakSurfacesByFixLcNumber();

  #pragma omp parallel for firstprivate(lL, lK, triNodes, triDofs, shapes) schedule(guided)
  for (unsigned int indTri = 0; indTri < elementCount; indTri++) {
    unsigned int fixLcNumber = tris.getFixLCNumber(indTri);

    // check if a weak surface exists by the FixLC number
    auto weakSurface = weakSurfaces.find(fixLcNumber);
    if (weakSurface == weakSurfaces.end()) {
      continue;
    }
    tris.loadNodes(indTri, triNodes.data());
    q.loadEquNodes(triNodes.data(), triNodes.data() + nodesPerElement, triDofs.data());

    // assemble local matrix for weak anchoring
    assembleLocalWeakAnchoringMatrix(indTri, lK, lL, triNodes, triDofs, shapes, q, geom, weakSurface->second, params.S0);
    // add to global matrix
    addToGlobalMatrix(lK, lL, q.getDofMap(), triNodes);
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
        L1{lc.L1()}, L2{lc.L2()}, L3{lc.L3()}, L4{lc.L4()}, L6{lc.L6()},
        deleps{lc.deleps()}
        {
  Log::info("Creating steady state solver for Q-tensor with isThreeElasticConstants={}, isChiral={}, isSymmetricMatrix={}", isThreeElasticConstants, isChiral, isSymmetricMatrix);
}

LCSolverResult SteadyStateLCSolver::solve(SolutionVector &q, const SolutionVector &v, const Geometry &geom, SimulationState &simulationState) {
  initialiseMatrixSystem(q, geom, 0);

  LCSolverParams params = {A, B, C, L1, L2, L3, L4, L6, deleps, simulationState.dt(), 0, lc.S0()}; // todo: this is replication

  Stopwatch assemblyStopwatch;
  Stopwatch solverStopwatch;

  assemblyStopwatch.start();
  assembleMatrixSystem(q, v, geom, params);
  assemblyStopwatch.stop();

  // X = K^-1 * r
  solverStopwatch.start();
  bool solverConverged = solveMatrixSystem(*K, *L, *X);
  solverStopwatch.stop();

  // q(m+1) = q(m) - dq(m)
  q.incrementFreeDofs(*X, -1.0); // decrement the result of the solver from the Q-tensor

  return {
    LCSolverType::STEADY_STATE,
    1,
    maxAbs(*X),
    solverConverged,
    false, // only a single iteration is always run
    {assemblyStopwatch.elapsedSeconds(), solverStopwatch.elapsedSeconds()}
  };
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

TimeSteppingLCSolver::TimeSteppingLCSolver(const LC &lc, const SolverSettings &solverSettings, double maxError, const Alignment &alignment, int maxNewtonIterations) :
ImplicitLCSolver(lc, solverSettings, alignment), maxError(maxError), maxNewtonIterations(maxNewtonIterations) {
  Log::info("Creating time stepping solver for Q-tensor, isTreeElasticConstants={}, isSymmetricMatrix={}, maxNewtonIterations={}",
            isThreeElasticConstants, isSymmetricMatrix, maxNewtonIterations);
}

bool isFixedNode(idx i) {
  return i == NOT_AN_INDEX;
}

void addToGlobalMassMatrix(const SpaMtrix::DenseMatrix &lK, const std::vector<unsigned int> &tetDofs, SpaMtrix::IRCMatrix &M) {
  const unsigned int npt = lK.getNumRows();
  const unsigned int numFreeNodes = M.getNumRows() / 5;

  for (idx i = 0; i < npt; i++) {
    const idx iDof = tetDofs[i];
    if (isFixedNode(iDof)) { continue; } // this row should not even exist in the matrix

    for (idx j = 0; j < npt; j++) {
      const idx jDof = tetDofs[j];
      if (isFixedNode(jDof)) { continue; }

      for (int dof = 0; dof < 5; dof ++) {
        M.sparse_add(iDof + dof * numFreeNodes, jDof + dof * numFreeNodes, lK(i, j));
      }
    }
  }
}

void assembleElementMassMatrix(SpaMtrix::DenseMatrix &lK, const std::vector<idx> &tetNodes, double tetDeterminant, const Geometry &geom, TetShapeFunction &shapes) {
  const idx npt = tetNodes.size();
  assert(npt == shapes.getNumPointsPerElement());
  assert(lK.getNumRows() == npt);
  assert(lK.getNumCols() == npt);

  lK.setAllValuesTo(0.0);
  std::vector<Vec3> coordinates(npt, Vec3());
  geom.getCoordinates().loadCoordinates(tetNodes.data(), tetNodes.data() + npt, coordinates.data());
  for (auto &coord : coordinates) {
    coord *= 1e-6;
  }

  for (; shapes.hasNextPoint(); shapes.nextPoint()) {
    shapes.initialiseElement(coordinates.data(), tetDeterminant);
    const double mul = shapes.getWeight() * tetDeterminant;
    for (idx i = 0; i < npt; i++) {
      for (idx j = 0; j < npt; j++) {
        lK(i, j) += mul * shapes.N(i) * shapes.N(j);
      }
    }
  }
}

void assembleGlobalMassMatrix(SpaMtrix::IRCMatrix &M, const Geometry &geom, const SolutionVector &q) {
  const unsigned int elementOrder = getElementOrder(geom.getTetrahedra().getElementType());
  TetShapeFunction shapes = TetShapeFunction(elementOrder);
  shapes.setIntegrationPoints(Keast4); // 4th order is enough for both linear and quadratic elements

  const Mesh tetMesh = geom.getTetrahedra();
  const unsigned int tetCount = tetMesh.getnElements();
  const unsigned int npt = tetMesh.getnNodes();
  SpaMtrix::DenseMatrix lK(npt, npt);
  std::vector<idx> tetNodes(npt, 0);
  std::vector<idx> tetDofs(npt, 0);

  for (idx tetIndex = 0; tetIndex < tetCount; tetIndex++) {
    if (tetMesh.getMaterialNumber(tetIndex) != MAT_DOMAIN1) {
      continue;
    }

    tetMesh.loadNodes(tetIndex, tetNodes.data());
    q.loadEquNodes(tetNodes.data(), tetNodes.data() + npt, tetDofs.data());

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
  LCSolverParams params = {lc.A(), lc.B(), lc.C(), lc.L1(), lc.L2(), lc.L3(), lc.L4(), lc.L6(), lc.deleps(), simulationState.dt(), lc.u1(), lc.S0()};
  double maxDq = 0., firstDq = 0.;
  int iter = 0;

  q.copyFreeDofsTo(*q1); // q1 is q from previous time step
  SpaMtrix::Vector r((*q1).getLength());
  SpaMtrix::Vector q2((*q1).getLength());

  // prediction based on dqdt from previous time step: q = q + dt * dqdt
  *dqdt *= params.dt;
  q.incrementFreeDofs(*dqdt);

  const double timeMultiplier = 2. * params.u1 / params.dt;
  bool solverConverged = true;

  Stopwatch assemblyStopwatch;
  Stopwatch solverStopwatch;

  do {
    q.copyFreeDofsTo(q2); // make sure q2 is updated with the latest q values from previous iteration

    assemblyStopwatch.start();
    assembleMatrixSystem(q, v, geom, params); // updates K and L using q2

    if (isFirstRun) {
      *f_prev = *L;
      isFirstRun = false;
    }
    modifySystemForTimeStepping(r, *K, *q1, q2, *f_prev, *L, *M, timeMultiplier);

    assemblyStopwatch.stop();

    // X = K^-1 * r
    solverStopwatch.start();
    solverConverged = solverConverged && solveMatrixSystem(*K, r, *X);
    solverStopwatch.stop();
    maxDq = maxAbs(*X);

    if (iter == 0) {
      firstDq = maxDq;
      Log::enableInfoNewline(false);
      Log::info("Newton iterations. dQ={:.4e}", maxDq);
    } else {
      Log::append_info(",{} {:.4e}", iter == maxNewtonIterations ? "|" : "", maxDq);
    }

    // q(m+1) = q(m) - dq(m)
    q.incrementFreeDofs(*X, -1); // decrement the result of the solver from the Q-tensor

    iter++;

  } while (maxDq > maxError && iter < maxNewtonIterations);
  Log::enableInfoNewline(true);

  if (iter >= maxNewtonIterations) {
    Log::warn("Newton iterations: max iterations {} reached", maxNewtonIterations);
  }

  // save RHS vector for next time step
  *f_prev = *L;

  // calculate dq/dt and save it for next time step
  q.copyFreeDofsTo(*dqdt);
  *dqdt -= *q1;
  *dqdt *= 1.0 / params.dt;
  double maxDqdt = maxAbs(*dqdt);
  Log::append_info(". dQ/dt={:.4e}", maxDqdt);

  return {
    LCSolverType::TIME_STEPPING,
    iter,
    firstDq,
    solverConverged,
    iter >= maxNewtonIterations,
    {assemblyStopwatch.elapsedSeconds(), solverStopwatch.elapsedSeconds()}
  };
}
// </editor-fold>