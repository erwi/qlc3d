#include <potential/potential-solver.h>
#include <electrodes.h>
#include <util/logging.h>
#include <memory>
#include <geometry.h>
#include <lc.h>
#include <solutionvector.h>
#include <solver-settings.h>
#include <geom/vec3.h>
#include <geom/coordinates.h>
#include <dofmap.h>
#include <thread>
#include <fe/gaussian-quadrature.h>
#include <spamtrix_luincpreconditioner.hpp>
#include <spamtrix_iterativesolvers.hpp>
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_densematrix.hpp>
#include <fe/fe-util.h>
#include <util/exception.h>

//<editor-fold desc="private">
PotentialSolver::PotentialSolver(Electrodes& electrodes, std::shared_ptr<LC> lc, std::shared_ptr<SolverSettings> solverSettings):
  electrodes(electrodes), lc(lc), solverSettings(solverSettings) {
  K = nullptr;
  L = nullptr;
  V = nullptr;

  S0 = lc->S0();
  eper_lc = lc->eps_per();
  deleps = lc->eps_par() - lc->eps_per();
  efe = 2.0 / (9 * S0) * (lc->e1() + 2 * lc->e3());
  efe2 = 4.0 / (9 * S0 * S0) * (lc->e1() - lc->e3());
}

PotentialSolver::~PotentialSolver() = default;

void setupSingleBlock(const Mesh &tetMesh,
                      const DofMap &dofMap,
                      SpaMtrix::MatrixMaker &mm){
  const idx nodesPerElement = tetMesh.getnNodes();

  std::vector<unsigned int> elemNodes(nodesPerElement, 0);

  for (idx it = 0; it < tetMesh.getnElements(); it++) {
    tetMesh.loadNodes(it, &elemNodes[0]);

    for (idx i = 0; i < nodesPerElement; ++i) {
      if (dofMap.isFixedNode(elemNodes[i])) {
        continue;
      }
      unsigned int iDof = dofMap.getDof(elemNodes[i]);
      mm.addNonZero(iDof, iDof);
      for (idx j = i + 1; j < nodesPerElement; ++j) {
        if (dofMap.isFixedNode(elemNodes[j])) {
          continue;
        }
        unsigned int jDof = dofMap.getDof(elemNodes[j]);

        mm.addNonZero(iDof, jDof);
        mm.addNonZero(jDof, iDof);
      }
    }
  }
}

bool PotentialSolver::isPotentialSolutionRequired(const SolutionVector &v) const {
  return v.getnFixed() > 0 || electrodes.isPotentialCalculationRequired();
  // todo: possibly also required if flexoelectric effect is present
}

void PotentialSolver::setUniformEField(SolutionVector &vOut, const Coordinates &coordinates) const {
  Vec3 E = electrodes.getElectricField();
  double Emag = E.norm();
  Vec3 Ehat = E.normalized();

  //1. Calculate centre of structure
  auto bounds = coordinates.findBoundingBox();
  Vec3 centre = bounds.centre();

  //2. Loop over each node and calculate distance to centre along the E-filed direction
  for (unsigned int i = 0 ; i < coordinates.size() ; i++) {
    Vec3 pos = coordinates.getPoint(i);
    Vec3 vec = pos - centre;

    // want distance along EField, i.e. dot product
    double dist = vec.dot(Ehat);
    // set potential value as distance*magnitude
    double value = dist * Emag;
    vOut.setValue(i, 0, value);
  }
}

void PotentialSolver::createPotentialMatrix(const Geometry &geom, const DofMap &dofMap) {
  if (!electrodes.isPotentialCalculationRequired()) {
    Log::info("Potential calculation not required. Creating empty matrix for Potential solver.");
    K = std::make_unique<SpaMtrix::IRCMatrix>();
  }

  const idx N = dofMap.getnFreeNodes();
  SpaMtrix::MatrixMaker mm(N, N);
  setupSingleBlock(geom.getTetrahedra(), dofMap, mm);
  idx nnz = mm.calcNumNonZeros();
  Log::info("Created sparse matrix of size {}x{}, with {} non-zeros for potential solver.", N, N, nnz);

  K = std::unique_ptr<SpaMtrix::IRCMatrix>(mm.newIRCMatrix());
  L = std::make_unique<SpaMtrix::Vector>(N);
  V = std::make_unique<SpaMtrix::Vector>(N);
}

void PotentialSolver::initialiseMatrixSystem(const SolutionVector &vOut, const Geometry &geom) {
  if (K == nullptr) {
    createPotentialMatrix(geom, vOut.getDofMap());
  }

  *K = 0.0; // reset matrix to zero
  *L = 0.0; // reset RHS vector to zero

  // set initial potential values guess from output vector
  auto &dofMap = vOut.getDofMap();
  for (idx i = 0; i < dofMap.getnDof(); i++) {
    if (dofMap.isFreeNode(i)) {
      auto iDof = dofMap.getDof(i);
      (*V)[iDof] = vOut.getValue(i);
    }
  }
}

void PotentialSolver::assembleMatrixSystem(const SolutionVector &v, const SolutionVector &q, const Geometry &geom) {
  assembleVolume(v, q, geom);
  assembleNeumann(v, q, geom); // dielectric material test pass if this is commented out
}

void PotentialSolver::addToGlobalMatrix(const SpaMtrix::DenseMatrix &lK, const std::vector<double> &lL,
                                        const SolutionVector &v,
                                        const std::vector<unsigned int> &tetNodes,
                                        const std::vector<unsigned int> &tetDofs) {
#ifndef NDEBUG
  assert(lK.getNumRows() == tetNodes.size());
  assert(lK.getNumCols() == tetNodes.size());
  assert(lL.size() == tetNodes.size());
  assert(tetNodes.size() == tetDofs.size());
#endif

  // System formed is K[i,j] * v[j] = L[i] , where we'll be solving for v.
  // Some values of v[j] are known, i.e. the nodes are fixed. No rows/columns exist for these nodes and the fixed
  // values are multiplied and subtracted from the RHS vector L[i].

  // check if any dof is fixed and load the nodal values if required.
  const unsigned int numTetNodes = tetNodes.size();
  bool anyFixed = std::any_of(&tetDofs[0], &tetDofs[numTetNodes], [this](auto val){ return this->isFixedNode(val); });
  std::vector<double> values(numTetNodes, 0.);
  if (anyFixed) {
    v.loadValues(&tetNodes[0], &tetNodes[numTetNodes], &values[0]);
  }

  for (idx  i = 0; i < numTetNodes; i++) {
    const idx iDof = tetDofs[i];
    if (!isFreeNode(iDof)) { // this row/column doesn't even exist in the system as the value is already known
      continue;
    }

    double fixedContribution = 0.;
    for (idx j = 0; j < numTetNodes; j++) {
      const idx jDof = tetDofs[j];
      if (isFreeNode(jDof)) { // both i and j are free dofs, add the matrix contribution
        double *val = K->getValuePtr(iDof, jDof);
        if (val == nullptr) {
          RUNTIME_ERROR(fmt::format("Value at row {} and column {} is not found in the matrix for potential solution", iDof, jDof));
        }
        #pragma omp atomic
        *val += lK(i, j);
      } else { // j'th node is a fixed value. Add its contribution to the RHS vector L[i]
        fixedContribution += lK(i, j) * values[j];
      }
    }

    #pragma omp atomic
    (*L)[iDof] += lL[i] - fixedContribution;
  }
}

void PotentialSolver::assembleVolume(const SolutionVector &v, const SolutionVector &q, const Geometry &geom) {
  const unsigned int elementCount = geom.getTetrahedra().getnElements();
  const unsigned int nodesPerTet = geom.getTetrahedra().getnNodes();
  unsigned int elementOrder = getElementOrder(geom.getTetrahedra().getElementType());
  SpaMtrix::DenseMatrix lK(nodesPerTet, nodesPerTet);
  std::vector<double> lL(nodesPerTet, 0.0);
  std::vector<unsigned int> tetNodes(nodesPerTet, 0);
  std::vector<unsigned int> tetDofs(nodesPerTet, 0);

  TetShapeFunction shapes(elementOrder);
  shapes.setIntegrationPoints(Keast4);

  #pragma omp parallel for default(none) shared(geom, v, q, lc, K, L, electrodes, elementCount, nodesPerTet) firstprivate(lK, lL, tetNodes, tetDofs, shapes) schedule(guided)
  for (idx elementIndex = 0; elementIndex < elementCount; elementIndex++) {
    geom.getTetrahedra().loadNodes(elementIndex, &tetNodes[0]);
    v.loadEquNodes(&tetNodes[0], &tetNodes[nodesPerTet], &tetDofs[0]);

    localKL(geom, lK, lL, elementIndex, q, *lc, electrodes, shapes);

    addToGlobalMatrix(lK, lL, v, tetNodes, tetDofs);
  }
}

void PotentialSolver::assembleNeumann(const SolutionVector &v, const SolutionVector &q, const Geometry &geom) {

  auto &triMesh = geom.getTriangles();
  auto &tetMesh = geom.getTetrahedra();
  const unsigned int triCount = geom.getTriangles().getnElements();

  const unsigned int nodesPerTet = tetMesh.getnNodes();
  const unsigned int nodesPerTri = triMesh.getnNodes();
  SpaMtrix::DenseMatrix lK(nodesPerTet, nodesPerTet);
  std::vector<double> lL(nodesPerTet, 0.0);
  std::vector<unsigned int> triNodes(nodesPerTri, 0);
  std::vector<unsigned int> tetNodes(nodesPerTet, 0);
  std::vector<unsigned int> tetDofs(nodesPerTet, 0);

  #pragma omp parallel for default(none) shared(geom, v, q, lc, K, L, electrodes, triMesh, tetMesh, triCount, nodesPerTet) firstprivate(lK, lL, triNodes, tetNodes, tetDofs) schedule(guided)
  for (unsigned int indTri = 0; indTri < triCount; indTri++) {
    if (triMesh.getMaterialNumber(indTri) != MAT_NEUMANN) {
      continue;
    }

    if(tetMesh.getElementType() != ElementType::LINEAR_TETRAHEDRON) {
      throw std::runtime_error("Quadratic tetrahedra not supported for Neumann boundary conditions");
    }


    const unsigned int indTet = triMesh.getConnectedVolume(indTri);
    if (indTet == NOT_AN_INDEX) {
      throw std::runtime_error(fmt::format("Neumann boundary triangle {} not connected to a volume element", indTri));
    }
    const unsigned int tetMaterial = tetMesh.getMaterialNumber(indTet);
    if (!isLCMaterial(tetMaterial)) {
      continue;
    }
    GaussianQuadratureTet<7> shapes = gaussQuadratureTetBoundaryIntegral4thOrder();
    triMesh.loadNodes(indTri, &triNodes[0]);
    tetMesh.loadNodes(indTet, &tetNodes[0]);
    reorderBoundaryTetNodes(tetNodes, triNodes);
    v.loadEquNodes(&tetNodes[0], &tetNodes[nodesPerTet], &tetDofs[0]);

    localKLNeumann(geom.getCoordinates(), lK, lL, tetNodes, q,
                   triMesh.getDeterminant(indTri),
                   tetMesh.getDeterminant(indTet),
                   triMesh.getSurfaceNormal(indTri),
                   shapes);
    addToGlobalMatrix(lK, lL, v, tetNodes, tetDofs);
  }
}

bool PotentialSolver::isFixedNode(idx i) {
  return i == NOT_AN_INDEX;
}

bool PotentialSolver::isFreeNode(idx i) {
  return i < NOT_AN_INDEX;
}

void PotentialSolver::localKL(const Geometry &geom,
                              SpaMtrix::DenseMatrix &lK,
                              std::vector<double> &lL,
                              unsigned int elementIndex,
                              const SolutionVector &q,
                              const LC &lc,
                              const Electrodes &electrodes,
                              TetShapeFunction &s) {
                              //GaussianQuadratureTet<11> &s) {
  const unsigned int nodesPerTet = lL.size();
  const Mesh &mesh = geom.getTetrahedra();
  const bool isLCElement = isLCMaterial(mesh.getMaterialNumber(elementIndex));
  double eper = 0;

  if (!isLCElement) { // otherwise dielectric
    idx ind_de = mesh.getDielectricNumber(elementIndex) - 1; // -1 for 0 indexing
    eper = electrodes.getDielectricPermittivity(ind_de);
  }

  const bool isFlexoElectric = isLCElement && (efe != 0 || efe2 != 0);

  lK.setAllValuesTo(0.0);
  lL.assign(lL.size(), 0.0);

  std::vector<unsigned int> elemNodeInds(nodesPerTet, 0);
  std::vector<Vec3> elemCoords(nodesPerTet, Vec3());
  geom.getTetrahedra().loadNodes(elementIndex, &elemNodeInds[0]);
  geom.getCoordinates().loadCoordinates(&elemNodeInds[0], &elemNodeInds[nodesPerTet], &elemCoords[0]);
  for (auto &elemCoord : elemCoords) { elemCoord *= 1e-6; }

  const double determinant = geom.getTetrahedra().getDeterminant(elementIndex);

  //s.initialiseElement(&elemCoords[0], determinant);

  std::vector<qlc3d::TTensor> qNodal(nodesPerTet, qlc3d::TTensor());
  std::vector<qlc3d::DielectricPermittivity> eNodal(nodesPerTet, qlc3d::DielectricPermittivity());
  if (isLCElement) {
    q.loadQtensorValues(&elemNodeInds[0], &elemNodeInds[nodesPerTet], &qNodal[0]);

    for (unsigned int i = 0; i < nodesPerTet; i++) {
      eNodal[i] = qlc3d::DielectricPermittivity::fromTTensor(qNodal[i], S0, deleps, eper_lc);
    }
  }

  // for each Gauss integration point
  for (; s.hasNextPoint(); s.nextPoint()) {
    s.initialiseElement(&elemCoords[0], determinant); // TODO: for linear elements this can be done only once outside the loop. move this inside nextPoint??
    double mul = s.getWeight() * determinant;

    double exx = 0, eyy = 0, ezz = 0, exy = 0, exz = 0, eyz = 0;
    if (isLCElement) {
      s.sampleAll(&eNodal[0], exx, eyy, ezz, exy, exz, eyz);
    } else {
      for (unsigned int i = 0; i < nodesPerTet; i++) { exx += s.N(i) * eper; }
      eyy = ezz = exx;
      exy = exz = eyz = 0.;
    }

    if (isFlexoElectric) {
      double q1, q2, q3, q4, q5;
      s.sampleQ(&qNodal[0], q1, q2, q3, q4, q5);
      double q1x, q2x, q3x, q4x, q5x, q1y, q2y, q3y, q4y, q5y, q1z, q2z, q3z, q4z, q5z;
      s.sampleQX(&qNodal[0], q1x, q2x, q3x, q4x, q5x);
      s.sampleQY(&qNodal[0], q1y, q2y, q3y, q4y, q5y);
      s.sampleQZ(&qNodal[0], q1z, q2z, q3z, q4z, q5z);

      for (unsigned int i = 0; i < nodesPerTet; i++) {
        double ShRx = s.Nx(i);
        double ShRy = s.Ny(i);
        double ShRz = s.Nz(i);
        // Flexoelectric polarisation terms, minus, since minus residual formed
        lL[i] -= -efe * mul * ((ShRx * (q2x + q3y + q5z) + ShRy * (q3x + q4z - q2y)
                                + ShRz * (q4y + q5x)) / rt2 + (2 * ShRz * q1z - ShRx * q1x - ShRy * q1y) / rt6);

        lL[i] -= -efe2 * mul * ((ShRx * q1 * q1x + ShRy * q1 * q1y + 4 * ShRz * q1 * q1z) / 6
                                + (ShRx * (q2 * (q2x + q3y + q5z) + q3 * (q3x - q2y + q4z) + q5 * (q4y + q5x))
                                   + ShRy * (q2 * (q2y - q3x - q4z) + q3 * (q2x + q3y + q5z) + q4 * (q4y + q5x))
                                   + ShRz * (q4 * (q4z - q2y + q3x) + q5 * (q3y + q2x + q5z))) / 2
                                + (ShRx * (-q1 * (q2x + q3y + q5z) - q2 * q1x - q3 * q1y + 2 * q5 * q1z)
                                   + ShRy * (q1 * (q2y - q3x - q4z) + q2 * q1y - q3 * q1x + 2 * q4 * q1z)
                                   + ShRz * (2 * q1 * (q4y + q5x) - q4 * q1y - q5 * q1x)) * rt3 / 6);
      }
    }

    for (unsigned int i = 0; i < nodesPerTet; i++) {
      for (unsigned int j = 0; j < nodesPerTet; j++) {
        lK(i, j) -= mul * (
                s.Nx(i) * s.Nx(j) * exx +
                s.Ny(i) * s.Ny(j) * eyy +
                s.Nz(i) * s.Nz(j) * ezz +
                (s.Nx(i) * s.Ny(j) + s.Ny(i) * s.Nx(j)) * exy +
                (s.Nx(i) * s.Nz(j) + s.Nz(i) * s.Nx(j)) * exz +
                (s.Ny(i) * s.Nz(j) + s.Nz(i) * s.Ny(j)) * eyz);
      }
    }
  }
}

void PotentialSolver::localKLNeumann(
        const Coordinates &coordinates,
        SpaMtrix::DenseMatrix &lK,
        std::vector<double> &lL,
        const std::vector<unsigned int> &tetNodes,
        const SolutionVector &q,
        double triDet,
        double tetDet,
        const Vec3 &n,
        GaussianQuadratureTet<7> &shapes) {
  using std::vector;
  const unsigned int nodesPerTet = lL.size();
  const bool isFlexoelectric = (efe != 0 || efe2 != 0);

  lK.setAllValuesTo(0.0);
  lL.assign(lL.size(), 0.0);

  vector<Vec3> tetCoords(nodesPerTet, Vec3());
  coordinates.loadCoordinates(&tetNodes[0], &tetNodes[nodesPerTet], &tetCoords[0]);
  for (auto &tetCoord : tetCoords) { tetCoord *= 1e-6; }

  shapes.initialiseElement(&tetCoords[0], tetDet);

  vector<qlc3d::TTensor> qNodal(nodesPerTet, qlc3d::TTensor());
  q.loadQtensorValues(&tetNodes[0], &tetNodes[nodesPerTet], &qNodal[0]);

  vector<qlc3d::DielectricPermittivity> eNodal(nodesPerTet, qlc3d::DielectricPermittivity());
  for (unsigned int i = 0; i < nodesPerTet; i++) {
    eNodal[i] = qlc3d::DielectricPermittivity::fromTTensor(qNodal[i], S0, deleps, eper_lc);
  }

  for (;shapes.hasNextPoint(); shapes.nextPoint()) {
    double mul = shapes.weight() * triDet;

    if (isFlexoelectric) {
      double q1, q2, q3, q4, q5;
      shapes.sampleQ(qNodal, q1, q2, q3, q4, q5);
      double q1x, q2x, q3x, q4x, q5x, q1y, q2y, q3y, q4y, q5y, q1z, q2z, q3z, q4z, q5z;
      shapes.sampleQX(qNodal, q1x, q2x, q3x, q4x, q5x);
      shapes.sampleQY(qNodal, q1y, q2y, q3y, q4y, q5y);
      shapes.sampleQZ(qNodal, q1z, q2z, q3z, q4z, q5z);

      for (unsigned int i = 0; i < nodesPerTet; i++) {
        const double Ni = shapes.N(i);

        // Flexoelectric polarisation terms, minus, since minus residual formed
        lL[i] -= -efe * mul * Ni * ((n.x() * (q2x + q3y + q5z) + n.y() * (q3x + q4z - q2y)
                                     + n.z() * (q4y + q5x)) / rt2 + (2 * n.z() * q1z - n.x() * q1x - n.y() * q1y) / rt6);

        lL[i] -= -efe2 * mul * Ni * ((n.x() * q1 * q1x + n.y() * q1 * q1y + 4 * n.z() * q1 * q1z) / 6
                                     + (n.x() * (q2 * (q2x + q3y + q5z) + q3 * (q3x - q2y + q4z) + q5 * (q4y + q5x))
                                        + n.y() * (q2 * (q2y - q3x - q4z) + q3 * (q2x + q3y + q5z) + q4 * (q4y + q5x))
                                        + n.z() * (q4 * (q4z - q2y + q3x) + q5 * (q3y + q2x + q5z))) / 2
                                     + (n.x() * (-q1 * (q2x + q3y + q5z) - q2 * q1x - q3 * q1y + 2 * q5 * q1z)
                                        + n.y() * (q1 * (q2y - q3x - q4z) + q2 * q1y - q3 * q1x + 2 * q4 * q1z)
                                        + n.z() * (2 * q1 * (q4y + q5x) - q4 * q1y - q5 * q1x)) * rt3 / 6);
      }
    }

    double exx, eyy, ezz, exy, exz, eyz;
    shapes.samplePermittivity(eNodal, exx, eyy, ezz, exy, exz, eyz);

    for (unsigned int i = 0; i < nodesPerTet; i++) {
      double Ni = shapes.N(i);
      for (unsigned int j = 0; j < nodesPerTet; j++) {
        const double Njx = shapes.Nx(j);
        const double Njy = shapes.Ny(j);
        const double Njz = shapes.Nz(j);

        // negative as surface normal is pointing in towards LC, but convention in FE equations is outward
        lK(i, j) -= mul * Ni * (((exx - 1) * Njx + exy * Njy + exz * Njz) * n.x()
                                + (exy * Njx + (eyy - 1) * Njy + eyz * Njz) * n.y()
                                + (exz * Njx + eyz * Njy + (ezz - 1) * Njz) * n.z());
      }
    }
  }
}

void PotentialSolver::solveMatrixSystem(SolutionVector &v) {
  idx size = K->getNumRows();
  idx maxiter     = solverSettings->getV_GMRES_Maxiter();
  idx restart     = solverSettings->getV_GMRES_Restart();
  maxiter = maxiter < size ? maxiter : size;
  restart = restart < maxiter ? restart : maxiter;
  double toler = solverSettings->getV_GMRES_Toler();
  SpaMtrix::LUIncPreconditioner LU(*K); // TODO: consider updating this less frequently, it's likely to change slowly
  SpaMtrix::IterativeSolvers solver(maxiter, restart, toler);
  if (!solver.gmres(*K, *V, *L, LU)) {
    Log::warn("GMRES did not converge in {} iterations when solving for potential. Tolerance achieved is {}.",
              solver.maxIter, solver.toler);
  }

  // Copy free values of solution to correct location in output vector
  auto &dofMap = v.getDofMap();
  for (idx i = 0 ; i < v.getnDoF() ; i++) {
    if (dofMap.isFixedNode(i)) {
      continue;
    }

    v.setValue(i, 0, (*V)[dofMap.getDof(i)]);
  }
}

double calculateConditionNumber(const SpaMtrix::IRCMatrix &m) {
  int rows = m.getNumRows();

  double minDiag = abs(m.getValue(0, 0));
  double maxDiag = abs(m.getValue(0, 0));
  for (int i = 0; i < rows; i++) {
    double v = abs(m.getValue(i, i));
    minDiag = std::min(minDiag, v);
    maxDiag = std::max(maxDiag, v);
  }
  return maxDiag / minDiag;
}

//</editor-fold>

//<editor-fold desc="public">
void PotentialSolver::solvePotential(SolutionVector &vOut,
                                     const SolutionVector &q,
                                     const Geometry &geom) {

  if (!isPotentialSolutionRequired(vOut)) {
    vOut.setValuesTo(0);
    if (electrodes.hasElectricField()) {
      setUniformEField(vOut, geom.getCoordinates());
    }
    return;
  }

  //auto elementType = geom.getTetrahedra().getElementType();
  //if (elementType != ElementType::LINEAR_TETRAHEDRON) {
  //  throw NotYetImplementedException("Potential solver only supports linear tetrahedral elements, got elementType " + toString(elementType));
  //}

  initialiseMatrixSystem(vOut, geom);

  assembleMatrixSystem(vOut, q, geom);

  solveMatrixSystem(vOut);
}

void PotentialSolver::onGeometryChanged() {
  K = nullptr;
  L = nullptr;
  V = nullptr;
}

//</editor-fold>