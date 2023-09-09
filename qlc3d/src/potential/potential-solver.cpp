#include <potential/potential-solver.h>
#include <electrodes.h>
#include <util/logging.h>
#include <spamtrix_matrixmaker.hpp>
#include <memory>
#include <geometry.h>
#include <solutionvector.h>


void setupSingleBlock(const Geometry &geom,
                      const SolutionVector &sol,
                      SpaMtrix::MatrixMaker &mm){

  const idx nodesPerElement = geom.t->getnNodes();
  vector<idx> eqn(nodesPerElement, 0);

  for (idx it = 0; it < geom.t->getnElements(); it++) {
    unsigned int nn[4];
    geom.getTetrahedra().loadNodes(it, nn);

    for (idx i = 0; i < nodesPerElement; ++i) {
      eqn[i] = sol.getEquNode(nn[i]);
    }

    for (idx i = 0; i < nodesPerElement; ++i) {
      if (eqn[i] == NOT_AN_INDEX) {
        continue;
      }

      mm.addNonZero(eqn[i],eqn[i]); // Diagonal
      for (idx j = i + 1; j < nodesPerElement; ++j) {
        if (eqn[j] == NOT_AN_INDEX) { // IGNORE FIXED NODES
          continue;
        }
        mm.addNonZero(eqn[i], eqn[j]);
        mm.addNonZero(eqn[j], eqn[i]);
      }
    }
  }
}


void PotentialSolver::createPotentialMatrix(const Geometry &geom, const SolutionVector &sol) {
  if (!electrodes->getCalcPot()) {
    Log::info("Potential calculation not required. Creating empty matrix for Potential solver.");
    K = std::make_unique<SpaMtrix::IRCMatrix>();
  }

  const idx N = sol.getnFreeNodes();
  SpaMtrix::MatrixMaker mm(N,N);
  setupSingleBlock(geom, sol, mm);
  idx nnz = mm.calcNumNonZeros();
  Log::info("Created sparse matrix of size {}x{}, with {} non-zeros for potential solver.", N, N, nnz);

  K = std::unique_ptr<SpaMtrix::IRCMatrix>(mm.newIRCMatrix());
}

void PotentialSolver::solvePotential(SolutionVector &vOut,
                                     const SolutionVector &q,
                                     const Geometry &geom) {

  if (this->K == nullptr) {
    createPotentialMatrix(geom, vOut);
  }

  calcpot3d(*this->K, vOut, q, *lc, geom, *solverSettings, *electrodes);
}

void PotentialSolver::onGeometryChanged() {
  K = nullptr;
}