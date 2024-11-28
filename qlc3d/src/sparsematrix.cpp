#include <sparsematrix.h>
#include <geometry.h>
#include <solutionvector.h>
#include <vector>

#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_matrixmaker.hpp>
#include <util/logging.h>
#include <dofmap.h>

using namespace std;

void setupSingleBlock(const Geometry &geom,
                      const DofMap &dofMap,
                      const idx &MatNum,
                      SpaMtrix::MatrixMaker &mm){
    // BOOK-KEEPING OF EQU-NODES
    auto &tets = geom.getTetrahedra();
    const idx npt = tets.getnNodes(); // 4 FOR LINEAR TETS
    vector<idx> eqn(npt,0);              // KEEPS EQU NODES FOR ELEMENT

    for (idx it = 0 ; it < tets.getnElements() ; it++){
        // MAKE SURE ONLY ELEMENTS WITH MATERIAL NUMBER MatNum ARE USED
        // IF MatNum HAS BEEN DEFINED
        if (MatNum){
            if ( MatNum != tets.getMaterialNumber(it)  ){
                continue;
            }
        }

        //idx* nn = geom.t->getPtrToElement(it); // SHORTCUT TO ELEMENT NODES
        unsigned int nn[4];
        geom.getTetrahedra().loadNodes(it, nn);

        // CONVERT TO EQU NODES FOR THIS ELEMENT
        for (idx i = 0 ; i < npt ;++i){
            eqn[i] = dofMap.getDof(nn[i]);
        }

        // ADD EQU NODES TO MATRIX
        for (idx i = 0 ; i < npt ; ++i){
            // IGNORE FIXED NODES
            if ( eqn[i] == DofMap::NOT_DOF) {
                continue;
            }

            mm.addNonZero(eqn[i],eqn[i]); // DIAGONAL i,i
            for (idx j = i+1 ; j < npt ; ++j ) {
                if (eqn[j] == DofMap::NOT_DOF) { // IGNORE FIXED NODES
                    continue;
                }
                mm.addNonZero(eqn[i], eqn[j]);
                mm.addNonZero(eqn[j], eqn[i]);
            }
        }
    }
}

SpaMtrix::IRCMatrix createPotentialMatrix(Geometry &geom,
                                          SolutionVector &sol,
                                          const int &MatNum,
                                          const Electrodes &electrodes) {
    if (!electrodes.isPotentialCalculationRequired()) {
        Log::info("Potential calculation not required. Creating empty matrix for Potential solver.");
        return SpaMtrix::IRCMatrix();
    }

    const idx N = sol.getnFreeNodes();
    SpaMtrix::MatrixMaker mm(N,N);
    setupSingleBlock(geom, sol.getDofMap(), MatNum, mm);
    idx nnz = mm.calcNumNonZeros();
    Log::info("Created sparse matrix of size {}x{}, with {} non-zeros for potential solver.", N, N, nnz);

    return mm.getIRCMatrix();
}

std::unique_ptr<SpaMtrix::IRCMatrix> createQMatrix(const Geometry &geom,
                                                   const SolutionVector &q,
                                                   const int& materialNumber) {
  SpaMtrix::MatrixMaker mm(q.getnFreeNodes(),q.getnFreeNodes());
  const idx N = q.getnFreeNodes() * 5;
  setupSingleBlock(geom, q.getDofMap(), materialNumber, mm);  // CREATE SPARSITY PATTERN FOR COMPONENT q1
  mm.expandBlocks(4);                     // EXPAND SPARSITY PATTERN FOR q2->q5 COMPONENTS
  idx nnz = mm.calcNumNonZeros();
  Log::info("Created sparse matrix of size {}x{}, with {} non-zeros for Q-tensor solver.", N, N, nnz);

  return std::unique_ptr<SpaMtrix::IRCMatrix>(mm.newIRCMatrix());
}

std::unique_ptr<SpaMtrix::IRCMatrix> createQMassMatrix(const Geometry &geom,
                                                       const SolutionVector &q,
                                                       const int& materialNumber) {
  SpaMtrix::MatrixMaker mm(q.getnFreeNodes(),q.getnFreeNodes());
  setupSingleBlock(geom, q.getDofMap(), materialNumber, mm); // Create sparsity pattern for q1 only
  mm.expandDiagonal(5);                       // Expand sparsity pattern for q2->q5 components along diagonal
  idx nnz = mm.calcNumNonZeros();
  SpaMtrix::IRCMatrix* massMatrix = mm.newIRCMatrix();

  Log::info("Created sparse mass-matrix of size {}x{}, with {} non-zeros for Q-tensor solver.",
            massMatrix->getNumRows(), massMatrix->getNumCols(), nnz);

  return std::unique_ptr<SpaMtrix::IRCMatrix>(massMatrix);
}

SpaMtrix::IRCMatrix createQMatrix(Geometry &geom,
                                  SolutionVector &q,
                                  const int &MatNum) {
    SpaMtrix::MatrixMaker mm(q.getnFreeNodes(),q.getnFreeNodes());
    const idx N = q.getnFreeNodes() * 5;
    setupSingleBlock(geom, q.getDofMap(), MatNum, mm);  // CREATE SPARSITY PATTERN FOR COMPONENT q1
    mm.expandBlocks(4);                     // EXPAND SPARSITY PATTERN FOR q2->q5 COMPONENTS
    idx nnz = mm.calcNumNonZeros();
    Log::info("Created sparse matrix of size {}x{}, with {} non-zeros for Q-tensor solver.", N, N, nnz);

    return mm.getIRCMatrix();
}
