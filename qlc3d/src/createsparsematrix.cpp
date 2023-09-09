//#include "qlc3d.h"
#include <time.h>
#include <geometry.h>
#include <material_numbers.h>
#include <solutionvector.h>
#include <solver-settings.h>
#include <list>
#include <vector>
#include <set>
#include <line.h>

#include <omp.h>
// SPAMTRIX INCLUDES
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_matrixmaker.hpp>
#include <util/logging.h>
using namespace std;

void setupSingleBlock(Geometry &geom,
                      SolutionVector &sol,
                      const idx &MatNum,
                      SpaMtrix::MatrixMaker &mm){
    // BOOK-KEEPING OF EQU-NODES
    const idx npt = geom.t->getnNodes(); // 4 FOR LINEAR TETS
    vector<idx> eqn(npt,0);              // KEEPS EQU NODES FOR ELEMENT

    for (idx it = 0 ; it < geom.t->getnElements() ; it++){
        // MAKE SURE ONLY ELEMENTS WITH MATERIAL NUMBER MatNum ARE USED
        // IF MatNum HAS BEEN DEFINED
        if (MatNum){
            if ( MatNum != geom.t->getMaterialNumber(it)  ){
                continue;
            }
        }

        //idx* nn = geom.t->getPtrToElement(it); // SHORTCUT TO ELEMENT NODES
        unsigned int nn[4];
        geom.getTetrahedra().loadNodes(it, nn);

        // CONVERT TO EQU NODES FOR THIS ELEMENT
        for (idx i = 0 ; i < npt ;++i){
            eqn[i] = sol.getEquNode( nn[i] );
        }

        // ADD EQU NODES TO MATRIX
        for (idx i = 0 ; i < npt ; ++i){
            // IGNORE FIXED NODES
            if ( eqn[i] == NOT_AN_INDEX ){
                continue;
            }

            mm.addNonZero(eqn[i],eqn[i]); // DIAGONAL i,i
            for (idx j = i+1 ; j < npt ; ++j ){
                if (eqn[j] == NOT_AN_INDEX){ // IGNORE FIXED NODES
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
    if (!electrodes.getCalcPot()) {
        Log::info("Potential calculation not required. Creating empty matrix for Potential solver.");
        return SpaMtrix::IRCMatrix();
    }

    const idx N = sol.getnFreeNodes();
    SpaMtrix::MatrixMaker mm(N,N);
    setupSingleBlock(geom, sol, MatNum, mm);
    idx nnz = mm.calcNumNonZeros();
    Log::info("Created sparse matrix of size {}x{}, with {} non-zeros for potential solver.", N, N, nnz);

    return mm.getIRCMatrix();
}

SpaMtrix::IRCMatrix createQMatrix(Geometry &geom,
                                  SolutionVector &q,
                                  const int &MatNum) {
    SpaMtrix::MatrixMaker mm(q.getnFreeNodes(),q.getnFreeNodes());
    const idx N = q.getnFreeNodes() * 5;
    setupSingleBlock(geom, q, MatNum, mm);  // CREATE SPARSITY PATTERN FOR COMPONENT q1
    mm.expandBlocks(4);                     // EXPAND SPARSITY PATTERN FOR q2->q5 COMPONENTS
    idx nnz = mm.calcNumNonZeros();
    Log::info("Created sparse matrix of size {}x{}, with {} non-zeros for Q-tensor solver.", N, N, nnz);

    return mm.getIRCMatrix();
}
