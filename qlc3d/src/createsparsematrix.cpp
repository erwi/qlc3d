//#include "qlc3d.h"
#include <time.h>
#include <geometry.h>
#include <material_numbers.h>
#include <solutionvector.h>
#include <settings.h>
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

// CREATES DANGLY FROM LINES
void create_dangly_matrix(vector <Line>& lines,
                          vector< list <unsigned int> >& dangly){
    /*! Creates a dangly sparse matrix of node pairs (i.e. lines)*/

    // FIND MAXIMUM NODE NUMBER IN lines
    size_t MaxNode = 0;
    for (size_t ind = 0 ; ind < lines.size() ; ind++ ){
        MaxNode = MaxNode > (size_t)lines[ind].L[1] ?  MaxNode: (size_t)lines[ind].L[1] ;
    }

    // PRE-ALLOCATE
    list <unsigned int> empty;
    dangly.clear();
    dangly.assign( MaxNode+1 , empty);

    // CREATE DIAGONAL 0->MaxNode
    for (size_t diag = 0 ; diag < MaxNode+1 ; diag++){
        dangly[diag].push_back( (unsigned int) diag);
    }

    // CREATE ACTUAL NNZ's
    for( size_t ind = 0 ; ind < lines.size() ; ind++){
        dangly[ lines[ind].L[0] ].push_back( lines[ind].L[1] );
        dangly[ lines[ind].L[1] ].push_back( lines[ind].L[0] );

    }

    // Remove repeated node indexes from columns. This step could be avoided if
    // sets was used instead of lists. However, sets (a tree data structure) tend
    // to (maybe) consume more memory...
#ifndef DEBUG
#pragma omp parallel for
#endif
    for( size_t ind = 0 ; ind < dangly.size() ; ind++){
        dangly[ind].sort();
        dangly[ind].unique();
    }
}

// CREATES DANGLY FROM GEOMETRY
void create_dangly_matrix(vector< list <idx> > & dangly,
                          const Geometry& geom,
                          const SolutionVector& sol,
                          const idx& MatNum){
    /* Creates a dangly sparse matrix of a Geometry, SolutionVector and material number */

    dangly.clear();
    list <idx> empty;
    dangly.assign( sol.getnFreeNodes() , empty ); // pre-allocate columns
    const idx npt = geom.t->getnNodes();
    vector <idx> eqn(npt, 0);
    eqn.resize(npt);

    for (idx it = 0 ; it < geom.t->getnElements() ; it++){ // LOOP OVER EACH ELEMENT
        if ( (!MatNum) ||    // if ignore material numebr OR
             ( MatNum == geom.t->getMaterialNumber(it))){// if correct material
            //idx* nn = geom.t->getPtrToElement( it );    // shrotcut to element node indexes
            idx nn[4];
            geom.getTetrahedra().loadNodes(it, nn);
            for (unsigned int i = 0 ; i < npt ; i++ ){   // GET EQU NODES FOR THIS ELEMENT
                eqn[i] = sol.getEquNode( nn[i] );
            }

            // LOOP OVER LOCAL NODES
            for (idx i = 0 ; i < npt ; i++){
                // DONT INSERT FIXED NODES
                if (eqn[i] != NOT_AN_INDEX ){
                    dangly[eqn[i]].push_back( eqn[i] ); // DIAGONAL ENTRY
                    for (unsigned int j = i+1 ; j < npt ; j++){
                        // IF NODE IS NOT FIXED, CAN INSERT OFF DIAGONAL TOO
                        if (eqn[j] != NOT_AN_INDEX){
                            dangly[eqn[i]].push_back( eqn[j] );
                            dangly[eqn[j]].push_back( eqn[i] );
                        }// end if j not fixed
                    }// end for node j
                }// end if i not fixed
            }// end for node i
        }// end if correct material
    }// end for tets



    // Remove repeated node indexes from columns. This step could be avoided if
    // sets was used instead of lists. However, sets (a tree data structure) tend
    // to (maybe) consume more memory...
#ifndef DEBUG
#   pragma omp parallel for
#endif
    for (unsigned int i = 0 ; i < dangly.size() ; i ++){
        dangly[i].sort(); // SORT
        dangly[i].unique(); // remove repetitions
    }
}

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
