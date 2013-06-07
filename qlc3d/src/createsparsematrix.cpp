//#include "qlc3d.h"
#include <time.h>
#include <geometry.h>
#include <material_numbers.h>
#include <solutionvector.h>
#include <sparsematrix.h>
#include <settings.h>
#include <list>
#include <vector>
#include <set>
#include <line.h>

#include <omp.h>
// SPAMTRIX INCLUDES
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_matrixmaker.hpp>
using namespace std;

void print_dangly_list( vector < list <unsigned int> > dl){
    // Debug function that prints a dangly list matrix
    for (unsigned int i = 0 ; i < dl.size() ; i ++){
        printf("column %i :" , i );
        list < unsigned int> :: iterator itr;
        for (itr = dl[i].begin() ; itr != dl[i].end() ; ++itr){
            printf("%u ", *itr);
        }
        printf("\n");
    }
}

void print_dangly_set( vector < set <unsigned int> > dl){
    // Debug function that prints a dangly set matrix
    for (unsigned int i = 0 ; i < dl.size() ; i ++){
        printf("column %i :" , i );
        set < unsigned int> :: iterator itr;
        for (itr = dl[i].begin() ; itr != dl[i].end() ; ++itr){
            printf("%u ", *itr);
        }
        printf("\n");
    }
}

// CREATES DANGLY FROM LINES
void create_dangly_matrix(vector <Line>& lines,
                          vector< list <unsigned int> >& dangly)
{
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
                            const idx& MatNum)
{
    /* Creates a dangly sparse matrix of a Geometry, SolutionVector and material number */

    dangly.clear();
    list <idx> empty;


    dangly.assign( sol.getnFreeNodes() , empty ); // pre-allocate columns

    const idx npt = geom.t->getnNodes();
    vector <idx> eqn(npt, 0);

    eqn.resize(npt);
    for (idx it = 0 ; it < geom.t->getnElements() ; it++) // LOOP OVER EACH ELEMENT
    {
        if ( (!MatNum) ||    // if ignore material numebr OR
             ( MatNum == geom.t->getMaterialNumber(it))  )// if correct material
        {
            idx* nn = geom.t->getPtrToElement( it );    // shrotcut to element node indexes

            for (unsigned int i = 0 ; i < npt ; i++ )   // GET EQU NODES FOR THIS ELEMENT
                eqn[i] = sol.getEquNode( nn[i] );


            // LOOP OVER LOCAL NODES
            for (idx i = 0 ; i < npt ; i++)
            {
                // DONT INSERT FIXED NODES
                if (eqn[i] != NOT_AN_INDEX )
                {
                    dangly[eqn[i]].push_back( eqn[i] ); // DIAGONAL ENTRY

                    for (unsigned int j = i+1 ; j < npt ; j++)
                    {
                        // IF NODE IS NOT FIXED, CAN INSERT OFF DIAGONAL TOO
                        if (eqn[j] != NOT_AN_INDEX)
                        {
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

void convert_sets_to_arrays( vector<list <unsigned int> > &ds,
                             const unsigned int dim,
                             SparseMatrix& K)
{
    /*! Converts dangly matrix linked lists to a proper column copressed sparse matrix. The
        Sparse matrix is expanded by the factor 'dim'. i.e., the row and column count are multiplied
        by it and the sparsity pattern is copied to fill the new rows/cols.*/

    // CALCULATE AMOUNT OF MEMRY NEEDED
    // find nnz
    unsigned int nnz = 0;
    unsigned int ncols_s = ds.size();      // starting number of columns, before expansion
    unsigned int ncols_f = ds.size() * dim; // final number of columns after expansion
    for (size_t i = 0 ; i < ncols_s ; i++)
        nnz+=ds[i].size();

    nnz = dim*dim * nnz;

    // Allocate memory for arrays
    size_t Isize = nnz*sizeof(int);
    size_t Jsize = (ncols_f+1) * sizeof(int);
    int *Ir = (int*) malloc( Isize );
    int *Jc = (int*) malloc( Jsize );
    if ( (Ir == NULL) || (Jc==NULL)){
        printf("\nerror while creating sparse matrix. out of memory - bye!");
        exit(1);
    }

    // fill in data from dangly.
    // loop over columns and rows expanding the matrix to multiple dimensions
    nnz = 0;
    list <unsigned int>:: const_iterator r;
    for (unsigned int i = 0 ; i < ncols_f ; i ++){ // loop over columns

        Jc[i] = nnz;
        unsigned int col = i % ds.size();

        //unsigned int n_row_nnzs =ds[col].size(); // number of nonzeros in this column before expansion
        for (unsigned int rc = 0; rc < dim ; rc++ ){ // loop over row dimensions
            for (r = ds[ col ].begin(); r!=ds[ col ].end(); ++r){ // loop over rows, dim 1
                Ir[nnz] = *r + ( rc*ncols_s );
                nnz++;
            }// end for rows dim 1
        }// end for dims
    }// end for cols
    Jc[ncols_f] = nnz;

    ds.clear();

    K.MakeSparseMatrix(ncols_f,ncols_f, nnz, Ir, Jc);

}





void setupSingleBlock(Geometry &geom,
                      SolutionVector &sol,
                      const idx &MatNum,
                      SpaMtrix::MatrixMaker &mm)
{
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

        idx* nn = geom.t->getPtrToElement(it); // SHORTCUT TO ELEMENT NODES

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
                                          const Electrodes &electrodes)
{
    // CHECK WHETHER POTENTIAL WILL NEED TO BE CALCULATED
    // IF YES MAKE MATRIX, IF NOT CREATE EMPTY MATRIX
    if ( electrodes.getCalcPot() ){
        const idx N = sol.getnFreeNodes();
        cout << "Matrix size : " << N <<"x"<<N; fflush(stdout);
        SpaMtrix::MatrixMaker mm(N,N);
        setupSingleBlock(geom, sol, MatNum, mm);
        idx nnz = mm.calcNumNonZeros();
        cout << " nnz = " << nnz << endl;
        return mm.getIRCMatrix();
    }
    else // EMPTY MATRIX
    {
        cout << "Matrix size : 0x0"<< endl;
        return SpaMtrix::IRCMatrix();// Kpot;
    }
}



    SpaMtrix::IRCMatrix createQMatrix(Geometry &geom,
                                      SolutionVector &q,
                                      const int &MatNum)
    {

        SpaMtrix::MatrixMaker mm(q.getnFreeNodes(),q.getnFreeNodes());

        const idx N = q.getnFreeNodes() * 5;
        cout << "Matrix Size : " << N <<"x" << N; fflush(stdout);
        setupSingleBlock(geom, q, MatNum, mm);  // CREATE SPARSITY PATTERN FOR COMPONENT q1
        mm.expandBlocks(4);                     // EXPAND SPARSITY PATTERN FOR q2->q5 COMPONENTS
        idx nnz = mm.calcNumNonZeros();
        cout << " nnz = " << nnz << endl;

        return mm.getIRCMatrix();
    }

    SparseMatrix* createSparseMatrix(vector <Line>& lines){
        /*! Crates a SparseMatrix object where non-zeros are determined by the node indexes in the input parameter lines*/


        vector <list <idx> > dangly;
        create_dangly_matrix( lines , dangly); // CREATE DANGLY SET MATRIX OF NODE PAIRS

        SparseMatrix* K = new SparseMatrix();
        convert_sets_to_arrays( dangly, 1, *K );
        return K;
    }
