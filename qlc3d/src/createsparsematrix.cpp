//#include "qlc3d.h"
#include <time.h>
#include <geometry.h>
#include <solutionvector.h>
#include <sparsematrix.h>
#include <list>
#include <vector>
#include <set>
#include <line.h>

// SPAMTRIX INCLUDES
#include <ircmatrix.h>
#include <matrixmaker.h>
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
                          Geometry& geom,
                          SolutionVector& sol,
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
                             SparseMatrix& K){
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



SparseMatrix* createSparseMatrix(Geometry& geom,
                                 SolutionVector& sol,
                                 const int& MatNum){
    /*! Creates a SparseMatrix object based on the input geometry and number of dimensions in sol (i.e. is
sol is for Q, matrix size is [5*np X 5*np], and for V [npXnp].
 \param geom = input geometry
 \param sol = solution vector for which the sparse matrix is generated
 \param MatNum = Material number of elements that are to be included in the sparse matrix.
    Default value of MatNum is 0, which results in a sparse matrix for all elements (use this for potential solution).
 \return Returns pointer to created sparse matrix.
*/

    vector <list <idx> > dangly_set;
    create_dangly_matrix( dangly_set,
                          geom,
                          sol ,
                          MatNum );

    SparseMatrix* K = new SparseMatrix();
    convert_sets_to_arrays(dangly_set, sol.getnDimensions(),*K);
    return K;
}

IRCMatrix createPotentialMatrix(Geometry &geom,
                                SolutionVector &sol,
                                const int &MatNum)
{
    // MATRIX SIZE
    const idx N = sol.getnFreeNodes();
    MatrixMaker mm(N,N);
    cout << "Matrix size : " << N <<"x"<<N; fflush(stdout);
    // BOOK-KEEPING OF EQU-NODES
    const idx npt = geom.t->getnNodes(); // 4 FOR LINEAR TETS
    vector<idx> eqn(npt,0);              // KEEPS EQU NODES FOR ELEMENT

    for (idx it = 0 ; it < geom.t->getnElements() ; it++)
    {
        idx* nn = geom.t->getPtrToElement(it); // SHORTCUT TO ELEMENT NODES

        // CONVERT TO EQU NODES FOR THIS ELEMENT
        for (idx i = 0 ; i < npt ;++i)
            eqn[i] = sol.getEquNode( nn[i] );

        // ADD EQU NODES TO MATRIX
        for (idx i = 0 ; i < npt ; ++i)
        {
            // IGNORE FIXED NODES
            if ( eqn[i] == NOT_AN_INDEX )
                continue;

            mm.addNonZero(eqn[i],eqn[i]); // DIAGONAL
            for (idx j = i+1 ; j < npt ; ++j )
            {
                if (eqn[j] == NOT_AN_INDEX)
                    continue;
                mm.addNonZero(eqn[i], eqn[j]);
                mm.addNonZero(eqn[j], eqn[i]);
            }
        }
    }
    //IRCMatrix K = mm.getIRCMatrix();
    idx nnz = mm.calcNumNonZeros();
    cout << " nnz = " << nnz << endl;
    return mm.getIRCMatrix();


}


IRCMatrix createQMatrix(Geometry &geom,
                        SolutionVector &q,
                        const int &MatNum)
{

    // MATRIX SIZE
    const idx N = q.getnFreeNodes() * 5;
    MatrixMaker mm(N,N);
    cout << "Matrix Size : " << N <<"x" << N; fflush(stdout);
    //const idx nFree = q.getnFreeNodes(); // NUMBER OF FREE NODES
    const idx npLC  = q.getnDoF();
    // BOOK-KEEPING OF EQU-NODES
    const idx npt = geom.t->getnNodes() ; // 4 FOR LINEAR TETS
    vector<idx> eqn(npt,0);              // KEEPS EQU NODES FOR ELEMENT

    for (idx it = 0 ; it < geom.t->getnElements() ; it++)
    {
        if (geom.t->getMaterialNumber(it) != MatNum)
            continue;

        idx* nn = geom.t->getPtrToElement(it); // SHORTCUT TO ELEMENT NODES

        // CONVERT TO EQU NODES FOR THIS ELEMENT
        for (idx i = 0 ; i < npt ;++i)
            eqn[i] = q.getEquNode( nn[i] );

        for (idx i = 0 ; i < npt ; ++i)
        {
            if (eqn[i] == NOT_AN_INDEX)
                continue;



            // ADD DIAGONAL ENTRIES
#pragma omp parallel for
            for (idx rb = 0 ; rb < 5 ; ++rb)        // FOR ROW BLOCK
            {
                const idx eqr = q.getEquNode( nn[i] + rb*npLC);
                for (idx cb = 0 ; cb < 5 ; ++cb)    // FOR COL BLOCK
                {
                    const idx eqc = q.getEquNode(nn[i] + cb*npLC);
                    mm.addNonZero(eqr, eqc);
                }
            }
            for (idx j = i+1 ; j < npt ; ++j)
            {
                if (eqn[j] == NOT_AN_INDEX)
                    continue;
//#pragma omp parallel for
                for (idx rb = 0 ; rb < 5 ; ++rb)        // FOR ROW BLOCK
                {
                    const idx eqr = q.getEquNode( nn[i] + rb*npLC);
                    for (idx cb = 0 ; cb < 5 ; ++cb)    // FOR COL BLOCK
                    {
                        const idx eqc = q.getEquNode(nn[j] + cb*npLC);
                        mm.addNonZero(eqr, eqc);
                        mm.addNonZero(eqc, eqr);
                    }
                }// end for row block rb
            }// end for j
        }// end for i



    }// end for elements
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
