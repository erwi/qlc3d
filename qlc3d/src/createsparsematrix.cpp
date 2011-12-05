//#include "qlc3d.h"
#include <time.h>
#include <geometry.h>
#include <solutionvector.h>
#include <sparsematrix.h>
#include <list>
#include <vector>
#include <set>
#include <line.h>
using namespace std;

void print_dangly_list( vector < list <unsigned int> > dl){
    // Debug function that prints a dangly list matrix
    for (unsigned int i = 0 ; i < dl.size() ; i ++){
        printf("column %i :" , i );
        list < unsigned int> :: iterator itr;
        for (itr = dl[i].begin() ; itr != dl[i].end() ; itr ++){
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
        for (itr = dl[i].begin() ; itr != dl[i].end() ; itr ++){
            printf("%u ", *itr);
        }
        printf("\n");
    }
}


void create_dangly_matrix(vector <Line>& lines, vector< list <unsigned int> >& dangly){
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
#       pragma omp parallel for
        for( size_t ind = 0 ; ind < dangly.size() ; ind++){
            dangly[ind].sort();
            dangly[ind].unique();
        }



}

void create_dangly_matrix(vector< list <unsigned int> > & dangly,
                            Geometry& geom,
                            SolutionVector& sol,
                            const int& MatNum)
{
    /* Creates a dangly sparse matrix of a Geometry, SolutionVector and material number */

    dangly.clear();
    list <unsigned int> empty;


    dangly.assign( sol.getnFreeNodes() , empty ); // pre-allocate columns

    const unsigned int npt = (unsigned int ) geom.t->getnNodes();
    //int* eqn = new int[ npt ]; // LOCAL EQU NODES
    vector<int> eqn(npt, 0);
    eqn.resize(npt);
    for (size_t it = 0 ; it < (size_t) geom.t->getnElements() ; it++) // LOOP OVER EACH ELEMENT
    {
        if ((!MatNum) ||    // if ignore material numebr OR
           ( MatNum == geom.t->getMaterialNumber(it))  ){// if correct material


            int* nn = geom.t->getPtrToElement( it );    // shrotcut to element node indexes

            for (unsigned int i = 0 ; i < npt ; i++ )   // GET EQU NODES FOR THIS ELEMENT
                eqn[i] = sol.getEquNode( nn[i] );


            // LOOP OVER LOCAL NODES
            for (unsigned int i = 0 ; i < npt ; i++){
                // ALWAYS INSERT DIAGONAL
                if (eqn[i] != SolutionVector::FIXED_NODE )
                {
                    dangly[eqn[i]].push_back( eqn[i] ); // DIAGONAL ENTRY

                for (unsigned int j = i+1 ; j < npt ; j++)
                {
                    // IF NODE IS NOT FIXED, CAN INSERT OFF DIAGONAL TOO
                    //if (    (!sol.getIsFixed( nn[i] ) ) ||
                    //        ( ! sol.getIsFixed( nn[j] ) ) )

                    if (eqn[j] != SolutionVector::FIXED_NODE)
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
#   pragma omp parallel for
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
            for (r = ds[ col ].begin(); r!=ds[ col ].end(); r ++){ // loop over rows, dim 1
                Ir[nnz] = *r + ( rc*ncols_s );
                nnz++;
            }// end for rows dim 1
        }// end for dims
    }// end for cols
    Jc[ncols_f] = nnz;

    ds.clear();

    K.MakeSparseMatrix(ncols_f,ncols_f, nnz, Ir, Jc);

 }



SparseMatrix* createSparseMatrix(Geometry& geom, SolutionVector& sol, const int& MatNum){
/*! Creates a SparseMatrix object based on the input geometry and number of dimensions in sol (i.e. is
sol is for Q, matrix size is [5*np X 5*np], and for V [npXnp].
 \param geom = input geometry
 \param sol = solution vector for which the sparse matrix is generated
 \param MatNum = Material number of elements that are to be included in the sparse matrix.
    Default value of MatNum is 0, which results in a sparse matrix for all elements (use this for potential solution).
 \return Returns pointer to created sparse matrix.
*/

    vector <list <unsigned int> > dangly_set;
    create_dangly_matrix( dangly_set, geom, sol ,MatNum );

    SparseMatrix* K = new SparseMatrix();
    convert_sets_to_arrays(dangly_set, sol.getnDimensions(),*K);
    return K;
}



SparseMatrix* createSparseMatrix(vector <Line>& lines){
/*! Crates a SparseMatrix object where non-zeros are determined by the node indexes in the input parameter lines*/

    vector <list <unsigned int> > dangly;
    create_dangly_matrix( lines , dangly); // CREATE DANGLY SET MATRIX OF NODE PAIRS

    SparseMatrix* K = new SparseMatrix();
    convert_sets_to_arrays( dangly, 1, *K );
    return K;
}
