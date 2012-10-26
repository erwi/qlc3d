//#include <sparsematrix.h>
#include <settings.h>
//#ifndef COMPLEX
//    #define COMPLEX std::complex<double>
//#endif

//#include <compcol_double.h>
//#include <cg.h>			// Solvers - conjugate gradient
//#include <gmres.h>		// Solvers - GMRES

//#include <icpre_double.h>
//#include <diagpre_double.h>
//#include <ilupre_double.h>
//#include MATRIX_H
#include <stdlib.h>
#include <iostream>

#include <ircmatrix.h>
#include <vector.h>
#include <iterativesolvers.h>
#include <luincpreconditioner.h>
#include <diagpreconditioner.h>



using namespace std;

// solves Kb=x, using preconditioned conjugate gradient method
void solve_pcg(IRCMatrix &K, Vector &b, Vector &x, Settings* settings )
{

    idx size = K.getNumRows();

	//for (int i  = 0 ; i <size ; i ++)
	//b[i] = 0;
	//printf("%f\n",b[i]);

    //fflush(stdout);
    // Create SparseLib++ data structures. Changes have been made to the original
    // SparseLib code in order to avoid copying data between K and SparseLib Matrix A.
    // This should allow for larger meshes in Win32 without running out of memory.
    // Also, this is a bit faster.

    // Create empty sparse matrix A and set array pointers from K.
    //CompCol_Mat_double A;
    //A.point_to(size, nnz, K->P, K->I , K->J );

    // Same with vectors - Create empty ones and then set pointers to data
    //VECTOR_double X;
    //VECTOR_double B;
    //B.point_to( b, A.dim(0) );
    //X.point_to( x, A.dim(0) );

    // PCG settings...
    /*
    int return_flag =10;
		int maxiter = settings->getQ_PCG_Maxiter();
		double toler =settings->getQ_PCG_Toler();
		//printf("Q-PCG maxiter = %i, toler = %f\n",maxiter,toler);
	// Solves with different preconditioners...
        if (settings->getQ_PCG_Preconditioner() == DIAG_PRECONDITIONER ){
			DiagPreconditioner_double D(A); // diagonal preconditioning, ~+3 times faster than cholesky
           return_flag = CG(A,X,B,D,maxiter,toler);
		}
        else if (settings->getQ_PCG_Preconditioner() == IC_PRECONDITIONER ){
			ICPreconditioner_double D(A);
			return_flag = CG(A,X,B,D,maxiter,toler);
		}
        else if (settings->getQ_PCG_Preconditioner() == ILU_PRECONDITIONER ){
			CompCol_ILUPreconditioner_double D(A); // compressed column format ILU
			return_flag = CG(A,X,B,D,maxiter,toler);
		}

	if (return_flag == 1) // if no convergence, print warning msg.
			printf("PCG did not converge in %i iterations \nTolerance achieved is %f\n",maxiter,toler);
*/
    //copy solution back to solution vector
#ifndef DEBUG
#pragma omp parallel for
#endif
    for (idx i = 0; i < size ; i++)
        x[i]*= -1;//-X(i);

}

void solve_gmres(IRCMatrix &K, Vector &b, Vector &x, Settings* settings ){

    idx nnz = K.getnnz();
    idx size = K.getNumRows();
	// Create SparseLib++ data structures
        //CompCol_Mat_double A;
        //A.point_to( size, nnz, K->P, K->I, K->J);
        //A = CompCol_Mat_double(size,size,nnz,K->P,K->I,K->J);
	//convert solution vector and RHS vector to SparseLib++
        //VECTOR_double X;// = VECTOR_double(x,A.dim(0));
        //VECTOR_double B;// = VECTOR_double(b,A.dim(0));
        //X.point_to(x, A.dim(0));
        //B.point_to(b, A.dim(0));
    //Vector X(x, size);
    //Vector B(b, size);

	// GMRES settings...
        idx return_flag =10;
        idx maxiter 	= settings->getQ_GMRES_Maxiter();
        idx restart 	= settings->getQ_GMRES_Restart();
        real toler      = settings->getQ_GMRES_Toler();

        Preconditioner *M;
        //LUIncPreconditioner LU(K);
        M = new DiagPreconditioner(K);

        //M = new LUIncPreconditioner(K);


        omp_set_num_threads(settings->getnThreads());
        if (!IterativeSolvers::gmres(K,x,b,*M,maxiter,restart, toler))
            printf("GMRES did not converge in %i iterations \nTolerance achieved is %f\n",maxiter,toler);

        delete M;
        //MATRIX_double H(restart+1, restart, 0.0);	// storage for upper Hessenberg H

		// Solves with different preconditioners...
        /*
        if (settings->getQ_GMRES_Preconditioner() == DIAG_PRECONDITIONER )
		{
			DiagPreconditioner_double D(A); // diagonal preconditioning, ~+3 times faster than cholesky
			return_flag = GMRES(A,X,B,D,H,restart,maxiter,toler);
		}
		else if (settings->getQ_GMRES_Preconditioner() == IC_PRECONDITIONER )
		{
			ICPreconditioner_double D(A);
			return_flag = GMRES(A,X,B,D,H,restart,maxiter,toler);
		}
		else if (settings->getQ_GMRES_Preconditioner() == ILU_PRECONDITIONER )
		{
			CompCol_ILUPreconditioner_double D(A); // compressed column format ILU
			return_flag = GMRES(A,X,B,D,H,restart,maxiter,toler);
		}

		if (return_flag == 1)
			printf("GMRES did not converge in %i iterations \nTolerance achieved is %f\n",maxiter,toler);
    */

	//copy solution back to solution vector
//#ifndef DEBUG
//    #pragma omp parallel for
//#endif
    //for (int i = 0; i < size ; i++)
    //    x[i] *= -1.0;//X(i);

}
