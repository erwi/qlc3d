#include <iostream>
#include <stdlib.h>
#include <stdlib.h>
#include "meshtest.h"
#include "compcol_double.h"

#include "cg.h"
#include "icpre_double.h"
#include "vecdefs.h"
#include "c:/sparselib_1_6/include/gmres.h"
#include "ilupre_double.h"
//#include "ilupre.h"
//#include "blas1.h"
#include MATRIX_H

using namespace std;



//void solve_pcg(double *P, int *I, int *J, double *b, double *x, int size)
void solve_gmres(SparseMatrix *Kred, double *b, double *x)
{

	int nnz = Kred->nnz; // number of nonzero entries
	int size = Kred->rows;
	
	CompCol_Mat_double A;
	double *P = Kred->P;
	int *I = Kred->I;
	int *J = Kred->J;
	A = CompCol_Mat_double(size,size,nnz,P,I,J);
		
	VECTOR_double X = VECTOR_double(A.dim(1),0);
	VECTOR_double B = VECTOR_double(A.dim(1),0);
	
	for (int i = 0; i < A.dim(1) ; i++)
	{
		X(i) = x[i];
		B(i) = b[i];
	}
	
	
	CompCol_ILUPreconditioner_double D(A);
	
	int restart = 32;
	MATRIX_double H(restart+1,restart,0.0);
	
	int return_flag =10;
	int maxiter =10000;
	
	double toler = 1e-6;
	printf("size(A) = [%i,%i],A = %f\n, size(B)=[%i],B=%f",A.dim(0),A.dim(1),A(0,0),B.size(),B(0));
	return_flag = GMRES(A,X,B,D,H,restart,maxiter,toler);

	if (return_flag == 1)
			printf("GMRES did not converge in %i iterations \nTolerance achieved is %f\n",maxiter,toler);
	else
			printf("GMRES converget at iteration %i, to accuracy %f\n",maxiter,toler);
	
	for (int i = 0; i < A.dim(1) ; i++)
		x[i] = X(i);
	
	
}


