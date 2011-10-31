

#include <iostream>
#include <stdlib.h>
#include <stdlib.h>
#include <vector>
#include <compcol_double.h>
#include <cg.h>
#include <icpre_double.h>
#include <sparsematrix.h>
//#include "meshtest.h"

using std::vector;

/*

void sparse_multiply(SparseMatrix *K, SolutionVector *v, double **L)
{
// multiplies sparse matrix by dirichlet nodes and stores result in L
	int ndof = v->nDoF;
	double *dL = (double*)malloc(ndof * sizeof(double)); //allocate memory for result
	memset(dL,0, ndof * sizeof(double));
	
	double *dirichlet = (double*)malloc(ndof * sizeof(double));
	memset(dirichlet,0,ndof * sizeof(double));
	
	int i,j;

	// construct dirichlet vector
	for (i = 0; i < v->nFixed; i++ )
		dirichlet[v->FixedNodes[i]] = v->FixedValues[i];
	for(i=0 ; i < ndof ; i++)
	{
		//printf("L[%i] = %f ",i,(*(L+i)));
		for (j = 0 ; j < ndof ; j++)
		{
		
		// L(i)  += K[i,j] * dirichlet[j];
			
		dL[i] +=   dirichlet[j] * K->sparse_get(i,j);
		
		}//end for j
        }// end for
	*L = dL;
}
*/

/*
void solver(SolutionVector *v,SparseMatrix *K)
{
  
  printf("solving...");
  
  // take care of fixed nodes
  
  //		1 . Multiply by dirichlet nodes
  double *L;
  sparse_multiply(K, v, &L);
 
  int i,j;
 
 int ii[16] = {0,1,2,3,0,1,2,3,0,1,2,3};
 int jj[5] = {0,4,8,8,12};
 double pp[12] = {11,21,31,41,12,22,32,42,14,24,34,44}; 
 
 SparseMatrix k;
 k.MakeSparseMatrix(4, 4, 12, &ii[0],&jj[0],&pp[0]);


// for (i = 0; i<4 ; i++)
//{//
//	for (j=0 ; j<4 ; j++)
//	{
//		k.sparse_set(i,j,pp[j*4 +i]);
//	}
//}

k.PrintMatrix();
 
 
vector<unsigned int> free; 
 //for (i = 1; i<2; i++)// K->rows ; i++)
 //free.push_back(i);
  free.push_back(0);
  free.push_back(3);
  //free.push_back(3);
printf("\nfree = ")  ;
 for (j = 0; j<free.size() ; j++)
  printf("%i ",*(free.begin()+j));

printf("\n P =");
for (j = 0; j<k.nnz ; j++)
	printf(" %1.0f",k.P[j]);
	printf("\n I =");
for (j = 0; j<k.nnz ; j++)
	printf(" %i",k.I[j]);
printf("\n J =");
for (j = 0; j<k.cols+1 ; j++)
	printf(" %i",k.J[j]);
		
	
 k.Reduce(&free);  
   
printf("\n P =");
for (j = 0; j<k.nnz ; j++)
	printf(" %1.0f",k.P[j]); 
	printf("\n");
for (j = 0; j<k.nnz ; j++)
	printf(" %i",k.I[j]);	
printf("\n J =");
for (j = 0; j<k.cols+1 ; j++)
	printf(" %i",k.J[j]);	
	
k.PrintMatrix();	
  //CompCol_Mat_double A;
  //double *P = (double*)malloc(3*sizeof(double));
  //int *I	= (int*)malloc(3*sizeof(int));
  //int *J    = (int*)malloc(4*sizeof(int));
  
  //P[0] = 1; P[1] = 1; P[2] = 1;
  //J[0] = 0; J[1] = 1; J[2] = 2; J[3] = 3;
  //I[0] = 0; I[1] = 1; I[2] = 2;
  
  //A = CompCol_Mat_double(v->nDoF,v->nDoF,nnz,P,I,J);
    
  //MV_Vector_double x =MV_Vector_double(3);
  //MV_Vector_double b =MV_Vector_double(3);
  //x(MV_VecIndex(0,2))=0;
  //b(0)=-7; b(1)=89; b(2)=-888;
  //A*x;
//  printf("b= [%f, %f, %f] done\n",b(0),b(1),b(2));
  
 //DiagPreconditioner_double D(A);
	//ICPreconditioner_double D(A);
	//int result =0;
	//int maxiter =100;
	//double toler = 1e-6;
	//result = CG(A,x,b,D,maxiter,toler);
 //printf("result= %i, x = [%f,%f,%f]\n",result,x(0),x(1),x(2));
 //return 0;
 
 printf("\nsolved \n");
}
*/
