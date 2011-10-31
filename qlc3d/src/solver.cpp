

#include <iostream>
#include <stdlib.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "qlc3d.h"


void sparse_multiply(SparseMatrix *K, SolutionVector *v, double **L)
{

// multiplies sparse matrix by dirichlet nodes and stores result in L
	int ndof = v->getnDoF();
	double *dL = (double*)malloc(ndof * sizeof(double)); //allocate memory for result
	memset(dL,0, ndof * sizeof(double));
	
	double *dirichlet = (double*)malloc(ndof * sizeof(double));
	memset(dirichlet,0,ndof * sizeof(double));
	
	int i,j;

	// construct dirichlet vector
	for (i = 0; i < v->getnFixed(); i++ )
		dirichlet[v->FixedNodes[i]] = v->FixedValues[i];
	for(i=0 ; i < ndof ; i++)
	{
		//printf("L[%i] = %f ",i,(*(L+i)));
		for (j = 0 ; j < ndof ; j++)
		{
		
		// L(i)  += K[i,j] * dirichlet[j];
			
		dL[i] +=  dirichlet[j] * K->sparse_get(i,j);
		//printf("diri=%f\n",dirichlet[j]);
		}//end for j
	}// end for i
	*L = dL;
	
}

void sparse_multiply2(SparseMatrix *K, SolutionVector *v, double **L)
{
	//K->PrintArrays();
	int r, c, row;
	int ndof = v->getnDoF();
	double *dL = (double*)malloc(ndof * sizeof(double)); //allocate memory for result
	memset(dL,0, ndof * sizeof(double));
	
	for (c = 0 ; c < K->rows ; c++)
		{
		//printf("col %i\n",c);
		for (r = K->J[c]; r < K->J[c+1] ; r++)
			{
			row = K->I[r];
			//printf("	row %i = %f\n",row, K->P[ row ]);
			dL[ row ] += K->P[ r ] * v->Values[c];
			}
		}
		*L = dL;
}

void solver(SolutionVector *v,SparseMatrix *K, double* L)
{	// does what the name says 
 
 // int i, j;
  //int nFree = v->nDoF - v->nFixed;
  //v->setValuesTo(0);
  //v->setToFixedValues();
   
// take care of fixed nodes
//1 . Multiply by dirichlet nodes
/*
	printf("RHS vector - ");
	double *L;
	L = (double *)malloc(v->getnDoF()*sizeof(double)); 
	sparse_multiply2(K, v, &L);
 
// 2. set fixed nodes and columns to 0 - except diagonal, which is set to 1	
//	 this could be optimised for speed by setting whole columns / rows, instead of searching individual
//	matrix values. 
	printf(" %i Fixed potential nodes\n", v->getnFixed());
	int row, column;
	
	for (i = 0; i< v->getnFixed() ; i++)
	{
		//printf("%i ",i);
		row = v->FixedNodes[i];
		for (j = 0; j < v->getnDoF() ; j++)
		{
			column = j ;
			if (row==column) // if diagonal
				K->sparse_reset(row,column,1.0);
			else
			{
				K->sparse_reset(row,column,0); 
				K->sparse_reset(column,row,0);
			}
		}
	}
	//K->PrintMatrix();
	  	
// Solve KL=v	(wth a minus sign hack - sort this out to make consistent) 
*/	
	//K->PrintMatrix();
	
	//for (int i = 0; i < K->rows ; i++)
	//	printf("L[%i] = %f\n",i,L[i]);
        //solve_pcg(K,L,v->Values);
	
        //v->setToFixedValues(); // set fixed nodes back to fixed values
	//v->PrintValues();
}
