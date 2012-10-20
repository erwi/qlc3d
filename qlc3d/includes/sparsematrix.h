#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <vector>
# include <algorithm>
# include <omp.h>
//# include "sparsematrixll.h"
# include <string>
#include "solutionvector.h"
using std::vector;

class SparseMatrix
{
private:
	vector <int> indFixedDiagonals;     // index to fixed diagonal nodes
	vector <int> indFixedOffDiagonals;  // index ro fixed off diagonal nodes
							
	
public:
	
	double *P;
	int *I, *J;
	int nnz, rows, cols;

	SparseMatrix();
	~SparseMatrix();
	
	void MakeSparseMatrix(int r, int c, int nonzero, int * i, int *j);
        //void MakeSparseMatrix(int r, int c, int nonzero, int * i, int *j,double *p);

	double sparse_get(int row, int column);
	void sparse_set(int const &row,int const &column,double val);
	void sparse_reset(int row,int column,double val);
	
	void sparse_add(const int &row, const int &column, const double &val) ;
	void sparse_add(int row, int column, double val, int Thread);
	
	int sparse_lfind_index(const int  &row, const int &column);
	int sparse_find_index(const int  &row, const int &column);
	
	void MultiplyArray(double *in, double *result);
	void ReleasePointers();
	void setAllValuesTo(double val); // set all values to val, e.g. for clearing matrix
	void setAllValuesTo(double val,int nThreads); // set all values, taking into accound multiple threads;
	void setFixedIndexes(SolutionVector* sv); //sets values of indexes to fixed nodes 
	void setFixed(); // sets fixed diagonals to 1 and off diagonals to 0
	
	void ExpandForOMP(int nThreads); // expands P array for slave threads
	void SumForOMP(int nThreads); // adds slave thread data to master thread data 
	
	
	void PrintInfo();
	void PrintMatrix();
	void PrintArrays();
	void SPY();
	void PrintMatlab(double *L, const char* fname);
	void PrintMatlab();
	void PrintDiagonal();
	void DetectZeroDiagonals();
	//void sparse_reset(int row, int column ,double val);
	//double sparse_get(int row, int column);

};

#endif

