#include<math.h>


#include <sparsematrix.h>
#define MB 1048576
SparseMatrix::SparseMatrix()
{
    nnz= 0;
    rows =0;
    cols =0;
    P = NULL;
    I = NULL;
    J = NULL;

    //SMLL = NULL;
}
void SparseMatrix::MakeSparseMatrix(int r, int c, int nonzero, int * i, int *j){
    /*! Sets SparseMatrix row and column indexes to input pointers (takes ownership of memory, with
  deallocation on destruction). Creates new array P for values. */

    nnz = nonzero;
    rows = r;
    cols = c;

    size_t Psize = nnz*sizeof(double);
    P = (double*)malloc(Psize);
    if (P == NULL){
        printf("\nerror - SparseMatrix::MakeSparseMAtrix\ncould not allocate %iMB for P - bye \n", (int) Psize / MB );
        exit(0);

    }

    memset(P,0,nnz*sizeof(double));

    I = i;
    J = j;

}
//void SparseMatrix::MakeSparseMatrix(int r, int c, int nonzero, int * i, int *j,double *p){
//
//        MakeSparseMatrix(r,c,nonzero,i,j);
//	// set precalculated vector p of values
//	int x;
//	#pragma omp parallel for
//	for (x = 0 ; x<nnz ; x++)
//		P[x] = p[x];
//}

SparseMatrix::~SparseMatrix(){//destructor
    if(P!=NULL)
        free(P);
    if(I!=NULL)
        free(I);
    if(J!=NULL)
        free(J);
}
double SparseMatrix::sparse_get(int row, int column)
{

    double krc=-100.0;
    int c1 = J[column];
    int c2 = J[column+1];
    int index  =-1;

    //mexPrintf("c1= %i, c2=%i, row= %i column = %i\n",c1,c2,row,column);
    for(int i=c1;i<c2;i++)
    {
        //mexPrintf("Ir[%i]= %i\n",i,Ir[i]);
        if(I[i]==row)
        {
            index = i;
            break;
        }
    }//end for

    if(index>=0)
    {
        krc = P[index];
        //printf("sparse matrix = %f\n",krc) ;
    }
    else
    {
        krc = 0;
        //printf("zero by default[%i,%i]\n",row,column);
    }//end if

    //mexPrintf("krc %f",krc*1e20);

    return krc;

}//end double sparse_value 
void SparseMatrix::sparse_add(const int &row, const int &column, const double &val) 
{

    int index = sparse_find_index(row,column);
    if (index>=0)
    {

#ifndef DEBUG
#pragma omp atomic
#endif
        P[index]+=val;

        //printf("%f\n",val);
    }
    else
    {
	printf("SparseMatrix::sparse_add - could not find [row,column]=[%i,%i]\n",row,column);
	exit(1);
    }
}//end double sparse_add
void SparseMatrix::sparse_add(int row, int column, double val, int Threadnum)
{
    int index = sparse_lfind_index(row,column);
    if (index>=0)
    {

        P[index+nnz*Threadnum]+=val;
    }
    else
    {
	printf("SparseMatrix::sparse_add - could not find [row,column]=[%i,%i]\n",row,column);
	exit(1);
    }



}

int SparseMatrix::sparse_lfind_index(const int &row, const int &column)
{ // linear search for matrix index - use binary search instead

    int c1 = J[column];
    int c2 = J[column+1];
    int index  =-1;
#ifdef DEBUG
    if ((rows==0) || (cols == 0) || (nnz==0) || (P==NULL) || (I==NULL) || (J==NULL))
    {
        printf("SparseMatrix::sparse_find_index - trying to access empty matrix\n");
        exit(1);
    }
    else if ((row>rows-1)||(column>cols-1))
    {
        printf("SparseMatrix::sparse_find_index - trying to access [%i,%i], size of matrix = [%i,%i]\n",row,column,rows,cols);
    }
#endif
    for(int i=c1;i<c2;i++)
    {

        if(I[i]==row)
        {
            index = i;
            break;
        }
    }//end for
    return index;







}
int SparseMatrix::sparse_find_index(const int &row, const int &column)
{

    if ((row>=rows) || (column >= cols))
    {
        printf("error - SparseMatrix::sparse_find_index - trying to acces row col [%i,%i], when matrix size is [%i x %i], bye!\n", row, column, rows,cols);
        exit(1);
	
    }

    int k1, k2,k;
    k = 0;
    k1 = J[column];
    k2 = J[column+1];
    while (k1<=k2)
    {
        k = (k1+k2)>>1;
	
        if (I[k]==row)
            break;
        else if(I[k] < row)
            k1 = k+1;
        else
            k2 = k-1;

	
    }
    if (I[k] != row)
    {
        k=-1; // entry not found return negative index

    }

    return k;




}

void SparseMatrix::sparse_set(int const &row, int const &column,double val)
{// sets -i.e. overwrites matrix value at row, column with val
    int index = sparse_find_index(row,column);
    if (index>=0)
    {

        P[index]=val;

    }
    else
    {
	printf("SparseMatrix::sparse_set - could not find [row,column]=[%i,%i]\n",row,column);
	exit(1);
    }
} // end void sparse set

void SparseMatrix::sparse_reset(int row,int column,double val)
{// sets -i.e. overwrites matrix value at row, column with val - if matrix entry found
    int index = sparse_find_index(row,column);
    if (index>=0)
    {
        //printf("setting [%i,%i]\n",row,column);
        P[index]=val;
    }
    else
    {
	
	//printf("[%i,%i]\n",row,column);
    }
} // end void sparse set


void SparseMatrix::MultiplyArray(double *in, double *result)
{//multiplies Matrix with array, stores result in result. array and result are assumed to be correct sizes...


    //K->PrintArrays();
    int r, c, row;
    //int ndof = v->nDoF;

    for (c = 0 ; c < rows ; c++)
    {
        //printf("col %i\n",c);
        for (r = J[c]; r < J[c+1] ; r++)
        {
            row = I[r];
            //printf("	row %i = %f\n",row, K->P[ row ]);
            result[ row ] += P[ r ] * in[c];
        }
    }

}
void SparseMatrix::ReleasePointers()
{// releases pointers to P, I and J - needed for SparseLib++ 
    if (P!=NULL)
        free(P);
    if (I!=NULL)
        free(I);
    if (J!=NULL)
        free(J);
    P = NULL;
    I = NULL;
    J = NULL;
    nnz  = 0;
    rows = 0;
    cols = 0;
}
void SparseMatrix::setAllValuesTo(double val)
{
    //printf("nnz = %i\n",nnz);
    int i;
#ifndef DEBUG
#pragma omp parallel for
#endif
    for (i = 0 ; i < nnz ; i++)
        P[i] = val;

}

void SparseMatrix::setFixedIndexes(SolutionVector* sv)
{ //sets values of indexes to fixed nodes 
    indFixedDiagonals.clear();
    indFixedOffDiagonals.clear();

    int row, column, tempind;


    for (idx i = 0; i < sv->getnFixed() * sv->getnDimensions() ; i++)  // loop rows
    {

        row = sv->FixedNodes[i];
        for (int j = 0; j < cols ; j++) // loop columns
        {
            column = j ;
            if (row==column) // if diagonal
            {
                tempind = sparse_find_index(row,column);
                if (tempind>=0)	indFixedDiagonals.push_back(tempind);
            }
            else
            {

                tempind = sparse_find_index(row,column);
                if (tempind>=0)	indFixedOffDiagonals.push_back(tempind);
                tempind = sparse_find_index(column,row);
                if (tempind>=0)	indFixedOffDiagonals.push_back(tempind);
            }
        }
    }


    //printf("nD = %i, nOD = %i\n", indFixedDiagonals.size(), indFixedOffDiagonals.size())	;
}// end setFixed Indexes
void SparseMatrix::setFixed()
{
    vector <int>::iterator itr;
    for (itr = indFixedDiagonals.begin(); itr != indFixedDiagonals.end() ; itr++)
        P[*itr] = 1;

    for (itr = indFixedOffDiagonals.begin(); itr != indFixedOffDiagonals.end() ; itr++)
        P[*itr] = 0;

}
void SparseMatrix::ExpandForOMP(int nThreads)
{// reallocates P to an expanded size of nThreads*sizeof(P)

    if (P == NULL)
    {
        printf("SparseMatrix::ExpandForMultiThreading(int) - error, Matrix has not been initialissed yet- bye!\n");
        exit(1);
    }
    else
    {
        if (nThreads > 1)
        {
            free(P);
            P = (double*) malloc(nnz*nThreads*sizeof(double));
            memset(P,0,nnz*nThreads*sizeof(double));
        }
    }

}
void SparseMatrix::SumForOMP(int nThreads)
{
    if (nThreads >1)
    {
        int i, th;
#ifndef DEBUG
#pragma omp for private(i,th)
#endif
        for (i = 0 ; i < nnz ; i++)
        {
            for ( th = 1 ; th < nThreads ; th ++)
            {
                P[i] += P[i+th*nnz]; 	// add slave data to master
                P[i+th*nnz] = 0;		// reset to zero
            }
        }
    }

}


void SparseMatrix::PrintInfo()
{
    printf("\n [rows,columns,nnz] = [%i,%i,%i]\n",rows,cols,nnz);

}

void SparseMatrix::PrintDiagonal()
{
    PrintInfo();
    for (int i = 0 ; i < rows ; i++)
    {
        printf("diagonal[%i] = %f\n",i,sparse_get(i,i) );
	
    }

}


void SparseMatrix::PrintMatrix()
{
    // debugging...
    printf("\nSparse Matrix is:\n");
    printf("[rows,cols,nnz] = [%i,%i,%i]\n",rows,cols,nnz);

    int i,j;

    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {
            printf("\t%0.0f ",sparse_get(i,j));
        }
        printf("\n\n");
    }

}

void SparseMatrix::PrintArrays()
{
    int i;
    printf("\n\nP=\t");
    for (i=0;i<nnz;i++)
	printf("%1.5f\t",P[i]); // print values
    printf("\nI=\t");
    for (i=0;i<nnz;i++)
	printf("%i\t",I[i]);	// print row indexes
    printf("\nJ=\t");
    for (i=0;i<cols+1;i++)
	printf("%i\t",J[i]);	// print column indexes

}

void SparseMatrix::SPY()
{
    printf("\nSparse Matrix is:\n");
    printf("[rows,cols,nnz] = [%i,%i,%i]\n",rows,cols,nnz);


    int i,j,exists;
    printf("\t");

    for (i=0;i<cols;i++)	// print colum numbers
        printf("%i",i);

    printf("\n");
    for (i=0;i<rows;i++)
    {
        printf("%i\t",i);		// printf row numbers
        for (j=0;j<cols;j++)
        {

            exists = sparse_lfind_index(i,j);
            //if ((i==1)&&(j==0))
            //{
            //	printf("\n[1,0] = %i\n",sparse_find_index(i,j));
            //}

            if (exists >= 0)
                printf("#");//,exists);
            else
                printf(".");
        }
        printf("\n");
    }

}//end void SPY

void SparseMatrix::PrintMatlab(double *L, const char* fname=NULL)
{// generates MATLAB m-file of sparse matrix
    FILE *fid = NULL;
    if (fname == NULL)
        fid = fopen("sparsematrix.m","wt");
    else
        fid = fopen(fname,"wt");
    if (fid)
    {
        //SPY();
        int i, j;
        //fprintf(fid,"K = sparse(");
        double *P;
        int *I, *J;
        P = (double*)malloc(rows*cols*sizeof(double));
        I = (int*)malloc(rows*cols*sizeof(double));
        J = (int*)malloc(cols*rows*sizeof(int));
        if ((P== NULL) || (I == NULL) || (J == NULL))
        {
            printf("error - SparseMatrix::PrintMatlab - Could not allocate memory");
            exit(1);
        }



        for (i = 0; i < rows; i++)
        {
            for (j = 0; j < cols ; j++)
            {
                P[i*cols + j] = sparse_get(i,j);
                I[i*cols + j] = i+1;
                J[i*cols + j] = j+1;
            }
        } // end for i
	
        fprintf(fid,"i = [");
        for (i = 0; i < rows*cols ; i++)
            fprintf(fid, "%i ",I[i]);
        fprintf(fid,"];\n");
	
        fprintf(fid,"j = [");
        for (i = 0; i < rows*cols ; i++)
            fprintf(fid, "%i ",J[i]);
        fprintf(fid,"];\n");

        fprintf(fid,"p = [");
        for (i = 0; i < rows*cols ; i++)
            fprintf(fid, "%e ",P[i]);
        fprintf(fid,"];\n");
	
        fprintf(fid,"sK = sparse(i,j,p);\n");
	
        fprintf(fid,"L = [");
        for( i = 0 ; i < rows ; i++)
            fprintf(fid, "%e ",L[i]);
        fprintf(fid,"];\n");

        fclose(fid);
        free(P);
        free(I);
        free(J);
    }
    else
    {
        printf("error - could not open file sparse.m for output\n");
    }




}

void SparseMatrix::PrintMatlab()
{// generates MATLAB m-file of sparse matrix
    FILE *fid = fopen("sparsematrix.m","wt");
    if (fid)
    {
        //SPY();
        int i, j;
        //fprintf(fid,"K = sparse(");
        double *P;
        int *I, *J;
        P = (double*)malloc(rows*cols*sizeof(double));
        I = (int*)malloc(rows*cols*sizeof(double));
        J = (int*)malloc(cols*rows*sizeof(int));
        if ((P== NULL) || (I == NULL) || (J == NULL))
        {
            printf("error - SparseMatrix::PrintMatlab - Could not allocate memory");
            exit(1);
        }



        for (i = 0; i < rows; i++)
        {
            for (j = 0; j < cols ; j++)
            {
                P[i*cols + j] = sparse_get(i,j);
                I[i*cols + j] = i+1;
                J[i*cols + j] = j+1;
            }
        } // end for i
	
        fprintf(fid,"i = [");
        for (i = 0; i < rows*cols ; i++)
            fprintf(fid, "%i ",I[i]);
        fprintf(fid,"];\n");
	
        fprintf(fid,"j = [");
        for (i = 0; i < rows*cols ; i++)
            fprintf(fid, "%i ",J[i]);
        fprintf(fid,"];\n");

        fprintf(fid,"p = [");
        for (i = 0; i < rows*cols ; i++)
            fprintf(fid, "%e ",P[i]);
        fprintf(fid,"];\n");
	
        fprintf(fid,"sK = sparse(i,j,p);\n");
	
        fclose(fid);
        free(P);
        free(I);
        free(J);
    }
    else
    {
        printf("error - could not open file sparse.m for output\n");
    }




}

void SparseMatrix::DetectZeroDiagonals()
{
    for (int i = 0 ; i < rows ; i ++)
    {
        if (sparse_get(i,i) == 0.0 )
        {
            printf("Zero diagonal in sparse matrix at [%i,%i]\n", i ,i);
            exit(1);
        }
	
    }

}
