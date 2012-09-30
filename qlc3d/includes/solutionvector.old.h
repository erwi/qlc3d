#include "mesh.h"

#include <stdio.h>
#include <vector>
using std::vector;
class SolutionVector
{

public:
int nDoF;
int nFixed;

int *FixedNodes;
double *FixedValues;
double *Values;
/*
void SetFixedNodes(vector<int> *Material, vector<double> *val , Mesh *e)
{
// set arrays FixedNodes to point to node numbers of Dirichlet nodes and FixedValues to contain the correponding Dirichlet values
	printf("setfixed");
	
	if (Material->size()!=val->size())
	{
		printf("\nsolutionvector::SetFixedNodes - number of fixed values must equal number of materials");
		exit(1);
	}
		nFixed = 0; // reset number of fixed nodes
	
	if (FixedNodes != NULL) free(FixedNodes);	// these may have been set earlier, and should be reset
	if (FixedValues!= NULL) free(FixedValues);
	printf("aaaaaa");
	

	vector<int> FxNodes (nDoF,-1); // start with all negative vector, 
	vector<int>::iterator i;
	vector<int>::iterator j;
	vector<int>::iterator t;
	vector<int> ind;
	
	//create index of fixed nodes, where location is node number and stored value material number
	int c=0;
	for (i=Material->begin(); i!=Material->end(); i++)
	{
		ind.clear();
		e->FindNodesMaterial(*i,&ind); // index to all nodes of this material
				
		//for (t = ind.begin() ; t!= ind.end(); t++)
		//	printf("node %i = material %i\n", *t,*i);
				
		for (j=ind.begin(); j!= ind.end(); j++)
			FxNodes[*j] = c; // store index to values vector
		
		c++;	// update index
	}
	
	// number of non-negative entries equals number of fixed nodes
	for (i = FxNodes.begin(); i != FxNodes.end(); i++)
		if (*i >= 0) nFixed++;
	
	//printf("number of fixed nodes %i\n",nFixed);
	// allocate memory for fixed nodesand values
	FixedNodes  = (int*)malloc(nFixed*sizeof(int));		
	FixedValues = (double*)malloc(nFixed*sizeof(double));
	
	// copy nodes and values into arrays
	c=0;
	int location=0;
	//printf("copying nodes and values into arrays");
	for (i = FxNodes.begin(); i != FxNodes.end(); i++)
	{			
		if (*i >= 0)
		{
			FixedNodes[c] = location;			//index to fixed node
			FixedValues[c] = *(val->begin()+*i);		// value of fixed node
			//printf("node %i =  [%i  %1.1f]\n",location, *i,*(val->begin()+*i));
			c++;
		}
		location ++;
	}
	
	// set fixed values
	for (c = 0 ; c < nFixed ; c++ )
			Values[FixedNodes[c]] = FixedValues[c];
	
	
	
}// end void SetFixedNodes
*/
void SetFixedNodes()//vector<int> *Material, vector<double> *val , Mesh *e)
{
	printf("aaaaaaaaaaaaaaaaaaa\n");
}
void SetValuesTo(double value)
{// sets all Values to value
	for (int i = 0 ; i < nDoF ; i ++)
		Values[i] = value;
}

void SetToFixedValues()
{// sets all fixed values to FixedValues
	for (int i = 0 ; i < nFixed; i++)
		Values[FixedNodes[i]] = FixedValues[i];
}

void PrintFixedNodes()
{//debugging
printf("number of fixed nodes is : %i\n",nFixed);	
int i;
	for (i=0; i<nFixed; i++)
		printf("fixed node %i, node number %i, value %1.3f\n",i,FixedNodes[i],FixedValues[i]);

}


void PrintValues()
{// prints all values
printf("\nnDoF = %i\n",nDoF);
int i;	
	for (i=0;i<nDoF;i++)
	 printf("Value[%i] = %f \n",i,Values[i]);
}// end void PrintValues



SolutionVector(int n)
{
	printf("solutionvector.h\n");
	nDoF = n;
	nFixed = 0;
	printf("n = %i\n",n);
	
	FixedNodes = NULL;   // intialise to NULL
    FixedValues = NULL;
	Values = NULL;

	
	
	if ((Values = (double*)malloc(nDoF*sizeof(double))) == NULL)
		{
		printf("SolutionVector::SolutionVector(int n) - could not allocate memory to hold values");
		exit(1);
		}
	memset(Values,'0',nDoF*sizeof(double));
		
}


~SolutionVector()
{
	
	if (FixedNodes != NULL)
		free(FixedNodes);
	if (FixedValues != NULL)
		free(FixedValues);
	if (Values != NULL)
		free(Values);
		
}


};// end class SolutionVector
