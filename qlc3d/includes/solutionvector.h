#ifndef SOLUTIONVECTOR_H
#define SOLUTIONVECTOR_H
//#include <stdio.h>
//#include <vector>
#include <omp.h>
#include <algorithm>
#include <alignment.h>
#include <electrodes.h>
#include <simu.h>
#include <mesh.h>
#include <geometry.h>
//#include "material_numbers.h"
using std::vector;
using std::list;
class SolutionVector
{


private:
	int nDoF;		// number number of degrees of freedom per dimension
	int nFixed;
	int nDimensions;// number of dofs per node. e.g potential = 1, Q = 5
	bool* IsFixed;
	int* Elim; 		// effective node ordering after taking into account periodic nodes that need not solving
	int* EquNodes;	// nodal equivalencies - i.e. periodic nodes
	
    int nFreeNodes; // number of independent degrees of freedom = nDoF - # number of nodes eliminated as periodic
	void setBooleanFixedNodeList(); // creates list of booleans (bool* IsFixed) for each node true=fixed node, false = free node
	
	
	void setCornerElim(	list <int>& corn0, // sets periodic equivalent nodes for 4 corners
					list <int>& corn1, 	// corn1[i] = corn0[i]
					list <int>& corn2,   // corn2[i] = corn0[i]
					list <int>& corn3,   // corn3[i] = corn0[i]
					int* Elim,
					const int& dim, // dimension along which corner extends, 0,1,2 -> x,y,z
					double* coords // node coordinates
					); 
	void setFaceElim( list <int>& face0, // face1[i] = face0[i]
				list <int>& face1,
				int* Elim,
				const int& norm, // face normal, 0,1,2 -> x,y,z
				double* coords); // pointer to node coordinates
					
public:
// DATA
	bool IsVector;
	int *FixedNodes;
	double *FixedValues;
	double *Values;
// END DATA


	~SolutionVector();
	SolutionVector();
	SolutionVector(int np);
	SolutionVector(int np, int dim);
	
	SolutionVector& operator=(const SolutionVector&);
	

	
	inline int getnDoF()const {return nDoF;} // npLC
	inline int getnFreeNodes()const {return nFreeNodes;} // nDof - periodic nodes
	inline int getnFixed()const{return nFixed;}
	inline int getnDimensions()const{return nDimensions;}
	inline int getEquNode(const int &n) const// returns equivalent node to n (for periodic surfaces etc.)
	{
		#ifdef DEBUG
			if (n>=nDoF*nDimensions)
			{printf("error - SolutionVector::getEquNode(int n) - j = %i is too big!! - bye\n",n);exit(1);}
		#endif
		if (Elim == NULL)
			return n;
		else
			return Elim[n];
	}
	inline double getValue(const int &n) const // gets the nth value
	{
		#ifdef DEBUG
			if ((n < 0 )  && (n >= getnDoF()) )
				{	printf("error - SolutionVector::getValue(int n) - when trying to access n = %i, bye!",n);	exit(1);}
		#endif
		return Values[n];
	}
	
	
	inline double getValue(const int &n , const int &dim) // gets the nth value of dimension dim;
	{
		#ifdef DEBUG
			if ((n < 0 )  && (n >= getnDoF()) && (dim < 0 ) && ( dim >= getnDimensions() ) )
			{printf("error - SolutionVector::getValue(int n, int dim) - when trying to access n = %i, dim = %i, bye!",n,dim);exit(1);}
		#endif
		return Values[n + dim*nDoF];
	}
	
	void setnDoF(int n);
	void setnFixed(int n);
	void Allocate(const unsigned int& np, const unsigned int& ndim = 1);
	void setnDimensions(int n);
	void setValuesTo(const double& value); // sets all values to value
	void setValuesTo(const double* values); // sets all values to those in array 'values'. length of 'values' must be correct
	void setValuesTo(const SolutionVector& other); // copies values from other SolutionVector
	void setValue(const unsigned int& n,const unsigned int& dim, const double& val);// sets nth value of dimension dim to val
	void setToFixedValues();
	void setFixedNodes(vector<int> *Material, vector<double> *val ,int *Elem,int *Mat,int nElem,int nNodes);
	void setFixedNodes(Alignment *alignment, int* e);
	void setFixedNodesQ(Alignment* alignment, Mesh* e);
	void setFixedNodesPot(Electrodes* electrodes, Mesh* e);
	void setFixedNodesPot(Electrodes* electrodes, Mesh* e, double CurrentTime); // fixed potential for switching
        void Resize(const unsigned int& n, const unsigned int& dim = 1); // resizes Values data, clears all data
        void ClearAll();

	void setPeriodicEquNodes(Geometry* geom); // use this for generating nodal periodic equivalency lists
	
	void ClearFixed(); // clears all fixed nodes and values
	void AddFixed(int mat, double val, Mesh* mesh); // adds fixed when all values are same in region mat . e.g. same potential on electrode
	void EnforceEquNodes(); // enforces periodicity

	void PrintFixedNodes();
	void PrintValues();
	void PrintElim();
	void PrintEquNodes();
	inline bool getIsFixed(const int &i)
	{
		return IsFixed[i];
	}
	bool test();


};

#endif


