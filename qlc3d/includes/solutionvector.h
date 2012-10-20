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
#include <globals.h>
#include <assert.h>
using std::vector;
using std::list;
using std::cout;
using std::endl;

namespace SolutionVectorNameSpace
{
class node{
public:
    int nodenum;
    int mat;
    node(const int& num, const int& mat):nodenum(num),mat(mat){}
    bool operator<(const node& other)const
    {
        return ( other.nodenum < nodenum );
    }
    bool operator==(const node& other)const
    {
        return (other.nodenum == nodenum);
    }
};
}

class SolutionVector
{
    static const double BIGNUM;

private:
    idx nDoF;		// number number of degrees of freedom per dimension
    idx nFixed;
    idx nDimensions;// number of dofs per node. e.g potential = 1, Q = 5
    idx* FixedNodeMaterial;  // material number of fixed nodes
    bool* IsFixed;
    idx* Elim; 		// effective node ordering after taking into account periodic nodes that need not solving NOT SAME AS PERIODIC EQUIVALENT NODE
    idx* EquNodes;	// nodal equivalencies - i.e. periodic nodes
    idx nFreeNodes; // number of independent degrees of freedom = nDoF - # number of nodes eliminated as periodic

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

    //static const int FIXED_NODE = -1;  // INDEX VALUE OF A FIXED NODE

    // DATA
    bool IsVector;
    idx *FixedNodes;        // INDEX TO EACH FIXED NODE
    double *FixedValues;    // HOLDS NODE VALUE FOR EACH FIXED NODE
    double *Values;
    // END DATA


    ~SolutionVector();
    SolutionVector();
    SolutionVector(idx np);
    SolutionVector(idx np, idx dim);

    SolutionVector& operator=(const SolutionVector&);

    inline idx getnDoF()const {return nDoF;} // npLC
    inline idx getnFreeNodes()const {return nFreeNodes;} // nDof - periodic nodes
    inline idx getnFixed()const{return nFixed;}
    inline idx getnDimensions()const{return nDimensions;}
    inline idx getEquNode(const idx n) const// returns equivalent node to n (for periodic surfaces etc.)
    {
#ifdef DEBUG
assert(n< getnDoF()*getnDimensions());
#endif
        if (Elim)
            return Elim[n];
        else
            return n;
    }

    inline double getValue(const idx n) const // gets the nth value
    {
#ifdef DEBUG
assert( n<getnDoF() );
#endif
        return Values[n];
    }


    inline double getValue(const idx n , const idx dim) const // gets the nth value of dimension dim;
    {
#ifdef DEBUG
        assert(n<getnDoF() );
        assert(dim<getnDimensions());
#endif
        return Values[n + dim*nDoF];
    }

    void setnDoF(idx n);
    void setnFixed(idx n);
    void Allocate(const idx np, const idx ndim = 1);
    void setnDimensions(idx n);
    void setValuesTo(const double& value); // sets all values to value
    void setValuesTo(const double* values); // sets all values to those in array 'values'. length of 'values' must be correct
    void setValuesTo(const SolutionVector& other); // copies values from other SolutionVector
    void setValue(const idx n,const idx dim, const double val);// sets nth value of dimension dim to val
    void setToFixedValues();
    void setFixedNodes(vector<int> *Material,
                       vector<double> *val ,
                       idx *Elem,
                       idx *Mat,
                       idx nElem,
                       idx nNodes);
    void setFixedNodes(Alignment *alignment, int* e);
    void setFixedNodesQ(Alignment* alignment, Mesh* e);
    void setFixedNodesPot(Electrodes* electrodes);
    void setFixedNodesPot(Electrodes& electrodes, Mesh* e, double CurrentTime); // fixed potential for switching

    void allocateFixedNodesArrays(Geometry& geom); // ALLOCATES FIXED NODE ARRAYS FOR POTENTIAL.

    void Resize(const unsigned int& n, const unsigned int& dim = 1); // resizes Values data, clears all data
    void ClearAll();

    void setPeriodicEquNodes(Geometry* geom); // use this for generating nodal periodic equivalency lists

    void ClearFixed(); // clears all fixed nodes and values
    void AddFixed(int mat, double val, Mesh* mesh); // adds fixed when all values are same in region mat . e.g. same potential on electrode
    //void changeFixedValue( int mat, double val,




    void EnforceEquNodes(const Geometry& geom); // enforces periodicity
    void PrintFixedNodes();
    void PrintValues();
    void PrintElim();
    void PrintEquNodes();
    void PrintIsFixed();
    inline bool getIsFixed(const size_t i)
    {
#ifdef DEBUG
assert(i < getnDoF()*getnDimensions() );
#endif
        return IsFixed[i];
    }
    bool test();


};

#endif


