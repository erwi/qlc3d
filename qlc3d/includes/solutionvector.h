#ifndef SOLUTIONVECTOR_H
#define SOLUTIONVECTOR_H
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

namespace qlc3d {
    class TTensor;
    class Director;
}

class SolutionVector {
    static const double BIGNUM;

private:
    idx nDoF;		// number of degrees of freedom per dimension
    idx nFixed;
    idx nDimensions; // number of dofs per node. e.g potential = 1, Q = 5
    std::vector<idx> fixedNodeMaterial;
    std::vector<idx> elim; // effective node ordering after taking into account periodic nodes that need not solving NOT SAME AS PERIODIC EQUIVALENT NODE
    std::vector<bool> isFixed;
    std::vector<idx> equNodes; // nodal equivalencies - i.e. periodic nodes TODO: looks like this might never be used!?
    idx nFreeNodes; // number of independent degrees of freedom = nDoF - # number of nodes eliminated as periodic
    std::vector<double> values;

    std::vector<idx> fixedNodes; // index to each fixed node
    std::vector<double> fixedValues; // value of each fixed node
    void setBooleanFixedNodeList(); // creates list of booleans (bool* IsFixed) for each node true=fixed node, false = free node
public:
    ~SolutionVector();
    SolutionVector();

    /**
     * Solution vector with total of values np * dim values.
     * @param np number of points
     * @param dim number of values/dimensions per point. E.g. 1 for scalar potential, or 5 for Q-tensor.
     */
    SolutionVector(idx np, idx dim);

    SolutionVector& operator=(const SolutionVector&);

    inline idx getnDoF()const {return nDoF;} // npLC
    inline idx getnFreeNodes()const {return nFreeNodes;} // nDof - periodic nodes
    inline idx getnFixed()const{return nFixed;}
    inline idx getnDimensions()const{return nDimensions;}
    inline idx getEquNode(const idx n) const {// returns equivalent node to n (for periodic surfaces etc.)
#ifdef DEBUG
assert(n< getnDoF()*getnDimensions());
#endif
        if (!elim.empty()) {
          return elim.at(n);
        } else {
          return n;
        }
    }

    [[nodiscard]] inline double getValue(const idx n) const { return values[n]; }
    [[nodiscard]] inline double getValue(const idx n , const idx dim) const { return values[n + dim*nDoF]; }

    void setnFixed(idx n);
    void setValuesTo(const double& value); // sets all values to value
    void setValuesTo(const SolutionVector& other); // copies values from other SolutionVector
    void setValue(const idx n,const idx dim, const double val);// sets nth value of dimension dim to val
    void setToFixedValues();
    void setFixedNodesQ(Alignment* alignment, Mesh* e);
    void setFixedNodesPot(Electrodes* electrodes);
    void allocateFixedNodesArrays(Geometry& geom); // ALLOCATES FIXED NODE ARRAYS FOR POTENTIAL.
    void Resize(const unsigned int& n, const unsigned int& dim = 1); // resizes Values data, clears all data
    void ClearAll();
    void setPeriodicEquNodes(Geometry* geom); // use this for generating nodal periodic equivalency lists
    void ClearFixed(); // clears all fixed nodes and values
    void EnforceEquNodes(const Geometry& geom); // enforces periodicity
    [[nodiscard]] inline bool getIsFixed(const size_t i) {
#ifdef DEBUG
assert(i < getnDoF()*getnDimensions() );
#endif
        return isFixed.at(i);
    }
    bool test();

    // Q-tensor related only
    //! set all of the 5 tensor values for te n'th DoF.
    void setValue(idx n, const qlc3d::TTensor &t);

    void setValue(idx i, const qlc3d::Director &d);

    [[nodiscard]] qlc3d::Director getDirector(idx i) const;
    [[nodiscard]] std::vector<qlc3d::Director> getDirector() const;
    [[nodiscard]] double& operator[](const unsigned int i) { return values[i]; }
    [[nodiscard]] int64_t hashCode();
};

#endif


