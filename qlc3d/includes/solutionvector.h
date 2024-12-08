#ifndef SOLUTIONVECTOR_H
#define SOLUTIONVECTOR_H
#include <omp.h>
#include <algorithm>
#include <memory>
#include <alignment.h>
#include <electrodes.h>
#include <simu.h>
#include <mesh.h>
#include <geometry.h>
#include <globals.h>
#include <cassert>
#include <dofmap.h>

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

private:
  /** Number of degrees of freedom per dimension, including fixed and periodic nodes. */
  idx nDoF;		// number of degrees of freedom per dimension
  idx nDimensions; // number of dofs per node. e.g potential = 1, Q = 5
  unsigned int numFixedNodes;
  /** Number of independent degrees of freedom per dimension. Fixed nodes an periodic nodes are not counted. */
    //idx nFreeNodes; // number of independent degrees of freedom = nDoF - # number of nodes eliminated as periodic
    std::vector<double> values;
    std::unique_ptr<DofMap> dofMap;

public:
    SolutionVector();

    /**
     * Solution vector with total of values np * dim values.
     * @param np number of points
     * @param dim number of values/dimensions per point. E.g. 1 for scalar potential, or 5 for Q-tensor.
     */
    SolutionVector(idx np, idx dim);

    SolutionVector& operator=(const SolutionVector&);

    /** Number of degrees of freedom per dimension, including fixed and periodic nodes. a.k.a npLC */
    [[nodiscard]] inline idx getnDoF()const {return nDoF;} // npLC
    /** Number of independent degrees of freedom per dimension. Fixed nodes an periodic nodes are not counted. */
    [[nodiscard]] inline idx getnFreeNodes()const {return dofMap->getnFreeNodes(); }//nFreeNodes;} // nDof - periodic nodes
    /** Number of fixed nodes per dimension. */
    [[nodiscard]] inline idx getnFixed() const { return numFixedNodes; }
    [[nodiscard]] inline idx getnDimensions() const { return nDimensions; }
    /** returns equivalent node to n (for periodic surfaces etc.) or max value (NOT_AN_INDEX) if node is fixed */
    [[nodiscard]] inline idx getEquNode(const idx n) const { return dofMap->getDof(n); }
    [[nodiscard]] const DofMap &getDofMap() const;
    void loadEquNodes(const idx *start, const idx *end, idx *equNodesOut) const;

    /** raw array access to values, ignores dimensions */
    [[nodiscard]] inline double getValue(const idx n) const { return values[n]; }
    /** Get the n'th value of the i'th dimension. Both n and i are 0 based, so 0 is first value */
    [[nodiscard]] inline double getValue(const idx n , const idx dim) const { return values[n + dim * nDoF]; }

    void setValuesTo(const double& value); // sets all values to value
    void setValuesTo(const SolutionVector& other); // copies values from other SolutionVector
    void setValue(const idx n,const idx dim, const double val);// sets nth value of dimension dim to val

    void initialiseLcBoundaries(Geometry &geom, const Alignment &alignment);
    void initialisePotentialBoundaries(const std::unordered_map<unsigned int, double> &potentialByElectrode, Geometry &geom);
    /** Set the values for potential at electrode nodes. */
    void setFixedPotentialValues(const Mesh &triangles, const std::unordered_map<unsigned int, double> &potentialByElectrode);

    /** DEPRECATED, this can probably be removed */
    void Resize(const unsigned int& n, const unsigned int& dim = 1); // resizes Values data, clears all data
    void ClearAll();

    // Q-tensor related only
    //! set all of the 5 tensor values for te n'th DoF.
    void setValue(idx n, const qlc3d::TTensor &t);

    void setValue(idx i, const qlc3d::Director &d);

    void loadValues(const idx *start, const idx *end, double *valuesOut) const;
    void loadQtensorValues(const idx *start, const idx *end, qlc3d::TTensor* tensorOut) const;

    [[nodiscard]] qlc3d::Director getDirector(idx i) const;
    [[nodiscard]] std::vector<qlc3d::Director> getDirector() const;
    [[nodiscard]] double& operator[](const unsigned int i) { return values[i]; }
    [[nodiscard]] int64_t hashCode();

    /** Increment the values of the free degrees of freedom by the given container with '[]' indexing operator. */
    template <typename T>
    void incrementFreeDofs(const T& array, const double factor = 1.0) {
      const idx numDofs = getnDoF();
      const idx numDims = getnDimensions();
      for (idx dim = 0; dim < numDims; dim++) {
        for (idx i = 0; i < nDoF; i++) {
          const idx index = i + dim * numDofs;
          const idx equIndex = getEquNode(index);

          if (equIndex != NOT_AN_INDEX) {
            values[index] += array[equIndex] * factor;
          }
        }
      }
    }

  /** Copy the values of the free degrees of freedom to the given container with '[]' indexing operator. */
  template <typename T>
  void copyFreeDofsTo(T& array) const {
    idx nDofs = getnDoF() * getnDimensions();
    for (idx i = 0; i < nDofs; i++) {
      const idx indDof = getEquNode(i);
      if (indDof == NOT_AN_INDEX) {
        continue;
      }
      array[indDof] = getValue(i);
    }
  }

};

#endif


