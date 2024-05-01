#ifndef MESH_H
#define MESH_H
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <memory>
#include <cstring>
#include <cmath>
#include <globals.h>
using std::vector;
using std::set;
using std::list;

class Coordinates;
class Vec3;

class Mesh {
private:
  const idx Dimension;  // number of dimensions of mesh - 2 for tris and 3 for tets
  const idx nNodes;     // number of nodes per element

  idx nElements;  //total number of elements
  std::vector<idx> nodes;
  std::vector<idx> materials;
  std::vector<double> determinants;
  std::vector<Vec3> surfaceNormals;
  std::vector<idx> connectedVolumes;
  double TotalSize;   // Total volume/area of the mesh

public:
  Mesh(unsigned int dimension, unsigned int nodesPerElement);

  ~Mesh();

  static std::shared_ptr<Mesh> triangleMesh() {
    return std::make_shared<Mesh>(2, 3);
  }

  static std::shared_ptr<Mesh> tetMesh() {
    return std::make_shared<Mesh>(3, 4);
  }

    inline idx getnElements() const {
      return nodes.size() / getnNodes();
    }
    /** number of nodes per element */
    inline idx getnNodes() const {
        return nNodes;
    }
    inline idx getDimension() const {
        return Dimension;   // number of dimesnions of mesh ( 3 / 2 for tets / tris)
    }
    idx getConnectedVolume(const idx e) const;  // returns index to connected volume element, or -1 if not connected to LC1
    inline idx getNode(const idx e, const idx n) const { // returns node n of element e
#ifdef DEBUG
        if ((e >= nElements) || (n >= nNodes)) {
            printf("error - Mesh::getNode(int,int) - index to node out of bounds, bye! ");
            printf("requested elem %u, node %u\n", e, n);
            exit(1);
        }
#endif
    return nodes[e * nNodes + n];
    }

    idx getMaterialNumber(const idx e) const;   // returns material number of element e
    idx getFixLCNumber(const idx e) const;  // gets alignment layer number, i.e. FixLC 1, 2, 3...
    idx getDielectricNumber(const idx e) const; // gets dielectric materials number, i.e. Dielectric 1, 2, 3 ....

    double getDeterminant(const idx i) const; // returns value of determinant of element i

    void setElementData(std::vector<unsigned int> &&nodes, std::vector<unsigned int> &&materials);

    void setConnectedVolume(Mesh *vol);     // sets indexes to connected LC volume elements

    /** copy all node values to the mesh */
    void setAllNodes(idx *nodes);
    void setSurfaceNormal(idx i, const Vec3 &normal);
    void setnElements(idx nnelem);                  // set numbero of elements TODO: should only depend on nodes vectors size

    /** writes all nodes of given element to nodesOut */
    void loadNodes(idx elementIndex, idx *nodesOut) const;
    void removeElements(std::set <idx> &index);         // removes elements in index. index must be sorted in ascending order
    void ClearMesh();                               // clears all data in mesh object
    void appendElements(const vector <idx> &nodeValues, const vector <idx> &materialValues);    //adds new elements and materials, assuming element types match existing elements (nodes/per element)

    // Creates list of all nodes belonging to elements of material mat
    void listNodesOfMaterial(std::vector <idx> &nodes, const idx mat) const;
    /** Deprecated use the other listFixLCSurfaceNodes(FixLC num) instead */
    void listFixLCSurfaces(std::vector <idx> &nodes, const idx FixLCNum) const; // list all nodes of given FixLC surface number (FixLCNum = 1,2,3...)
    /**
     * Return a vector of all surface nodes by given FIXLC number.
     */
    std::set<idx> listFixLCSurfaceNodes(const idx FixLCNum) const;

    bool containsCoordinate(idx elem, const Coordinates& coordinates, const Vec3 p) const; // checks whether point p is within element elem
    void CompleteNodesSet(const idx elem, std::vector<idx> &nodes) const; // completes nodes vector with those from element, if nodes is empty returns all elements
    void calculateDeterminants3D(const Coordinates &coords); // calculates determinants of all elements
    void calculateSurfaceNormals(const Coordinates &coords, Mesh *tets = NULL);
    void CopyMesh(Mesh *rhs);   // makes this a copy of Mesh* rhs - why does operator= overloading not work???
    void ScaleDeterminants(const double &s);  // scales all determinants by s, e.g. to go to microns
    void calcLocCoords(const idx elem, const Coordinates &coordinates, const Vec3 &targetPoint, double localCoordinates[4]) const; // calculates 4 local coordinates of coordinate cord in element elem
    [[nodiscard]] Vec3 elementCentroid(unsigned int i, const Coordinates &coordinates) const;
    [[nodiscard]] Vec3 getSurfaceNormal(unsigned int i) const;
    void gen_p_to_elem(vector<set <idx> > &p_to_elem) const; // generates index from points to mesh elements
    // index number of non-existent neighbours elements equals total number of elements
    // i.e. 1 too large to use as an index
};

#endif



