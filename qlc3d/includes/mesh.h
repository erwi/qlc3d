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
#include <functional>
#include <unordered_set>
#include <fmt/format.h>

using std::vector;
using std::set;
using std::list;

class Coordinates;
class Vec3;

enum class ElementType {
  UNKNOWN = 0,
  LINEAR_TRIANGLE = 1,
  LINEAR_TETRAHEDRON = 2,
  QUADRATIC_TRIANGLE = 3,
  QUADRATIC_TETRAHEDRON = 4,
};

[[nodiscard]] inline std::string toString(ElementType elementType) {
  return elementType == ElementType::LINEAR_TRIANGLE ? "LINEAR_TRIANGLE" :
         elementType == ElementType::LINEAR_TETRAHEDRON ? "LINEAR_TETRAHEDRON" :
         elementType == ElementType::QUADRATIC_TRIANGLE ? "QUADRATIC_TRIANGLE" :
         elementType == ElementType::QUADRATIC_TETRAHEDRON ? "QUADRATIC_TETRAHEDRON" :
         "UNKNOWN";
}

template <>
class fmt::formatter<ElementType> {
public:
  constexpr auto parse (format_parse_context& ctx) { return ctx.begin(); }
  template <typename Context>
  constexpr auto format (ElementType const& t, Context& ctx) const {
    return format_to(ctx.out(), "{}", toString(t));
  }
};

[[nodiscard]] inline unsigned int getNodesPerElement(ElementType elementType) {
  switch (elementType) {
    case ElementType::LINEAR_TRIANGLE:
      return 3;
    case ElementType::QUADRATIC_TRIANGLE:
      return 6;
    case ElementType::LINEAR_TETRAHEDRON:
      return 4;
    case ElementType::QUADRATIC_TETRAHEDRON:
      return 10;
    default:
      return 0; // unknown
  }
}



class Mesh {
private:
  const idx Dimension;  // number of dimensions of mesh - 2 for tris and 3 for tets

  ElementType elementType_;

  idx nElements;  //total number of elements
  std::vector<idx> nodes;
  std::vector<idx> materials;
  std::vector<double> determinants;
  std::vector<Vec3> surfaceNormals;
  std::vector<idx> connectedVolumes;
  double TotalSize;   // Total volume/area of the mesh

public:
  Mesh(unsigned int dimension, ElementType elementType);

  ~Mesh();

  static std::shared_ptr<Mesh> triangleMesh() {
    return std::make_shared<Mesh>(2, ElementType::UNKNOWN);
  }

  static std::shared_ptr<Mesh> tetMesh() {
    return std::make_shared<Mesh>(3, ElementType::UNKNOWN);
  }

    inline idx getnElements() const {
      return nodes.size() / getnNodes();
    }
    /** Element order - 1 for linear elements, 2 for quadratic elements */
    //[[nodiscard]] inline idx getElementOrder() const { return elementOrder; }
    [[nodiscard]] inline ElementType getElementType() const { return elementType_; }
    /** number of nodes per element */
    [[nodiscard]] inline unsigned int getnNodes() const { return getNodesPerElement(elementType_); }
    /** number of dimensions of mesh - 2 for tris and 3 for tets */
    [[nodiscard]] inline idx getDimension() const { return Dimension; }
    idx getConnectedVolume(const idx e) const;  // returns index to connected volume element, or -1 if not connected to LC1
    inline idx getNode(const idx e, const idx n) const { // returns node n of element e
#ifdef DEBUG
        if ((e >= nElements) || (n >= nNodes)) {
            printf("error - Mesh::getNode(int,int) - index to node out of bounds, bye! ");
            printf("requested elem %u, node %u\n", e, n);
            exit(1);
        }
#endif
    return nodes[e * getnNodes() + n];
    }

    [[nodiscard]] idx getMaterialNumber(idx e) const;   // returns material number of element e
    /**
     * gets alignment layer number, i.e. FixLC 1, 2, 3... for the indexed triangle element.
     * If the element is not an alignment layer, returns 0.
     */
    [[nodiscard]] idx getFixLCNumber(idx e) const;  // gets alignment layer number, i.e. FixLC 1, 2, 3...
    [[nodiscard]] idx getDielectricNumber(idx e) const; // gets dielectric materials number, i.e. Dielectric 1, 2, 3 ....
    [[nodiscard]] double getDeterminant(idx i) const; // returns value of determinant of element i

    void setElementData(ElementType elementType, std::vector<unsigned int> &&nodes, std::vector<unsigned int> &&materials);

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
    void listNodesOfMaterial(std::vector <idx> &nodes, idx mat) const;
    /** Deprecated use the other listFixLCSurfaceNodes(FixLC num) instead */
    //void listFixLCSurfaces(std::vector <idx> &nodes, idx FixLCNum) const; // list all nodes of given FixLC surface number (FixLCNum = 1,2,3...)
    /**
     * Return a vector of all surface nodes by given FIXLC number.
     */
    std::unordered_set<idx> listFixLCSurfaceNodes(const idx FixLCNum) const;
    std::unordered_set<idx> findElectrodeSurfaceNodes(idx electrodeNumber) const;
    [[nodiscard]] std::vector<unsigned int> findElementsWhere(std::function<bool(unsigned int)> &predicate) const;

    /**
     * Find set of all node indexes from elements where predicate evaluates to true
     * The input argument for predicate is index to element.
     */
    [[nodiscard]] std::set<unsigned int> findNodesWhere(std::function<bool(unsigned int)> &predicate) const;

    bool containsCoordinate(idx elem, const Coordinates& coordinates, const Vec3 p) const; // checks whether point p is within element elem
    void CompleteNodesSet(const idx elem, std::vector<idx> &nodes) const; // completes nodes vector with those from element, if nodes is empty returns all elements
    void calculateDeterminants3D(const Coordinates &coords); // calculates determinants of all elements
    void calculateSurfaceNormals(const Coordinates &coords, Mesh *tets = NULL);
    void CopyMesh(Mesh *rhs);   // makes this a copy of Mesh* rhs - why does operator= overloading not work???
    void ScaleDeterminants(const double &s);  // scales all determinants by s, e.g. to go to microns
    void calcLocCoords(const idx elem, const Coordinates &coordinates, const Vec3 &targetPoint, double localCoordinates[4]) const; // calculates 4 local coordinates of coordinate cord in element elem
    [[nodiscard]] Vec3 elementCentroid(unsigned int i, const Coordinates &coordinates) const;
    [[nodiscard]] Vec3 getSurfaceNormal(unsigned int i) const;
    [[nodiscard]] bool hasSurfaceNormals() const { return !surfaceNormals.empty(); }
    void gen_p_to_elem(vector<set <idx> > &p_to_elem) const; // generates index from points to mesh elements
    // index number of non-existent neighbours elements equals total number of elements
    // i.e. 1 too large to use as an index
};

#endif



