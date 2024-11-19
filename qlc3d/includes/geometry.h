#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <mesh.h>
#include <material_numbers.h>
#include <cfloat> // MAXIMUM DOUBLE VALUES NEEDED IN COORDINATE COMPARISON
#include <vector>
#include <list>
#include <memory>
#include <iostream>
#include <algorithm>
#include <limits>
#include <alignment.h>
#include <regulargrid.h>
#include <globals.h>
#include <geom/aabox.h>

class RegularGrid; // forward declaration
class Coordinates; // forward declaration
class Vec3;
class PeriodicNodesMapping;

using namespace std;
class Geometry {
private:
  std::shared_ptr<Coordinates> coordinates_;
  /** Tetrahedra mesh */
  std::shared_ptr<Mesh> t;
  /** Triangles mesh */
  std::shared_ptr<Mesh> e;

  unsigned int npLC;          // number of LC nodes
  std::vector<Vec3> nodeNormals;
  AABox boundingBox;
  RegularGrid *regularGrid;

  std::unique_ptr<PeriodicNodesMapping> periodicNodesMapping;

  void updateMaxNodeNumbers(); // Updates MaxNodeNumbers for surface and tet meshes after a node renumbering

public:
    // UNFORTUNATE HACKERY... SPECIAL ERROR INDEX VALUE FOR AN UNSIGNED INDEX THAT WAS NOT FOUND
    static const unsigned int NOT_AN_INDEX;// = std::numeric_limits<unsigned int>::max();

    Geometry();
    ~Geometry();
    /** for testing purposes */
    void setCoordinates(const std::shared_ptr<Coordinates>& coordinates);

    /** for testing purposes */
    void setTetrahedra(const std::shared_ptr<Mesh>& tetrahedra);

    /** for testing purposes */
    void setTriangles(const std::shared_ptr<Mesh>& triangles);

    void setMeshData(unsigned int elementOrder, const std::shared_ptr<Coordinates> &coordinates,
                     std::vector<unsigned int> &&tetNodes, std::vector<unsigned int> &&tetMaterials,
                     std::vector<unsigned int> &&triNodes, std::vector<unsigned int> &&triMaterials);

    void addCoordinates(const vector<double> &coords);   // adds new coordinates to end of existing ones
    void calculateNodeNormals();      // calculates surface node normals
    void setnpLC(int n);        // set number of LC nodes
    void ReorderDielectricNodes();  // reorder nodes so that dielectric material nodes are last

    void ClearGeometry();       // clears all data for geometry

    void genIndToTetsByCoords(vector <unsigned int> &returnIndex,   // Return index
                              const Coordinates &targetCoordinates,    // coordinates to search
                              const bool &terminateOnError = true, // terminate app. if coord not found
                              const bool &requireLCElement = false); // only LC elements are considered

    bool brute_force_search(unsigned int &ind,  // found index, return value
                            const Vec3 &crd,   // target coordinate to find
                            const bool &terminateOnError = true, // whether application terminates if coord is not found
                            const bool &requireLCElement = false // only LC elements are considered
                           );
    /** finds which tet contains target point */
    size_t recursive_neighbour_search(const Vec3 &targetPoint,
                                      const vector< set < unsigned int> > &p_to_t,
                                      const size_t &currentTetIndex,
                                      std::set<size_t> &tetHistory,
                                      const bool &requireLCElement = false   // only LC elements are considered
                                     );

    void makeRegularGrid(const size_t &nx, // GENERATES REGULAR GRID LOOKUP INDEXES AND WEIGHTS
                         const size_t &ny,
                         const size_t &nz);

    void setTo(Geometry *geom);                     // makes this = geom

    [[nodiscard]] unsigned int getnp() const;
    [[nodiscard]] unsigned int getnpLC() const { return npLC; }
    [[nodiscard]] const Coordinates& getCoordinates() const;
    [[nodiscard]] const AABox& getBoundingBox() const { return boundingBox; }
    [[nodiscard]] const std::vector<Vec3>& getNodeNormals() const;
    [[nodiscard]] Vec3 getNodeNormal(unsigned int i) const;
    [[nodiscard]] const Mesh& getTetrahedra() const { return *t; }
    [[nodiscard]] Mesh& getTetrahedra() { return const_cast<Mesh&>(*t); }
    [[nodiscard]] const Mesh& getTriangles() const { return *e; }
    [[nodiscard]] Mesh& getTriangles() { return const_cast<Mesh&>(*e); }
    [[nodiscard]] RegularGrid* getRegularGrid() const { return regularGrid; }

    [[nodiscard]] PeriodicNodesMapping& createPeriodicNodesMapping();
    /** releases memory for periodic nodes mapping */
    void clearPeriodicNodesMapping();
};
#endif

