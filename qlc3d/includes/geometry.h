#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <mesh.h>
#include <material_numbers.h>
#include <float.h> // MAXIMUM DOUBLE VALUES NEEDED IN COORDINATE COMPARISON
#include <vector>
#include <list>
#include <memory>
#include <iostream>
#include <algorithm>
#include <limits>
#include <alignment.h>
#include <regulargrid.h>
#include <globals.h>
#define EPS 1e-7
class RegularGrid; // forward declaration
class Coordinates; // forward declaration
class Vec3;

using namespace std;
class Geometry {
private:
    std::shared_ptr<Coordinates> coordinates_;

    unsigned int npLC;          // number of LC nodes
    std::vector<Vec3> nodeNormals;

    double Xmin;
    double Xmax;
    double Ymin;
    double Ymax;
    double Zmin;
    double Zmax;
    bool left_right_is_periodic;
    bool front_back_is_periodic;
    bool top_bottom_is_periodic;
    vector < list <int> > peri_equ_nodes;
    vector<size_t> periNodes_;
    void setEdgePeriNodes(list <size_t> &edge0,
                          list <size_t> &edge1,
                          const int &dim);   // edge direction 0,1,2 -> x,y,z
    void setFacePeriNodes(list <size_t> &face0,
                          list <size_t> &face1,
                          const int &norm);    // face normal 0,1,2 -> x,y,z
    void updateMaxNodeNumbers(); // Updates MaxNodeNumbers for surface and tet meshes after a node renumbering
public:
    // UNFORTUNATE HACKERY... SPECIAL ERROR INDEX VALUE FOR AN UNSIGNED INDEX THAT WAS NOT FOUND
    static const unsigned int NOT_AN_INDEX;// = std::numeric_limits<unsigned int>::max();
    std::shared_ptr<Mesh> t;                        // volume mesh
    std::shared_ptr<Mesh> e;                        // surface mesh
    RegularGrid *regularGrid;
    Geometry();
    ~Geometry();
    void setCoordinates(const std::shared_ptr<Coordinates>& coordinates);

    void setMeshData(const std::shared_ptr<Coordinates> &coordinates,
                     std::vector<unsigned int> &&tetNodes, std::vector<unsigned int> &&tetMaterials,
                     std::vector<unsigned int> &&triNodes, std::vector<unsigned int> &&triMaterials);

    void addCoordinates(const vector<double> &coords);   // adds new coordinates to end of existing ones
    void calculateNodeNormals();      // calculates surface node normals
    void setnp(int n);
    void setnpLC(int n);        // set number of LC nodes
    void ReorderDielectricNodes();  // reorder nodes so that dielectric material nodes are last
    void makePeriEquNodes();    // generates periodic equivalent nodes index
    void ClearGeometry();       // clears all data for geometry
    bool getleft_right_is_periodic() const;
    bool getfront_back_is_periodic() const;
    bool gettop_bottom_is_periodic() const;

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
    void checkForPeriodicGeometry();    // detects type of periodicity of the strucuture

    size_t getPeriodicEquNode(const size_t &i) const { // RETURNS INDEX TO NODE PERIODIC TO i
        if (i < periNodes_.size())
            return periNodes_[i];
        else
            return i;
    }

    unsigned int getnp() const;

    unsigned int getnpLC() const  {
        return npLC;
    }

    [[nodiscard]] const Coordinates& getCoordinates() const {
      if (!coordinates_) {
        throw std::runtime_error("Coordinates not set");
      }
      return *coordinates_;
    }

    double getpX(int i)const;   // return node coordinates at node i
    double getpY(int i)const;
    double getpZ(int i)const;
    double getXmin();
    double getXmax();
    double getYmin();
    double getYmax();
    double getZmin();
    double getZmax();
    double getAbsXDist(int i , double x);   // gets absolute distance between x-coord of node i and x
    double getAbsYDist(int i , double y);   //
    double getAbsZDist(int i , double z);   //
    double getAbsDistSqr(const unsigned int i , const double *const coord) const;
    //void genIndWeakSurfaces(Alignment &alignment);  // generates index to weak surface elements
    //
    // NodeNormal methods
    [[nodiscard]] const std::vector<Vec3>& getNodeNormals() const;
    Vec3 getNodeNormal(unsigned int i) const;
    void countNodeReferences(vector <int> &refc, Mesh &mesh); // counts the number of times each node is used in mesh. DEBUG

    [[nodiscard]] const Mesh& getTetrahedra() const { return *t; }
    [[nodiscard]] Mesh& getTetrahedra() { return const_cast<Mesh&>(*t); }
    [[nodiscard]] const Mesh& getTriangles() const { return *e; }
};
#endif

