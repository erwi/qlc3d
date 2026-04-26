#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <mesh/mesh.h>
#include <vector>
#include <memory>
#include <geom/aabox.h>

class Coordinates;
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

  unsigned int npLC;
  std::vector<Vec3> nodeNormals;
  AABox boundingBox;

  std::unique_ptr<PeriodicNodesMapping> periodicNodesMapping;

  void updateMaxNodeNumbers();

public:
    static const unsigned int NOT_AN_INDEX;

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

    void addCoordinates(const vector<double> &coords);
    void calculateNodeNormals();
    void setnpLC(int n);
    void ReorderDielectricNodes();

    void ClearGeometry();

    void genIndToTetsByCoords(vector <unsigned int> &returnIndex,
                              const Coordinates &targetCoordinates,
                              const bool &terminateOnError = true,
                              const bool &requireLCElement = false);

    bool brute_force_search(unsigned int &ind,
                            const Vec3 &crd,
                            const bool &terminateOnError = true,
                            const bool &requireLCElement = false
                           );

    void setTo(Geometry *geom);

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

    [[nodiscard]] PeriodicNodesMapping& createPeriodicNodesMapping();
    /** releases memory for periodic nodes mapping */
    void clearPeriodicNodesMapping();
};
#endif

