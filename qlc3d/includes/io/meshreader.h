#ifndef PROJECT_QLC3D_MESHREADER_H
#define PROJECT_QLC3D_MESHREADER_H
#include <string>
#include <filesystem>
#include <utility>
#include <globals.h>
#include <vector>
#include <geom/vec3.h>

class GmshFileData;

enum class MeshFormat {
    UNKNOWN_FORMAT,
    GID_MESH,
    GMSH_ASCII
};

class MeshFormatInspector {
    static bool isGidMesh(std::ifstream &fin) ;
    static bool isGmshAsciiMesh(std::ifstream &fin) ;
public:
    static MeshFormat inspectFormat(const std::filesystem::path &fileName) ;
};

class RawMeshData {

public:
  std::vector<Vec3> points;
  std::vector<idx> tetNodes;
  std::vector<idx> tetMaterials;
  std::vector<idx> triNodes;
  std::vector<idx> triMaterials;
  RawMeshData(std::vector<Vec3> points,
              std::vector<idx> tetNodes, std::vector<idx> tetMaterials,
              std::vector<idx> triNodes, std::vector<idx> triMaterials);
};

/**
 * Read mesh data from provided mesh file path. The format of the file is deduced by inspecting it's contents.
 */
class MeshReader {
    static void readGmsMesh(const std::filesystem::path &fileName,
                            std::vector<Vec3> &pointsOut,
                            std::vector<idx> &tetNodes,
                            std::vector<idx> &tetMaterials,
                            std::vector<idx> &triNodes,
                            std::vector<idx> &triMaterials);

    static void copyGmshTriangleData(const GmshFileData &data, idx** trisOut, idx** triMaterialsOut, idx *numTrisOut);

public:
    static RawMeshData readMesh(const std::filesystem::path &fileName);
};

void ReadGiDMesh3D(const std::filesystem::path &fileName,
                   std::vector<Vec3> &pointsOut,
                   std::vector<idx> &tetNodes,
                   std::vector<idx> &tetMaterials,
                   std::vector<idx> &triNodes,
                    std::vector<idx> &triMaterials);
#endif //PROJECT_QLC3D_MESHREADER_H
