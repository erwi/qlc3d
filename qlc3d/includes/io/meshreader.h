#ifndef PROJECT_QLC3D_MESHREADER_H
#define PROJECT_QLC3D_MESHREADER_H
#include <string>
#include <filesystem>
#include <globals.h>

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

/**
 * Read mesh data from provided mesh file path. The format of the file is deduced by inspecting it's contents.
 */
class MeshReader {
    static void readGmsMesh(const std::filesystem::path &fileName,
                            double **pointsOut,
                            idx *numPointsOut,
                            idx **tetsOut,
                            idx *numTetsOut,
                            idx **trisOut,
                            idx *numTrisOut,
                            idx **tetMaterials,
                            idx **triMaterials);

    static void copyGmshCoordinateData(const GmshFileData &data, double **pointsOut, idx *numPointsOut);
    static void copyGmshTriangleData(const GmshFileData &data, idx** trisOut, idx** triMaterialsOut, idx *numTrisOut);
    static void copyGmshTetData(const GmshFileData &data, idx **tetsOut, idx ** tetMaterialsOut, idx *numTetsOut);
public:
    static void readMesh(const std::filesystem::path &fileName,
                         double **pointsOut, // TODO: get rid of these raw arrays and return a mesh object
                         idx *numPointsOut,
                         idx **tetsOut,
                         idx *numTetsOut,
                         idx **trisOut,
                         idx *numTrisOut,
                         idx **tetMaterialsOut,
                         idx **triMaterialsOut);
};

void ReadGiDMesh3D(const std::filesystem::path &fileName,
                   double **p,
                   idx *np,
                   idx **t,
                   idx *nt,
                   idx **e,
                   idx *ne,
                   idx **matt,
                   idx **mate);
#endif //PROJECT_QLC3D_MESHREADER_H
