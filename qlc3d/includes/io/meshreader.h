#ifndef PROJECT_QLC3D_MESHREADER_H
#define PROJECT_QLC3D_MESHREADER_H
#include <string>

enum class MeshFormat {
    UNKNOWN_FORMAT,
    GID_MESH,
    GMSH_ASCII
};

class MeshFormatInspector {
    static bool isGidMesh(std::ifstream &fin) ;
    static bool isGmshAsciiMesh(std::ifstream &fin) ;
public:
    static MeshFormat inspectFormat(const std::string &fileName) ;
};


#endif //PROJECT_QLC3D_MESHREADER_H
