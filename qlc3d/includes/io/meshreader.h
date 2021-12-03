#ifndef PROJECT_QLC3D_MESHREADER_H
#define PROJECT_QLC3D_MESHREADER_H
#include <string>
#include <globals.h>

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

void ReadGiDMesh3D(const std::string &fileName,
                   double **p,
                   idx *np,
                   idx **t,
                   idx *nt,
                   idx **e,
                   idx *ne,
                   idx **matt,
                   idx **mate);

void readBinaryMesh(std::string filename ,  // same as above
                    double *&p,
                    idx *&t, idx *&tmat,
                    idx *&e, idx *&emat,
                    idx *np, idx *nt, idx *ne);

#endif //PROJECT_QLC3D_MESHREADER_H
