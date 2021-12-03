#include <io/meshreader.h>
#include <filesystem>
#include <fstream>

using namespace std;

/**
 * Checks up to 256 bytes of the file to try to check whether it is a text or binary file.
 * The input stream is re-wound to the start of file before the function returns.
 *
 * @param fin file input stream.
 * @return true if the file contains text, false if binary data
 */
bool isTextFile(std::ifstream &fin) {
    char someBytes[100];
    streamsize numRead = fin.readsome(someBytes, sizeof (someBytes));

    for (auto i = 0 ; i < numRead ; ++i) {
        if (!std::isprint(someBytes[0])) {
            return false;
        }
    }

    // rewind to start of file.
    fin.seekg(0, ios::beg);
    return true;
}

bool MeshFormatInspector::isGidMesh(std::ifstream &fin) {
    string line;
    getline(fin, line);
    fin.seekg(0, ios::beg);
    // GiD meshes start with a line like:
    // MESH    dimension 3 ElemType Tetrahedra  Nnode 4
    // Note the suspicious white-space. We split the comparison so that it's ignored
    return ((line.find("MESH") == 0) &&
            (line.find("dimension 3 ElemType Tetrahedra") != string::npos));
}

bool MeshFormatInspector::isGmshAsciiMesh(std::ifstream &fin) {
    string line;
    getline(fin, line);
    fin.seekg(0, ios::beg);

    // expect a valid Gmsh ascii mesh to start with the line
    // $MeshFormat
    return line.find("$MeshFormat") == 0;
}

MeshFormat MeshFormatInspector::inspectFormat(const std::string &fileName) {
    if (!filesystem::exists(fileName)) {
        throw invalid_argument("no such mesh file " + fileName);
    }

    ifstream fin;
    fin.open(fileName);

    if (!isTextFile(fin)) {
        throw runtime_error("detected unsupported, binary, file " + fileName);
    }

    if (isGidMesh(fin)) {
        return MeshFormat::GID_MESH;
    }
    else if (isGmshAsciiMesh(fin)) {
        return MeshFormat::GMSH_ASCII;
    }
    else {
        return MeshFormat::UNKNOWN_FORMAT;
    }
}
