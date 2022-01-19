#include <io/meshreader.h>
#include <io/gmsh-read.h>
#include <filesystem>
#include <util/logging.h>
#include <util/exception.h>

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

void MeshReader::copyGmshCoordinateData(const GmshFileData &data, double **pointsOut, idx *numPointsOut) {
    const shared_ptr<SectionNodes> &nodesData = data.getNodes();
    *numPointsOut = nodesData->_numNodes;

    auto coordinates = nodesData->_coordinates;
    *pointsOut = (double*) malloc(3 * (*numPointsOut) * sizeof(double));
    std::copy(coordinates.begin(), coordinates.end(), *pointsOut);
}

void MeshReader::copyGmshTriangleData(const GmshFileData &data, idx **trisOut, idx **triMaterialsOut, idx *numTrisOut) {
    auto &elements = data.getElements();
    *numTrisOut = elements->_numTriangles;

    // element indices
    *trisOut = (idx*) malloc(3 * (*numTrisOut) * sizeof(idx));
    std::copy(elements->_triangleIndices.begin(), elements->_triangleIndices.end(), *trisOut);

    // material numbers
    GmshPhysicalNamesMapper mapper(data.getPhysicalNames(), data.getEntities(), data.getElements());
    auto materialNumbers = mapper.mapTriangleNamesToMaterialNumbers();
    *triMaterialsOut = (idx*) malloc(*numTrisOut * sizeof(idx));
    std::copy(materialNumbers.begin(), materialNumbers.end(), *triMaterialsOut);
}

void MeshReader::copyGmshTetData(const GmshFileData &data, idx **tetsOut, idx **tetMaterialsOut, idx *numTetsOut) {
    auto &elements = data.getElements();
    *numTetsOut = elements->_numTetrahedra;

    // element indices
    *tetsOut = (idx*) malloc(4 * (*numTetsOut) * sizeof (idx));
    std::copy(elements->_tetrahedraIndices.begin(), elements->_tetrahedraIndices.end(), *tetsOut);

    // material numbers
    GmshPhysicalNamesMapper mapper(data.getPhysicalNames(), data.getEntities(), data.getElements());
    auto materialNumbers = mapper.maptTetrahedraNamesToMaterialNumbers();
    *tetMaterialsOut = (idx*) malloc(*numTetsOut * sizeof(idx));
    std::copy(materialNumbers.begin(), materialNumbers.end(), *tetMaterialsOut);
}

void MeshReader::readGmsMesh(const std::string &fileName, double **pointsOut, idx *numPointsOut, idx **tetsOut, idx *numTetsOut, idx **trisOut, idx *numTrisOut,
                             idx **tetMaterials, idx **triMaterials) {
    GmshFileReader reader;
    auto meshData = reader.readGmsh(fileName);

    // copy gmsh data into output arrays
    copyGmshCoordinateData(*meshData, pointsOut, numPointsOut);
    copyGmshTriangleData(*meshData, trisOut, triMaterials, numTrisOut);
    copyGmshTetData(*meshData, tetsOut, tetMaterials, numTetsOut);
}

void MeshReader::readMesh(const std::string &fileName,
                          double **pointsOut,
                          idx *numPointsOut,
                          idx **tetsOut,
                          idx *numTetsOut,
                          idx **trisOut,
                          idx *numTrisOut,
                          idx **tetMaterialsOut,
                          idx **triMaterialsOut) {
    MeshFormat format = MeshFormatInspector::inspectFormat(fileName);

    switch (format) {
        case MeshFormat::GID_MESH:
            Log::info("Reading GiD mesh from {}.", fileName);
            ReadGiDMesh3D(fileName, pointsOut, numPointsOut, tetsOut, numTetsOut, trisOut, numTrisOut, tetMaterialsOut, triMaterialsOut);
            break;
        case MeshFormat::GMSH_ASCII:
            Log::info("Reading Gmsh mesh from {}", fileName);
            readGmsMesh(fileName, pointsOut, numPointsOut, tetsOut, numTetsOut, trisOut, numTrisOut, tetMaterialsOut, triMaterialsOut);
            break;
        case MeshFormat::UNKNOWN_FORMAT:
            RUNTIME_ERROR(fmt::format("Could not determine mesh format of file {}.", fileName))
        default:
            RUNTIME_ERROR(fmt::format("Unhandled mesh format - did not read mesh from {}.", fileName));
    }
}