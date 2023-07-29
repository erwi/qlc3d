#include <io/meshreader.h>
#include <io/gmsh-read.h>
#include <filesystem>
#include <util/logging.h>
#include <util/exception.h>
#include <geom/vec3.h>

using namespace std;

RawMeshData::RawMeshData(std::vector<Vec3> points,
                         std::vector<idx> tetNodes, std::vector<idx> tetMaterials,
                         std::vector<idx> triNodes, std::vector<idx> triMaterials)
  : points(std::move(points)),
  tetNodes(std::move(tetNodes)), tetMaterials(std::move(tetMaterials)),
  triNodes(std::move(triNodes)), triMaterials(std::move(triMaterials)) { }

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

MeshFormat MeshFormatInspector::inspectFormat(const std::filesystem::path &fileName) {
    if (!filesystem::exists(fileName)) {
        throw invalid_argument("no such mesh file " + fileName.string());
    }

    ifstream fin;
    fin.open(fileName);

    if (!isTextFile(fin)) {
        throw runtime_error("detected unsupported, binary, file " + fileName.string());
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

void MeshReader::readGmsMesh(const std::filesystem::path &fileName,
                             std::vector<Vec3> &pointsOut,
                             std::vector<idx> &tetNodes,
                             std::vector<idx> &tetMaterials,
                             std::vector<idx> &triNodes,
                             std::vector<idx> &triMaterials) {

    GmshFileReader reader;
    auto meshData = reader.readGmsh(fileName);

    pointsOut = meshData->getNodes()->_coordinates;
    tetNodes = meshData->getElements()->_tetrahedraIndices;
    triNodes = meshData->getElements()->_triangleIndices;

    // converts gmsh material names to material numbers
    GmshPhysicalNamesMapper mapper(meshData->getPhysicalNames(), meshData->getEntities(), meshData->getElements());
    tetMaterials = mapper.maptTetrahedraNamesToMaterialNumbers();
    triMaterials = mapper.mapTriangleNamesToMaterialNumbers();
}

RawMeshData MeshReader::readMesh(const std::filesystem::path &fileName) {
    MeshFormat format = MeshFormatInspector::inspectFormat(fileName);

    std::vector<Vec3> points;
    std::vector<idx> tetNodes;
    std::vector<idx> tetMaterials;
    std::vector<idx> triNodes;
    std::vector<idx> triMaterials;

    switch (format) {
        case MeshFormat::GID_MESH:
            Log::info("Reading GiD mesh from {}.", fileName.string());
            ReadGiDMesh3D(fileName, points, tetNodes, tetMaterials, triNodes, triMaterials);
            break;
        case MeshFormat::GMSH_ASCII:
            Log::info("Reading Gmsh mesh from {}", fileName.string());
            MeshReader::readGmsMesh(fileName, points, tetNodes, tetMaterials, triNodes, triMaterials);
            break;
        case MeshFormat::UNKNOWN_FORMAT:
            RUNTIME_ERROR(fmt::format("Could not determine mesh format of file {}.", fileName))
        default:
            RUNTIME_ERROR(fmt::format("Unhandled mesh format - did not read mesh from {}.", fileName));
    }

    return {points, tetNodes, tetMaterials, triNodes, triMaterials};
}