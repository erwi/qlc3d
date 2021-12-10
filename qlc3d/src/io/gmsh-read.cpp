#include <io/gmsh-read.h>
#include <material_numbers.h>
#include <cassert>
#include <filesystem>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>


using namespace std;

void GmshFileReader::split(const string &str, const string &pattern, vector<string> &destination) const {
    size_t pos = 0;
    destination.clear();

    while (pos < str.size()) {
        size_t end = str.find(pattern, pos);

        end = min(end, str.size());
        destination.push_back(str.substr(pos, end - pos));

        pos = end + 1;
    }
}

void GmshFileReader::log(const string &str) const {
    cout << str << endl;
}

/**
 * Convert Gmsh physical names into qlc3d material numbers
 * @param physicalNamesByTag
 * @return map
 */
std::unordered_map<size_t, int> GmshPhysicalNamesMapper::mapPhysicalNamesToNumbers(
        const std::unordered_map<size_t, std::string> &physicalNamesByTag) {
    // qlc3d material number by gmsh physical name tag
    unordered_map<size_t, int> materialNumberByTag;

    for (auto &entry : physicalNamesByTag) {
        size_t tag = entry.first;
        string name = entry.second;
        int materialNumber = qlc3d::MATERIAL_NUMBER_BY_NAME.at(name);
        materialNumberByTag.emplace(tag, materialNumber);
    }
    return materialNumberByTag;
}

std::vector<unsigned int> GmshPhysicalNamesMapper::mapTriangleNamesToMaterialNumbers() {
    auto materialNumberByTag = mapPhysicalNamesToNumbers(_physicalNames->_physicalNames);

    vector<unsigned int> triangleMaterialsOut(_elements->_numTriangles, MAT_INVALID);
    for (size_t i = 0; i < _elements->_numTriangles; i++) {
        // the surface should have one or maybe two physical names. E.g "periodic" or "fixlc1" and "electrode1"
        int surfaceTagNumber = _elements->_triangleEntityTags[i];
        SurfaceTag st = _entities->_surfaceTags.at(surfaceTagNumber);
        const auto &physicalTags = st._physicalTags;

        if (physicalTags.empty() || physicalTags.size() > 2) {
            throw runtime_error("surface " + to_string(surfaceTagNumber) + " should have 1 or 2 physical names" +
            ", but found " + to_string(physicalTags.size()));
        }

        // if multiple physical names exist, sum them together when converting to qlc3d material number
        int materialNumber = 0;
        for (auto &p : physicalTags) {
            materialNumber += materialNumberByTag.at(p);
        }

        triangleMaterialsOut[i] = materialNumber;
    }

    return triangleMaterialsOut;
}

std::vector<unsigned int> GmshPhysicalNamesMapper::maptTetrahedraNamesToMaterialNumbers() {
    auto materialNumberByTag = GmshPhysicalNamesMapper::mapPhysicalNamesToNumbers(_physicalNames->_physicalNames);

    vector<unsigned int> tetrahedraMaterialsOut(_elements->_numTetrahedra, MAT_INVALID);
    for (size_t i = 0; i < _elements->_numTetrahedra; i++) {
        // the volume should have exactly one physical name, "domain1" or "dielectricX"
        int volumeTagNumber = _elements->_tetrahedraEntityTags[i];
        VolumeTag vt = _entities->_volumeTags.at(volumeTagNumber);
        const auto &physicalTags = vt._physicalTags;

        if (physicalTags.size() != 1) {
            throw runtime_error("volume " + to_string(volumeTagNumber) + " should have 1 physical name, but found "
                + to_string(physicalTags.size()));
        }

        int materialNumber = materialNumberByTag.at(physicalTags[0]);
        tetrahedraMaterialsOut[i] = materialNumber;
    }

    return tetrahedraMaterialsOut;
}

std::runtime_error GmshFileReader::fileReadException(const std::string &message) {
    return runtime_error("error reading line " + to_string(_lineNumber) + " of file " + _fileName + " : " + message);
}

std::unique_ptr<SectionMeshFormat> GmshFileReader::readMeshFormat() {
    log("reading mesh format");
    string line;
    readLine(line);

    vector<string> splits;
    split(line, " ", splits);

    double fileVersion = stod(splits[0]);
    int fileType = stoi(splits[1]);
    int dataSize = stoi(splits[2]);

    readLine(line);
    if (line.find("$EndMeshFormat")) {
        throw fileReadException("expected $EndMeshFormat, but found "  + line);
    }

    return std::make_unique<SectionMeshFormat>(fileVersion, fileType, dataSize);
}

std::unique_ptr<SectionPhysicalNames> GmshFileReader::readPhysicalNames() {
    log("reading physical names");
    string line;
    vector<string> splits;

    readLine(line);
    size_t numNames = stoul(line);

    std::unordered_map<size_t, string> materials;

    while (readLine(line)) {
        if (line.find("$EndPhysicalNames") == 0) {
            return make_unique<SectionPhysicalNames>(numNames, materials);
        }
        split(line, " ", splits);
        // line contains:
        // physical-dimension physical-tag "physical-name"
        size_t physicalTag = stoul(splits[1]);
        // remove leading and trailing double quotes
        string physicalName = splits[2];
        size_t start = physicalName.find_first_of('\"') + 1;
        size_t end = physicalName.find_last_of('\"');
        physicalName = physicalName.substr(start, end - start);

        // convert the physical name to all-lowercase
        std::transform(physicalName.begin(), physicalName.end(), physicalName.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        if (qlc3d::MATERIAL_NUMBER_BY_NAME.count(physicalName) == 0) {
            throw fileReadException("physical name = \"" + physicalName + "\" not recognised");
        }

        materials[physicalTag] = physicalName;
    }

    throw fileReadException("expected $EndPhysicalNames tag");
}

std::unique_ptr<SectionEntities> GmshFileReader::readEntities() {
    // Entities are points, lines, faces, and volumes of the geometry. They
    // are different from nodes and elements of the mesh.
    log("reading entities");
    string line;
    vector<string> splits;

    // first line contains: number of points, lines/curves, surfaces, volumes
    readLine(line);
    split(line, " ", splits);
    size_t numPoints = stoul(splits[0]);
    size_t numCurves = stoul(splits[1]);
    size_t numSurfaces = stoul(splits[2]);
    size_t numVolumes = stoul(splits[3]);

    // skip over points and lines definitions as they are not needed for anything
    for (size_t i = 0; i < numPoints + numCurves; i++) {
        readLine(line);
    }

    // read surfaces
    std::unordered_map<int, SurfaceTag> surfaceTags;
    for (size_t i = 0; i < numSurfaces; i++) {
        readLine(line);
        split(line, " ", splits);
        // the line contains:
        // tag, minX, minY, minZ, maxX, maxY, maxZ, numPhysicalTags, physicalTags..., numBoundingSurfaces, surfaceTags...
        // We're only interested in the 'tag' and physicalTags
        int surfaceTag = stoi(splits[0]);
        size_t numPhysicalTags = stoul(splits[7]);
        vector<int> physicalTags;
        for (size_t j = 8; j < 8 + numPhysicalTags; j++) {
            physicalTags.push_back(stoi(splits[j]));
        }
        surfaceTags.emplace(surfaceTag, SurfaceTag(surfaceTag, physicalTags));
    }

    // read volumes
    std::unordered_map<int, VolumeTag> volumeTags;
    for (size_t i = 0; i < numVolumes; i++) {
        readLine(line);
        split(line, " ", splits);
        // the line contains:
        // volumeTag(int) minX(double) minY(double) minZ(double)
        //    maxX(double) maxY(double) maxZ(double)
        //    numPhysicalTags(size_t) physicalTags(int) ...
        //    numBoundngSurfaces(size_t) surfaceTags(int) ...
        // We're only interested in the 'volumeTag' and 'physicalTags'
        int volumeTag = stoi(splits[0]);
        size_t numPhysicalTags = stoul(splits[7]);
        vector<int> physicalTags;
        for (size_t j = 8; j < 8 + numPhysicalTags; j++) {
            physicalTags.push_back(stoi(splits[j]));
        }
        volumeTags.emplace(volumeTag, VolumeTag(volumeTag, physicalTags));
    }

    // the next line should indicate end of entities
    readLine(line);
    if (line.find("$EndEntities") == string::npos) {
        throw fileReadException("expected $EndEntities tag");
    }
    return std::make_unique<SectionEntities>(numPoints, numCurves, numSurfaces, numVolumes,
                                             std::move(surfaceTags),
                                             std::move(volumeTags));
}

std::unique_ptr<SectionNodes> GmshFileReader::readNodes() {
    string line;
    vector<string> splits;

    readLine(line);
    split(line, " ", splits);

    size_t numEntityBlocks = stoul(splits[0]);
    size_t numNodes = stoul(splits[1]);
    size_t minNodeTag = stoul(splits[2]);
    size_t maxNodeTag = stoul(splits[3]);
    log("reading " + to_string(numNodes) + " nodes");

    if (minNodeTag != 1) {
        throw fileReadException("expected minNodeTag = 1 , got " + to_string(minNodeTag));
    }

    if (maxNodeTag != numNodes) {
        throw fileReadException("expected maxNodeTag = " + to_string(numNodes) + ", got " + to_string(maxNodeTag));
    }

    // read the actual coordinate values
    vector<double> coordinates(3 * numNodes, std::numeric_limits<double>::quiet_NaN()); // initialise all coordinates as NaN

    // read until end of $Nodes section
    size_t counter = 0;
    while(readLine(line)) {
        if (line.find("$EndNodes") == 0) {
            assert(coordinates.size() / 3 == numNodes);
            return make_unique<SectionNodes>(numEntityBlocks, numNodes, minNodeTag, maxNodeTag, std::move(coordinates));
        }

        split(line, " ", splits);
        // if line contains :
        //  1 number, it indicates the next coordinate index
        //  3 numbers, it contains the coordinate xyz value
        //  4 numbers, it contains [entityDim(int) entityTag(int) parametric(int; 0 or 1) numNodesInBlock(size_t)],

        if (splits.size() == 3) {
            coordinates[counter++] = stod(splits[0]);
            coordinates[counter++] = stod(splits[1]);
            coordinates[counter++] = stod(splits[2]);
        }
    }

    // no end of Nodes section tag found
    throw fileReadException("expected $EndNodes tag");
}

std::unique_ptr<SectionElements> GmshFileReader::readElements() {
    log("reading elements");
    string line;
    vector<string> splits;

    readLine(line);
    split(line, " ", splits);

    // size_t numEntityBlocks = stoul(splits[0]);
    size_t numElements = stoul(splits[1]);
    size_t minElementTag = stoul(splits[2]);
    size_t maxElementTag = stoul(splits[3]);

    if (minElementTag != 1) {
        throw fileReadException("expected minElementTag = 1, got " + to_string(minElementTag));
    }
    if (maxElementTag != numElements) {
        throw fileReadException("expected maxElementTag = " + to_string(numElements) + ", got " + to_string(maxElementTag));
    }

    // Read the rest of the $Elements section
    vector<size_t> triangles;   // node indices to triangles
    vector<size_t> tetrahedra;  // node indices to tetrahedra
    vector<int> triangleTags;   // surface tag used in geometry creation
    vector<int> tetrahedraTags; // volume tag used geometry creation

    while (readLine(line)) {
        if (line.find("$EndElements") == 0) {
            size_t numTriangles = triangles.size() / 3;
            size_t numTetrahedra = tetrahedra.size() / 4;
            log("triangles count = " + to_string(numTriangles) + " tetrahedra count = " + to_string(numTetrahedra));

            return make_unique<SectionElements>(
                    triangles.size() / 3,
                    tetrahedra.size() / 4,
                    std::move(triangles),
                    std::move(triangleTags),
                    std::move(tetrahedra),
                    std::move(tetrahedraTags));
        }

        split(line, " ", splits);
        if (splits.size() != 4) {
            throw fileReadException("expected element block descriptor with 4 values, got " + to_string(splits.size()) + " values");
        }

        // int entityDim = stoi(splits[0]); // 1, 2, or 3
        int entityTag = stoi(splits[1]); // surface or volume number of the geometry, not mesh.
        int elementType = stoi(splits[2]); // 2 = triangle, 4 = tetrahedron
        size_t numElementsInBlock = stoul(splits[3]);

        // read one block of elements. All elements in the block should be of same type and "material"
        for (size_t i = 0; i < numElementsInBlock; i++) {
            readLine(line);
            if (elementType == SectionElements::ELEMENT_TYPE_TRIANGLE_3_NODES) {
                split(line, " ", splits);
                if (splits.size() != 4) {
                    throw fileReadException("expected 4 values when reading triangle element, got " + to_string(splits.size()));
                }
                triangles.push_back(stoul(splits[1]) - 1);
                triangles.push_back(stoul(splits[2]) - 1);
                triangles.push_back(stoul(splits[3]) - 1);
                triangleTags.push_back(entityTag);
            } else if (elementType == SectionElements::ELEMENT_TYPE_TETRAHEDRON_4_NODES) {
                split(line, " ", splits);
                if (splits.size() != 5) {
                    throw fileReadException("expected 5 values when reading tetrahedron element, got " + to_string(splits.size()));
                }
                tetrahedra.push_back(stoul(splits[1]) - 1);
                tetrahedra.push_back(stoul(splits[2]) - 1);
                tetrahedra.push_back(stoul(splits[3]) - 1);
                tetrahedra.push_back(stoul(splits[4]) - 1);
                tetrahedraTags.push_back(entityTag);
            } else {
                // do nothing, we don't care about this element type
            }
        }
    }

    // no end of Elements section tag found
    throw fileReadException("expected $EndElements tag");
}

bool GmshFileReader::readLine(std::string &line) {
    _lineNumber++;
    return (bool) getline(_fin, line);
}

std::shared_ptr<GmshFileData> GmshFileReader::readGmsh(const string &fileName) {
    _fileName = fileName;
    _lineNumber = 0;

    if (!std::filesystem::exists(fileName)) {
        throw std::invalid_argument("no such mesh file " + fileName);
    }
    _fin.open(fileName);

    string line;
    auto data = make_shared<GmshFileData>();
    while(readLine(line)) {
        if (line.find("$MeshFormat") == 0) {
            data->setMeshFormat(readMeshFormat());
        }

        if (line.find("$Nodes") == 0) {
            data->setNodes(readNodes());
        }

        if (line.find("$Elements") == 0) {
            data->setElements(readElements());
        }

        if (line.find("$PhysicalNames") == 0) {
            data->setPhysicalNames(readPhysicalNames());
        }

        if (line.find("$Entities") == 0) {
            data->setEntities(readEntities());
        }
    }

    // Check that all required sections were loaded
    if (data->getMeshFormat() == nullptr) {
        throw fileReadException("No $MeshFormat section found");
    }

    if (data->getPhysicalNames() == nullptr) {
        throw fileReadException("No $PhysicalNames section found");
    }

    if (data->getNodes() == nullptr) {
        throw fileReadException("No $Nodes section found");
    }

    if (data->getElements() == nullptr) {
        throw fileReadException("No $Elements section found");
    }

    return data;
}
