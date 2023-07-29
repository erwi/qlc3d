#ifndef PROJECT_QLC3D_GMSH_READ_H
#define PROJECT_QLC3D_GMSH_READ_H
#include <string>
#include <fstream>
#include <filesystem>
#include <utility>
#include <memory>
#include <vector>
#include <unordered_map>
#include "util/logging.h"
#include "geom/vec3.h"

class SectionMeshFormat {
public:
    const double _versionNumber;
    const int _fileType;
    /** dataSize = sizeof(<floating point number used>) */
    const int _dataSize;

    SectionMeshFormat(double versionNumber, int fileType, int dataSize) :
            _versionNumber{versionNumber}, _fileType{fileType}, _dataSize{dataSize} {}
};

struct SectionPhysicalNames {
    /** number of defined physical names */
    const size_t _numNames;
    /** human readable names by name id tag */
    const std::unordered_map<size_t, std::string> _physicalNames;
    SectionPhysicalNames(size_t numNames, std::unordered_map<size_t, std::string> &physicalNames) :
            _numNames {numNames},
            _physicalNames {physicalNames}{}
};

// "Tag" lines within the Entities section. The tags for points and lines exist, but are not declared here
// as we ignore them and just skip over them when reading the file. We're only interested in surfaces and
// volumes. Furthermore, the tag lines in the file contain additional information (e.g. x,y,z coordinates etc.)
// but these too are skipped and only the information we need is included here (mainly material numbers).
// NOTE: the tags ("ID") are only unique within the entity type. That is, a surface may have same tag number
// as a volume.

struct EntityTag {
    /** the "ID" of this entity */
    const int _tag;
    /** the physical material tags of this entity; */
    const std::vector<int> _physicalTags;
    EntityTag(int tag, std::vector<int> physicalTags):
        _tag {tag}, _physicalTags{std::move(physicalTags)} { }
};

struct SurfaceTag : EntityTag {
    SurfaceTag(int tag, std::vector<int> physicalTags) : EntityTag {tag, std::move(physicalTags)} { }
};
struct VolumeTag : EntityTag {
    VolumeTag(int tag, std::vector<int> physicalTags) : EntityTag {tag, std::move(physicalTags)} { }
};

struct SectionEntities {
    /** number of points in the structure, not mesh points */
    const size_t _numPoints;
    /** number of curves/ines in the structure */
    const size_t _numCurves;
    /** number of surfaces in the structure */
    const size_t _numSurfaces;
    /** number of volumes in the structure */
    const size_t _numVolumes;

    const std::unordered_map<int, SurfaceTag> _surfaceTags;
    const std::unordered_map<int, VolumeTag> _volumeTags;

    SectionEntities(size_t numPoints, size_t numCurves, size_t numSurfaces, size_t numVolumes,
                    std::unordered_map<int, SurfaceTag>&& surfaceTags,
                    std::unordered_map<int, VolumeTag>&& volumeTags):
                    _numPoints { numPoints }, _numCurves { numCurves }, _numSurfaces { numSurfaces }, _numVolumes { numVolumes },
                    _surfaceTags { surfaceTags }, _volumeTags { volumeTags } { }

    /** returns the physical tag numbers for the given surface tag */
    std::vector<int> getPhysicalTagsForSurface(int surfaceTag) const;
    /** returns the physical tag numbers for the given volume tag */
    std::vector<int> getPhysicalTagsForVolume(int volumeTag) const;
};

struct SectionNodes {
    const size_t _numEntityBlocks;
    const size_t _numNodes;
    const size_t _minNodeTag;
    const size_t _maxNodeTag;
    const std::vector<Vec3> _coordinates;

    SectionNodes(size_t numEntityBlocks,
                 size_t numNodes,
                 size_t minNodeTag,
                 size_t maxNodeTag,
                 std::vector<Vec3> && coordinates) :
            _numEntityBlocks{numEntityBlocks}, _numNodes{numNodes}, _minNodeTag{minNodeTag},
            _maxNodeTag{maxNodeTag}, _coordinates{coordinates} { }
};

struct SectionElements {
    static const int ELEMENT_TYPE_TRIANGLE_3_NODES = 2;
    static const int ELEMENT_TYPE_TETRAHEDRON_4_NODES = 4;

    const unsigned int _numTriangles;
    const unsigned int _numTetrahedra;

    const std::vector<unsigned int> _triangleIndices;
    /** id tag to geometry surface */
    const std::vector<unsigned int> _triangleEntityTags;

    const std::vector<unsigned int> _tetrahedraIndices;
    /** id tag to geometry volume */
    const std::vector<unsigned int> _tetrahedraEntityTags;

    SectionElements(unsigned int numTriangles,
                    unsigned int numTetrahedra,
                    std::vector<unsigned int> &&triangleIndices,
                    std::vector<unsigned int> &&triangleEntityTags,
                    std::vector<unsigned int> &&tetrahedraIndices,
                    std::vector<unsigned int> &&tetrahedraEntityTags) :
        _numTriangles{numTriangles}, _numTetrahedra{numTetrahedra},
        _triangleIndices { triangleIndices },
        _triangleEntityTags { triangleEntityTags },
        _tetrahedraIndices { tetrahedraIndices },
        _tetrahedraEntityTags { tetrahedraEntityTags } { }
};

class GmshFileData {
    std::shared_ptr<SectionMeshFormat> _meshFormat;
    std::shared_ptr<SectionPhysicalNames> _physicalNames;
    std::shared_ptr<SectionEntities> _entities;
    std::shared_ptr<SectionNodes> _nodes;
    std::shared_ptr<SectionElements> _elements;

public:
    void setMeshFormat(std::shared_ptr<SectionMeshFormat> meshFormat) { _meshFormat = std::move(meshFormat); }
    void setPhysicalNames(std::shared_ptr<SectionPhysicalNames> &&physicalNames) { _physicalNames = std::move(physicalNames); }
    void setEntities(std::shared_ptr<SectionEntities> &&entities) { _entities = entities; }
    void setNodes(std::shared_ptr<SectionNodes> &&nodes) { _nodes = std::move(nodes); }
    void setElements(std::shared_ptr<SectionElements> &&elements) { _elements = std::move(elements); }
    void validateMeshData() const;

    [[nodiscard]] const std::shared_ptr<SectionMeshFormat> &getMeshFormat() const { return _meshFormat; }
    [[nodiscard]] const std::shared_ptr<SectionPhysicalNames> &getPhysicalNames() const { return _physicalNames; }
    [[nodiscard]] const std::shared_ptr<SectionEntities> &getEntities() const { return _entities; }
    [[nodiscard]] const std::shared_ptr<SectionNodes> &getNodes() const { return _nodes; }
    [[nodiscard]] const std::shared_ptr<SectionElements> &getElements() const { return _elements; }
};

class GmshPhysicalNamesMapper {
    const std::shared_ptr<SectionPhysicalNames> _physicalNames;
    const std::shared_ptr<SectionEntities> _entities;
    const std::shared_ptr<SectionElements> _elements;

    static std::unordered_map<size_t, int> mapPhysicalNamesToNumbers(const std::unordered_map<size_t, std::string> &physicalNamesByTag);
public:
    GmshPhysicalNamesMapper(const std::shared_ptr<SectionPhysicalNames> physicalNames,
                            const std::shared_ptr<SectionEntities> entities,
                            const std::shared_ptr<SectionElements> elements) :
        _physicalNames { physicalNames }, _entities { entities }, _elements { elements } { }

    /**
     * Convert Gmsh "PhysicalNames" strings to qlc3d material number arrays
     * @return vector of same size as number of triangles. The values in the vector correspond to the
     * material number of the triangle element with same ordinal.
     */
    std::vector<unsigned int> mapTriangleNamesToMaterialNumbers();

    /**
     * Convert Gmsh "PhysicalNames" strings to qlc3d material number arrays
     * @return vector of same size as number of tetrahedra. The values in the vector correspond to the
     * material number of the tetrahedron element with same ordinal.
     */
    std::vector<unsigned int> maptTetrahedraNamesToMaterialNumbers();
};

/**
 * Reads the raw data from a Gmsh (ASCII) mesh file. Not all available data, supported by the file format is read,
 * only whatever is needed by qlc3d, i.e. nodes, triangles, tetrahedra, and material definitions.
 *
 * See https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format for more details about the Gmsh file format.
 */
class GmshFileReader {
    std::filesystem::path _fileName;
    size_t _lineNumber = 0;
    std::ifstream _fin;

    void split(const std::string &in, const std::string &pattern, std::vector<std::string> &splitsOut) const;

    /** Read one line of text into input arg line. Return true on success, false if failed (e.g. EOF) */
    bool readLine(std::string &line);
    std::unique_ptr<SectionMeshFormat> readMeshFormat();
    std::unique_ptr<SectionPhysicalNames> readPhysicalNames();
    std::unique_ptr<SectionEntities> readEntities();
    std::unique_ptr<SectionNodes> readNodes();
    std::unique_ptr<SectionElements> readElements();

public:
    std::shared_ptr<GmshFileData> readGmsh(const std::filesystem::path &fileName);
};

#endif //PROJECT_QLC3D_GMSH_READ_H
