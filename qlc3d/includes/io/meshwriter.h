#ifndef PROJECT_QLC3D_MESHWRITER_H
#define PROJECT_QLC3D_MESHWRITER_H

#include <filesystem>

class Coordinates;
class Mesh;

/**
 * Writes a mesh in the GiD/LCView text mesh format.
 *
 * @param coordinates Coordinates for the mesh nodes.
 * @param tetrahedra Volume tetrahedra to write.
 * @param triangles Boundary triangles to write.
 * @param fileName Output file path.
 */
void writeMeshFile(const Coordinates &coordinates,
                   const Mesh &tetrahedra,
                   const Mesh &triangles,
                   const std::filesystem::path &fileName);

#endif //PROJECT_QLC3D_MESHWRITER_H
