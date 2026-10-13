#include <io/meshwriter.h>

#include <cassert>
#include <cstdio>
#include <string>
#include <vector>

#include <geom/coordinates.h>
#include <mesh/element-split-convert.h>
#include <mesh/mesh.h>
#include <util/exception.h>

void writeLinearTetrahedra(FILE *fid, const Mesh &tetrahedra) {
  assert(tetrahedra.getElementType() == ElementType::LINEAR_TETRAHEDRON);
  for (unsigned i = 0; i < tetrahedra.getnElements(); i++) {
    fprintf(fid, "%i\t%i\t%i\t%i\t%i\t%i\n",
            i + 1, tetrahedra.getNode(i, 0) + 1, tetrahedra.getNode(i, 1) + 1,
            tetrahedra.getNode(i, 2) + 1, tetrahedra.getNode(i, 3) + 1,
            tetrahedra.getMaterialNumber(i));
  }
}

void writeQuadraticTetrahedra(FILE *fid, const Mesh &tetrahedra, const Coordinates &coordinates) {
  assert(tetrahedra.getElementType() == ElementType::QUADRATIC_TETRAHEDRON);

  std::vector<unsigned int> tetrahedronNodes(tetrahedra.getnNodes());
  for (unsigned i = 0; i < tetrahedra.getnElements(); i++) {
    tetrahedra.loadNodes(i, tetrahedronNodes.data());
    const std::vector<std::vector<unsigned int>> splitTetrahedra =
        splitQuadraticTetrahedronToLinear(tetrahedronNodes, coordinates);

    const unsigned int material = tetrahedra.getMaterialNumber(i);
    for (unsigned j = 0; j < splitTetrahedra.size(); j++) {
      fprintf(fid, "%i\t%i\t%i\t%i\t%i\t%i\n", 8 * i + j + 1,
              splitTetrahedra[j][0] + 1, splitTetrahedra[j][1] + 1,
              splitTetrahedra[j][2] + 1, splitTetrahedra[j][3] + 1, material);
    }
  }
}

void writeLinearTriangles(FILE *fid, const Mesh &triangles) {
  assert(triangles.getElementType() == ElementType::LINEAR_TRIANGLE);
  for (unsigned i = 0; i < triangles.getnElements(); i++) {
    fprintf(fid, "%i\t%i\t%i\t%i\t%i\n", i + 1, triangles.getNode(i, 0) + 1,
            triangles.getNode(i, 1) + 1, triangles.getNode(i, 2) + 1,
            triangles.getMaterialNumber(i));
  }
}

void writeQuadraticTriangles(FILE *fid, const Mesh &triangles) {
  assert(triangles.getElementType() == ElementType::QUADRATIC_TRIANGLE);
  std::vector<unsigned int> triangleNodes(triangles.getnNodes());
  for (unsigned i = 0; i < triangles.getnElements(); i++) {
    triangles.loadNodes(i, triangleNodes.data());
    const auto splitTriangles = splitQuadraticTriangleToLinear(triangleNodes);

    const unsigned int material = triangles.getMaterialNumber(i);
    for (unsigned j = 0; j < splitTriangles.size(); j++) {
      fprintf(fid, "%i\t%i\t%i\t%i\t%i\n", 4 * i + j + 1,
              splitTriangles[j][0] + 1, splitTriangles[j][1] + 1,
              splitTriangles[j][2] + 1, material);
    }
  }
}

void writeMeshFile(const Coordinates &coordinates,
                   const Mesh &tetrahedra,
                   const Mesh &triangles,
                   const std::filesystem::path &fileName) {
  const ElementType elementType = tetrahedra.getElementType();
  FILE *fid = fopen(fileName.string().c_str(), "w");
  const unsigned int numberOfPoints = coordinates.size();
  if (fid != nullptr) {
    fputs("MESH    dimension 3 ElemType Tetrahedra  Nnode 4\nCoordinates\n", fid);
    for (unsigned int i = 0; i < numberOfPoints; i++) {
      const auto &point = coordinates.getPoint(i);
      fprintf(fid, "%i\t%f\t%f\t%f\n", i + 1, point.x(), point.y(), point.z());
    }
    fprintf(fid, "end coordinates\n\nElements\n");

    if (elementType == ElementType::LINEAR_TETRAHEDRON) {
      writeLinearTetrahedra(fid, tetrahedra);
    } else if (elementType == ElementType::QUADRATIC_TETRAHEDRON) {
      writeQuadraticTetrahedra(fid, tetrahedra, coordinates);
    } else {
      RUNTIME_ERROR("Only linear and quadratic tetrahedra are supported, got element with " +
                    std::to_string(tetrahedra.getnNodes()) + " nodes.");
    }

    fprintf(fid, "end elements\nMESH    dimension 3 ElemType Triangle  Nnode 3\nCoordinates\nend coordinates\n\nElements\n");
    if (elementType == ElementType::LINEAR_TETRAHEDRON) {
      writeLinearTriangles(fid, triangles);
    } else if (elementType == ElementType::QUADRATIC_TETRAHEDRON) {
      writeQuadraticTriangles(fid, triangles);
    } else {
      RUNTIME_ERROR("Only linear and quadratic tetrahedra are supported, got element with " +
                    std::to_string(tetrahedra.getnNodes()) + " nodes.");
    }

    fprintf(fid, "end elements\n");
    fclose(fid);
  } else {
    RUNTIME_ERROR("Could not open file for output mesh: " + fileName.string());
  }
}
