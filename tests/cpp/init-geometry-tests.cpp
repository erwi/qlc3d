#include <catch.h>
#include <inits.h>
#include <geometry.h>
#include <io/meshreader.h>
#include <geom/coordinates.h>
#include <electrodes.h>
#include <alignment.h>
#include <material_numbers.h>
#include <test-util.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <unordered_set>

namespace {

std::string readTextFile(const std::filesystem::path &path) {
  std::ifstream fin(path);
  REQUIRE(fin.is_open());
  std::ostringstream buffer;
  buffer << fin.rdbuf();
  return buffer.str();
}

std::string replaceOnce(std::string text, const std::string &from, const std::string &to) {
  const auto pos = text.find(from);
  REQUIRE(pos != std::string::npos);
  text.replace(pos, from.size(), to);
  return text;
}

std::filesystem::path quadraticMeshResource() {
  return {TestUtil::RESOURCE_SMALL_CUBE_QUADRATIC_GMSH_MESH};
}

} // namespace


TEST_CASE("Reorder quadratic element node order") {
  // GMSH quadratic tet element node ordering is different from the usual used in literature.
  // The difference is that the last two nodes are swapped. When reading the GMSH mesh and preparing
  // the geometry, we should swap the last two nodes if they are in the GMSH order.

  Geometry geom;
  Electrodes electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));

  prepareGeometry(geom,
                  TestUtil::RESOURCE_SMALL_CUBE_QUADRATIC_GMSH_MESH,
                  electrodes,
                  alignment);

  Mesh &tets = geom.getTetrahedra();
  const Coordinates &coords = geom.getCoordinates();
  REQUIRE(tets.getnElements() == 24);

  // check node order in each tetrahedra
  for (unsigned int i = 0; i < tets.getnElements(); i++) {
    unsigned int tetNodes[10] = {};
    Vec3 coordNodes[10] = {};

    tets.loadNodes(i, tetNodes);
    coords.loadCoordinates(tetNodes, tetNodes + 10, coordNodes);

    // assume that the first 4 nodes are the same as in linear tetrahedra, i.e. corner nodes
    // Check the mid-point nodes positions are as expected w.r.t. to the corner nodes.
    // Not GMSH node ordering.

    // Node 4 should be the mid-point of edge 0-1
    Vec3 expectedPos = (coordNodes[0] + coordNodes[1]) * 0.5;
    REQUIRE(coordNodes[4].equals(expectedPos, 1e-9));

    // Node 5 should be the mid-point of edge 1-2
    expectedPos = (coordNodes[1] + coordNodes[2]) * 0.5;
    REQUIRE(coordNodes[5].equals(expectedPos, 1e-9));

    // Node 6 should be the mid-point of edge 0-2
    expectedPos = (coordNodes[2] + coordNodes[0]) * 0.5;
    REQUIRE(coordNodes[6].equals(expectedPos, 1e-9));

    // Node 7 should be the mid-point of edge 0-3
    expectedPos = (coordNodes[0] + coordNodes[3]) * 0.5;
    REQUIRE(coordNodes[7].equals(expectedPos, 1e-9));

    // Node 8 should be the mid-point of edge 1-3
    expectedPos = (coordNodes[1] + coordNodes[3]) * 0.5;
    REQUIRE(coordNodes[8].equals(expectedPos, 1e-9));

    // Node 9 should be the mid-point of edge 2-3
    expectedPos = (coordNodes[2] + coordNodes[3]) * 0.5;
    REQUIRE(coordNodes[9].equals(expectedPos, 1e-9));
  }
}

TEST_CASE("Quadratic tetra mid-edge nodes are snapped when within tolerance") {
  Geometry geom;
  Electrodes electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));

  auto meshText = readTextFile(quadraticMeshResource());
  meshText = replaceOnce(meshText, "0.4999999999986718 0 0", "0.5004 0 0");

  auto meshFile = TestUtil::TemporaryFile::withContents(meshText);
  RawMeshData meshData = MeshReader::readMesh(meshFile.name());
  const unsigned int midNodeIndex = meshData.tetNodes[4];
  const unsigned int cornerAIndex = meshData.tetNodes[0];
  const unsigned int cornerBIndex = meshData.tetNodes[1];
  const Vec3 expectedMidpoint = (meshData.points[cornerAIndex] + meshData.points[cornerBIndex]) * 0.5;

  prepareGeometry(geom, meshFile.name(), electrodes, alignment);

  REQUIRE(geom.getCoordinates().getPoint(midNodeIndex).equals(expectedMidpoint, 1e-12));
}

TEST_CASE("Quadratic tetra mid-edge nodes beyond tolerance fail validation") {
  Geometry geom;
  Electrodes electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));

  auto meshText = readTextFile(quadraticMeshResource());
  meshText = replaceOnce(meshText, "0.4999999999986718 0 0", "0.503 0 0");

  auto meshFile = TestUtil::TemporaryFile::withContents(meshText);

  REQUIRE_THROWS_AS(prepareGeometry(geom, meshFile.name(), electrodes, alignment), std::runtime_error);
}

TEST_CASE("Dielectric node reordering keeps LC nodes first") {
  Geometry geom;
  Electrodes electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));

  prepareGeometry(geom,
                  TestUtil::RESOURCE_UNIT_CUBE_DIELECTRIC_NEUMAN_QUADRATIC_GMSH_MESH,
                  electrodes,
                  alignment);

  const Mesh &tets = geom.getTetrahedra();
  REQUIRE(tets.getElementType() == ElementType::QUADRATIC_TETRAHEDRON);
  REQUIRE(geom.getnpLC() > 0);
  REQUIRE(geom.getnpLC() < geom.getnp());

  std::unordered_set<idx> lcNodes;
  std::unordered_set<idx> dielectricOnlyNodes;

  for (idx i = 0; i < tets.getnElements(); ++i) {
    idx nodeIndices[10] = {};
    tets.loadNodes(i, nodeIndices);

    if (tets.getMaterialNumber(i) <= MAT_DOMAIN7) {
      for (idx node : nodeIndices) {
        lcNodes.insert(node);
      }
    }
  }

  for (idx i = 0; i < tets.getnElements(); ++i) {
    if (tets.getMaterialNumber(i) <= MAT_DOMAIN7) {
      continue;
    }

    idx nodeIndices[10] = {};
    tets.loadNodes(i, nodeIndices);
    for (idx node : nodeIndices) {
      if (lcNodes.find(node) == lcNodes.end()) {
        dielectricOnlyNodes.insert(node);
      }
    }
  }

  REQUIRE_FALSE(lcNodes.empty());
  REQUIRE_FALSE(dielectricOnlyNodes.empty());
  REQUIRE(geom.getnpLC() == lcNodes.size());

  for (idx node : lcNodes) {
    REQUIRE(node < geom.getnpLC());
  }

  for (idx node : dielectricOnlyNodes) {
    REQUIRE(node >= geom.getnpLC());
  }
}

