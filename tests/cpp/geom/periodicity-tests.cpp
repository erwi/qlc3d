#include <catch.h>
#include <io/gmsh-read.h>
#include <test-util.h>
#include <io/meshreader.h>
#include <mesh.h>
#include <geom/coordinates.h>
#include <geom/periodicity.h>

TEST_CASE("Periodicity mapping") {

  auto meshData = MeshReader::readMesh(TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH);

  Mesh triMesh(2, 3);
  triMesh.setElementData(std::move(meshData.triNodes), std::move(meshData.triMaterials));

  Mesh tetMesh(3, 4);
  tetMesh.setElementData(std::move(meshData.tetNodes), std::move(meshData.tetMaterials));




  Coordinates coordinates(std::move(meshData.points));
  triMesh.setConnectedVolume(&tetMesh);

  triMesh.calculateSurfaceNormals(coordinates, &tetMesh);

  auto periodicityType = PeriodicityType(triMesh);

  PeriodicNodesMapping pnm(triMesh, coordinates, periodicityType);

  pnm.initialisePeriodicNodes(triMesh, coordinates);


}