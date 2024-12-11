#include <catch.h>
#include <dofmap.h>

TEST_CASE("Initialise DofMap") {
  DofMap dofMap(10, 5);

  for (unsigned int i = 0; i < 10; i++) {
    for (unsigned int j = 0; j < 5; j++) {
      REQUIRE(dofMap.getDof(i, j) == DofMap::NOT_DOF);
    }
  }
}

TEST_CASE("Calculate DofMap when there are no fixed or periodic nodes") {
  // ARRANGE
  unsigned int numNodes = 3;
  unsigned int numDimensions = 2;

  DofMap dofMap(numNodes, numDimensions);
  std::unordered_set<unsigned int> fixedNodes;
  std::vector<unsigned int> periodicNodesMapping = {0, 1, 2};

  // ACT
  dofMap.calculateMapping(fixedNodes, periodicNodesMapping);

  // ASSERT
  REQUIRE(dofMap.getnFreeNodes() == numNodes); // num free nodes per dimension
  REQUIRE(dofMap.getnDof() == numNodes);
  REQUIRE(dofMap.getnDimensions() == numDimensions);

  REQUIRE(dofMap.getDof(0, 0) == 0);
  REQUIRE(dofMap.getDof(1, 0) == 1);
  REQUIRE(dofMap.getDof(2, 0) == 2);

  REQUIRE(dofMap.getDof(0, 1) == 3);
  REQUIRE(dofMap.getDof(1, 1) == 4);
  REQUIRE(dofMap.getDof(2, 1) == 5);
}

TEST_CASE("Calculate DofMap when periodic nodes exist, but no fixed nodes") {
  // ARRANGE
  unsigned int numNodes = 3;
  unsigned int numDimensions = 2;

  DofMap dofMap(numNodes, numDimensions);
  std::unordered_set<unsigned int> fixedNodes;
  std::vector<unsigned int> periodicNodesMapping = {0, 0, 2};

  // ACT
  dofMap.calculateMapping(fixedNodes, periodicNodesMapping);

  // ASSERT
  REQUIRE(dofMap.getnFreeNodes() == 2); // num free nodes per dimension
  REQUIRE(dofMap.getnDimensions() == numDimensions);
  REQUIRE(dofMap.getnDof() == numNodes);

  REQUIRE(dofMap.getDof(0, 0) == 0);
  REQUIRE(dofMap.getDof(1, 0) == 0);
  REQUIRE(dofMap.getDof(2, 0) == 1);

  REQUIRE(dofMap.getDof(0, 1) == 2);
  REQUIRE(dofMap.getDof(1, 1) == 2);
  REQUIRE(dofMap.getDof(2, 1) == 3);
}

TEST_CASE("Calculate DofMap when fixed nodes exist, but no periodic nodes") {
  // ARRANGE
  unsigned int numNodes = 3;
  unsigned int numDimensions = 2;

  DofMap dofMap(numNodes, numDimensions);
  std::unordered_set<unsigned int> fixedNodes = {0, 1};
  std::vector<unsigned int> periodicNodesMapping = {0, 1, 2};

  // ACT
  dofMap.calculateMapping(fixedNodes, periodicNodesMapping);

  // ASSERT
  REQUIRE(dofMap.getnFreeNodes() == 1); // num free nodes per dimension
  REQUIRE(dofMap.getnDimensions() == numDimensions);
  REQUIRE(dofMap.getnDof() == numNodes);

  REQUIRE(dofMap.getDof(0, 0) == DofMap::NOT_DOF);
  REQUIRE(dofMap.getDof(1, 0) == DofMap::NOT_DOF);
  REQUIRE(dofMap.getDof(2, 0) == 0);

  REQUIRE(dofMap.getDof(0, 1) == DofMap::NOT_DOF);
  REQUIRE(dofMap.getDof(1, 1) == DofMap::NOT_DOF);
  REQUIRE(dofMap.getDof(2, 1) == 1);
}

TEST_CASE("Calculate DofMap when both fixed and periodic nodes exist") {
  // ARRANGE
  unsigned int numNodes = 4;
  unsigned int numDimensions = 2;

  DofMap dofMap(numNodes, numDimensions);
  std::unordered_set<unsigned int> fixedNodes = {1};
  std::vector<unsigned int> periodicNodesMapping = {0, 1, 0, 3};

  // ACT
  dofMap.calculateMapping(fixedNodes, periodicNodesMapping);

  // ASSERT
  REQUIRE(dofMap.getnFreeNodes() == 2);
  REQUIRE(dofMap.getnDimensions() == numDimensions);
  REQUIRE(dofMap.getnDof() == numNodes);

  REQUIRE(dofMap.getDof(0, 0) == 0);
  REQUIRE(dofMap.getDof(1, 0) == DofMap::NOT_DOF);
  REQUIRE(dofMap.getDof(2, 0) == 0);
  REQUIRE(dofMap.getDof(3, 0) == 1);

  REQUIRE(dofMap.getDof(0, 1) == 2);
  REQUIRE(dofMap.getDof(1, 1) == DofMap::NOT_DOF);
  REQUIRE(dofMap.getDof(2, 1) == 2);
  REQUIRE(dofMap.getDof(3, 1) == 3);
}
