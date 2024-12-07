#include <fixednodes.h>
#include <mesh.h>
#include <util/exception.h>
#include <util/logging.h>
#include <material_numbers.h>

void FixedNodes::setFixedNodesPot(const Mesh &triangles, const std::unordered_map<unsigned int, double> &potentialByElectrode) {
  fixedValueByNodeIndex.clear();

  // iterate over each potentialByElectrode
  for (const auto& [electrodeNumber, potential] : potentialByElectrode) {
    Log::info("Setting Electrode{} nodes to potential {}", electrodeNumber, potential);

    std::unordered_set<unsigned int> nodeIndices = triangles.findElectrodeSurfaceNodes(electrodeNumber);
    for (auto &i :nodeIndices) {
      fixedValueByNodeIndex[i] = potential;
    }
  }

  Log::info("{} fixed nodes set", fixedValueByNodeIndex.size());
}

void FixedNodes::setFixedNodesAndValues(const std::unordered_map<unsigned int, double> &fixedNodesAndValues) {
  fixedValueByNodeIndex = fixedNodesAndValues;
}