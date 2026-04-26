#include <geom/tet-mesh-search.h>
#include <mesh/mesh.h>
#include <geom/coordinates.h>
#include <geom/aabox.h>
#include <geom/vec3.h>
#include <material_numbers.h>
#include <util/exception.h>
#include <util/logging.h>
#include <globals.h>

#include <algorithm>
#include <cfloat>
#include <limits>
#include <set>
#include <vector>

const unsigned int TetMeshSearch::NOT_AN_INDEX = std::numeric_limits<unsigned int>::max();

TetMeshSearch::TetMeshSearch(const Mesh& tets, const Coordinates& coords, const AABox& bounds)
  : tets_(tets), coords_(coords), bounds_(bounds) {}

bool TetMeshSearch::bruteForceSearch(unsigned int& ind,
                                      const Vec3& crd,
                                      bool terminateOnError,
                                      bool requireLCElement) const {
  for (idx i = 0; i < tets_.getnElements(); i++) {
    bool found = tets_.containsCoordinate(i, coords_, crd);
    if (found) {
      if (requireLCElement) {
        if (tets_.getMaterialNumber(i) <= MAT_DOMAIN7) {
          ind = i;
          return true;
        }
      } else {
        ind = i;
        return true;
      }
    }
  }

  if (terminateOnError) {
    RUNTIME_ERROR(fmt::format("Brute force search could not find coordinate at ({})", crd));
  }
  return false;
}

size_t TetMeshSearch::recursiveNeighbourSearch(const Vec3& target,
                                                const std::vector<std::set<unsigned int>>& p_to_t,
                                                size_t currentTet,
                                                std::set<size_t>& history,
                                                bool requireLCElement) const {
  bool found = tets_.containsCoordinate(currentTet, coords_, target);
  if (found) {
    if (!requireLCElement) {
      return currentTet;
    } else {
      if (tets_.getMaterialNumber(currentTet) <= MAT_DOMAIN7) {
        return currentTet;
      }
    }
  }

  history.insert(currentTet);
  // Recursion limit to avoid stack overflow on large grids / concave meshes
  static const idx RECURSION_LIMIT = 1000;
  if ((idx)history.size() > RECURSION_LIMIT) return NOT_AN_INDEX;

  // Collect all neighbour elements by visiting every node of the current element
  std::vector<unsigned int> neighs;
  for (idx i = 0; i < tets_.getnNodes(); i++) {
    int n = tets_.getNode(currentTet, i);
    neighs.insert(neighs.end(), p_to_t[n].begin(), p_to_t[n].end());
  }
  std::sort(neighs.begin(), neighs.end());
  auto itr = std::unique(neighs.begin(), neighs.end());
  neighs.erase(itr, neighs.end());

  // Build distance-to-centroid list for all neighbours
  std::vector<double> dists;
  dists.reserve(neighs.size());
  for (size_t i = 0; i < neighs.size(); i++) {
    Vec3 centroid = tets_.elementCentroid(neighs[i], coords_);
    dists.push_back(centroid.distanceSquared(target));
  }

  size_t indn = std::min_element(dists.begin(), dists.end()) - dists.begin();

  while (dists[indn] < DBL_MAX) {
    size_t indt = neighs[indn];
    if (history.find(indt) == history.end()) {
      size_t indexFound = recursiveNeighbourSearch(target, p_to_t, indt, history, requireLCElement);
      if (indexFound != NOT_AN_INDEX) {
        return indexFound;
      }
    }
    dists[indn] = DBL_MAX;
    indn = std::min_element(dists.begin(), dists.end()) - dists.begin();
  }

  return NOT_AN_INDEX;
}

void TetMeshSearch::genIndToTetsByCoords(std::vector<unsigned int>& returnIndex,
                                          const Coordinates& targets,
                                          bool terminateOnError,
                                          bool requireLCElement) const {
  if (tets_.getnElements() == 0) {
    RUNTIME_ERROR("No tetrahedra elements defined");
  }

  returnIndex.clear();
  unsigned int nt = (unsigned int)tets_.getnElements();
  returnIndex.assign(targets.size(), nt);

  std::vector<std::set<unsigned int>> p_to_t;
  tets_.gen_p_to_elem(p_to_t);

  // Find most central tet as the starting point for all searches
  Vec3 structureCentroid = bounds_.centre();
  std::set<size_t> initHistory;
  unsigned int midTet = (unsigned int)recursiveNeighbourSearch(structureCentroid, p_to_t, 0, initHistory);
  if (midTet == NOT_AN_INDEX) {
    midTet = 0; // Fall back to element 0 if centre not found (concave or hollow mesh)
  }

  for (unsigned int n = 0; n < targets.size(); n++) {
    Vec3 targetPoint = targets.getPoint(n);
    std::set<size_t> searchHistory;

    size_t t0 = recursiveNeighbourSearch(targetPoint, p_to_t, midTet, searchHistory, requireLCElement);
    if (t0 != NOT_AN_INDEX) {
      returnIndex[n] = (unsigned int)t0;
    } else {
      unsigned int tetIndex = 0;
      if (bruteForceSearch(tetIndex, targetPoint, terminateOnError, requireLCElement)) {
        returnIndex[n] = tetIndex;
      } else {
        returnIndex[n] = NOT_AN_INDEX;
        Log::info("Could not find regular grid point {} at ({}) in volume mesh. "
                  "Assuming it is outside the mesh and continuing.", n, targetPoint);
      }
    }
  }
}

