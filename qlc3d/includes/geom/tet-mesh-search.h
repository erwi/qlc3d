#pragma once
#include <vector>
#include <set>

class Mesh;
class Coordinates;
class AABox;
class Vec3;

/**
 * Spatial search service for tetrahedral meshes.
 *
 * Given a tetrahedral mesh, coordinates, and bounding box, provides methods to
 * locate which tetrahedron contains a given point. Works correctly for both
 * LINEAR_TETRAHEDRON (TET4) and QUADRATIC_TETRAHEDRON (TET10) element types,
 * since only the four corner-node coordinates are used for geometry operations.
 */
class TetMeshSearch {
public:
  /** Sentinel value indicating that no tetrahedron was found for a query point. */
  static const unsigned int NOT_AN_INDEX;

  /**
   * @param tets   The tetrahedral mesh to search.
   * @param coords Node coordinates for the mesh.
   * @param bounds Precomputed axis-aligned bounding box of the mesh.
   */
  TetMeshSearch(const Mesh& tets, const Coordinates& coords, const AABox& bounds);

  /**
   * For each point in @p targets find the index of the tetrahedron that contains it,
   * storing results in @p returnIndex. Points not found within any tetrahedron are
   * stored as @p NOT_AN_INDEX.
   *
   * @param returnIndex      Output: element index per target point, or NOT_AN_INDEX.
   * @param targets          Query coordinates.
   * @param terminateOnError Throw if a point cannot be located in any element.
   * @param requireLCElement Only LC-material elements (material <= MAT_DOMAIN7) may be returned.
   */
  void genIndToTetsByCoords(std::vector<unsigned int>& returnIndex,
                             const Coordinates& targets,
                             bool terminateOnError = true,
                             bool requireLCElement = false) const;

  /**
   * Single-point brute-force linear search over all elements.
   * @param ind              Output: element index when found.
   * @param crd              The coordinate to locate.
   * @param terminateOnError Throw if the point is not found.
   * @param requireLCElement Only LC-material elements may be returned.
   * @return true if the point was found; false otherwise (only when terminateOnError == false).
   */
  bool bruteForceSearch(unsigned int& ind,
                        const Vec3& crd,
                        bool terminateOnError = true,
                        bool requireLCElement = false) const;

private:
  const Mesh& tets_;
  const Coordinates& coords_;
  const AABox& bounds_;

  /**
   * Greedy nearest-neighbour walk to find the tetrahedron containing @p target.
   * Starts from @p currentTet and follows neighbours whose centroid is closest
   * to @p target.
   *
   * @param target         The coordinate to locate.
   * @param p_to_t         Node-to-element adjacency list.
   * @param currentTet     Index of the element to start from.
   * @param history        Set of already-visited element indices (prevents infinite loops).
   * @param requireLCElement Only LC-material elements may be returned.
   * @return Index of the containing element, or NOT_AN_INDEX if not found.
   */
  size_t recursiveNeighbourSearch(const Vec3& target,
                                   const std::vector<std::set<unsigned int>>& p_to_t,
                                   size_t currentTet,
                                   std::set<size_t>& history,
                                   bool requireLCElement = false) const;
};

