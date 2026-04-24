#ifndef QLC3D_REFINEMENT_TET_CLASSIFIER_H
#define QLC3D_REFINEMENT_TET_CLASSIFIER_H

#include <globals.h>
#include <line.h>
#include <geometry.h>
#include <refinement.h>
#include <vector>
#include <set>

namespace qlc3d::refinement {

/**
 * @brief Result of a full refinement classification pass.
 *
 * i_tet[i] contains the refinement type code for element i:
 *   0           — not refined
 *   GREEN1_TET  — one bisected edge  → 2 child tets
 *   GREEN2_TET  — two bisected edges → 3 or 4 child tets
 *   GREEN3_TET  — three bisected edges forming a face loop → 4 child tets
 *   RED_TET     — all six edges bisected → 8 child tets
 *
 * lines is the deduplicated, sorted list of all bisected edges.
 * t_to_l[i] maps element i to the set of indices into lines that it shares.
 */
struct ClassificationResult {
    std::vector<idx>            i_tet;   ///< Refinement type per tet
    std::vector<idx>            i_tri;   ///< Refinement type per triangle
    std::vector<Line>           lines;   ///< Deduplicated bisected edges
    std::vector<std::set<idx>>  t_to_l;  ///< Tet → bisected-edge index
    std::vector<std::set<idx>>  e_to_l;  ///< Triangle → bisected-edge index
};

/**
 * @brief Collect all edges of the red tetrahedra into a deduplicated sorted list.
 *
 * @param tets   The tetrahedral mesh.
 * @param i_tet  Refinement type per element; only RED_TET entries are processed.
 * @return Sorted, unique vector of Lines covering all red-tet edges.
 */
std::vector<Line> collectRedTetEdges(const Mesh& tets,
                                     const std::vector<idx>& i_tet);

/**
 * @brief Count bisected edges per element and build the element→edge index.
 *
 * @param elements     Mesh (tets or triangles).
 * @param p_to_elem    Node-to-element index (built from elements).
 * @param lines        Bisected edge list.
 * @param i_elem       Per-element bisected-edge count (output, resized to nElements).
 * @param m_to_l       Per-element set of indices into lines (output, resized to nElements).
 */
void countEdgesPerElement(const Mesh& elements,
                          const std::vector<std::set<idx>>& p_to_elem,
                          const std::vector<Line>& lines,
                          std::vector<idx>& i_elem,
                          std::vector<std::set<idx>>& m_to_l);

/**
 * @brief Map per-element edge counts to refinement type codes.
 *
 * Values ≥ 4 are promoted to RED_TET (6). Value 0 elements are unchanged.
 *
 * @param i_tet  In/out: updated in-place to refinement type codes.
 * @return       Counts of each refinement type.
 */
Num_Ref_Tet assignRefinementTypes(std::vector<idx>& i_tet);

/**
 * @brief Resolve the green-3 vs. red ambiguity for tets with 3 bisected edges.
 *
 * A tet with 3 bisected edges and 3 unique nodes is green-3 (closed face loop).
 * If it has 4 unique nodes it is actually red and is promoted to RED_TET.
 *
 * @param i_tet   In/out: mis-classified green-3 entries promoted to RED_TET.
 * @param lines   Bisected edge list.
 * @param t_to_l  Element→edge index.
 */
void resolveGreen3RedAmbiguity(std::vector<idx>& i_tet,
                               const std::vector<Line>& lines,
                               const std::vector<std::set<idx>>& t_to_l);

/**
 * @brief Classify all tetrahedra given an initial set of red-tet indices.
 *
 * Runs the iterative expansion loop until no new red tets are produced,
 * then classifies surface triangles. Periodic boundary expansion is
 * delegated to the existing expand_periodic_boundaries helper.
 *
 * @param geom     The geometry (tets and triangles read; coordinates used
 *                 only by the periodic expansion helper).
 * @param i_tet_in Initial classification vector: RED_TET for selected tets,
 *                 0 elsewhere. Must be the same size as the number of tets.
 * @return         Fully expanded ClassificationResult.
 */
ClassificationResult classifyRefinement(Geometry& geom,
                                        std::vector<idx> i_tet_in);

} // namespace qlc3d::refinement

#endif // QLC3D_REFINEMENT_TET_CLASSIFIER_H


