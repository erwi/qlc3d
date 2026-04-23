#ifndef QLC3D_REFINEMENT_TET_SPLITTER_H
#define QLC3D_REFINEMENT_TET_SPLITTER_H

#include <globals.h>
#include <array>
#include <vector>

namespace qlc3d::refinement {

/**
 * @brief Result returned by a pure refinement splitter.
 *
 * The connectivity is stored in a flat array: 4 node indices per tetrahedron
 * and 3 node indices per triangle. The materials vector contains one entry per
 * child element.
 */
struct SplitResult {
    std::vector<::idx> nodes;
    std::vector<::idx> materials;
};

/**
 * @brief Resolved node ordering for a green-2 tetrahedron split.
 *
 * nA and nB are the endpoints of the first bisected edge, nC and nD are the
 * endpoints of the second bisected edge, and isSharedNode indicates whether
 * the two edges share one endpoint.
 */
struct Green2Nodes {
    ::idx nA;
    ::idx nB;
    ::idx nC;
    ::idx nD;
    bool isSharedNode;
};

/**
 * @brief Resolved node ordering for a green-3 tetrahedron split.
 *
 * nA, nB, and nC are the three nodes on the bisected face, while nD is the
 * apex node opposite that face.
 */
struct Green3Nodes {
    ::idx nA;
    ::idx nB;
    ::idx nC;
    ::idx nD;
};

/**
 * @brief Resolve the node ordering for a green-2 tetrahedron split.
 *
 * @param tetNodes The four corner nodes of the tetrahedron.
 * @param lines The two bisected edges in the tetrahedron.
 * @return The resolved node ordering and whether the two bisected edges share
 *         a node.
 */
Green2Nodes resolveGreen2Nodes(const std::array<::idx, 4> &tetNodes,
                               const std::array<std::array<::idx, 2>, 2> &lines);

/**
 * @brief Resolve the node ordering for a green-3 tetrahedron split.
 *
 * @param tetNodes The four corner nodes of the tetrahedron.
 * @param lines The three bisected edges that form the split face.
 * @return The three face nodes in sorted order and the apex node.
 */
Green3Nodes resolveGreen3Nodes(const std::array<::idx, 4> &tetNodes,
                               const std::array<std::array<::idx, 2>, 3> &lines);

/** Split a tet with 1 bisected edge into 2 child tets. */
SplitResult splitGreen1Tet(::idx nA, ::idx nB, ::idx nC, ::idx nD,
                           ::idx nAB,
                           ::idx material);

/** Split a tet with 2 bisected edges (disjoint: AB and CD) into 4 child tets. */
SplitResult splitGreen2aTet(::idx nA, ::idx nB, ::idx nC, ::idx nD,
                            ::idx nAB, ::idx nCD,
                            ::idx material);

/** Split a tet with 2 bisected edges sharing node A (AB and AC) into 3 child tets. */
SplitResult splitGreen2bTet(::idx nA, ::idx nB, ::idx nC, ::idx nD,
                            ::idx nAB, ::idx nAC,
                            ::idx material);

/** Split a tet with 3 bisected edges forming a loop (AB, AC, BC) into 4 child tets. */
SplitResult splitGreen3Tet(::idx nA, ::idx nB, ::idx nC, ::idx nD,
                           ::idx nAB, ::idx nAC, ::idx nBC,
                           ::idx material);

/** Split a red tet (all 6 edges bisected) into 8 child tets. */
SplitResult splitRedTet(::idx nA, ::idx nB, ::idx nC, ::idx nD,
                        ::idx nAB, ::idx nAC, ::idx nAD,
                        ::idx nBC, ::idx nBD, ::idx nCD,
                        ::idx material);

/** Split a triangle with 1 bisected edge into 2 child triangles. */
SplitResult splitTri1(::idx nA, ::idx nB, ::idx nC, ::idx nAB, ::idx material);

/** Split a triangle with 2 bisected edges sharing node A into 3 child triangles. */
SplitResult splitTri2(::idx nA, ::idx nB, ::idx nC, ::idx nAB, ::idx nAC, ::idx material);

/** Split a triangle with all 3 edges bisected into 4 child triangles. */
SplitResult splitTri3(::idx nA, ::idx nB, ::idx nC, ::idx nAB, ::idx nAC, ::idx nBC, ::idx material);

} // namespace qlc3d::refinement

#endif // QLC3D_REFINEMENT_TET_SPLITTER_H


