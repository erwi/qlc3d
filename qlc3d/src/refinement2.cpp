#include <refinement.h>
#include <algorithm>
#include <array>
#include <geometry.h>
#include <line.h>
#include <set>
#include <stdexcept>
#include <vector>
#include <qlc3d.h>
#include <util/exception.h>
#include <util/logging.h>
#include <geom/coordinates.h>
#include <refinement/midpoint-node-lookup.h>
#include <refinement/tet-splitter.h>

using qlc3d::refinement::MapMidpointNodeLookup;
using qlc3d::refinement::SplitResult;
using qlc3d::refinement::Green2Nodes;
using qlc3d::refinement::Green3Nodes;
using qlc3d::refinement::resolveGreen2Nodes;
using qlc3d::refinement::resolveGreen3Nodes;
using qlc3d::refinement::splitGreen1Tet;
using qlc3d::refinement::splitGreen2aTet;
using qlc3d::refinement::splitGreen2bTet;
using qlc3d::refinement::splitGreen3Tet;
using qlc3d::refinement::splitRedTet;
using qlc3d::refinement::splitTri1;
using qlc3d::refinement::splitTri2;
using qlc3d::refinement::splitTri3;

namespace {

using TetNodes = std::array<idx, 4>;
using TriNodes = std::array<idx, 3>;

MapMidpointNodeLookup createNodeNumbersMatrix(Geometry &geom, const std::vector<Line> &lines) {
    /*!Creates midpoint-node lookup for new nodes*/

    MapMidpointNodeLookup lookup;
    const idx nold = geom.getnp();
    idx i = 0;
    for (const auto &l : lines) {
        lookup.registerEdge(l.L[0], l.L[1], nold + i);
        ++i;
    }
    return lookup;
}

std::array<std::array<idx, 2>, 1> collectLineNodes(const std::vector<Line> &lines, const std::set<idx> &lineIndices) {
    if (lineIndices.size() != 1) {
        RUNTIME_ERROR("Expected exactly one bisected line.");
    }

    std::array<std::array<idx, 2>, 1> edges{};
    edges[0] = {lines[*lineIndices.begin()].L[0], lines[*lineIndices.begin()].L[1]};
    return edges;
}

template <std::size_t N>
std::array<std::array<idx, 2>, N> collectLineNodes(const std::vector<Line> &lines, const std::set<idx> &lineIndices) {
    if (lineIndices.size() != N) {
        RUNTIME_ERROR("Unexpected number of bisected lines for refinement case.");
    }

    std::array<std::array<idx, 2>, N> edges{};
    std::size_t i = 0;
    for (const auto lineIndex : lineIndices) {
        edges[i++] = {lines[lineIndex].L[0], lines[lineIndex].L[1]};
    }
    return edges;
}

TetNodes getTetNodes(const Geometry &geom, idx elem) {
    const auto &tets = geom.getTetrahedra();
    return TetNodes{
        tets.getNode(elem, 0),
        tets.getNode(elem, 1),
        tets.getNode(elem, 2),
        tets.getNode(elem, 3)
    };
}

TriNodes getTriNodes(const Geometry &geom, idx elem) {
    const auto &tris = geom.getTriangles();
    return TriNodes{
        tris.getNode(elem, 0),
        tris.getNode(elem, 1),
        tris.getNode(elem, 2)
    };
}

TriNodes resolveTri2Nodes(const std::array<std::array<idx, 2>, 2> &lines) {
    std::array<idx, 4> nodes{
        lines[0][0], lines[0][1], lines[1][0], lines[1][1]
    };
    std::sort(nodes.begin(), nodes.end());

    auto rep = std::adjacent_find(nodes.begin(), nodes.end());
    if (rep == nodes.end()) {
        RUNTIME_ERROR("Unable to find shared node for triangle green-2 refinement.");
    }

    const idx sharedNode = *rep;
    std::array<idx, 2> remaining{};
    std::size_t remainingCount = 0;
    for (const auto node : nodes) {
        if (node != sharedNode) {
            remaining[remainingCount++] = node;
        }
    }

    if (remainingCount != 2) {
        RUNTIME_ERROR("Unable to resolve triangle green-2 non-shared nodes.");
    }

    std::sort(remaining.begin(), remaining.end());
    return TriNodes{sharedNode, remaining[0], remaining[1]};
}

template <typename Container>
void appendSplitResult(Container &newNodes, std::vector<idx> &newMaterials, const SplitResult &splitResult) {
    newNodes.insert(newNodes.end(), splitResult.nodes.begin(), splitResult.nodes.end());
    newMaterials.insert(newMaterials.end(), splitResult.materials.begin(), splitResult.materials.end());
}

} // namespace

void create_new_coordinates(Geometry& geom, const std::vector<Line>& lines, std::vector<Vec3>& new_p) {
    new_p.clear();
    new_p.reserve(lines.size() * 3); // reserve space for 3 coordinates per new node
    auto &coords = geom.getCoordinates();
    // LOOP OVER LINES, CALCULATE MID-POINT LOCATION AND ADD TO NEW COORDINATES
    for (unsigned int i = 0; i < lines.size(); i++) {
        auto p1 = coords.getPoint(lines[i].L[0]);
        auto p2 = coords.getPoint(lines[i].L[1]);
        new_p.push_back(Vec3::mean(p1, p2));
    }
}

void count_refinement_types(int& nred, int& ngreen1, int& ngreen2, int& ngreen3, std::vector<unsigned int>& i_tet) {
    /*! Counts the number of tets of each refinement type */
    nred = 0;
    ngreen1 = 0;
    ngreen2 = 0;
    ngreen3 = 0;
    for (size_t i = 0; i < i_tet.size(); i++) {
        if (i_tet[i] == RED_TET) nred++;
        if (i_tet[i] == GREEN1_TET) ngreen1++;
        if (i_tet[i] == GREEN2_TET) ngreen2++;
        if (i_tet[i] == GREEN3_TET) ngreen3++;
    }
}

/**
 * Makes sure the node ordering of tets in splitResult is such that the tet Jacobian determinant is positive
 * @param coords
 * @param splitResult
 */
void repairTetNodeOrder(const Coordinates &coords, SplitResult &splitResult) {
    const size_t numTets = splitResult.materials.size();

    for (size_t i = 0; i < numTets; i++) {
        auto p0 = coords.getPoint(splitResult.nodes[i * 4]);
        auto p1 = coords.getPoint(splitResult.nodes[i * 4 + 1]);
        auto p2 = coords.getPoint(splitResult.nodes[i * 4 + 2]);
        auto p3 = coords.getPoint(splitResult.nodes[i * 4 + 3]);
        if (det3D(p0, p1, p2, p3) < 0) {
            // swap last two nodes to change sign of determinant
            std::swap(splitResult.nodes[i * 4 + 2], splitResult.nodes[i * 4 + 3]);
        }
    }

}

void create_new_elements(Geometry& geom,
                         std::vector<idx>& i_tet,
                         std::vector<idx>& i_tri,
                         const std::vector<Line>& lines,
                         const std::vector<std::set<idx> >& t_to_l,
                         const std::vector<std::set<idx> >& e_to_l,
                         std::vector<Vec3>& new_p,
                         std::vector<idx>& new_t,
                         std::vector<idx>& new_mat_t,
                         std::vector<idx>& new_e,
                         std::vector<idx>& new_mat_e) {
    // CREATES NEW NODES, ELEMENTS (TETS + TRIS ) AND MATERIAL VALUES ARRAYS.

    MapMidpointNodeLookup nnumbers = createNodeNumbersMatrix(geom, lines);
    create_new_coordinates(geom, lines, new_p);
    geom.appendCoordinates(new_p);
    auto coords = geom.getCoordinates();

    auto &tets = geom.getTetrahedra();
    auto &tris = geom.getTriangles();

    for (idx i = 0; i < static_cast<idx>(i_tet.size()); i++) {
        const auto material = tets.getMaterialNumber(i);
        const auto tetNodes = getTetNodes(geom, i);

        switch (i_tet[i]) {
            case 0:
                break;
            case GREEN1_TET: {
                const auto lineIndex = *t_to_l[i].begin();
                const idx nAB = nnumbers.lookup(lines[lineIndex].L[0], lines[lineIndex].L[1]);
                auto res = splitGreen1Tet(tetNodes[0], tetNodes[1], tetNodes[2], tetNodes[3], nAB, material);
                repairTetNodeOrder(coords, res);
                appendSplitResult(new_t, new_mat_t, res);
                break;
            }
            case GREEN2_TET: {
                const auto edgeNodes = collectLineNodes<2>(lines, t_to_l[i]);
                const auto nodes = resolveGreen2Nodes(tetNodes, edgeNodes);
                if (nodes.isSharedNode) {
                    auto res = splitGreen2bTet(nodes.nA, nodes.nB, nodes.nC, nodes.nD,
                                              nnumbers.lookup(nodes.nA, nodes.nB),
                                              nnumbers.lookup(nodes.nA, nodes.nC),
                                              material);
                    repairTetNodeOrder(coords, res);
                    appendSplitResult(new_t, new_mat_t,res);
                } else {
                    auto res = splitGreen2aTet(nodes.nA, nodes.nB, nodes.nC, nodes.nD,
                                              nnumbers.lookup(nodes.nA, nodes.nB),
                                              nnumbers.lookup(nodes.nC, nodes.nD),
                                              material);
                    repairTetNodeOrder(coords, res);
                    appendSplitResult(new_t, new_mat_t, res);
                }
                break;
            }
            case GREEN3_TET: {
                const auto edgeNodes = collectLineNodes<3>(lines, t_to_l[i]);
                const auto nodes = resolveGreen3Nodes(tetNodes, edgeNodes);
                auto res = splitGreen3Tet(nodes.nA, nodes.nB, nodes.nC, nodes.nD,
                                                 nnumbers.lookup(nodes.nA, nodes.nB),
                                                 nnumbers.lookup(nodes.nA, nodes.nC),
                                                 nnumbers.lookup(nodes.nB, nodes.nC),
                                                 material);
                repairTetNodeOrder(coords, res);
                appendSplitResult(new_t, new_mat_t, res);
                break;
            }
            case RED_TET: {
                auto res = splitRedTet(tetNodes[0], tetNodes[1], tetNodes[2], tetNodes[3],
                                              nnumbers.lookup(tetNodes[0], tetNodes[1]),
                                              nnumbers.lookup(tetNodes[0], tetNodes[2]),
                                              nnumbers.lookup(tetNodes[0], tetNodes[3]),
                                              nnumbers.lookup(tetNodes[1], tetNodes[2]),
                                              nnumbers.lookup(tetNodes[1], tetNodes[3]),
                                              nnumbers.lookup(tetNodes[2], tetNodes[3]),
                                              material);
                repairTetNodeOrder(coords, res);
                appendSplitResult(new_t, new_mat_t, res);
                break;
            }
            default:
                Log::error("Unknown tetrahedron refinement type {} for element {}.", i_tet[i], i);
                RUNTIME_ERROR("Unexpected tetrahedron refinement type in create_new_elements.");
        }
    }

    for (idx i = 0; i < static_cast<idx>(i_tri.size()); i++) {
        const auto material = tris.getMaterialNumber(i);

        if (i_tri[i] == 0) {
            continue;
        }

        if (i_tri[i] == 1) {
            const auto lineIndex = *e_to_l[i].begin();
            std::vector<idx> triNodesVec{lines[lineIndex].L[0], lines[lineIndex].L[1]};
            tris.CompleteNodesSet(i, triNodesVec);
            if (triNodesVec.size() != 3) {
                Log::error("Triangle {} expected 3 nodes after completion but got {}.", i, triNodesVec.size());
                RUNTIME_ERROR("Unable to resolve triangle green-1 node ordering.");
            }
            appendSplitResult(new_e, new_mat_e,
                              splitTri1(triNodesVec[0], triNodesVec[1], triNodesVec[2],
                                        nnumbers.lookup(triNodesVec[0], triNodesVec[1]),
                                        material));
        } else if (i_tri[i] == 2) {
            const auto edgeNodes = collectLineNodes<2>(lines, e_to_l[i]);
            const auto nodes = resolveTri2Nodes(edgeNodes);
            appendSplitResult(new_e, new_mat_e,
                              splitTri2(nodes[0], nodes[1], nodes[2],
                                        nnumbers.lookup(nodes[0], nodes[1]),
                                        nnumbers.lookup(nodes[0], nodes[2]),
                                        material));
        } else if (i_tri[i] == 3) {
            const auto triNodes = getTriNodes(geom, i);
            const auto edgeNodes = collectLineNodes<3>(lines, e_to_l[i]);
            (void) edgeNodes;
            appendSplitResult(new_e, new_mat_e,
                              splitTri3(triNodes[0], triNodes[1], triNodes[2],
                                        nnumbers.lookup(triNodes[0], triNodes[1]),
                                        nnumbers.lookup(triNodes[0], triNodes[2]),
                                        nnumbers.lookup(triNodes[1], triNodes[2]),
                                        material));
        } else {
            Log::error("Unknown triangle refinement type {} for element {}.", i_tri[i], i);
            RUNTIME_ERROR("Unexpected triangle refinement type in create_new_elements.");
        }
    }
}


