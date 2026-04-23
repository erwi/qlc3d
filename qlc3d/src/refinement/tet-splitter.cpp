#include <refinement/tet-splitter.h>

#include <algorithm>
#include <initializer_list>
#include <set>
#include <stdexcept>

namespace qlc3d::refinement {
namespace {

::idx findTetApex(const std::array<::idx, 4> &tetNodes, const std::set<::idx> &faceNodes) {
    for (const auto node : tetNodes) {
        if (!faceNodes.count(node)) {
            return node;
        }
    }
    throw std::runtime_error("Unable to resolve tetrahedron apex node.");
}

std::set<::idx> edgeNodeSet(const std::array<std::array<::idx, 2>, 2> &lines) {
    std::set<::idx> nodes;
    for (const auto &line : lines) {
        nodes.insert(line[0]);
        nodes.insert(line[1]);
    }
    return nodes;
}

std::set<::idx> edgeNodeSet(const std::array<std::array<::idx, 2>, 3> &lines) {
    std::set<::idx> nodes;
    for (const auto &line : lines) {
        nodes.insert(line[0]);
        nodes.insert(line[1]);
    }
    return nodes;
}

SplitResult makeTetSplit(std::initializer_list<::idx> nodes, ::idx material, std::size_t childCount) {
    SplitResult result;
    result.nodes.assign(nodes.begin(), nodes.end());
    result.materials.assign(childCount, material);
    return result;
}

} // namespace

Green2Nodes resolveGreen2Nodes(const std::array<::idx, 4> &tetNodes,
                               const std::array<std::array<::idx, 2>, 2> &lines) {
    const auto uniqueNodes = edgeNodeSet(lines);
    if (uniqueNodes.size() == 4) {
        return Green2Nodes{
            lines[0][0], lines[0][1], lines[1][0], lines[1][1], false
        };
    }

    if (uniqueNodes.size() != 3) {
        throw std::runtime_error("green2 refinement must have either 3 or 4 unique line nodes.");
    }

    std::array<::idx, 4> nodes{
        lines[0][0], lines[0][1], lines[1][0], lines[1][1]
    };
    std::sort(nodes.begin(), nodes.end());

    auto rep = std::adjacent_find(nodes.begin(), nodes.end());
    if (rep == nodes.end()) {
        throw std::runtime_error("Unable to find shared node for green2 refinement.");
    }

    const ::idx sharedNode = *rep;
    std::array<::idx, 2> remaining{};
    std::size_t remainingCount = 0;
    for (const auto node : nodes) {
        if (node != sharedNode) {
            remaining[remainingCount++] = node;
        }
    }

    if (remainingCount != 2) {
        throw std::runtime_error("Unable to resolve non-shared nodes for green2 refinement.");
    }

    std::sort(remaining.begin(), remaining.end());
    const std::set<::idx> faceNodes{sharedNode, remaining[0], remaining[1]};
    return Green2Nodes{sharedNode, remaining[0], remaining[1], findTetApex(tetNodes, faceNodes), true};
}

Green3Nodes resolveGreen3Nodes(const std::array<::idx, 4> &tetNodes,
                               const std::array<std::array<::idx, 2>, 3> &lines) {
    const auto faceNodes = edgeNodeSet(lines);
    if (faceNodes.size() != 3) {
        throw std::runtime_error("green3 refinement must have exactly 3 unique face nodes.");
    }

    std::array<::idx, 3> nodes{};
    std::size_t i = 0;
    for (const auto node : faceNodes) {
        nodes[i++] = node;
    }
    if (i != 3) {
        throw std::runtime_error("Unable to resolve green3 face nodes.");
    }

    return Green3Nodes{nodes[0], nodes[1], nodes[2], findTetApex(tetNodes, faceNodes)};
}

SplitResult splitGreen1Tet(::idx nA, ::idx nB, ::idx nC, ::idx nD,
                           ::idx nAB,
                           ::idx material) {
    return makeTetSplit({nA, nC, nD, nAB,
                         nAB, nC, nD, nB}, material, 2);
}

SplitResult splitGreen2aTet(::idx nA, ::idx nB, ::idx nC, ::idx nD,
                            ::idx nAB, ::idx nCD,
                            ::idx material) {
    return makeTetSplit({nA, nAB, nCD, nC,
                         nA, nAB, nD, nCD,
                         nAB, nB, nCD, nC,
                         nAB, nB, nD, nCD}, material, 4);
}

SplitResult splitGreen2bTet(::idx nA, ::idx nB, ::idx nC, ::idx nD,
                            ::idx nAB, ::idx nAC,
                            ::idx material) {
    return makeTetSplit({nA, nAB, nD, nAC,
                         nAC, nAB, nD, nC,
                         nAB, nB, nD, nC}, material, 3);
}

SplitResult splitGreen3Tet(::idx nA, ::idx nB, ::idx nC, ::idx nD,
                           ::idx nAB, ::idx nAC, ::idx nBC,
                           ::idx material) {
    return makeTetSplit({nA, nD, nAB, nAC,
                         nB, nD, nAB, nBC,
                         nC, nD, nAC, nBC,
                         nD, nAB, nAC, nBC}, material, 4);
}

SplitResult splitRedTet(::idx nA, ::idx nB, ::idx nC, ::idx nD,
                        ::idx nAB, ::idx nAC, ::idx nAD,
                        ::idx nBC, ::idx nBD, ::idx nCD,
                        ::idx material) {
    return makeTetSplit({nA, nAB, nAC, nAD,
                         nB, nBC, nAB, nBD,
                         nC, nAC, nBC, nCD,
                         nD, nAD, nCD, nBD,
                         nAB, nAC, nAD, nBD,
                         nAB, nAC, nBD, nBC,
                         nAC, nAD, nBD, nCD,
                         nAC, nBC, nCD, nBD}, material, 8);
}

SplitResult splitTri1(::idx nA, ::idx nB, ::idx nC, ::idx nAB, ::idx material) {
    return makeTetSplit({nA, nAB, nC,
                         nB, nC, nAB}, material, 2);
}

SplitResult splitTri2(::idx nA, ::idx nB, ::idx nC, ::idx nAB, ::idx nAC, ::idx material) {
    return makeTetSplit({nA, nAB, nAC,
                         nAC, nAB, nC,
                         nAB, nB, nC}, material, 3);
}

SplitResult splitTri3(::idx nA, ::idx nB, ::idx nC, ::idx nAB, ::idx nAC, ::idx nBC, ::idx material) {
    return makeTetSplit({nA, nAB, nAC,
                         nAB, nB, nBC,
                         nAC, nBC, nC,
                         nAC, nAB, nBC}, material, 4);
}

} // namespace qlc3d::refinement



