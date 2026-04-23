#include <catch.h>
#include <refinement/tet-splitter.h>

#include <array>
#include <algorithm>
#include <set>
#include <vector>

using qlc3d::refinement::SplitResult;
using qlc3d::refinement::Green2Nodes;
using qlc3d::refinement::Green3Nodes;
using qlc3d::refinement::splitGreen1Tet;
using qlc3d::refinement::splitGreen2aTet;
using qlc3d::refinement::splitGreen2bTet;
using qlc3d::refinement::splitGreen3Tet;
using qlc3d::refinement::splitRedTet;
using qlc3d::refinement::splitTri1;
using qlc3d::refinement::splitTri2;
using qlc3d::refinement::splitTri3;
using qlc3d::refinement::resolveGreen2Nodes;
using qlc3d::refinement::resolveGreen3Nodes;

namespace {

void requireMaterialsPreserved(const SplitResult &result, idx material, std::size_t expectedCount) {
    REQUIRE(result.materials.size() == expectedCount);
    REQUIRE(std::all_of(result.materials.begin(), result.materials.end(), [material](idx m) {
        return m == material;
    }));
}

void requireNodesEqual(const SplitResult &result, const std::vector<idx> &expectedNodes) {
    REQUIRE(result.nodes == expectedNodes);
}

void requireCornersAppear(const SplitResult &result, std::initializer_list<idx> corners) {
    std::set<idx> nodeSet(result.nodes.begin(), result.nodes.end());
    for (idx corner : corners) {
        REQUIRE(nodeSet.count(corner) > 0);
    }
}

std::size_t countOccurrences(const std::vector<idx> &nodes, idx needle) {
    return static_cast<std::size_t>(std::count(nodes.begin(), nodes.end(), needle));
}

std::size_t countEndpointOccurrences(const std::array<std::array<idx, 2>, 3> &lines, idx needle) {
    std::size_t count = 0;
    for (const auto &line : lines) {
        count += static_cast<std::size_t>(line[0] == needle);
        count += static_cast<std::size_t>(line[1] == needle);
    }
    return count;
}

} // namespace

TEST_CASE("resolveGreen2Nodes: two disjoint edges produce isSharedNode == false") {
    const std::array<idx, 4> tetNodes{10, 11, 12, 13};
    const std::array<std::array<idx, 2>, 2> lines{{{10, 11}, {12, 13}}};

    const Green2Nodes nodes = resolveGreen2Nodes(tetNodes, lines);

    REQUIRE_FALSE(nodes.isSharedNode);
    REQUIRE(nodes.nA == 10);
    REQUIRE(nodes.nB == 11);
    REQUIRE(nodes.nC == 12);
    REQUIRE(nodes.nD == 13);
}

TEST_CASE("resolveGreen2Nodes: two edges sharing a node produce isSharedNode == true") {
    const std::array<idx, 4> tetNodes{10, 11, 12, 13};
    const std::array<std::array<idx, 2>, 2> lines{{{10, 11}, {10, 12}}};

    const Green2Nodes nodes = resolveGreen2Nodes(tetNodes, lines);

    REQUIRE(nodes.isSharedNode);
    REQUIRE(nodes.nA == 10);
    REQUIRE(nodes.nB == 11);
    REQUIRE(nodes.nC == 12);
    REQUIRE(nodes.nD == 13);
}

TEST_CASE("resolveGreen3Nodes: returns 3 unique loop nodes and 1 apex") {
    const std::array<idx, 4> tetNodes{10, 11, 12, 13};
    const std::array<std::array<idx, 2>, 3> lines{{{11, 12}, {10, 12}, {10, 11}}};

    const Green3Nodes nodes = resolveGreen3Nodes(tetNodes, lines);

    REQUIRE(nodes.nA == 10);
    REQUIRE(nodes.nB == 11);
    REQUIRE(nodes.nC == 12);
    REQUIRE(nodes.nD == 13);
    REQUIRE(countEndpointOccurrences(lines, nodes.nA) == 2);
    REQUIRE(countEndpointOccurrences(lines, nodes.nB) == 2);
    REQUIRE(countEndpointOccurrences(lines, nodes.nC) == 2);
}

TEST_CASE("splitGreen1Tet produces two child tets") {
    const idx material = 7;
    auto result = splitGreen1Tet(10, 11, 12, 13, 101, material);

    REQUIRE(result.nodes.size() == 8);
    requireMaterialsPreserved(result, material, 2);
    requireCornersAppear(result, {10, 11, 12, 13});
    REQUIRE(countOccurrences(result.nodes, 101) == 2);
    requireNodesEqual(result, {10, 12, 13, 101,
                               101, 12, 13, 11});
}

TEST_CASE("splitGreen2aTet produces four child tets") {
    const idx material = 9;
    auto result = splitGreen2aTet(10, 11, 12, 13, 101, 102, material);

    REQUIRE(result.nodes.size() == 16);
    requireMaterialsPreserved(result, material, 4);
    requireCornersAppear(result, {10, 11, 12, 13});
    REQUIRE(countOccurrences(result.nodes, 101) == 4);
    REQUIRE(countOccurrences(result.nodes, 102) == 4);
    requireNodesEqual(result, {10, 101, 102, 12,
                               10, 101, 13, 102,
                               101, 11, 102, 12,
                               101, 11, 13, 102});
}

TEST_CASE("splitGreen2bTet produces three child tets") {
    const idx material = 5;
    auto result = splitGreen2bTet(10, 11, 12, 13, 101, 102, material);

    REQUIRE(result.nodes.size() == 12);
    requireMaterialsPreserved(result, material, 3);
    requireCornersAppear(result, {10, 11, 12, 13});
    REQUIRE(countOccurrences(result.nodes, 101) == 3);
    REQUIRE(countOccurrences(result.nodes, 102) == 2);
    requireNodesEqual(result, {10, 101, 13, 102,
                               102, 101, 13, 12,
                               101, 11, 13, 12});
}

TEST_CASE("splitGreen3Tet produces four child tets") {
    const idx material = 8;
    auto result = splitGreen3Tet(10, 11, 12, 13, 101, 102, 103, material);

    REQUIRE(result.nodes.size() == 16);
    requireMaterialsPreserved(result, material, 4);
    requireCornersAppear(result, {10, 11, 12, 13});
    REQUIRE(countOccurrences(result.nodes, 13) == 4);
    REQUIRE(countOccurrences(result.nodes, 101) == 3);
    REQUIRE(countOccurrences(result.nodes, 102) == 3);
    REQUIRE(countOccurrences(result.nodes, 103) == 3);
    requireNodesEqual(result, {10, 13, 101, 102,
                               11, 13, 101, 103,
                               12, 13, 102, 103,
                               13, 101, 102, 103});
}

TEST_CASE("splitRedTet produces eight child tets") {
    const idx material = 11;
    auto result = splitRedTet(10, 11, 12, 13, 101, 102, 103, 104, 105, 106, material);

    REQUIRE(result.nodes.size() == 32);
    requireMaterialsPreserved(result, material, 8);
    requireCornersAppear(result, {10, 11, 12, 13});
    REQUIRE(countOccurrences(result.nodes, 102) == 6);
    REQUIRE(countOccurrences(result.nodes, 104) == 4);
    REQUIRE(countOccurrences(result.nodes, 105) == 6);
    REQUIRE(countOccurrences(result.nodes, 106) == 4);
    requireNodesEqual(result, {10, 101, 102, 103,
                               11, 104, 101, 105,
                               12, 102, 104, 106,
                               13, 103, 106, 105,
                               101, 102, 103, 105,
                               101, 102, 105, 104,
                               102, 103, 105, 106,
                               102, 104, 106, 105});
}

TEST_CASE("splitTri1 produces two child triangles") {
    const idx material = 3;
    auto result = splitTri1(1, 2, 3, 101, material);

    REQUIRE(result.nodes.size() == 6);
    requireMaterialsPreserved(result, material, 2);
    REQUIRE(countOccurrences(result.nodes, 101) == 2);
    requireNodesEqual(result, {1, 101, 3,
                               2, 3, 101});
}

TEST_CASE("splitTri2 produces three child triangles") {
    const idx material = 4;
    auto result = splitTri2(1, 2, 3, 101, 102, material);

    REQUIRE(result.nodes.size() == 9);
    requireMaterialsPreserved(result, material, 3);
    REQUIRE(countOccurrences(result.nodes, 101) == 3);
    REQUIRE(countOccurrences(result.nodes, 102) == 2);
    requireNodesEqual(result, {1, 101, 102,
                               102, 101, 3,
                               101, 2, 3});
}

TEST_CASE("splitTri3 produces four child triangles") {
    const idx material = 6;
    auto result = splitTri3(1, 2, 3, 101, 102, 103, material);

    REQUIRE(result.nodes.size() == 12);
    requireMaterialsPreserved(result, material, 4);
    REQUIRE(countOccurrences(result.nodes, 101) == 3);
    REQUIRE(countOccurrences(result.nodes, 102) == 3);
    REQUIRE(countOccurrences(result.nodes, 103) == 3);
    requireNodesEqual(result, {1, 101, 102,
                               101, 2, 103,
                               102, 103, 3,
                               102, 101, 103});
}




