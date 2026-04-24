#include <catch.h>
#include <refinement/tet-classifier.h>
#include <geometry.h>
#include <inits.h>
#include <electrodes.h>
#include <alignment.h>
#include <refinement.h>
#include <test-util.h>
#include <algorithm>
#include <vector>
#include <set>

using namespace qlc3d::refinement;

// Helper: build a tet Mesh from a flat list of node indices (4 per tet) and materials
static std::shared_ptr<Mesh> makeTetMesh(const std::vector<idx>& nodes,
                                         const std::vector<idx>& mats) {
    auto m = Mesh::tetMesh();
    m->setElementData(ElementType::LINEAR_TETRAHEDRON,
                      std::vector<unsigned int>(nodes.begin(), nodes.end()),
                      std::vector<unsigned int>(mats.begin(), mats.end()));
    return m;
}

// ===== Phase 1: collectRedTetEdges =====

TEST_CASE("collectRedTetEdges: no red tets returns empty list") {
    auto mesh = makeTetMesh({0,1,2,3}, {1});
    std::vector<idx> i_tet = {0}; // not red
    auto lines = collectRedTetEdges(*mesh, i_tet);
    REQUIRE(lines.empty());
}

TEST_CASE("collectRedTetEdges: single red tet produces 6 edges") {
    auto mesh = makeTetMesh({0,1,2,3}, {1});
    std::vector<idx> i_tet = {RED_TET};
    auto lines = collectRedTetEdges(*mesh, i_tet);
    REQUIRE(lines.size() == 6);
}

TEST_CASE("collectRedTetEdges: two adjacent red tets sharing a face deduplicate shared edges") {
    // tet0: nodes 0,1,2,3; tet1: nodes 0,1,2,4
    // shared face 0-1-2 contributes edges 0-1, 0-2, 1-2 (3 shared edges)
    // total unique: 6 + 6 - 3 = 9
    auto mesh = makeTetMesh({0,1,2,3, 0,1,2,4}, {1,1});
    std::vector<idx> i_tet = {RED_TET, RED_TET};
    auto lines = collectRedTetEdges(*mesh, i_tet);
    REQUIRE(lines.size() == 9);
}

TEST_CASE("collectRedTetEdges: two fully disjoint red tets produce 12 edges") {
    auto mesh = makeTetMesh({0,1,2,3, 4,5,6,7}, {1,1});
    std::vector<idx> i_tet = {RED_TET, RED_TET};
    auto lines = collectRedTetEdges(*mesh, i_tet);
    REQUIRE(lines.size() == 12);
}

TEST_CASE("collectRedTetEdges: non-red tets are ignored") {
    auto mesh = makeTetMesh({0,1,2,3, 4,5,6,7}, {1,1});
    std::vector<idx> i_tet = {GREEN1_TET, GREEN2_TET};
    auto lines = collectRedTetEdges(*mesh, i_tet);
    REQUIRE(lines.empty());
}

TEST_CASE("collectRedTetEdges: result is sorted") {
    auto mesh = makeTetMesh({3,2,1,0, 4,5,6,7}, {1,1});
    std::vector<idx> i_tet = {RED_TET, RED_TET};
    auto lines = collectRedTetEdges(*mesh, i_tet);
    REQUIRE(std::is_sorted(lines.begin(), lines.end()));
}

// ===== Phase 2: countEdgesPerElement =====

TEST_CASE("countEdgesPerElement: empty lines produces all-zero counts") {
    auto mesh = makeTetMesh({0,1,2,3}, {1});
    std::vector<std::set<idx>> p_to_elem;
    mesh->gen_p_to_elem(p_to_elem);
    std::vector<Line> lines;
    std::vector<idx> i_elem;
    std::vector<std::set<idx>> m_to_l;
    countEdgesPerElement(*mesh, p_to_elem, lines, i_elem, m_to_l);
    REQUIRE(i_elem.size() == 1);
    REQUIRE(i_elem[0] == 0);
}

TEST_CASE("countEdgesPerElement: single line shared by two tets increments both") {
    // two tets sharing edge 0-1
    auto mesh = makeTetMesh({0,1,2,3, 0,1,4,5}, {1,1});
    std::vector<std::set<idx>> p_to_elem;
    mesh->gen_p_to_elem(p_to_elem);
    std::vector<Line> lines = {Line(0,1)};
    std::vector<idx> i_elem;
    std::vector<std::set<idx>> m_to_l;
    countEdgesPerElement(*mesh, p_to_elem, lines, i_elem, m_to_l);
    REQUIRE(i_elem.size() == 2);
    REQUIRE(i_elem[0] == 1);
    REQUIRE(i_elem[1] == 1);
}

TEST_CASE("countEdgesPerElement: edge not in any element produces no increments") {
    auto mesh = makeTetMesh({0,1,2,3}, {1});
    std::vector<std::set<idx>> p_to_elem;
    mesh->gen_p_to_elem(p_to_elem);
    // Node 99 doesn't exist in mesh
    std::vector<Line> lines = {Line(0,99)};
    std::vector<idx> i_elem;
    std::vector<std::set<idx>> m_to_l;
    countEdgesPerElement(*mesh, p_to_elem, lines, i_elem, m_to_l);
    REQUIRE(i_elem[0] == 0);
}

TEST_CASE("countEdgesPerElement: m_to_l contains correct line index") {
    auto mesh = makeTetMesh({0,1,2,3}, {1});
    std::vector<std::set<idx>> p_to_elem;
    mesh->gen_p_to_elem(p_to_elem);
    std::vector<Line> lines = {Line(0,1)};
    std::vector<idx> i_elem;
    std::vector<std::set<idx>> m_to_l;
    countEdgesPerElement(*mesh, p_to_elem, lines, i_elem, m_to_l);
    REQUIRE(m_to_l[0].count(0) == 1);
}

TEST_CASE("countEdgesPerElement: single tet with all 6 of its own edges scores 6") {
    auto mesh = makeTetMesh({0,1,2,3}, {1});
    std::vector<std::set<idx>> p_to_elem;
    mesh->gen_p_to_elem(p_to_elem);
    std::vector<Line> lines = {Line(0,1),Line(0,2),Line(0,3),Line(1,2),Line(1,3),Line(2,3)};
    std::vector<idx> i_elem;
    std::vector<std::set<idx>> m_to_l;
    countEdgesPerElement(*mesh, p_to_elem, lines, i_elem, m_to_l);
    REQUIRE(i_elem[0] == 6);
}

TEST_CASE("countEdgesPerElement: output vectors are correctly sized") {
    auto mesh = makeTetMesh({0,1,2,3, 4,5,6,7, 8,9,10,11}, {1,1,1});
    std::vector<std::set<idx>> p_to_elem;
    mesh->gen_p_to_elem(p_to_elem);
    std::vector<Line> lines = {Line(0,1)};
    std::vector<idx> i_elem;
    std::vector<std::set<idx>> m_to_l;
    countEdgesPerElement(*mesh, p_to_elem, lines, i_elem, m_to_l);
    REQUIRE(i_elem.size() == 3);
    REQUIRE(m_to_l.size() == 3);
}

// ===== Phase 3: assignRefinementTypes =====

TEST_CASE("assignRefinementTypes: all zeros stays all zeros") {
    std::vector<idx> i_tet = {0, 0, 0};
    auto nrt = assignRefinementTypes(i_tet);
    REQUIRE(nrt.red == 0);
    REQUIRE(nrt.green1 == 0);
    REQUIRE(nrt.green2 == 0);
    REQUIRE(nrt.green3 == 0);
    for (auto v : i_tet) REQUIRE(v == 0);
}

TEST_CASE("assignRefinementTypes: counts each type correctly") {
    std::vector<idx> i_tet = {GREEN1_TET, GREEN1_TET, GREEN2_TET, GREEN3_TET, RED_TET, RED_TET, RED_TET};
    auto nrt = assignRefinementTypes(i_tet);
    REQUIRE(nrt.green1 == 2);
    REQUIRE(nrt.green2 == 1);
    REQUIRE(nrt.green3 == 1);
    REQUIRE(nrt.red == 3);
}

TEST_CASE("assignRefinementTypes: value 6 stays 6") {
    std::vector<idx> i_tet = {6};
    auto nrt = assignRefinementTypes(i_tet);
    REQUIRE(i_tet[0] == RED_TET);
    REQUIRE(nrt.red == 1);
}

TEST_CASE("assignRefinementTypes: returned counts match updated i_tet") {
    std::vector<idx> i_tet = {0, GREEN1_TET, RED_TET, GREEN3_TET};
    auto nrt = assignRefinementTypes(i_tet);
    idx countRed = std::count(i_tet.begin(), i_tet.end(), (idx)RED_TET);
    idx countG1 = std::count(i_tet.begin(), i_tet.end(), (idx)GREEN1_TET);
    idx countG3 = std::count(i_tet.begin(), i_tet.end(), (idx)GREEN3_TET);
    REQUIRE((idx)nrt.red == countRed);
    REQUIRE((idx)nrt.green1 == countG1);
    REQUIRE((idx)nrt.green3 == countG3);
}

// ===== Phase 3: resolveGreen3RedAmbiguity =====

TEST_CASE("resolveGreen3RedAmbiguity: 3-unique-node green3 stays green3") {
    // 3 lines forming a closed face loop: 0-1, 1-2, 0-2 -> nodes {0,1,2}
    std::vector<Line> lines = {Line(0,1), Line(0,2), Line(1,2)};
    std::vector<idx> i_tet = {GREEN3_TET};
    std::vector<std::set<idx>> t_to_l = {{0,1,2}};
    resolveGreen3RedAmbiguity(i_tet, lines, t_to_l);
    REQUIRE(i_tet[0] == GREEN3_TET);
}

TEST_CASE("resolveGreen3RedAmbiguity: 4-unique-node green3 promoted to red") {
    // 3 lines with 4 unique nodes: 0-1, 2-3, 0-2 -> nodes {0,1,2,3}
    std::vector<Line> lines = {Line(0,1), Line(0,2), Line(2,3)};
    std::vector<idx> i_tet = {GREEN3_TET};
    std::vector<std::set<idx>> t_to_l = {{0,1,2}};
    resolveGreen3RedAmbiguity(i_tet, lines, t_to_l);
    REQUIRE(i_tet[0] == RED_TET);
}

TEST_CASE("resolveGreen3RedAmbiguity: non-green3 elements are not touched") {
    std::vector<Line> lines = {Line(0,1)};
    std::vector<idx> i_tet = {GREEN1_TET, RED_TET};
    std::vector<std::set<idx>> t_to_l = {{0}, {0}};
    resolveGreen3RedAmbiguity(i_tet, lines, t_to_l);
    REQUIRE(i_tet[0] == GREEN1_TET);
    REQUIRE(i_tet[1] == RED_TET);
}

TEST_CASE("resolveGreen3RedAmbiguity: throws on unexpected node count") {
    // 3 lines all sharing the same two endpoints -> 2 unique nodes
    std::vector<Line> lines = {Line(0,1), Line(0,1), Line(0,1)};
    std::vector<idx> i_tet = {GREEN3_TET};
    std::vector<std::set<idx>> t_to_l = {{0,1,2}};
    REQUIRE_THROWS(resolveGreen3RedAmbiguity(i_tet, lines, t_to_l));
}

// ===== Phase 4: classifyRefinement integration tests =====

namespace {
    std::filesystem::path meshResourcePath(const char* filename) {
        return std::filesystem::absolute(std::filesystem::path(__FILE__))
            .parent_path().parent_path().parent_path() / "resources" / filename;
    }

    Alignment makeAlignmentForClassifier() {
        Alignment alignment;
        alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
        alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
        return alignment;
    }
    Electrodes makeElectrodesForClassifier() {
        std::vector<std::shared_ptr<Electrode>> ev;
        ev.emplace_back(std::make_shared<Electrode>(1, std::vector<double>{0}, std::vector<double>{0}));
        ev.emplace_back(std::make_shared<Electrode>(2, std::vector<double>{0}, std::vector<double>{0}));
        return Electrodes::withElectrodePotentials(ev);
    }
    void loadSmallCube(Geometry& geom) {
        Alignment alignment = makeAlignmentForClassifier();
        Electrodes electrodes = makeElectrodesForClassifier();
        prepareGeometry(geom, meshResourcePath("gmsh-small-cube.msh"), electrodes, alignment, {1, 1, 1});
    }
}

TEST_CASE("classifyRefinement: no red tets returns all-zero classification") {
    Geometry geom; loadSmallCube(geom);
    idx nTets = geom.getTetrahedra().getnElements();
    std::vector<idx> i_tet(nTets, 0);
    auto result = classifyRefinement(geom, i_tet);
    REQUIRE(result.i_tet.size() == (size_t)nTets);
    REQUIRE(std::all_of(result.i_tet.begin(), result.i_tet.end(), [](idx v){ return v == 0; }));
    REQUIRE(result.lines.empty());
}

TEST_CASE("classifyRefinement: i_tet size equals number of tets in mesh") {
    Geometry geom; loadSmallCube(geom);
    idx nTets = geom.getTetrahedra().getnElements();
    std::vector<idx> i_tet(nTets, 0);
    i_tet[0] = RED_TET;
    auto result = classifyRefinement(geom, i_tet);
    REQUIRE(result.i_tet.size() == (size_t)nTets);
}

TEST_CASE("classifyRefinement: i_tri size equals number of triangles in mesh") {
    Geometry geom; loadSmallCube(geom);
    idx nTets = geom.getTetrahedra().getnElements();
    idx nTris = geom.getTriangles().getnElements();
    std::vector<idx> i_tet(nTets, 0);
    i_tet[0] = RED_TET;
    auto result = classifyRefinement(geom, i_tet);
    REQUIRE(result.i_tri.size() == (size_t)nTris);
}

TEST_CASE("classifyRefinement: all lines in result reference valid node indices") {
    Geometry geom; loadSmallCube(geom);
    idx nTets = geom.getTetrahedra().getnElements();
    idx nNodes = geom.getnp();
    std::vector<idx> i_tet(nTets, 0);
    i_tet[0] = RED_TET;
    auto result = classifyRefinement(geom, i_tet);
    for (const auto& line : result.lines) {
        REQUIRE(line.L[0] < nNodes);
        REQUIRE(line.L[1] < nNodes);
        REQUIRE(line.L[0] < line.L[1]);
    }
}

TEST_CASE("classifyRefinement: t_to_l indices are within bounds of lines") {
    Geometry geom; loadSmallCube(geom);
    idx nTets = geom.getTetrahedra().getnElements();
    std::vector<idx> i_tet(nTets, 0);
    i_tet[0] = RED_TET;
    auto result = classifyRefinement(geom, i_tet);
    idx nLines = (idx) result.lines.size();
    for (const auto& s : result.t_to_l) {
        for (idx li : s) {
            REQUIRE(li < nLines);
        }
    }
}

TEST_CASE("classifyRefinement: type codes are in valid set after classification") {
    Geometry geom; loadSmallCube(geom);
    idx nTets = geom.getTetrahedra().getnElements();
    std::vector<idx> i_tet(nTets, 0);
    i_tet[0] = RED_TET;
    auto result = classifyRefinement(geom, i_tet);
    for (idx v : result.i_tet) {
        REQUIRE((v == 0 || v == GREEN1_TET || v == GREEN2_TET || v == GREEN3_TET || v == RED_TET));
    }
}

