#include <catch.h>
#include <refinement/element-selector.h>
#include <refinement/refinement-spec.h>
#include <geometry.h>
#include <solutionvector.h>
#include <inits.h>
#include <electrodes.h>
#include <alignment.h>
#include <refinement.h>
#include <test-util.h>
#include <geom/coordinates.h>

using namespace qlc3d::refinement;

// ── Helpers ──────────────────────────────────────────────────────────────────

namespace {

Alignment makeAlignment() {
    Alignment alignment;
    alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
    alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
    return alignment;
}

Electrodes makeElectrodes() {
    std::vector<std::shared_ptr<Electrode>> ev;
    ev.emplace_back(std::make_shared<Electrode>(1, std::vector<double>{0}, std::vector<double>{0}));
    ev.emplace_back(std::make_shared<Electrode>(2, std::vector<double>{0}, std::vector<double>{0}));
    return Electrodes::withElectrodePotentials(ev);
}

void loadSmallCube(Geometry& geom) {
    Alignment alignment = makeAlignment();
    Electrodes electrodes = makeElectrodes();
    prepareGeometry(geom, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH, electrodes, alignment, {1, 1, 1});
}

} // namespace

// ── fromSpec: factory dispatch ────────────────────────────────────────────────

TEST_CASE("ElementSelector::fromSpec returns ChangeSelector for Change type") {
    auto spec = RefinementSpec::makeExplicit("change", 1, 0, {0.05}, {}, {}, {});
    auto sel = ElementSelector::fromSpec(*spec);
    REQUIRE(sel != nullptr);
    REQUIRE(sel->needsQTensor() == true);
}

TEST_CASE("ElementSelector::fromSpec returns SphereSelector for Sphere type") {
    auto spec = RefinementSpec::makeExplicit("sphere", 1, 0, {0.3}, {0.5}, {0.5}, {0.5});
    auto sel = ElementSelector::fromSpec(*spec);
    REQUIRE(sel != nullptr);
    REQUIRE(sel->needsQTensor() == false);
}

TEST_CASE("ElementSelector::fromSpec returns BoxSelector for Box type") {
    auto spec = RefinementSpec::makeExplicit("box", 1, 0, {}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0});
    auto sel = ElementSelector::fromSpec(*spec);
    REQUIRE(sel != nullptr);
    REQUIRE(sel->needsQTensor() == false);
}

// ── SphereSelector ────────────────────────────────────────────────────────────

TEST_CASE("SphereSelector: sphere covering entire mesh marks all LC tets") {
    // The small cube mesh fits in [0,1]^3; a sphere of radius 2 centred at
    // (0.5, 0.5, 0.5) encloses all nodes.
    auto spec = RefinementSpec::makeExplicit("sphere", 1, 0, {2.0}, {0.5}, {0.5}, {0.5});
    auto sel = ElementSelector::makeSphere(*spec);

    Geometry geom; loadSmallCube(geom);
    idx numTets = geom.getTetrahedra().getnElements();
    std::vector<idx> i_tet(numTets, 0);

    sel->selectTets(geom, nullptr, 0, i_tet);

    long redCount = std::count(i_tet.begin(), i_tet.end(), static_cast<idx>(RED_TET));
    REQUIRE(redCount > 0);
    // All LC tets should be red when the sphere covers everything
    long lcTets = 0;
    for (idx i = 0; i < numTets; ++i) {
        if (geom.getTetrahedra().getMaterialNumber(i) <= MAT_DOMAIN7) ++lcTets;
    }
    REQUIRE(redCount == lcTets);
}

TEST_CASE("SphereSelector: tiny sphere outside mesh marks no tets") {
    // Sphere is outside the [0,1]^3 domain
    auto spec = RefinementSpec::makeExplicit("sphere", 1, 0, {0.1}, {5.0}, {5.0}, {5.0});
    auto sel = ElementSelector::makeSphere(*spec);

    Geometry geom; loadSmallCube(geom);
    idx numTets = geom.getTetrahedra().getnElements();
    std::vector<idx> i_tet(numTets, 0);

    sel->selectTets(geom, nullptr, 0, i_tet);

    long redCount = std::count(i_tet.begin(), i_tet.end(), static_cast<idx>(RED_TET));
    REQUIRE(redCount == 0);
}

// ── BoxSelector ───────────────────────────────────────────────────────────────

TEST_CASE("BoxSelector: box covering entire mesh marks all LC tets") {
    // Box covers the full [0,1]^3 domain
    auto spec = RefinementSpec::makeExplicit("box", 1, 0, {}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0});
    auto sel = ElementSelector::makeBox(*spec);

    Geometry geom; loadSmallCube(geom);
    idx numTets = geom.getTetrahedra().getnElements();
    std::vector<idx> i_tet(numTets, 0);

    sel->selectTets(geom, nullptr, 0, i_tet);

    long redCount = std::count(i_tet.begin(), i_tet.end(), static_cast<idx>(RED_TET));
    long lcTets = 0;
    for (idx i = 0; i < numTets; ++i) {
        if (geom.getTetrahedra().getMaterialNumber(i) <= MAT_DOMAIN7) ++lcTets;
    }
    REQUIRE(redCount == lcTets);
}

TEST_CASE("BoxSelector: box outside mesh marks no tets") {
    auto spec = RefinementSpec::makeExplicit("box", 1, 0, {}, {5.0, 6.0}, {5.0, 6.0}, {5.0, 6.0});
    auto sel = ElementSelector::makeBox(*spec);

    Geometry geom; loadSmallCube(geom);
    idx numTets = geom.getTetrahedra().getnElements();
    std::vector<idx> i_tet(numTets, 0);

    sel->selectTets(geom, nullptr, 0, i_tet);

    long redCount = std::count(i_tet.begin(), i_tet.end(), static_cast<idx>(RED_TET));
    REQUIRE(redCount == 0);
}

TEST_CASE("BoxSelector: partial box selects corner-adjacent tets but not all") {
    // Select only nodes at x=0 face (x in [0, 0]). The mesh has 8 corner nodes,
    // 4 of which are on the x=0 face. Tets touching those nodes should be selected
    // but tets on the x=1 face only should not.
    // Use a thin slab at x=0 with a small positive tolerance.
    auto spec = RefinementSpec::makeExplicit("box", 1, 0, {}, {0.0, 0.01}, {0.0, 1.0}, {0.0, 1.0});
    auto sel = ElementSelector::makeBox(*spec);

    Geometry geom; loadSmallCube(geom);
    idx numTets = geom.getTetrahedra().getnElements();
    std::vector<idx> i_tet(numTets, 0);

    sel->selectTets(geom, nullptr, 0, i_tet);

    long redCount = std::count(i_tet.begin(), i_tet.end(), static_cast<idx>(RED_TET));
    long lcTets = 0;
    for (idx i = 0; i < numTets; ++i) {
        if (geom.getTetrahedra().getMaterialNumber(i) <= MAT_DOMAIN7) ++lcTets;
    }
    REQUIRE(redCount > 0);
    REQUIRE(redCount < lcTets);
}

// ── ChangeSelector ────────────────────────────────────────────────────────────

TEST_CASE("ChangeSelector: throws when q is null") {
    auto spec = RefinementSpec::makeExplicit("change", 1, 0, {0.01}, {}, {}, {});
    auto sel = ElementSelector::makeChange(*spec);

    Geometry geom; loadSmallCube(geom);
    std::vector<idx> i_tet(geom.getTetrahedra().getnElements(), 0);

    REQUIRE_THROWS_AS(sel->selectTets(geom, nullptr, 0, i_tet), std::invalid_argument);
}

TEST_CASE("ChangeSelector: uniform q marks no tets") {
    // All Q values the same → maxDeltaQ == 0 < any positive threshold
    auto spec = RefinementSpec::makeExplicit("change", 1, 0, {0.01}, {}, {}, {});
    auto sel = ElementSelector::makeChange(*spec);

    Geometry geom; loadSmallCube(geom);
    idx npLC = geom.getnpLC();
    SolutionVector q(npLC, 5);
    for (idx n = 0; n < npLC; ++n)
        for (int d = 0; d < 5; ++d)
            q.setValue(n, d, 0.5);

    std::vector<idx> i_tet(geom.getTetrahedra().getnElements(), 0);
    sel->selectTets(geom, &q, 0, i_tet);

    long redCount = std::count(i_tet.begin(), i_tet.end(), static_cast<idx>(RED_TET));
    REQUIRE(redCount == 0);
}

TEST_CASE("ChangeSelector: large variation marks tets") {
    // Very low threshold so that any variation triggers selection
    auto spec = RefinementSpec::makeExplicit("change", 1, 0, {1e-6}, {}, {}, {});
    auto sel = ElementSelector::makeChange(*spec);

    Geometry geom; loadSmallCube(geom);
    idx npLC = geom.getnpLC();
    SolutionVector q(npLC, 5);
    // Assign q[0] linearly varying across nodes so neighbouring nodes differ
    auto& coords = geom.getCoordinates();
    for (idx n = 0; n < npLC; ++n) {
        auto p = coords.getPoint(n);
        q.setValue(n, 0, p.x() + p.y() + p.z()); // varies from ~0 to ~3
        for (int d = 1; d < 5; ++d) q.setValue(n, d, 0.5);
    }

    std::vector<idx> i_tet(geom.getTetrahedra().getnElements(), 0);
    sel->selectTets(geom, &q, 0, i_tet);

    long redCount = std::count(i_tet.begin(), i_tet.end(), static_cast<idx>(RED_TET));
    REQUIRE(redCount > 0);
}

TEST_CASE("ChangeSelector: very high threshold marks no tets") {
    auto spec = RefinementSpec::makeExplicit("change", 1, 0, {1e10}, {}, {}, {});
    auto sel = ElementSelector::makeChange(*spec);

    Geometry geom; loadSmallCube(geom);
    idx npLC = geom.getnpLC();
    SolutionVector q(npLC, 5);
    auto& coords = geom.getCoordinates();
    for (idx n = 0; n < npLC; ++n) {
        auto p = coords.getPoint(n);
        q.setValue(n, 0, p.x());
        for (int d = 1; d < 5; ++d) q.setValue(n, d, 0.5);
    }

    std::vector<idx> i_tet(geom.getTetrahedra().getnElements(), 0);
    sel->selectTets(geom, &q, 0, i_tet);

    long redCount = std::count(i_tet.begin(), i_tet.end(), static_cast<idx>(RED_TET));
    REQUIRE(redCount == 0);
}




