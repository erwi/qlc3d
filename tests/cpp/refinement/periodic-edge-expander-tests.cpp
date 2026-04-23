#include <catch.h>
#include <refinement/periodic-edge-expander.h>
#include <geometry.h>
#include <inits.h>
#include <electrodes.h>
#include <alignment.h>
#include <test-util.h>
#include <line.h>
#include <geom/periodicity.h>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

namespace {

Alignment makeAlignment() {
    Alignment a;
    a.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
    a.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
    return a;
}

Electrodes makeElectrodes() {
    std::vector<std::shared_ptr<Electrode>> ev;
    ev.emplace_back(std::make_shared<Electrode>(1, std::vector<double>{0}, std::vector<double>{0}));
    ev.emplace_back(std::make_shared<Electrode>(2, std::vector<double>{0}, std::vector<double>{0}));
    return Electrodes::withElectrodePotentials(ev);
}

void loadFrontBackPeriodicGeom(Geometry& geom) {
    // pseudo-2d-neumann.msh: 1×0.1×1 with front/back periodic boundaries
    auto alignment = makeAlignment();
    auto electrodes = makeElectrodes();
    prepareGeometry(geom, TestUtil::RESOURCE_PSEUDO_2D_NEUMANN_GMSH_MESH, electrodes, alignment, {1, 1, 1});
}

void loadNonPeriodicGeom(Geometry& geom) {
    // unit-cube-neumann.msh has left/right/front/back Neumann surfaces and no periodic boundaries
    auto alignment = makeAlignment();
    auto electrodes = makeElectrodes();
    prepareGeometry(geom, TestUtil::RESOURCE_UNIT_CUBE_NEUMANN_GMSH_MESH, electrodes, alignment, {1, 1, 1});
}

} // anonymous namespace

// ===========================================================================
// WP3 – Catalogue tests
// ===========================================================================

TEST_CASE("PeriodicEdgeExpander: non-periodic geometry produces empty catalogues") {
    Geometry geom;
    loadNonPeriodicGeom(geom);

    PeriodicEdgeExpander expander(geom);

    REQUIRE(expander.frontEdges().empty());
    REQUIRE(expander.backEdges().empty());
    REQUIRE(expander.leftEdges().empty());
    REQUIRE(expander.rightEdges().empty());
    REQUIRE(expander.topEdges().empty());
    REQUIRE(expander.bottomEdges().empty());
}

TEST_CASE("PeriodicEdgeExpander: front/back mesh populates front and back catalogues") {
    Geometry geom;
    loadFrontBackPeriodicGeom(geom);

    PeriodicityType pt(geom.getTriangles());
    REQUIRE(pt.isFrontBackPeriodic());
    REQUIRE_FALSE(pt.isLeftRightPeriodic());

    PeriodicEdgeExpander expander(geom);

    // Both front and back edge lists must be non-empty
    REQUIRE_FALSE(expander.frontEdges().empty());
    REQUIRE_FALSE(expander.backEdges().empty());

    // Left/right/top/bottom must remain empty for a front/back-only mesh
    REQUIRE(expander.leftEdges().empty());
    REQUIRE(expander.rightEdges().empty());
}

TEST_CASE("PeriodicEdgeExpander: catalogue edges are deduplicated") {
    Geometry geom;
    loadFrontBackPeriodicGeom(geom);

    PeriodicEdgeExpander expander(geom);

    // Check that there are no duplicate entries in front and back catalogues
    auto checkUnique = [](const std::vector<Line>& v) {
        std::vector<Line> copy = v;
        std::sort(copy.begin(), copy.end());
        auto it = std::unique(copy.begin(), copy.end());
        return it == copy.end();
    };

    REQUIRE(checkUnique(expander.frontEdges()));
    REQUIRE(checkUnique(expander.backEdges()));
}

// ===========================================================================
// WP4 – findTranslation tests
// ===========================================================================
// findTranslation is private/static so we test its effect indirectly through
// expand(). We also expose specific scenarios through expand on a carefully
// controlled catalogue.
//
// For direct testing we build a small geometry and call expand() with
// hand-crafted edges.

TEST_CASE("PeriodicEdgeExpander: expand on non-periodic geometry is a no-op") {
    Geometry geom;
    loadNonPeriodicGeom(geom);

    PeriodicEdgeExpander expander(geom);

    // No edges at all – should stay empty
    std::vector<Line> lines;
    expander.expand(lines);
    REQUIRE(lines.empty());
}

TEST_CASE("PeriodicEdgeExpander: expand adds back mirror for a front edge") {
    Geometry geom;
    loadFrontBackPeriodicGeom(geom);

    PeriodicEdgeExpander expander(geom);

    // Take the first front edge as our bisected edge
    REQUIRE_FALSE(expander.frontEdges().empty());
    Line frontEdge = expander.frontEdges().front();

    std::vector<Line> lines = {frontEdge};
    expander.expand(lines);

    // After expansion the list must contain more than just the front edge
    REQUIRE(lines.size() > 1);

    // The mirror (back edge) must be in the back catalogue
    const auto& backCat = expander.backEdges();
    bool foundInBack = false;
    for (const Line& l : lines) {
        if (l == frontEdge) continue;
        if (std::find(backCat.begin(), backCat.end(), l) != backCat.end()) {
            foundInBack = true;
            break;
        }
    }
    REQUIRE(foundInBack);
}

TEST_CASE("PeriodicEdgeExpander: expand on already-present mirror produces no duplicates") {
    Geometry geom;
    loadFrontBackPeriodicGeom(geom);

    PeriodicEdgeExpander expander(geom);

    REQUIRE_FALSE(expander.frontEdges().empty());
    Line frontEdge = expander.frontEdges().front();

    std::vector<Line> lines = {frontEdge};
    expander.expand(lines);
    std::size_t sizeAfterFirst = lines.size();

    // Call expand again – should not add more edges
    expander.expand(lines);
    REQUIRE(lines.size() == sizeAfterFirst);
}

TEST_CASE("PeriodicEdgeExpander: expand all front catalogue edges succeeds") {
    Geometry geom;
    loadFrontBackPeriodicGeom(geom);

    PeriodicEdgeExpander expander(geom);

    // Use a copy of all front edges as input to expand
    std::vector<Line> lines = expander.frontEdges();
    REQUIRE_NOTHROW(expander.expand(lines));

    // Result must include at least the back mirrors too
    const auto& backCat = expander.backEdges();
    for (const Line& bl : backCat) {
        REQUIRE(std::find(lines.begin(), lines.end(), bl) != lines.end());
    }
}

TEST_CASE("PeriodicEdgeExpander: expand all back catalogue edges succeeds") {
    Geometry geom;
    loadFrontBackPeriodicGeom(geom);

    PeriodicEdgeExpander expander(geom);

    std::vector<Line> lines = expander.backEdges();
    REQUIRE_NOTHROW(expander.expand(lines));

    const auto& frontCat = expander.frontEdges();
    for (const Line& fl : frontCat) {
        REQUIRE(std::find(lines.begin(), lines.end(), fl) != lines.end());
    }
}


