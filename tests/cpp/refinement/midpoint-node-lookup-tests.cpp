#include <catch.h>
#include <refinement/midpoint-node-lookup.h>

using qlc3d::refinement::MapMidpointNodeLookup;

TEST_CASE("lookup returns correct node for a registered edge") {
    MapMidpointNodeLookup lookup;
    lookup.registerEdge(2, 5, 12);

    REQUIRE(lookup.lookup(2, 5) == 12);
}

TEST_CASE("lookup returns same node for reversed edge") {
    MapMidpointNodeLookup lookup;
    lookup.registerEdge(2, 5, 12);

    REQUIRE(lookup.lookup(5, 2) == 12);
}

TEST_CASE("lookup throws for an unknown edge") {
    MapMidpointNodeLookup lookup;
    lookup.registerEdge(2, 5, 12);

    REQUIRE_THROWS_AS(lookup.lookup(1, 2), std::out_of_range);
}


