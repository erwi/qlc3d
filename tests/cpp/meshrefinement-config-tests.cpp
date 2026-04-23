#include <catch.h>
#include <meshrefinement.h>

TEST_CASE("RefinementConfig validation") {

    SECTION("Change type with non-empty values constructs without throwing") {
        REQUIRE_NOTHROW(RefinementConfig("Change", {}, {}, {0.1}, {}, {}, {}));
    }

    SECTION("Change type with empty values throws") {
        REQUIRE_THROWS_AS(RefinementConfig("Change", {}, {}, {}, {}, {}, {}), std::invalid_argument);
    }

    SECTION("Sphere type with valid values, x, y, z constructs without throwing") {
        REQUIRE_NOTHROW(RefinementConfig("Sphere", {}, {}, {0.5}, {0.0}, {0.0}, {0.0}));
    }

    SECTION("Sphere type missing values throws") {
        REQUIRE_THROWS_AS(RefinementConfig("Sphere", {}, {}, {}, {0.0}, {0.0}, {0.0}), std::invalid_argument);
    }

    SECTION("Sphere type missing x throws") {
        REQUIRE_THROWS_AS(RefinementConfig("Sphere", {}, {}, {0.5}, {}, {0.0}, {0.0}), std::invalid_argument);
    }

    SECTION("Sphere type missing y throws") {
        REQUIRE_THROWS_AS(RefinementConfig("Sphere", {}, {}, {0.5}, {0.0}, {}, {0.0}), std::invalid_argument);
    }

    SECTION("Sphere type missing z throws") {
        REQUIRE_THROWS_AS(RefinementConfig("Sphere", {}, {}, {0.5}, {0.0}, {0.0}, {}), std::invalid_argument);
    }

    SECTION("Box type with valid paired x, y, z constructs without throwing") {
        REQUIRE_NOTHROW(RefinementConfig("Box", {}, {}, {}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}));
    }

    SECTION("Box type with odd-count coordinates throws") {
        REQUIRE_THROWS_AS(RefinementConfig("Box", {}, {}, {}, {0.0, 1.0, 2.0}, {0.0, 1.0, 2.0}, {0.0, 1.0, 2.0}), std::invalid_argument);
    }

    SECTION("Box type with mismatched x/y/z sizes throws") {
        REQUIRE_THROWS_AS(RefinementConfig("Box", {}, {}, {}, {0.0, 1.0}, {0.0, 1.0, 2.0, 3.0}, {0.0, 1.0}), std::invalid_argument);
    }

    SECTION("Unknown type string throws") {
        REQUIRE_THROWS_AS(RefinementConfig("Unknown", {}, {}, {0.1}, {}, {}, {}), std::invalid_argument);
    }

    SECTION("occursPeriodically returns true when both iterations and times are empty") {
        RefinementConfig config("Change", {}, {}, {0.1}, {}, {}, {});
        REQUIRE(config.occursPeriodically() == true);
    }

    SECTION("occursPeriodically returns false when iterations is non-empty") {
        RefinementConfig config("Change", {1}, {}, {0.1}, {}, {}, {});
        REQUIRE(config.occursPeriodically() == false);
    }

    SECTION("occursPeriodically returns false when times is non-empty") {
        RefinementConfig config("Change", {}, {0.1e-9}, {0.1}, {}, {}, {});
        REQUIRE(config.occursPeriodically() == false);
    }
}

