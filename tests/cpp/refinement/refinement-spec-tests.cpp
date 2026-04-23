#include <catch.h>
#include <refinement/refinement-spec.h>
#include <meshrefinement.h>

TEST_CASE("RefinementSpec factory methods") {

    SECTION("makeExplicit with iteration=5, time=0 -> isPeriodic is false, getIteration()==5") {
        auto spec = RefinementSpec::makeExplicit("Sphere", 5, 0, {0.5}, {0.0}, {0.0}, {0.0});
        REQUIRE_FALSE(spec->isPeriodic());
        REQUIRE(spec->getIteration() == 5);
    }

    SECTION("makeExplicit with iteration=0, time=1e-9 -> isPeriodic is false, getTime()==1e-9") {
        auto spec = RefinementSpec::makeExplicit("Sphere", 0, 1e-9, {0.5}, {0.0}, {0.0}, {0.0});
        REQUIRE_FALSE(spec->isPeriodic());
        REQUIRE(spec->getTime() == Approx(1e-9));
    }

    SECTION("makeExplicit with both iteration>0 and time>0 throws") {
        REQUIRE_THROWS_AS(RefinementSpec::makeExplicit("Sphere", 5, 1e-9, {0.5}, {0.0}, {0.0}, {0.0}), std::invalid_argument);
    }

    SECTION("makeExplicit with both <= 0 throws") {
        REQUIRE_THROWS_AS(RefinementSpec::makeExplicit("Sphere", 0, 0, {0.5}, {0.0}, {0.0}, {0.0}), std::invalid_argument);
    }

    SECTION("makePeriodic -> isPeriodic is true") {
        auto spec = RefinementSpec::makePeriodic("Sphere", {0.5}, {0.0}, {0.0}, {0.0});
        REQUIRE(spec->isPeriodic());
        REQUIRE(spec->getIteration() == 0);
        REQUIRE(spec->getTime() == 0.0);
    }

    SECTION("getRefIter for Change type equals values.size()") {
        auto spec = RefinementSpec::makeExplicit("Change", 1, 0, {0.1, 0.2, 0.3}, {}, {}, {});
        REQUIRE(spec->getRefIter() == 3);
    }

    SECTION("getRefIter for Box type equals x.size() / 2") {
        auto spec = RefinementSpec::makeExplicit("Box", 1, 0, {}, {0.0, 1.0, 0.0, 1.0}, {0.0, 1.0, 0.0, 1.0}, {0.0, 1.0, 0.0, 1.0});
        REQUIRE(spec->getRefIter() == 2);
    }

    SECTION("Validation errors are enforced - Change with empty values throws") {
        REQUIRE_THROWS_AS(RefinementSpec::makeExplicit("Change", 1, 0, {}, {}, {}, {}), std::invalid_argument);
    }

    SECTION("Validation errors are enforced - Box with odd coords throws") {
        REQUIRE_THROWS_AS(RefinementSpec::makeExplicit("Box", 1, 0, {}, {0.0, 1.0, 2.0}, {0.0, 1.0, 2.0}, {0.0, 1.0, 2.0}), std::invalid_argument);
    }

    SECTION("Validation errors are enforced - unknown type throws") {
        REQUIRE_THROWS_AS(RefinementSpec::makeExplicit("Unknown", 1, 0, {0.1}, {}, {}, {}), std::invalid_argument);
    }
}

TEST_CASE("RefinementConfig::toSpecs") {

    SECTION("Box config with no iters/times -> periodic spec") {
        RefinementConfig config("Box", {}, {}, {}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0});
        auto specs = config.toSpecs();
        REQUIRE(specs.size() == 1);
        REQUIRE(specs[0]->isPeriodic());
        REQUIRE(specs[0]->getType() == RefinementSpec::Type::Box);
    }

    SECTION("Change config with iteration -> explicit iteration spec") {
        RefinementConfig config("Change", {3}, {}, {0.05}, {}, {}, {});
        auto specs = config.toSpecs();
        REQUIRE(specs.size() == 1);
        REQUIRE_FALSE(specs[0]->isPeriodic());
        REQUIRE(specs[0]->getIteration() == 3);
        REQUIRE(specs[0]->getType() == RefinementSpec::Type::Change);
    }

    SECTION("Config with multiple iterations creates multiple specs") {
        RefinementConfig config("Change", {1, 2, 5}, {}, {0.05}, {}, {}, {});
        auto specs = config.toSpecs();
        REQUIRE(specs.size() == 3);
        REQUIRE(specs[0]->getIteration() == 1);
        REQUIRE(specs[1]->getIteration() == 2);
        REQUIRE(specs[2]->getIteration() == 5);
    }

    SECTION("Config with time creates time-based spec") {
        RefinementConfig config("Sphere", {}, {2e-9}, {0.5}, {0.0}, {0.0}, {0.0});
        auto specs = config.toSpecs();
        REQUIRE(specs.size() == 1);
        REQUIRE_FALSE(specs[0]->isPeriodic());
        REQUIRE(specs[0]->getTime() == Approx(2e-9));
    }
}

