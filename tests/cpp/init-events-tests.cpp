#include <catch.h>
#include <inits.h>
#include <eventlist.h>
#include <meshrefinement.h>

using namespace std;

TEST_CASE("Initialise mesh refinement events") {

    EventList eventList;
    REQUIRE(eventList.eventsInQueue() == false);
    REQUIRE(eventList.eventsInQueue() == false);

    MeshRefinement refinement;

    SECTION("periodically occurring mesh refinement") {
        // ARRANGE
        refinement.setRepRefIter(13);
        refinement.setRepRefTime(2e-9);
        vector<RefinementConfig> refConfig = {
            RefinementConfig("Sphere", {}, {}, {1}, {0}, {0}, {0})
        };

        refinement.setRefinementConfig(std::move(refConfig));

        // ACT
        createMeshRefinementEvents(refinement, eventList);

        // ASSERT
        REQUIRE(eventList.getPeriodicRefinementIteration() == 13);
        REQUIRE(eventList.getPeriodicRefinementTime() == 2e-9);

        // no events with explicit time or iteration should exist
        REQUIRE(eventList.getNumTimeEvents() == 0);
        REQUIRE(eventList.getNumIterationEvents() == 0);
    }

    SECTION("explicit time and iteration mesh refinement") {
        // ARRANGE: refinement occurring at iterations 1, 2, and time 0.1s
        vector<RefinementConfig> refConfig = {
            RefinementConfig("Box", {1, 2}, {0.1}, {}, {0, 1}, {0, 1}, {0, 1})
        };

        refinement.setRefinementConfig(std::move(refConfig));

        // ACT
        createMeshRefinementEvents(refinement, eventList);

        // ASSERT
        // no periodically occurring event should exist
        REQUIRE(eventList.getPeriodicRefinementIteration() == 0);
        REQUIRE(eventList.getPeriodicRefinementTime() == 0.);

        // two iteration events and one time event should exist
        REQUIRE(eventList.getNumIterationEvents() == 2);
        REQUIRE(eventList.getNumTimeEvents() == 1);
    }
}

