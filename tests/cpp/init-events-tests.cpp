#include <catch.h>
#include <inits.h>
#include <eventlist.h>
#include <meshrefinement.h>
#include "simulation-state.h"

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

    SECTION("create electrode switching events") {
      std::vector<std::shared_ptr<Electrode>> el;
      auto electrode1 = std::shared_ptr<Electrode>(new Electrode(1, {0.2}, {2}));
      auto electrode2 = std::shared_ptr<Electrode>(new Electrode(2, {0.4}, {4}));
      el.push_back(electrode1);
      el.push_back(electrode2);

      Electrodes electrodes(el, nullptr);

      createElectrodeSwitchingEvents(electrodes, eventList);

      REQUIRE(eventList.eventsInQueue());

      REQUIRE(4 == eventList.getNumTimeEvents()); // 3 for electrode 1, 1 fore electrode 2 implicit ones at time 0

      double nextEventTime = eventList.timeUntilNextEvent(0.);

      REQUIRE(0 == nextEventTime);

      SimulationState simulationState;
      simulationState.currentTime(nextEventTime);

      Event* e1 = eventList.getCurrentEvent(simulationState);
      Event* e2 = eventList.getCurrentEvent(simulationState);
      REQUIRE((e1->time == 0.0 && e2->time == 0.0));
      delete e1;
      delete e2;

      nextEventTime = eventList.timeUntilNextEvent(simulationState.currentTime());
      REQUIRE(0.2 == nextEventTime);
      simulationState.currentTime(0.2);
      Event* e3 = eventList.getCurrentEvent(simulationState);
      REQUIRE(e3->time == 0.2);
      delete e3;

      nextEventTime = eventList.timeUntilNextEvent(simulationState.currentTime());
      REQUIRE(0.2 == nextEventTime);
      simulationState.currentTime(0.4);
      Event* e4 = eventList.getCurrentEvent(simulationState);
      REQUIRE(e4->time == 0.4);
      delete e4;

      REQUIRE_FALSE(eventList.eventsInQueue());
    }
}

