#include <catch.h>
#include <inits.h>
#include <eventlist.h>
#include <meshrefinement.h>
#include "simulation-state.h"
#include <refinement/refinement-spec.h>

using namespace std;

TEST_CASE("Initialise mesh refinement events") {

    EventList eventList;
    REQUIRE(eventList.eventsInQueue() == false);
    REQUIRE(eventList.eventsInQueue() == false);

    MeshRefinement refinement;

    SECTION("periodically occurring mesh refinement") {
        // ARRANGE
        refinement.setRepRefIter(13);
        vector<RefinementConfig> refConfig = {
            RefinementConfig("Sphere", {}, {}, {1}, {0}, {0}, {0})
        };

        refinement.setRefinementConfig(std::move(refConfig));

        // ACT
        createMeshRefinementEvents(refinement, eventList);

        // ASSERT
        REQUIRE(eventList.getPeriodicRefinementIteration() == 13);

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

      Electrodes electrodes = Electrodes::withElectrodePotentials(el);

      createElectrodeSwitchingEvents(electrodes, eventList);

      REQUIRE(eventList.eventsInQueue());

      REQUIRE(4 == eventList.getNumTimeEvents()); // 3 for electrode 1, 1 fore electrode 2 implicit ones at time 0

      double nextEventTime = eventList.timeUntilNextEvent(0.);

      REQUIRE(0 == nextEventTime);

      SimulationState simulationState;
    simulationState.setCurrentTime(nextEventTime);

      Event* e1 = eventList.getCurrentEvent(simulationState);
      Event* e2 = eventList.getCurrentEvent(simulationState);
      REQUIRE((e1->time == 0.0 && e2->time == 0.0));
      delete e1;
      delete e2;

      nextEventTime = eventList.timeUntilNextEvent(simulationState.currentTime().getTime());
      REQUIRE(0.2 == nextEventTime);
    simulationState.setCurrentTime(0.2);
      Event* e3 = eventList.getCurrentEvent(simulationState);
      REQUIRE(e3->time == 0.2);
      delete e3;

      nextEventTime = eventList.timeUntilNextEvent(simulationState.currentTime().getTime());
      REQUIRE(0.2 == nextEventTime);
    simulationState.setCurrentTime(0.4);
      Event* e4 = eventList.getCurrentEvent(simulationState);
      REQUIRE(e4->time == 0.4);
      delete e4;

      REQUIRE_FALSE(eventList.eventsInQueue());
    }
}

TEST_CASE("createMeshRefinementEvents inserts typed RefinementSpec events") {

    SECTION("Change config at iteration 7 -> iteration event carries RefinementSpec with getIteration()==7 and type Change") {
        EventList eventList;
        MeshRefinement refinement;
        vector<RefinementConfig> refConfig = {
            RefinementConfig("Change", {7}, {}, {0.1}, {}, {}, {})
        };
        refinement.setRefinementConfig(std::move(refConfig));

        createMeshRefinementEvents(refinement, eventList);

        REQUIRE(eventList.getNumIterationEvents() == 1);

        SimulationState state;
        state.currentIteration(7);
        Event* e = eventList.getCurrentEvent(state);
        REQUIRE(e->getEventType() == EVENT_REFINEMENT);
        REQUIRE(e->getRefinementSpec() != nullptr);
        REQUIRE(e->getRefinementSpec()->getIteration() == 7);
        REQUIRE(e->getRefinementSpec()->getType() == RefinementSpec::Type::Change);
        delete e;
    }

    SECTION("Sphere config at time 2e-9 -> time event carries RefinementSpec with getTime()==2e-9 and isPeriodic()==false") {
        EventList eventList;
        MeshRefinement refinement;
        vector<RefinementConfig> refConfig = {
            RefinementConfig("Sphere", {}, {2e-9}, {0.5}, {0.0}, {0.0}, {0.0})
        };
        refinement.setRefinementConfig(std::move(refConfig));

        createMeshRefinementEvents(refinement, eventList);

        REQUIRE(eventList.getNumTimeEvents() == 1);

        SimulationState state;
        state.setCurrentTime(2e-9);
        Event* e = eventList.getCurrentEvent(state);
        REQUIRE(e->getEventType() == EVENT_REFINEMENT);
        REQUIRE(e->getRefinementSpec() != nullptr);
        REQUIRE(e->getRefinementSpec()->getTime() == Approx(2e-9));
        REQUIRE_FALSE(e->getRefinementSpec()->isPeriodic());
        delete e;
    }

    SECTION("Invalid RefinementConfig (empty values for Change type) throws before any event is inserted") {
        EventList eventList;
        MeshRefinement refinement;
        vector<RefinementConfig> refConfig;
        // Constructing an invalid RefinementConfig should throw at construction time
        REQUIRE_THROWS_AS(
            refConfig.push_back(RefinementConfig("Change", {1}, {}, {}, {}, {}, {})),
            std::invalid_argument
        );
        // No events should have been inserted
        REQUIRE(eventList.getNumIterationEvents() == 0);
    }
}
