#include <catch.h>
#include <eventlist.h>
#include <refinement/refinement-spec.h>
#include <simulation-state.h>

TEST_CASE("Event with RefinementSpec payload") {

    SECTION("makeRefinement(spec, iter) has correct type and occurrence") {
        auto spec = RefinementSpec::makeExplicit("Sphere", 3, 0, {0.5}, {0.0}, {0.0}, {0.0});
        Event e = Event::makeRefinement(std::move(spec), 3u);
        REQUIRE(e.getEventType() == EVENT_REFINEMENT);
        REQUIRE(e.getEventOccurrence() == EVENT_ITERATION);
        REQUIRE(e.getRefinementSpec() != nullptr);
        REQUIRE(e.iteration == 3);
    }

    SECTION("makeRefinement(spec, time) has correct occurrence and time") {
        auto spec = RefinementSpec::makeExplicit("Sphere", 0, 1e-9, {0.5}, {0.0}, {0.0}, {0.0});
        Event e = Event::makeRefinement(std::move(spec), 1e-9);
        REQUIRE(e.getEventOccurrence() == EVENT_TIME);
        REQUIRE(e.getRefinementSpec() != nullptr);
        REQUIRE(e.getRefinementSpec()->getTime() == Approx(1e-9));
    }

    SECTION("Event constructed with makeRefinement is destroyed without crash") {
        {
            auto spec = RefinementSpec::makeExplicit("Sphere", 3, 0, {0.5}, {0.0}, {0.0}, {0.0});
            Event e = Event::makeRefinement(std::move(spec), 3u);
            (void)e;
        } // destructor called here, no crash expected
        REQUIRE(true);
    }

    SECTION("Moving an Event transfers RefinementSpec ownership") {
        auto spec = RefinementSpec::makeExplicit("Sphere", 3, 0, {0.5}, {0.0}, {0.0}, {0.0});
        Event src = Event::makeRefinement(std::move(spec), 3u);
        const RefinementSpec* oldPtr = src.getRefinementSpec();
        REQUIRE(oldPtr != nullptr);

        Event dst = std::move(src);
        REQUIRE(src.getRefinementSpec() == nullptr);  // ownership transferred
        REQUIRE(dst.getRefinementSpec() == oldPtr);   // same address
    }

    SECTION("EventList::addRepRefInfo accepts makeRefinement event; getNumPeriodicRefinementObjects()==1") {
        EventList eventList;
        auto spec = RefinementSpec::makePeriodic("Sphere", {0.5}, {0.0}, {0.0}, {0.0});
        Event* e = new Event(Event::makeRefinement(std::move(spec), 0.0));
        eventList.addRepRefInfo(e);
        REQUIRE(eventList.getNumPeriodicRefinementObjects() == 1);
    }

    SECTION("insertIterEvent with makeRefinement event is retrieved by getCurrentEvent at correct iteration") {
        EventList eventList;
        auto spec = RefinementSpec::makeExplicit("Sphere", 7, 0, {0.5}, {0.0}, {0.0}, {0.0});
        Event* e = new Event(Event::makeRefinement(std::move(spec), 7u));
        eventList.insertIterEvent(e);

        SimulationState state;
        state.currentIteration(7);
        REQUIRE(eventList.eventOccursNow(state));

        Event* retrieved = eventList.getCurrentEvent(state);
        REQUIRE(retrieved->getEventType() == EVENT_REFINEMENT);
        REQUIRE(retrieved->getRefinementSpec() != nullptr);
        REQUIRE(retrieved->getRefinementSpec()->getIteration() == 7);
        delete retrieved;
    }

    SECTION("EventList destructor with unconsumed smart-pointer events does not crash") {
        // Create EventList with an unconsumed event and let it go out of scope
        {
            EventList eventList;
            auto spec = RefinementSpec::makeExplicit("Sphere", 5, 0, {0.5}, {0.0}, {0.0}, {0.0});
            Event* e = new Event(Event::makeRefinement(std::move(spec), 5u));
            eventList.insertIterEvent(e);
        } // EventList destructor called; event should be cleaned up
        REQUIRE(true);
    }
}


