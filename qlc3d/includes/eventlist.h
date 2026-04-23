#ifndef EVENTLIST_H
#define EVENTLIST_H

#include <stdio.h>
#include <memory>
#include <vector>
#include <list>
#include <electrodes.h>
#include <limits>
#include <refinement/refinement-spec.h>

class Simu;
class SimulationState;
class SimulationTime;
class Reader;
using namespace std;

enum EventType {
  EVENT_SWITCHING,
  EVENT_SAVE,
  EVENT_REFINEMENT,
  SAVE_PERIODICALLY_BY_TIME//
  //SAVE_PERIODICALLY_BY_ITERATION,
};
enum EventOccurrence {
    EVENT_TIME,         // EVENT OCCURRENCE DETERMINED BY TIME
    EVENT_ITERATION     // EVENT OCCURRENCE DETERMINED BY ITERATION NUMBER
};

class Event {
    /*! Event class, represents electrode switching, result saving etc. */
public:
    const EventType eventType_;
    const EventOccurrence eventOccurrence_;
    const double time;
    const unsigned int iteration;
private:
    void *eventData_;    // EVENT DATA IS AN OPTIONAL CUSTOM STRUCT/OBJECT WITH ADDITIONAL DATA ABOUT THE EVENT (used for switching)
    std::unique_ptr<RefinementSpec> refinementSpec_; // typed payload for refinement events
public:
    Event(const EventType &et, const double &time, void *const ed = NULL):
        eventType_(et),
        eventOccurrence_(EVENT_TIME),
        time(time),
        iteration(0),
        eventData_(ed) {
    }
    Event(const EventType &et, const unsigned int iter, void *const ed = NULL):
        eventType_(et),
        eventOccurrence_(EVENT_ITERATION),
        time(0),
        iteration(iter),
        eventData_(ed) {
    }

    // Move constructor to allow moving the unique_ptr payload
    Event(Event &&other) noexcept :
        eventType_(other.eventType_),
        eventOccurrence_(other.eventOccurrence_),
        time(other.time),
        iteration(other.iteration),
        eventData_(other.eventData_),
        refinementSpec_(std::move(other.refinementSpec_)) {
        other.eventData_ = nullptr;
    }

    ~Event() noexcept(false);

    /**
     * Create a typed refinement event that fires at the given iteration.
     * @param spec  Owned RefinementSpec payload.
     * @param iteration Iteration at which to fire.
     * @return New Event owning the spec.
     */
    static Event makeRefinement(std::unique_ptr<RefinementSpec> spec, unsigned int iteration);

    /**
     * Create a typed refinement event that fires at the given simulation time.
     * @param spec  Owned RefinementSpec payload.
     * @param time  Simulation time at which to fire.
     * @return New Event owning the spec.
     */
    static Event makeRefinement(std::unique_ptr<RefinementSpec> spec, double time);

    EventType getEventType()const {
        return eventType_;
    }
    EventOccurrence getEventOccurrence() const {
        return eventOccurrence_;
    }
    void *getEventDataPtr() {
        return eventData_;
    }

    void setEventDataPtr(void *ed) {
        eventData_ = ed;
    }

    /**
     * @return Pointer to the owned RefinementSpec, or nullptr if not a typed refinement event.
     */
    [[nodiscard]] const RefinementSpec* getRefinementSpec() const { return refinementSpec_.get(); }

    std::string toString() const;
    static const char *getEventString(const EventType);     // DEBUG, RETURN EVENT TYPE STRING
};

class EventList {
    /*!The EventList class manages the individual envents*/
private:
    // SPECIAL CONSTANTS
    static constexpr double NO_TIME_EVENTS = std::numeric_limits<double>::max();    // INDICATES THAT NO TIME EVENTS EXIST
    static constexpr size_t NO_ITER_EVENTS = std::numeric_limits<size_t>::max();    // INDICATES THAT NO ITER EVENT EXIST
    double nextTimeEvent_;          // TIME OF NEXT TIMED EVENT
    size_t nextIterEvent_;          // ITERATION NUMBER OF NEXT ITER EVENT
    list <Event *> timeEvents_;     // Collection od time/iteration events
    list <Event *> iterationEvents_;

    /** periodically (by iteration) occurring mesh refinement events */
    list <Event *> repRefinements_;
    size_t saveIter_;               // SAVE PERIOD IN ITERATIONS
    double saveTime_;               // SAVE PERIOD IN SECONDS
    size_t saveTimeCount_;          // KEEPS COUNT OF PROCESSED REOCCURRING SAVE ITERS SO FAR
    // REPEATED REFINEMENT EVENTS
    size_t repRefIter_ = 0;             // REFINEMENT PERIOD IN ITERATIONS 0 -> NEVER
    void prependReoccurringIterEvent(Event *iEvent); // ADDS REOCCURRING TIME EVENT TO FRONT OF QUEUE


    [[nodiscard]] bool isPeriodicSaveIteration(unsigned int iter) {
      return saveIter_ > 0 && iter % saveIter_ == 0;
    }

public:
    EventList();
    ~EventList();
    bool eventsInQueue() const; // return true if further events exist in queue
    /** Whether an event occurs now */
    bool eventOccursNow(const SimulationState &simulationState) const;
    void insertTimeEvent(Event *tEvent);
    void insertIterEvent(Event *iEvent);
    void setSaveIter(const size_t &si) {
        saveIter_ = si;
    }
    void setSaveTime(const double &st);
    size_t getSaveIter() const { return saveIter_; }
    double getSaveTime() const { return saveTime_; }
    void setRepRefIter(const size_t rri) { repRefIter_ = rri; }

    [[nodiscard]] size_t getPeriodicRefinementIteration() const { return repRefIter_; }
    [[nodiscard]] size_t getNumPeriodicRefinementObjects() const { return repRefinements_.size(); }

    /** total number of time events queued up */
    [[nodiscard]] size_t getNumTimeEvents() const { return timeEvents_.size(); }

    /** total number of iteration event queued up */
    [[nodiscard]] size_t getNumIterationEvents() const { return iterationEvents_.size(); }

    /** add periodically reoccurring refinement event */
    void addRepRefInfo(Event *repRefEvent);

    Event *getCurrentEvent(const SimulationState &simulationState);   // removes current event from queue and returns a copy of it
    [[nodiscard]] double timeUntilNextEvent(const double &currentTime) const;
    [[nodiscard]] SimulationTime nextEventTime() const;

    //! Adds reoccurring event to queue
    void manageReoccurringEvents(int currentIteration, const SimulationTime &currentTime, double timeStep);
    std::string toString() const;
};
#endif




