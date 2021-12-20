#ifndef EVENTLIST_H
#define EVENTLIST_H

#include <stdio.h>
#include <vector>
#include <list>
#include <electrodes.h>
#include <limits>

// Forward declarations
class Simu;
class SimulationState;
class Reader;
class RefInfo;
using namespace std;
// EVENT CAN OCCUR:
//      AT PRE-DEFINED TIMES
//      AT PRE-DEFINED ITERATIONS
//      WITH REPEATING TIME-INTERVALLS (EVERY Nth SECOND)
//      WITH REPEATING ITERATION-INTERVALLS (EVERY Nth ITERATION)
enum EventType {    EVENT_SWITCHING,    // SWITCH ELECTRODE
                    EVENT_SAVE,         // SAVE RESULT
                    EVENT_REFINEMENT,   // MESH REFINEMENT
                    EVENT_INVALID
               };
enum EventOccurrence {
    EVENT_TIME,         // EVENT OCCUREENCE DETERMINED BY TIME
    EVENT_ITERATION     // EVENT OCCURRENVE DETERMINED BY ITERATION NUMBER
};

class Event {
    /*! Event class, represents electrode switching, result saving etc. */
public:
    const EventType eventType_;
    const EventOccurrence eventOccurrence_;
    const double time;
    const unsigned int iteration;
private:
    void *eventData_;    // EVENT DATA IS AN OPTIONAL CUSTOM STRUCT/OBJECT WITH ADDITIONAL DATA ABOUT THE EVENT
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
    ~Event();

    /** Create new unmanaged periodic mesh refinement event. It will need to be manually deleted!*/
    static Event* ofPeriodicMeshRefinement(RefInfo *refInfo);

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
    static const char *getEventString(const EventType);     // DEBUG, RETURN EVENT TYPE STRING
};

class EventList {
    /*!The EventList class manages the individual envents*/
private:
    // SPECIAL CONSTANTS
    const double NO_TIME_EVENTS;    // INDICATES THAT NO TIME EVENTS EXIST
    const size_t NO_ITER_EVENTS;    // INDINCATES THAT NO ITER EVENT EXIST
    double nextTimeEvent_;          // TIME OF NEXT TIMED EVENT
    size_t nextIterEvent_;          // ITERATION NUMBER OF NEXT ITER EVENT
    list <Event *> timeEvents_;     // Collection od time/iteration events
    list <Event *> iterationEvents_;

    /** periodically (by iteration and/or time) occurring mesh refinement events */
    list <Event *> repRefinements_;
    size_t saveIter_;               // SAVE PERIOD IN ITERATIONS
    double saveTime_;               // SAVE PERIOD IN SECONDS
    size_t saveTimeCount_;          // KEEPS COUNT OF PROCESSED REOCCURRING SAVE ITERS SO FAR
    // REPEATED REFINMENT EVENTS
    size_t repRefIter_ = 0;             // REFINEMENT PERIOD IN ITERATIONS 0 -> NEVER
    double repRefTime_ = 0;             // REFINEMENT PERIOD IN SECONDS    0 -> NEVER
    size_t repRefTimeCount_;        // KEEPS COUNT OF PROCESSED ROCCURRING REFINEMENT EVENT SO FAR
    void prependReoccurringIterEvent(Event *iEvent); // ADDS REOCCURRING TIME EVENT TO FRONT OF QUEUE
public:
    EventList();
    ~EventList();
    bool eventsInQueue() const; // return true if further events exist in queue
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
    void setRepRefTime(const double rrt) { repRefTime_ = rrt; }

    [[nodiscard]] size_t getPeriodicRefinementIteration() const { return repRefIter_; }
    [[nodiscard]] double getPeriodicRefinementTime() const { return repRefTime_; }
    [[nodiscard]] size_t getNumPeriodicRefinementObjects() const { return repRefinements_.size(); }

    /** total number of time events queued up */
    [[nodiscard]] size_t getNumTimeEvents() const { return timeEvents_.size(); }

    /** total number of iteration event queued up */
    [[nodiscard]] size_t getNumIterationEvents() const { return iterationEvents_.size(); }

    /** add periodically reoccurring refinement event */
    void addRepRefInfo(Event *repRefEvent);

    Event *getCurrentEvent(const SimulationState &simulationState);   // removes current event from queue and returns a copy of it
    double timeUntilNextEvent(const double &currentTime) const;

    //! Adds reoccurring event to queue
    void manageReoccurringEvents(int currentIteration, double currentTime, double timeStep);
};
#endif




