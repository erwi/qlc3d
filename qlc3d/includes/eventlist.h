#ifndef EVENTLIST_H
#define EVENTLIST_H

#include <stdio.h>
#include <vector>
#include <list>
#include <electrodes.h>
#include <simu.h>
#include <limits>
using namespace std;


// EVENT CAN OCCUR:
//      AT PRE-DEFINED TIMES
//      AT PRE-DEFINED ITERATIONS
//      WITH REPEATING TIME-INTERVALLS (EVERY Nth SECOND)
//      WITH REPEATING ITERATION-INTERVALLS (EVERY Nth ITERATION)

enum EventType {    EVENT_SWITCHING,    // SWITCH ELECTRODE
                    EVENT_SAVE,         // SAVE RESULT
                    EVENT_REFIENEMENT,  // MESH REFINEMENT
                    EVENT_INVALID};     // BAD/ERROR EVENT
// BASE EVENT CLASS
class Event
{
    const EventType eventType_;
    void *eventData_;    // EVENT DATA IS AN OPTIONAL CUSTOM STRUCT/OBJECT WITH ADDITIONAL DATA ABOUT THE EVENT
public:
    Event(const EventType& et, void* const ed = NULL):
        eventType_(et),
        eventData_(ed)
    { }
    ~Event(){}
    EventType getEventType()const {return eventType_;}
    void* getEventDataPtr()
    {
        return eventData_;
    }

    void setEventDataPtr(void* ed)
    {
        eventData_ = ed;
    }
};

class TimeEvent:public Event
{
    double time_;
    TimeEvent():Event(EVENT_INVALID){} // PRIVATE CONSTRUCTOR FORCES USE OF OTHER CONSTRUCTORS

public:
    TimeEvent( const EventType& et, const double& time);
    ~TimeEvent() {}
    double getEventTime() const  {return time_;}
    bool operator <(const TimeEvent& other) const;
    bool occursNow(const Simu &simu) const;
};


class IterEvent:public Event
{
    size_t iteration_;
    IterEvent():Event(EVENT_INVALID){} // PRIVATE CONSTRUCTOR FORCES USE OF OTHER CONSTRUCTORS
public:
    IterEvent( const EventType& et, const size_t& iter);
    ~IterEvent() {};
    size_t getEventIteration()const {return iteration_;}
    bool operator <(const IterEvent& other) const;
    bool occursNow(const Simu& simu) const;
};

/*
class SwitchingEvent:public TimeEvent
{
    const double potential_;
    const size_t electrodeNumber_;
public:
    SwitchingEvent(const double& time,
                   const double& potential,
                   const size_t& electrodeNumber );
    ~SwitchingEvent(){};
    size_t getElectrodeNumber()const{return electrodeNumber_;}
    double getElectrodePotential()const{return potential_;}
};
*/





class EventList
{
    private:

        // SPECIAL CONSTANTS
        const double NO_TIME_EVENTS;    // INDICATES THAT NO TIME EVENTS EXIST
        const size_t NO_ITER_EVENTS;    // INDINCATES THAT NO ITER EVENT EXIST

        double nextTimeEvent_;  // TIME OF NEXT TIMED EVENT
        size_t nextIterEvent_;  // ITERATION NUMBER OF NEXT ITER EVENT

        list <TimeEvent*> timeEvents_;
        list <IterEvent*> iterationEvents_;
	
        // REOCCURRING EVENTS
        size_t saveIter_;    // SAVE PERIOD IN ITERATIONS

        double saveTime_;    // SAVE PERIOD IN SECONDS
        size_t saveTimeCount_; // KEEPS COUNT OF PROCESSED REOCCURRING SAVE ITERS SO FAR
        //TimeEvent saveTime_;    // SAVE PERIOD IN SECONDS

        void prependReoccurringIterEvent(IterEvent* iEvent);

    public:
        EventList();
        bool eventsInQueue() const; // return true if further events exist in queue
        bool eventOccursNow(const Simu& simu) const;

        void insertTimeEvent(TimeEvent* tEvent);
        void setSaveIter(const size_t& si){saveIter_ = si;}
        void setSaveTime(const double& st);//{saveTime_ = st;}

        Event* getCurrentEvent(const Simu& simu);    // removes current event from queue and returns a copy of it
        double timeUntilNextEvent(const Simu& simu) const;

        void manageReoccurringEvents(Simu& simu); // PERIODICALLY ADDS REOCCURING EVENTS TO EVENT QUEQUE

        void printEventList() const;

        size_t getSaveIter() const {return saveIter_;}
        double getSaveTime() const {return saveTime_;}
};





#endif




