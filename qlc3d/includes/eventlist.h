#ifndef EVENTLIST_H
#define EVENTLIST_H

#include <stdio.h>
#include <vector>
#include <list>
#include <electrodes.h>
#include <simu.h>
#include <limits>
#include <refinfo.h>
using namespace std;


// EVENT CAN OCCUR:
//      AT PRE-DEFINED TIMES
//      AT PRE-DEFINED ITERATIONS
//      WITH REPEATING TIME-INTERVALLS (EVERY Nth SECOND)
//      WITH REPEATING ITERATION-INTERVALLS (EVERY Nth ITERATION)

enum EventType {    EVENT_SWITCHING,    // SWITCH ELECTRODE
                    EVENT_SAVE,         // SAVE RESULT
                    EVENT_REFINEMENT,  // MESH REFINEMENT
                    EVENT_INVALID};     // BAD/ERROR EVENT

enum EventOccurrence{
    EVENT_TIME,         // EVENT OCCUREENCE DETERMINED BY TIME
    EVENT_ITERATION     // EVENT OCCURRENVE DETERMINED BY ITERATION NUMBER
};

// BASE EVENT CLASS
class Event
{
public:
    const EventType eventType_;
    const EventOccurrence eventOccurrence_;
    const double time;
    const unsigned int iteration;
private:
    void *eventData_;    // EVENT DATA IS AN OPTIONAL CUSTOM STRUCT/OBJECT WITH ADDITIONAL DATA ABOUT THE EVENT
public:



    Event(const EventType& et, const double& time, void* const ed = NULL):
        eventType_(et),
        eventOccurrence_(EVENT_TIME),
        time(time),
        iteration(0),
        eventData_(ed)
    { }

    Event(const EventType& et, const unsigned int iter, void* const ed = NULL):
        eventType_(et),
        eventOccurrence_(EVENT_ITERATION),
        time(0),
        iteration(iter),
        eventData_(ed)
    {}


    ~Event();
    EventType getEventType()const {return eventType_;}
    EventOccurrence getEventOccurrence() const {return eventOccurrence_;}
    void* getEventDataPtr()
    {
        return eventData_;
    }

    void setEventDataPtr(void* ed)
    {
        eventData_ = ed;
    }
    static const char* getEventString( const EventType );   // DEBUG, RETURN EVENT TYPE STRING
};

/*
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
*/

/*
class IterEvent:public Event
{
    size_t iteration_;
    IterEvent():Event(EVENT_INVALID){} // PRIVATE CONSTRUCTOR FORCES USE OF OTHER CONSTRUCTORS
public:
    IterEvent( const EventType& et, const size_t& iter);
    ~IterEvent() {}
    size_t getEventIteration()const {return iteration_;}
    bool operator <(const IterEvent& other) const;
    bool occursNow(const Simu& simu) const;

};
*/

// A CLASS FOR TAKING CARE OF PROPERLY DEALLOCATING
// ALL TYPES OF EVENTS
/*
class EventExterminator{

public:
    static void exterminateEvent(Event *&e)
    {
    // EVENTS CAN BE EITHER TIME-EVENTS OR ITEREVENTS

        // TRY TO CAST TO ITER EVENT, DESTROY AND EXIT FUNCTION
        IterEvent* ie = dynamic_cast <IterEvent*> (e);
        if (ie)
        {
            delete ie;
            e = NULL;
            return;
        }

        // TRY TO CAST TO TIME EVENT, DESTROY AND EXIT FUNCTION
        TimeEvent* te = dynamic_cast <TimeEvent*> (e);
        if (te)
        {
            delete te;
            e = NULL;
            return;
        }

        // IF NEITHER CAST DIDN'T WORK, SOMETHING IS REALLY WRONG
        printf("error in %s, could not destroy event - bye\n", __func__);
        exit(1);

    }
};
*/


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




//
//  EVENTLIST MANAGES EVENTS.
//
class EventList
{
    private:

        // SPECIAL CONSTANTS
        const double NO_TIME_EVENTS;    // INDICATES THAT NO TIME EVENTS EXIST
        const size_t NO_ITER_EVENTS;    // INDINCATES THAT NO ITER EVENT EXIST

        double nextTimeEvent_;  // TIME OF NEXT TIMED EVENT
        size_t nextIterEvent_;  // ITERATION NUMBER OF NEXT ITER EVENT

        list <Event*> timeEvents_;
        list <Event*> iterationEvents_;
        list <Event*> repRefinements_;  // ALL OF THESE WILL BE EXECUTED SIMULTANEOUSLY

        // REOCCURRING SAVE EVENTS
        size_t saveIter_;    // SAVE PERIOD IN ITERATIONS
        double saveTime_;    // SAVE PERIOD IN SECONDS
        size_t saveTimeCount_; // KEEPS COUNT OF PROCESSED REOCCURRING SAVE ITERS SO FAR

        // REPEATED REFINMENT EVENTS
        size_t repRefIter_;        // REFINEMENT PERIOD IN ITERATIONS 0 -> NEVER
        double repRefTime_;        // REFINEMENT PERIOD IN SECONDS    0 -> NEVER
        size_t repRefTimeCount_;   // KEEPS COUNT OF PROCESSED ROCCURRING REFINEMENT EVENT SO FAR


        void prependReoccurringIterEvent(Event* iEvent); // ADDS REOCCURRING TIME EVENT TO FRONT OF QUEUE

    public:
        EventList();
        ~EventList();
        bool eventsInQueue() const; // return true if further events exist in queue
        bool eventOccursNow(const Simu& simu) const;

        void insertTimeEvent(Event* tEvent);
        void insertIterEvent(Event* iEvent);
        void setSaveIter(const size_t& si){saveIter_ = si;}
        void setSaveTime(const double& st);//{saveTime_ = st;}
        size_t getSaveIter() const {return saveIter_;}
        double getSaveTime() const {return saveTime_;}

        void setRepRefIter( const size_t rri ) {repRefIter_ = rri;}
        void setRepRefTime( const double rrt ) {repRefTime_ = rrt;}
        void addRepRefInfo( Event* repRefEvent); // ADDS REPEATING EVENT TO repRefinements LIST


        Event* getCurrentEvent(const Simu& simu);    // removes current event from queue and returns a copy of it
        double timeUntilNextEvent(const Simu& simu) const;

        void manageReoccurringEvents(Simu& simu); // PERIODICALLY ADDS REOCCURING EVENTS TO EVENT QUEQUE

        void printEventList() const;


};





#endif




