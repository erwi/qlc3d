#include <eventlist.h>
#include <math.h>



//-------- TimeEvent DEFINITINS -------------
TimeEvent::TimeEvent(const EventType &et, const double &time):
    Event(et),
    time_(time)
{

}
bool TimeEvent::operator <(const TimeEvent& other) const
{
    return ( this->getEventTime() < other.getEventTime() );
}

bool compare(const TimeEvent* first, const TimeEvent* second)
{
    // COMPARISON OF POINTERS TO TIME EVENTS
    return ( first->getEventTime() < second->getEventTime() );
}


//-------- IterEvent DEFINITINS ------------
IterEvent::IterEvent(const EventType &et, const size_t &iter):
    Event(et),
    iteration_(iter)
{

}

//---------SwitchingEvent DEFINITIONS ---------------
SwitchingEvent::SwitchingEvent(const double &time,
                               const double &potential,
                               const size_t &electrodeNumber):
    TimeEvent(EVENT_SWITCHING, time),
    potential_(potential),
    electrodeNumber_(electrodeNumber)
{

}




bool IterEvent::operator <( const IterEvent& other) const
{
    return ( getEventIteration() < other.getEventIteration() );
}

// -------------EventList DEFINITIONS --------------

EventList::EventList():
    NO_TIME_EVENTS(std::numeric_limits<double>::max() ),
    NO_ITER_EVENTS(std::numeric_limits<size_t>::max() ),
    nextTimeEvent_(NO_TIME_EVENTS),
    nextIterEvent_(NO_ITER_EVENTS),
    saveIter_(0),
    saveTime_(0),
    saveTimeCount_(0)
{


}

bool EventList::eventsInQueue() const
{
// RETURNS BOOLEAN INDICATING WHETHER ANY EVENTS EXIST IN QUEUE
    return ( ( nextIterEvent_ < NO_ITER_EVENTS ) || // ITER EVENTS EXIST?
             ( nextTimeEvent_ < NO_TIME_EVENTS )    // TIME EVENTS EXIST?
                );

}

bool EventList::eventOccursNow(const Simu &simu) const
{
// RETURNS TRUE IF AN EVENT OCCURS RIGHT NOW!!

    // OTHER EXPLICITLY SPECIFIED EVENTS
    return ((nextIterEvent_==(size_t)simu.getCurrentIteration()) || // IF EVENT SCHEDULED FOR CURRENT ITERATION
            (nextTimeEvent_==simu.getCurrentTime() ) ) ;    // IF EVENT SCHEDULED FOR CURRENT TIME

}

Event* EventList::getCurrentEvent(const Simu &simu)
{
// EVENTS EXIST AS OBJECTS ON HEAP. ONLY POINTERS TO THESE
// ARE STORED IN THE EVENT LIST.
// THIS FUNCTION RETURNS A POINTER TO AN EVENT OBJECT
// AND REMOVES THE CORRESPONDING POINTER FROM QUEUE.
// ACTUAL DELETION OF OBJECT MUST BE PERFORMED ELSEWHERE


// FIRST CHECK TIME EVENTS
    if ( nextTimeEvent_ == simu.getCurrentTime() )
    {
        if ( ! timeEvents_.size() )
        {
            printf("error in %s, no time events in queue - bye!\n", __func__ );
            exit(1);
        }


        TimeEvent* timeEve( ( *timeEvents_.begin() ) );   // MAKE A POINTER TO OBJECT ON HEAP
        Event* eve = static_cast<Event*>( timeEve );      // CAST TO EVENT
        timeEvents_.erase( timeEvents_.begin() );         // REMOVE EVENT POINTER FROM QUEUE FRONT

        // UPDATE TIME FOR NEXT TIME EVENT MONITOR
        if ( timeEvents_.size() )// THIS NEED MODIFYING TO ALLOW FOR REPEATING EVENTS
        {
            nextTimeEvent_ = (*timeEvents_.begin() )->getEventTime();
        }
        else
        {
            nextTimeEvent_ = NO_TIME_EVENTS;    // INDICATE END OF QUEUE
        }

        return eve;
    }
// IF NO CURRENT TIME EVENTS LEFT, CHECK ITER EVENTS
    else
    if ( nextIterEvent_ == (size_t) simu.getCurrentIteration() )
    {
        if ( ! iterationEvents_.size() )
        {
            printf("error in %s, no iteration event is queue - bye !\n", __func__);
            exit(1);
        }

        const IterEvent* iterEve ( (*iterationEvents_.begin() ) ); // MAKE A POINTER TO OBJECT ON HEAP
        const Event* eve = static_cast<const Event*>( iterEve );   // CAST TO EVENT
        iterationEvents_.erase(iterationEvents_.begin());          // REMOVE EVENT POINTER FROM QUEUE FRONT

// UPDATE NEXT ITER EVENT MONITOR
        if ( iterationEvents_.size() )  // THIS NEED MODIFYING TO ALLOW FOR REPEATING EVENTS
        {
            nextIterEvent_ = ( *iterationEvents_.begin() )->getEventIteration();
        }
        else
        {
            nextIterEvent_ = NO_ITER_EVENTS;
        }

        return const_cast<Event*>(eve);
    }
    else
    {
        // ERROR. SOMETHING HAS GONE WRONG. NO CURRENT EVENT
        printf("error in %s. No current events - bye !\n", __func__);
        exit(1);

    }

}// end getCurrent Event


void EventList::setSaveTime(const double &st)
{
// CHECK AND SET REOCCURING SAVE EVENT TIME
    if ( st >= 0 )
    {
        saveTime_ = st;
    }
    else
    {
        printf("error in %s. Bad SaveTime = %e - bye!\n", __func__ , st);
        exit(1);
    }
}


void EventList::insertTimeEvent(TimeEvent *tEvent)
{
// ADDS TIMED EVENT TO EVENT QUEUE
    timeEvents_.push_back( tEvent );
    timeEvents_.sort( compare );    // USE COMPARISON OF *POINTERS* TO TIME EVENT (OTHERWISE SORTING BY ADDRESS VALUE)

// KEEP NEXT TIME EVENT TRACKER UP TO DATE
    if (tEvent->getEventTime() < nextTimeEvent_ )
        nextTimeEvent_ = tEvent->getEventTime();
}

void EventList::manageReoccurringEvents(Simu& simu)
{
// ADDS REOCCURRING EVENTS TO FRONT OF QUEUE IF NEEDED.
// NEW EVENTS WILL BE SCEDULED FOR NEXT ITERATION/TIME STEP.
// NOTE THAT THIS SHOULD BE CALLED AFTER ALL EVENTS FOR CURRENT
// ITERATION/TIMESTEP HAVE BEEN PROCESSED.

// ADD ITER EVENTS
    size_t nextIter = (size_t) simu.getCurrentIteration() + 1;
    double nextTime = simu.getCurrentTime() + simu.getdt(); // TIME AT NEXT ITERATION


// REOCCURRING SAVE ITERATION
    if ( ( saveIter_ > 0 ) &&
         ( nextIter % saveIter_ ) == 0 )
    {
        prependReoccurringIterEvent( new IterEvent(EVENT_SAVE ,  nextIter) );
    }

// REOCCURRING SAVE TIME
    if (saveTime_ > 0 )
    {

        // PREDICTED NUM SAVE EVENTS BY NEXT TIME STEP, IF dt IS NOT ADJUSTED
        size_t num_next = static_cast<size_t> (nextTime / saveTime_ );

        // IF A REOCCURRING SAVE EVENT BETWEEN NOW AND NEXT TIME STEP
        if (num_next > saveTimeCount_ )
        {
            // FIND EXACT TIME OF NEXT SAVE EVENT
            double t_next = (saveTimeCount_ + 1 ) * saveTime_;

            // THIS SHOULD NEVER HAPPEN! FORCE REOCCURRING SAVE TIME DUE TO FLOATING POINT ACCURACY PROBLEMS
            if ( t_next <= simu.getCurrentTime() )
            {
                t_next =  simu.getCurrentTime() + simu.getMindt();
                printf("Warning in %s, floating point numerical accuracy in determining reoccurring save time\n", __func__);
            }

            this->insertTimeEvent( new TimeEvent(EVENT_SAVE, t_next) );
            saveTimeCount_++;
        }
    }

// REOCCURRING REFINEMENT ITERATION

}

void EventList::prependReoccurringIterEvent(IterEvent *iEvent)
{
// ADDS NEW EVENT AT FRONT OF QUEUE

    // MAKE SURE ADDED EVENT REALLY BELONGS TO FRONT OF QUEUE
    if (nextIterEvent_ < iEvent->getEventIteration() )
    {
        printf("error in %s. Trying to add event with iteration number %u to front of queue when first existing event has number %u - bye!\n",
               __func__,
               iEvent->getEventIteration(),
               nextIterEvent_ );
        exit(1);
    }

    iterationEvents_.push_front( iEvent );
    nextIterEvent_ = iEvent->getEventIteration();
}


double EventList::timeUntilNextEvent(const Simu &simu) const
{
// RETURNS AMOUNT OF TIME UNTIL THE OCCURRENCE OF NEXT
// TIME EVENT, OR THE MAGIC NUMBER "NO_TIME_EVENTS" IF
// NO MORE TIME EVENTS EXIST IN QUEUE

    if (nextTimeEvent_ != NO_TIME_EVENTS)
        return nextTimeEvent_ - simu.getCurrentTime();
    else
        return NO_TIME_EVENTS;
}


void EventList::printEventList()const
{
// DEBUG PRINTING

    printf(" %u TimeEvents:\n", timeEvents_.size() );
    std::list<TimeEvent*>::const_iterator tei = timeEvents_.begin();
    for (tei = timeEvents_.begin(); tei != timeEvents_.end() ; tei++)
    {
        EventType et = (*tei)->getEventType();
        switch( et )
        {
        case(EVENT_SWITCHING):
            printf("EVENT_SWITCHING ");
            break;
        case(EVENT_SAVE):
            printf("EVENT_SAVE ");
            break;
        case(EVENT_REFIENEMENT):
            printf("EVENT_REFINEMENT ");
            break;
        default:
            printf("error, unknown envent type in %s - bye\n", __func__);
            exit(1);
        }

        printf("at t = %e\n", (*tei)->getEventTime() );

    }


    // ADD LOOP OVER ITER EVENTS


}


