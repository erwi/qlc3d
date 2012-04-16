#include <eventlist.h>
#include <math.h>
#include <refinfo.h>
const char* Event::getEventString(const EventType e)
{
    switch(e)
    {
    case(EVENT_SAVE):
        return "EVENT_SAVE";
    case(EVENT_SWITCHING):
        return "EVENT_SWITCHING";
    case(EVENT_REFINEMENT):
        return "EVENT_REFINEMENT";
    default:
        return "EVENT_INVALID";
    }
}

Event::~Event()
{
    // IF THIS EVENT CONTAINS A DATA POINTER, DELETE IT
    if (eventData_)
    {
        switch(eventType_)
        {
        case(EVENT_REFINEMENT):
        {
            RefInfo* ri = static_cast<RefInfo*>(eventData_);
            delete ri;
            eventData_ = NULL;
            break;
        }
        case(EVENT_SWITCHING):
        {
            SwitchingInstance* si = static_cast<SwitchingInstance*>(eventData_);
            delete si;
            eventData_ = NULL;
            break;
        }
        default:    // THIS SHOULD NEVER HAPPEN
            printf("something went wrong in %s - bye!", __func__);
            exit(1);
        }
    }
}

bool compare_timeptr(const Event* first, const Event* second)
{
    // COMPARISON OF POINTERS TO TIME EVENTS
    return ( first->time < second->time );
}

bool compare_iterptr(const Event* first, const Event* second)
{
    return ( first->iteration < second->iteration );
}

/*
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


*/
/*
//-------- IterEvent DEFINITINS ------------
IterEvent::IterEvent(const EventType &et, const size_t &iter):
    Event(et),
    iteration_(iter)
{

}




bool IterEvent::operator <( const IterEvent& other) const
{
    return ( getEventIteration() < other.getEventIteration() );
}
*/
// -------------EventList DEFINITIONS --------------

EventList::EventList():
    NO_TIME_EVENTS(std::numeric_limits<double>::max() ),
    NO_ITER_EVENTS(std::numeric_limits<size_t>::max() ),
    nextTimeEvent_(NO_TIME_EVENTS),
    nextIterEvent_(NO_ITER_EVENTS),
    saveIter_(0),
    saveTime_(0),
    saveTimeCount_(0),
    repRefIter_(0),
    repRefTime_(0),
    repRefTimeCount_(0)
{


}

EventList::~EventList()
{
    // ITERATE THROUGH LISTS DELETING REMAINING EVENT OBJECT POINTERS
    list<Event*> :: iterator itr;

    if (!timeEvents_.empty())
    {
#ifdef DEBUG
        printf("unhandled time Events \n");
#endif

        for (itr = timeEvents_.begin() ; itr != timeEvents_.end(); itr++)
            delete *itr;
    }

    if (!iterationEvents_.empty())
    {
#ifdef DEBUG
        printf("unhandled iteration Events \n");
#endif
        for (itr = iterationEvents_.begin(); itr!= iterationEvents_.end(); itr++)
            delete *itr;
    }

    if (!repRefinements_.empty() )
    {
#ifdef DEBUG
        printf("unhandled repRef Events \n");
#endif



        for (itr = repRefinements_.begin(); itr!= repRefinements_.end(); itr++)
            delete *itr;
    }


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
        if ( timeEvents_.empty() )
        {
            printf("error in %s, no time events in queue - bye!\n", __func__ );
            exit(1);
        }

        Event* eve =  *timeEvents_.begin();   // MAKE A POINTER TO OBJECT ON HEAP
        timeEvents_.erase( timeEvents_.begin() );         // REMOVE EVENT POINTER FROM QUEUE FRONT

        // UPDATE TIME FOR NEXT TIME EVENT MONITOR
        if ( !timeEvents_.empty() )// THIS NEED MODIFYING TO ALLOW FOR REPEATING EVENTS
        {
            nextTimeEvent_ = (*timeEvents_.begin() )->time;
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
        if ( iterationEvents_.empty() )
        {
            printf("error in %s, no iteration event is queue - bye !\n", __func__);
            exit(1);
        }

        Event* eve = *iterationEvents_.begin(); // MAKE A POINTER TO OBJECT ON HEAP
        iterationEvents_.erase(iterationEvents_.begin());          // REMOVE EVENT POINTER FROM QUEUE FRONT

// UPDATE NEXT ITER EVENT MONITOR
        if ( !iterationEvents_.empty() )  // THIS NEED MODIFYING TO ALLOW FOR REPEATING EVENTS
        {
            nextIterEvent_ = ( *iterationEvents_.begin() )->iteration;
        }
        else
        {
            nextIterEvent_ = NO_ITER_EVENTS;
        }

        return eve;
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


void EventList::insertTimeEvent(Event *tEvent)
{
// ADDS TIMED EVENT TO CORRECT POSITION IN EVENT QUEUE

    // MAKE SURE OCCURRENCE IS TIME
    if (tEvent->getEventOccurrence() != EVENT_TIME )
    {
        printf("error in %s, expected 'time' event - bye!\n", __func__);
        exit(1);
    }

    timeEvents_.push_back( tEvent );
    timeEvents_.sort( compare_timeptr );    // USE COMPARISON OF *POINTERS* TO TIME EVENT (OTHERWISE SORTING BY ADDRESS VALUE)

// KEEP NEXT TIME EVENT TRACKER UP TO DATE
    if (tEvent->time < nextTimeEvent_ )
        nextTimeEvent_ = tEvent->time;
}

void EventList::insertIterEvent(Event *iEvent)
{
// ADDS ITERATION EVENT TO CORRECT POSITION IN EVENT QUEUE

    // MAKE SURE OCCURENCE IN ITERATIONS
    if (iEvent->getEventOccurrence() != EVENT_ITERATION)
    {
        printf("error in %s, expected 'itertion' event - bye!\n", __func__);
        exit(1);
    }

    iterationEvents_.push_back( iEvent );
    iterationEvents_.sort(compare_iterptr);

// KEEP NEXT ITERATION EVENT TRACKER UP TO DATE
    if (iEvent->iteration < nextIterEvent_ )
        nextIterEvent_ = iEvent->iteration;
}

void EventList::addRepRefInfo(Event* e)
{
    this->repRefinements_.push_back( e );
}


void EventList::manageReoccurringEvents(Simu& simu)
{
// ADDS REOCCURRING EVENTS TO FRONT OF QUEUE IF NEEDED.
// NEW EVENTS WILL BE SCEDULED FOR NEXT ITERATION/TIME STEP.
// NOTE THAT THIS MUST BE CALLED AFTER ALL EVENTS FOR CURRENT
// ITERATION/TIMESTEP HAVE BEEN PROCESSED.

// ADD ITER EVENTS
    unsigned int nextIter = (unsigned int) simu.getCurrentIteration() + 1;
    double nextTime = simu.getCurrentTime() + simu.getdt(); // TIME AT NEXT ITERATION


// REOCCURRING SAVE ITERATION
    if ( ( saveIter_ > 0 ) &&
         ( nextIter % saveIter_ ) == 0 )
    {
        prependReoccurringIterEvent( new Event(EVENT_SAVE ,  nextIter) );
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

            this->insertTimeEvent( new Event(EVENT_SAVE, t_next) );
            saveTimeCount_++;
        }
    }

// REOCCURRING REFINEMENT ITERATION
    if ( (repRefIter_ > 0 ) &&
         (nextIter % repRefIter_) == 0 )
    {
        //
        // ALL REPEATING REFINEMENT DEFINITIONS NEED TO BE COPIED TO FRONT OF
        // EVENT QUEUE TO BE PROCESSED NEXT TIME STEP
        //
        std::list<Event*> ::iterator itr = repRefinements_.begin();
        for ( ;itr != repRefinements_.end() ; ++itr)
        {
            // GET PTR TO REPEATING REFINFO OBJECT
            RefInfo* repref = static_cast<RefInfo*> ( (*itr)->getEventDataPtr() );
            // CREATE COPY OF REPEATING REFINFO OBJECT
            RefInfo* cpy_repref = new RefInfo( *repref );
            Event* ev = new Event(EVENT_REFINEMENT, nextIter );         // EVENT FOR NEXT ITERATOR
            ev->setEventDataPtr( cpy_repref );
            this->prependReoccurringIterEvent( ev ); // ADD TO FRONT OF QUEUE
        }

    }
// REOCCURRING REFINEMENT TIME
    if  (repRefTime_>0)
    {
        printf("RepRefTime has not been implemented yet in %s\n", __func__ );
        exit(1);
    }



}

void EventList::prependReoccurringIterEvent(Event *iEvent)
{
// ADDS NEW EVENT AT FRONT OF QUEUE

    //MAKE SURE EVENT OCCURRENCE IS ITERATIONS
    if ( iEvent->getEventOccurrence() != EVENT_ITERATION )
    {
        printf("error in %s, event occurrence is not 'iterations' - bye!\n",__func__);
        exit(1);
    }

    // MAKE SURE ADDED EVENT REALLY BELONGS TO FRONT OF QUEUE
    if (nextIterEvent_ < iEvent->iteration )
    {
        printf("error in %s. Trying to add event with iteration number %u to front of queue when first existing event has number %u - bye!\n",
               __func__,
               iEvent->iteration,
               nextIterEvent_ );
        exit(1);
    }

    iterationEvents_.push_front( iEvent );
    nextIterEvent_ = iEvent->iteration;
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
    std::list<Event*>::const_iterator tei = timeEvents_.begin();
    for (tei = timeEvents_.begin(); tei != timeEvents_.end() ; ++tei)
    {
        EventType et = (*tei)->getEventType();
        double time = (*tei)->time;
        printf("Event type %s, at t = %e\n",Event::getEventString(et), time );
    }


    // ADD LOOP OVER ITER EVENTS
    printf(" %u IterEvent:\n", iterationEvents_.size() );
    std::list<Event*>::const_iterator iei = iterationEvents_.begin();
    for(; iei != iterationEvents_.end() ; ++iei)
    {
        EventType et = (*iei)->getEventType();
        unsigned int iter = (*iei)->iteration;
        printf("Event type %s at iteration %u\n", Event::getEventString(et), iter);
    }
}


