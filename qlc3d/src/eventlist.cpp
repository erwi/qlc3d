#include <eventlist.h>
#include <cmath>
#include <string>
#include <meshrefinement.h>
#include <simulation-state.h>
#include <util/logging.h>
#include <util/exception.h>

const char *Event::getEventString(const EventType e) {
  switch (e) {
    case (EVENT_SAVE):
      return "EVENT_SAVE";
    case (EVENT_SWITCHING):
      return "EVENT_SWITCHING";
    case (EVENT_REFINEMENT):
      return "EVENT_REFINEMENT";
    case (SAVE_PERIODICALLY_BY_TIME):
      return "SAVE_PERIODICALLY_BY_TIME";
    default:
      RUNTIME_ERROR(fmt::format("Unknown event type {}", e));
  }
}

Event::~Event() noexcept(false) {
    // IF THIS EVENT CONTAINS A DATA POINTER, DELETE IT
    if (eventData_) {
        switch (eventType_) {
        case (EVENT_REFINEMENT): {
            RefInfo *ri = static_cast<RefInfo *>(eventData_);
            delete ri;
            eventData_ = NULL;
            break;
        }
        case (EVENT_SWITCHING): {
            SwitchingInstance *si = static_cast<SwitchingInstance *>(eventData_);
            delete si;
            eventData_ = NULL;
            break;
        }
        default:    // THIS SHOULD NEVER HAPPEN
            RUNTIME_ERROR("Unhandled event type.");
        }
    }
}

Event * Event::ofPeriodicMeshRefinement(RefInfo *refInfo) {
    return new Event(EVENT_REFINEMENT, 0., (void*) refInfo);
}

bool compare_timeptr(const Event *first, const Event *second) {
    // COMPARISON OF POINTERS TO TIME EVENTS
    return (first->time < second->time);
}

bool compare_iterptr(const Event *first, const Event *second) {
    return (first->iteration < second->iteration);
}

std::string Event::toString() const {
    std::string eventString = fmt::format("Event type: {}, ", getEventString(eventType_));
    if (eventOccurrence_ == EVENT_TIME) {
        eventString += fmt::format("Time: {}", time);
    } else {
        eventString += fmt::format("Iteration: {}", iteration);
    }
    return eventString;
}

// -------------EventList DEFINITIONS --------------
EventList::EventList():
    nextTimeEvent_(NO_TIME_EVENTS),
    nextIterEvent_(NO_ITER_EVENTS),
    saveIter_(0),
    saveTime_(0),
    saveTimeCount_(0),
    repRefIter_(0),
    repRefTime_(0),
    repRefTimeCount_(0) {
}

EventList::~EventList() {
    // ITERATE THROUGH LISTS DELETING REMAINING EVENT OBJECT POINTERS
    list<Event *> :: iterator itr;
    if (!timeEvents_.empty()) {
      Log::warn("There are {} unhandled time events in the queue.", timeEvents_.size());
      auto eventString = this->toString();
      Log::warn("EventList: {}", eventString);
      for (itr = timeEvents_.begin(); itr != timeEvents_.end(); itr++) {
          delete *itr;
      }
    }
    if (!iterationEvents_.empty()) {
        Log::warn("There are {} unhandled iteration events in the queue.", iterationEvents_.size());
        for (itr = iterationEvents_.begin(); itr != iterationEvents_.end(); itr++) {
            delete *itr;
        }
    }
    if (!repRefinements_.empty()) {
        Log::warn("There are {} unhandled repeating refinement event in the queue.", repRefinements_.size());
        for (itr = repRefinements_.begin(); itr != repRefinements_.end(); itr++) {
            delete *itr;
        }
    }
}

bool EventList::eventsInQueue() const {
// RETURNS BOOLEAN INDICATING WHETHER ANY EVENTS EXIST IN QUEUE
    return ((nextIterEvent_ < NO_ITER_EVENTS) ||    // ITER EVENTS EXIST?
            (nextTimeEvent_ < NO_TIME_EVENTS)      // TIME EVENTS EXIST?
           );
}

bool EventList::eventOccursNow(const SimulationState &simulationState) const {
    return (nextIterEvent_ == (size_t)simulationState.currentIteration()) || // IF EVENT SCHEDULED FOR CURRENT ITERATION
            simulationState.currentTime().equals(nextTimeEvent_) ;    // IF EVENT SCHEDULED FOR CURRENT TIME
}

Event *EventList::getCurrentEvent(const SimulationState &simulationState) {
// EVENTS EXIST AS OBJECTS ON HEAP. ONLY POINTERS TO THESE
// ARE STORED IN THE EVENT LIST.
// THIS FUNCTION RETURNS A POINTER TO AN EVENT OBJECT
// AND REMOVES THE CORRESPONDING POINTER FROM QUEUE.
// ACTUAL DELETION OF OBJECT MUST BE PERFORMED ELSEWHERE
// FIRST CHECK TIME EVENTS
    if (simulationState.currentTime().equals(nextTimeEvent_)) {
        if (timeEvents_.empty()) {
            RUNTIME_ERROR("Time events queue is empty.");
        }
        Event *eve =  *timeEvents_.begin();   // MAKE A POINTER TO OBJECT ON HEAP
        timeEvents_.erase(timeEvents_.begin());           // REMOVE EVENT POINTER FROM QUEUE FRONT
        // UPDATE TIME FOR NEXT TIME EVENT MONITOR
        if (!timeEvents_.empty()) { // THIS NEED MODIFYING TO ALLOW FOR REPEATING EVENTS
            nextTimeEvent_ = (*timeEvents_.begin())->time;
        } else {
            nextTimeEvent_ = NO_TIME_EVENTS;    // INDICATE END OF QUEUE
        }
        return eve;
    }
// IF NO CURRENT TIME EVENTS LEFT, CHECK ITER EVENTS
    else if (nextIterEvent_ == (size_t) simulationState.currentIteration()) {
        if (iterationEvents_.empty()) {
            RUNTIME_ERROR("Iteration events queue is empty.");
        }
        Event *eve = *iterationEvents_.begin(); // MAKE A POINTER TO OBJECT ON HEAP
        iterationEvents_.erase(iterationEvents_.begin());          // REMOVE EVENT POINTER FROM QUEUE FRONT
// UPDATE NEXT ITER EVENT MONITOR
        if (!iterationEvents_.empty()) {  // THIS NEED MODIFYING TO ALLOW FOR REPEATING EVENTS
            nextIterEvent_ = (*iterationEvents_.begin())->iteration;
        } else {
            nextIterEvent_ = NO_ITER_EVENTS;
        }
        return eve;
    } else {
        // ERROR. SOMETHING HAS GONE WRONG. NO CURRENT EVENT
        RUNTIME_ERROR("No event available currently.");
    }
}

void EventList::setSaveTime(const double &st) {
// CHECK AND SET REOCCURING SAVE EVENT TIME
    if (st >= 0) {
        saveTime_ = st;
    } else {
        RUNTIME_ERROR(fmt::format("Negative save time {}", st));
    }
}

void EventList::insertTimeEvent(Event *tEvent) {
// ADDS TIMED EVENT TO CORRECT POSITION IN EVENT QUEUE
    // MAKE SURE OCCURRENCE IS TIME
    if (tEvent->getEventOccurrence() != EVENT_TIME) {
        RUNTIME_ERROR("Time event expected.");
    }
    timeEvents_.push_back(tEvent);
    timeEvents_.sort(compare_timeptr);      // USE COMPARISON OF *POINTERS* TO TIME EVENT (OTHERWISE SORTING BY ADDRESS VALUE)
// KEEP NEXT TIME EVENT TRACKER UP TO DATE
    if (tEvent->time < nextTimeEvent_)
        nextTimeEvent_ = tEvent->time;
}

void EventList::insertIterEvent(Event *iEvent) {
// ADDS ITERATION EVENT TO CORRECT POSITION IN EVENT QUEUE
    // MAKE SURE OCCURENCE IN ITERATIONS
    if (iEvent->getEventOccurrence() != EVENT_ITERATION) {
        RUNTIME_ERROR("Iteration event expected.");
    }
    iterationEvents_.push_back(iEvent);
    iterationEvents_.sort(compare_iterptr);
// KEEP NEXT ITERATION EVENT TRACKER UP TO DATE
    if (iEvent->iteration < nextIterEvent_)
        nextIterEvent_ = iEvent->iteration;
}

void EventList::addRepRefInfo(Event *e) {
    this->repRefinements_.push_back(e);
}


void EventList::manageReoccurringEvents(int currentIteration, const SimulationTime &currentTime, double timeStep) {
    unsigned int nextIter = currentIteration+ 1;
    //double nextTime = currentTime + timeStep; // TIME AT NEXT ITERATION
    // REOCCURRING SAVE ITERATION
    if (isPeriodicSaveIteration(nextIter)) {
        prependReoccurringIterEvent(new Event(EVENT_SAVE ,  nextIter));
    }
    // REOCCURRING SAVE TIME
    if (saveTime_ > 0) {

        // Find next instance in time when periodic save event should occur
        SimulationTime nextSaveTime((saveTimeCount_ + 1) * saveTime_);

        auto firstSaveEvent = std::find_if(timeEvents_.begin(), timeEvents_.end(),
                                        [](Event *e) { return e->getEventType() == SAVE_PERIODICALLY_BY_TIME; });

        if (firstSaveEvent == timeEvents_.end() || (*firstSaveEvent)->time > nextSaveTime.getTime()) {

          Log::info("Adding save event {} to occur at time {}, currentTime={}", saveTimeCount_, nextSaveTime.getTime(), currentTime);
          insertTimeEvent(new Event(SAVE_PERIODICALLY_BY_TIME, nextSaveTime.getTime()));
          saveTimeCount_++;
        }

        /*
        // PREDICTED NUM SAVE EVENTS BY NEXT TIME STEP, IF dt IS NOT ADJUSTED
        size_t num_next = static_cast<size_t>(nextTime / saveTime_);
        // IF A REOCCURRING SAVE EVENT BETWEEN NOW AND NEXT TIME STEP
        if (num_next > saveTimeCount_) {
            // FIND EXACT TIME OF NEXT SAVE EVENT
            double t_next = (saveTimeCount_ + 1) * saveTime_;
            this->insertTimeEvent(new Event(EVENT_SAVE, t_next));
            saveTimeCount_++;
        }
         */
    }


    // REOCCURRING REFINEMENT ITERATION
    if ((repRefIter_ > 0) &&
            (nextIter % repRefIter_) == 0) {
        //
        // ALL REPEATING REFINEMENT DEFINITIONS NEED TO BE COPIED TO FRONT OF
        // EVENT QUEUE TO BE PROCESSED NEXT TIME STEP
        //
        std::list<Event *> ::iterator itr = repRefinements_.begin();
        for (; itr != repRefinements_.end() ; ++itr) {
            // GET PTR TO REPEATING REFINFO OBJECT
            RefInfo *repref = static_cast<RefInfo *>((*itr)->getEventDataPtr());
            // CREATE COPY OF REPEATING REFINFO OBJECT
            RefInfo *cpy_repref = new RefInfo(*repref);
            Event *ev = new Event(EVENT_REFINEMENT, nextIter);          // EVENT FOR NEXT ITERATOR
            ev->setEventDataPtr(cpy_repref);
            this->prependReoccurringIterEvent(ev);   // ADD TO FRONT OF QUEUE
        }
    }
    // REOCCURRING REFINEMENT TIME
    if (repRefTime_ > 0) {
        RUNTIME_ERROR("\"RepRefTime\" has not been implemented yet.");
    }
}

void EventList::prependReoccurringIterEvent(Event *iEvent) {
    // ADDS NEW EVENT AT FRONT OF QUEUE
    //MAKE SURE EVENT OCCURRENCE IS ITERATIONS
    if (iEvent->getEventOccurrence() != EVENT_ITERATION) {
        RUNTIME_ERROR("Expected event occurrence to be iteration");
    }
    // MAKE SURE ADDED EVENT REALLY BELONGS TO FRONT OF QUEUE
    if (nextIterEvent_ < iEvent->iteration) {
        RUNTIME_ERROR(fmt::format("Trying to add event with iteration number {} to front of queue "
                                  "when first event in queue has number {}. ",
                                  (unsigned int) iEvent->iteration, (unsigned int) nextIterEvent_));
    }
    iterationEvents_.push_front(iEvent);
    nextIterEvent_ = iEvent->iteration;
}

double EventList::timeUntilNextEvent(const double &currentTime) const {
    // RETURNS AMOUNT OF TIME UNTIL THE OCCURRENCE OF NEXT
    // TIME EVENT, OR THE MAGIC NUMBER "NO_TIME_EVENTS" IF
    // NO MORE TIME EVENTS EXIST IN QUEUE
    if (nextTimeEvent_ != NO_TIME_EVENTS) {
        return nextTimeEvent_ - currentTime;
    } else {
        return NO_TIME_EVENTS;
    }
}


SimulationTime EventList::nextEventTime() const {
    return SimulationTime(nextTimeEvent_);
}


std::string EventList::toString() const {
  // append all events to a string
  std::string eventString;
  eventString += "{{";
  eventString += "nextTimeEvent_=" + std::to_string(nextTimeEvent_) + ", nextIterEvent_=" + std::to_string(nextIterEvent_) + ",";
  eventString += "timeEvents_=[";
  for (auto &event : timeEvents_) {
    eventString += event->toString();
  }
  eventString += "], iterationEvents_=[";
  for (auto &event : iterationEvents_) {
    eventString += event->toString();
  }
  // add repeating refinements
  eventString += "], repRefinements_=[";
  for (auto &event : repRefinements_) {
    eventString += event->toString();
  }
  // add periodically occurring events' periods
  eventString += "], saveIter_=" + std::to_string(saveIter_) + ", saveTime_=" + std::to_string(saveTime_) + ", saveTimeCount_=" + std::to_string(saveTimeCount_) + ",";
  eventString += "repRefIter_=" + std::to_string(repRefIter_) + ", repRefTime_=" + std::to_string(repRefTime_) + ", repRefTimeCount_=" + std::to_string(repRefTimeCount_);

  eventString += "}}";
  return eventString;
}