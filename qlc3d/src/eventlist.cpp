#include "../includes/eventlist.h"
#include <math.h>

Event::Event(double t, int n)
{
	Time = t;
	Number = n;
}

EventList::EventList()
{

}
bool Event::operator<(Event other)
{
	if (Time < other.Time)
		return true;
	else
		return false;
}


void EventList::setElectrodeEvents(Electrodes* electrodes)
{
    for (int i =  0; i < electrodes->nElectrodes ; i ++){
	for (int j = 0 ; j < electrodes->E[i]->getnTimes() ; j++){
	    InsertEvent(electrodes->E[i]->Time[j], EVENT_SWITCHING);
	}
    }
}
void EventList::InsertEvent(double time, int number)
{
	events.push_back( Event(time,number) );
	events.sort();
	
	/*
	list <Event>::iterator i;
	i = events.begin();
	if ( events.size() > 0 ) // find correcet location if not empty list
	{
		while ( (i->Time < time) && (i != events.end() ))
			i++;
	}
	
	events.insert(i,Event(time,number));
	*/

}
void EventList::printEventList()
{
	list<Event>::iterator i;

	printf("%i events \n",(int)events.size());
	
	for ( i = events.begin() ; i != events.end() ; i ++)
	{
		printf("T, E = [%e, %i]\n",  i->Time, i->Number);
	}

}
void EventList::RemoveFirstEvent()
{
	
		events.pop_front();
	
}
double EventList::getNextEventTime()
{
	if (getLength() >0)
		return(events.begin()->Time);//1e-3 * events.begin()->Time );
	else
		return -1;
}
int EventList::getNextEventNumber()
{
	return(events.begin()->Number);
}
int EventList::getLength()
{
	return events.size();

}

