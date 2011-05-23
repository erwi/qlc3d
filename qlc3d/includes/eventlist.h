#ifndef EVENTLIST_H
#define EVENTLIST_H

#define EVENT_SWITCHING 0



#include <stdio.h>
#include <vector>
#include <list>
#include "electrodes.h"

using namespace std;

class Event
{
	public:
		double Time;
		int Number;
		Event(double t, int n);
		
		bool operator<(Event other);
		
};

class EventList
{
	private:
		list <Event> events;

		void InsertEvent(double time, int number);
	
	
	public:
		EventList();
		void setElectrodeEvents(Electrodes* electrodes);
		void PrintEventList();
		double getNextEventTime();
		int    getNextEventNumber();
		int    getLength();
		void RemoveFirstEvent();
};
#endif




