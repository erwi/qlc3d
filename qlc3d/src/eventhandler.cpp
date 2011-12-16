#include <eventhandler.h>
#include <filesysfun.h>
#include <qlc3d.h>


void handleElectrodeSwitching(Event* currentEvent,
                              Electrodes& electr,
                              SolutionVector& v,
                              Simu& simu )
{
// SWITCHES ELECTRODE

    // CONVERT GENERIC EVENT TO A SWITCHING EVENT
    SwitchingEvent *se = static_cast<SwitchingEvent*>(currentEvent);

    if(!se)
    {
        printf("error in %s, NULL pointer received\n", __func__);
        exit(1);
    }

// IF TIME-STEPPING REDUCE STEP SIZE TO MINIMUM
    if ( simu.getdt()!= 0.0 )
        simu.setdt( simu.getMindt() );


    // FIND WHICH ELECTRODE IS SWITCHING, AND TO WHAT VALUE
    size_t electrodeNum = se->getElectrodeNumber();
    double potential    = se->getElectrodePotential();
    delete se;
    // SET THE NEW ELECTRODE VALUE
    electr.setElectrodePotential( electrodeNum, potential);

    // SET POTENTIAL BOUNDARY CONDITIONS FOR ALL ELECTRODES
    v.setFixedNodesPot(&electr );
    v.setToFixedValues();
}

void handleResultOutput(Simu& simu,
                        LC& lc,
                        Geometry& geom,
                        SolutionVector& v,
                        SolutionVector& q)
{
// TAKES CARE OF CALLING APPROPRIATE RESULT OUTPUT FUNCTION


    FilesysFun::setCurrentDirectory( simu.getSaveDir() ); // GOTO OUPUT DIR
    double* director(NULL);
    //double* vReg(NULL);


    if ( simu.getSaveFormat() & Simu::LCview )
    {
        printf("LCview \n");
        WriteResults::WriteResult(&simu, &lc, &geom, &v, &q);
    }
    if ( simu.getSaveFormat() & Simu::RegularVTK )
    {
        printf("VTK GRID\n");

        std::stringstream ss;
        std::string filename;
        ss << "regular"<<simu.getCurrentIteration() << ".vtk";
        ss >> filename;
        if (!director)
            director = tensortovector( q.Values, geom.getnpLC() );

        RegularGrid& rGrid = *geom.regularGrid;
        rGrid.writeVTKGrid( filename.c_str() ,
                                        v.Values,
                                        director,
                                        geom.getnpLC() );


    }


// CLEANUP AFTER ALL SAVING HAS BEEN DONE
    if (director) delete [] director;
    FilesysFun::setCurrentDirectory( simu.getCurrentDir() ); // GOTO EXECUTION DIR

}



void handleInitialEvents(EventList& evel,      // EVENT LIST
                  Electrodes& electr,   // ELECTRODES WITH POTENTIALS AND TIMING
                  Simu& simu,           // VARIOUS SIMU SETTINGS
                  SolutionVector& v,    // POTENTIAL SOLUTION
                  Geometry& geom,       // CURRENT MESH
                  MeshRefinement& ref,  // MESH REFINEMENT INFO
                  SparseMatrix* Kpot,   // POTENTIAL CALCULATION MATRIX
                  SolutionVector& q,    // Q-TENSOR
                  LC& lc,               // MATERIAL PARAMS.
                  Settings& settings)   // SPARSE SOLVER SETTINGS
{
// THIS IS ONLY CALLED BEFORE SIMULATION STARTS, DOES NOT
// NEED TO BE AS GENERAL AS handleEvents.
// TAKES CARE OF:
//      PRE-REFINEMENT
//      CALCULATING INITIAL POTENTIALS
//      OUTPUT result_initial FILE

// IF NEEDS PRE-REFINEMENT. DO IT FIRST


    while ( evel.eventOccursNow(simu) )
    {
        Event* currentEvent = evel.getCurrentEvent( simu ); // removes event from queue to be processed
        EventType et = currentEvent->getEventType();
        // REMOVE EVENT FROM LIST AND GET ITS TYPE
        ///EventType et = evel.popCurrentEvent( simu );

        // DEPENDING ON EVENT TYPE, DO STUFF
        switch (et)
        {
        case(EVENT_SAVE): // SAVE RESULTS
            break;
        case(EVENT_SWITCHING):  // SWITCH ELECTRODES
            handleElectrodeSwitching(currentEvent, electr, v, simu );
            break;
        case(EVENT_REFIENEMENT): // REFINE MESH - DON'T
            break;
        default:
            printf("error in %s, unknown event type - bye !\n", __func__);
            exit(1);
        }
    }

// ALWAYS CALCULATE INITIAL POTENTIAL
    calcpot3d( Kpot,
               &v,
               &q,
               &lc,
               geom,
               &settings,
               &electr);

// WRITE INITIAL RESULT FILE
    FilesysFun::setCurrentDirectory( simu.getSaveDir() );
    WriteResults::WriteResult(&simu, &lc, &geom, &v, &q );
    FilesysFun::setCurrentDirectory( simu.getCurrentDir() );

// ADD REOCCURRING EVENTS
    evel.manageReoccurringEvents(simu);

}

void handleEvents(EventList& evel,      // EVENT LIST
                  Electrodes& electr,   // ELECTRODES WITH POTENTIALS AND TIMING
                  Simu& simu,           // VARIOUS SIMU SETTINGS
                  SolutionVector& v,    // POTENTIAL SOLUTION
                  Geometry& geom,       // CURRENT MESH
                  MeshRefinement& ref,  // MESH REFINEMENT INFO
                  SparseMatrix* Kpot,   // POTENTIAL CALCULATION MATRIX
                  SolutionVector& q,    // Q-TENSOR
                  LC& lc,               // MATERIAL PARAMS.
                  Settings& settings)   // SPARSE SOLVER SETTINGS
{

// LEAVE IF NO EVENTS LEFT IN QUEUE
    if ( !evel.eventsInQueue() )    // event queue is empty
    {
        evel.manageReoccurringEvents( simu );
        return ;
    }


// EVENTS ARE ORDERED BY TIME/ITERATION NUMBER,
// BUT NOT ACCIRDING TO TYPE. HOWEVER, DIFFERENT
// EVENT TYPES SHOULD ALSO BE EXECUTED IN PARTICULAR ORDER
// (e.g. UPDATE POTENTIAL VALUES BEFORE SAVING NEW RESULT FILE)
// USE FOLLOWING FLAGS TO DETERMINE THIS
    bool recalculatePotential = false;
    bool saveResult = false;
    bool refineMesh =false;


// CHECK WHICH EVENTS ARE OCCURRING *NOW* AND SET CORRESPONFING
// FLAGS + OTHER PRE-EVENT PROCESSING
    while ( evel.eventOccursNow(simu) )
    {

        Event* currentEvent = evel.getCurrentEvent( simu ); // removes event from queue to be processed
        EventType et = currentEvent->getEventType();
        // REMOVE EVENT FROM LIST AND GET ITS TYPE
        ///EventType et = evel.popCurrentEvent( simu );

        // DEPENDING ON EVENT TYPE, DO STUFF
        switch (et)
        {
        case(EVENT_SAVE): // SAVE RESULTS

            saveResult = true;
            break;
        case(EVENT_SWITCHING):  // SWITCH ELECTRODES
            handleElectrodeSwitching(currentEvent, electr, v, simu );
            recalculatePotential = true;
            break;
        case(EVENT_REFIENEMENT): // REFINE MESH
            refineMesh = true;
            recalculatePotential = true;
            saveResult = true;
            break;
        default:
            printf("error in %s, unknown event type - bye !\n", __func__);
            exit(1);
        }
    }

// ADDS REOCCURRING EVENTS TO QUEUEU FOR NEXT ITERATION
    evel.manageReoccurringEvents(simu);

// IF MESH REFINEMENT
    if (refineMesh) {}

// IF ELECTRODE POTENTIALS HAVE CHANGED, POTENTIALS MUST BE RECALCULATED
    if ( recalculatePotential )
        calcpot3d( Kpot,
                   &v,
                   &q,
                   &lc,
                   geom,
                   &settings,
                   &electr);



// IF RESULT OUTPUT IS NEEDED - THIS SHOUL BE DONE LAST
    if (saveResult)
    {
        handleResultOutput(simu, lc, geom, v, q);
    }


    // IF TIME-STEPPING, REDUCE dt IF IT IS LARGER THAN
    // TIME UNTIL NEXT EVENT
    double tNext = evel.timeUntilNextEvent( simu );
    if ( tNext < simu.getdt() )
        simu.dt = tNext;
    else if (tNext < 0 )
    {
        printf("error , event missed - bye!\n");
        exit(1);
    }


}//end void HandleEvents


