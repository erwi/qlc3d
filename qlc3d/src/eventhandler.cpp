#include <eventhandler.h>
#include <filesysfun.h>
#include <qlc3d.h>


void handleElectrodeSwitching(Event* currentEvent,
                              Electrodes& electr,
                              SolutionVector& v,
                              Simu& simu )
{
// SWITCHES ELECTRODE


    // GET SWITCHING EVENT DATA
    SwitchingInstance* si = static_cast<SwitchingInstance*> ( currentEvent->getEventDataPtr() );
    if(!si)
    {
        printf("error in %s, NULL pointer received\n", __func__);
        exit(1);
    }

// IF TIME-STEPPING REDUCE STEP SIZE TO MINIMUM
    if ( simu.getdt()!= 0.0 )
        simu.setdt( simu.getMindt() );

    // IF SWITCHING INSTANCE IS A FLAG FOR UNIFORM ELECTRIC FIELD, CAN EXIT
    if (  si->electrodeNumber == SwitchingInstance::UNIFORM_E_FIELD )
        return;

    // SET THE NEW ELECTRODE VALUE. FIRST CHECK THAT
    electr.setElectrodePotential( si->electrodeNumber, si->potential);

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
       // printf("LCview \n");
        WriteResults::WriteLCViewResult(&simu, &lc, &geom, &v, &q);
    }
    if ( simu.getSaveFormat() & Simu::RegularVTK )
    {
        //printf("VTK GRID\n");

        std::stringstream ss;
        std::string filename;
        ss << "regularvtk"<<simu.getCurrentIteration() << ".vtk";
        ss >> filename;
        if (!director)
            director = tensortovector( q.Values, geom.getnpLC() );

        RegularGrid& rGrid = *geom.regularGrid;
        rGrid.writeVTKGrid( filename.c_str() ,
                            v.Values,
                            director,
                            geom.getnpLC() );


    }
    if ( simu.getSaveFormat() & Simu::RegularVecMat)
    {
        //printf("REGULAR GRID VECTORS MATLAB\n");
        std::stringstream ss;
        std::string filename;
        ss << "regularvec"<<simu.getCurrentIteration() << ".m";
        ss >> filename;
        if (!director)
            director = tensortovector(q.Values, geom.getnpLC() );

        RegularGrid& rGrid =  *geom.regularGrid;
        rGrid.writeVecMat( filename.c_str() ,       // WRITE REGULAR GRID RESULT FILE
                           v.Values,
                           director,
                           geom.getnpLC(),
                           simu.getCurrentTime());
    }
    if ( simu.getSaveFormat() & Simu::DirStackZ)
    {
        cout << "REGULAR DIR STACKZ" << endl;
        std::stringstream ss;
        std::string filename;
        ss << "dirstackz"<<simu.getCurrentIteration() <<".csv";
        ss >> filename;
        if (!director)
            director = tensortovector(q.Values, geom.getnpLC() );
        RegularGrid &rGrid = *geom.regularGrid;
        rGrid.writeDirStackZ(filename.c_str(),
                             director,
                             geom.getnpLC(),
                             simu.getCurrentTime() );

    }

// CLEANUP AFTER ALL SAVING HAS BEEN DONE
    if (director) delete [] director;
    FilesysFun::setCurrentDirectory( simu.getCurrentDir() ); // GOTO EXECUTION DIR

}



void handleInitialEvents(EventList& evel,          // EVENT LIST
                         Electrodes& electrodes,   // ELECTRODES WITH POTENTIALS AND TIMING
                         Alignment& alignment,     // ANCHORING CONDITIONS
                         Simu& simu,               // VARIOUS SIMU SETTINGS
                         Geometries& geometries,   // POINTERS TO CURRENT MESHES
                         SolutionVectors& solutionvectors, // PTRS TO SOLUTIONS
                         LC& lc,                   // MATERIAL PARAMS.
                         Settings& settings,       // SPARSE SOLVER SETTINGS
                         SpaMtrix::IRCMatrix& Kpot,
                         SpaMtrix::IRCMatrix& Kq
                         )
{
// THIS IS ONLY CALLED BEFORE SIMULATION STARTS, DOES NOT
// NEED TO BE AS GENERAL AS handleEvents.
// TAKES CARE OF:
//      PRE-REFINEMENT
//      CALCULATING INITIAL POTENTIALS
//      OUTPUT result_initial FILE

// IF NEEDS PRE-REFINEMENT. DO IT FIRST
    bool refineMesh = false;

    std::list<Event*> refEvents;    // REFINEMENT EVENTS EXECUTED TOGETHER
    while ( evel.eventOccursNow(simu) )
    {
        Event* currentEvent = evel.getCurrentEvent( simu ); // removes event from queue to be processed
        EventType et = currentEvent->getEventType();
        // REMOVE EVENT FROM LIST AND GET ITS TYPE
        ///EventType et = evel.popCurrentEvent( simu );

        // DEPENDING ON EVENT TYPE, DO STUFF
        switch (et)
        {
        case(EVENT_SAVE): // INITIAL RESULT IS ALWAYS WRITTEN. SEE BELOW
            delete currentEvent;
            break;
        case(EVENT_SWITCHING):  // SWITCH ELECTRODES
            handleElectrodeSwitching(currentEvent,
                                     electrodes,
                                     *solutionvectors.v,
                                     simu );
            delete currentEvent; // NOT NEEDED ANYMORE
            break;
        case(EVENT_REFINEMENT): // REFINE MESH
            refEvents.push_back( currentEvent );
            refineMesh = true;
            break;
        default:
            printf("error in %s, unknown event type - bye !\n", __func__);
            exit(1);
        }
    }
    if (refineMesh){

        handlePreRefinement(refEvents,
                            geometries,
                            solutionvectors,
                            simu,
                            alignment,
                            electrodes,
                            lc,
                            Kpot,
                            Kq); // defined in refinementhandler.cpp

    }
// ALWAYS CALCULATE INITIAL POTENTIAL
    calcpot3d( Kpot,
               solutionvectors.v,
               solutionvectors.q,
               &lc,
               *geometries.geom,
               &settings,
               &electrodes);

// WRITE INITIAL RESULT FILE. ALWAYS!

    handleResultOutput( simu,
                        lc,
                        *geometries.geom,
                        *solutionvectors.v,
                        *solutionvectors.q);

// ADD REOCCURRING EVENTS
    evel.manageReoccurringEvents(simu);

}

void reduceTimeStep(Simu& simu, EventList& evel)
{
// REDUCES TIME STEP SIZE IF NECESSARY, SO THAT NEXT ITERATION
// COINCIDES WITH NEXT TIME EVENT

    if (simu.getdt() == 0 ) return; // ONLY NEEDED WHEN TIME-STEPPING

    // FIND TIME UNTIL NEXT EVENT
    double tNext = evel.timeUntilNextEvent( simu );
    if (tNext < 0 )
    {
        printf("error in %s, event missed - bye %es.!\n", __func__,tNext);
        evel.printEventList();
        exit(1);
    }

    if ( tNext < simu.getdt() )
        simu.setdtForced( tNext );

}

void handleEvents(EventList& evel,          // EVENT LIST
                  Electrodes& electrodes,   // ELECTRODES WITH POTENTIALS AND TIMING
                  Alignment& alignment,     // ANCHORING CONDITIONS
                  Simu& simu,               // VARIOUS SIMU SETTINGS
                  Geometries& geometries,   // POINTERS TO CURRENT MESHES
                  SolutionVectors& solutionvectors, // PTRS TO SOLUTIONS
                  LC& lc,                   // MATERIAL PARAMS.
                  Settings& settings,       // SPARSE SOLVER SETTINGS
                  SpaMtrix::IRCMatrix &Kpot,
                  SpaMtrix::IRCMatrix &Kq
                  )
{


// LEAVE IF NO EVENTS LEFT IN QUEUE
    if ( !evel.eventsInQueue() )    // event queue is empty
    {
        evel.manageReoccurringEvents( simu );
        reduceTimeStep(simu, evel);
        return ;
    }


// EVENTS ARE ORDERED BY TIME/ITERATION NUMBER,
// BUT NOT ACCORDING TO TYPE. HOWEVER, DIFFERENT
// EVENT TYPES SHOULD ALSO BE EXECUTED IN PARTICULAR ORDER
// (e.g. UPDATE POTENTIAL VALUES BEFORE SAVING NEW RESULT FILE)
// USE FOLLOWING FLAGS TO DETERMINE THIS
    bool recalculatePotential = false;
    bool saveResult = false;
    bool refineMesh =false;


// CHECK WHICH EVENTS ARE OCCURRING *NOW* AND SET CORRESPONFING
// FLAGS + OTHER PRE-EVENT PROCESSING

    std::list<Event*> refEvents; // STORES REF-EVENTS THAT NEED TO BE EXECUTED

    while ( evel.eventOccursNow(simu) )
    {
        // REMOVE EVENT FROM LIST AND GET ITS TYPE
        Event* currentEvent = evel.getCurrentEvent( simu ); // removes event from queue to be processed
        EventType et = currentEvent->getEventType();

        // DEPENDING ON EVENT TYPE, DO STUFF
        switch (et)
        {
        case(EVENT_SAVE): // SAVE RESULTS
            saveResult = true;
            delete currentEvent; // NOT NEEDED ANYMORE
            break;
        case(EVENT_SWITCHING):  // SWITCH ELECTRODES
            handleElectrodeSwitching(currentEvent, electrodes, *solutionvectors.v, simu );
            delete currentEvent; // NOT NEEDED ANYMORE
            recalculatePotential = true;

            if ( (evel.getSaveIter() > 1) || (evel.getSaveTime()>0) ) // OUTPUT RESULT ON SWITCHING ITERATION
                saveResult = true;

            break;
        case(EVENT_REFINEMENT): // REFINE MESH
            refineMesh = true;
            recalculatePotential = true;
            saveResult = true;
            refEvents.push_back( currentEvent );
            break;
        default:
            printf("error in %s, unknown event type - bye !\n", __func__);
            exit(1);
        }
    }

// ADDS REOCCURRING EVENTS TO QUEUE FOR NEXT ITERATION
    evel.manageReoccurringEvents(simu);

// IF MESH REFINEMENT
    if (refineMesh)
    {
       //*
        handleMeshRefinement(refEvents,
                             geometries,
                             solutionvectors,
                             simu,
                             alignment,
                             electrodes,
                             lc,
                             Kpot,
                             Kq); // defined in refinementhandler.cpp
    //*/
}

// IF ELECTRODE POTENTIALS HAVE CHANGED, POTENTIALS MUST BE RECALCULATED
    if ( recalculatePotential )
        calcpot3d( Kpot,
                   solutionvectors.v,
                   solutionvectors.q,
                   &lc,
                   *geometries.geom,
                   &settings,
                   &electrodes);



// IF RESULT OUTPUT IS NEEDED - THIS SHOUL BE DONE LAST
    if (saveResult)
    {
        handleResultOutput(simu,
                           lc,
                           *geometries.geom,
                           *solutionvectors.v,
                           *solutionvectors.q);
    }


    // IF TIME-STEPPING, REDUCE dt IF IT IS LARGER THAN
    // TIME UNTIL NEXT EVENT
    reduceTimeStep( simu, evel );

}//end void HandleEvents


