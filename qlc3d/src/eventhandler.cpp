#include <eventhandler.h>
#include <filesysfun.h>
#include <qlc3d.h>
#include <simulation-state.h>
#include <filesystem>
#include <util/logging.h>
#include <io/result-output.h>

void handleElectrodeSwitching(Event *currentEvent,
                              Electrodes &electr,
                              SolutionVector &v,
                              Simu &simu,
                              SimulationState &simulationState) {
    // SWITCHES ELECTRODE
    // GET SWITCHING EVENT DATA
    SwitchingInstance *si = static_cast<SwitchingInstance *> ( currentEvent->getEventDataPtr());
    assert(si != nullptr);

    // If time stepping, reduce step size to minimum
    if (simu.simulationMode() == TimeStepping) {
        simulationState.dt(simu.getMindt());
    }

    // IF SWITCHING INSTANCE IS A FLAG FOR UNIFORM ELECTRIC FIELD, CAN EXIT
    if (si->electrodeNumber == SwitchingInstance::UNIFORM_E_FIELD) {
        return;
    }

    // SET THE NEW ELECTRODE VALUE. FIRST CHECK THAT
    electr.setElectrodePotential(si->electrodeNumber, si->potential);

    // SET POTENTIAL BOUNDARY CONDITIONS FOR ALL ELECTRODES
    v.setFixedNodesPot(&electr);
    v.setToFixedValues();
}

/*
void handleResultOutput(SimulationState &simulationState,
                        Simu &simu,
                        double S0,
                        Geometry &geom,
                        SolutionVector &v,
                        SolutionVector &q) {
    int currentIteration = simulationState.currentIteration();
    double currentTime = simulationState.currentTime();
    std::string currentDirectory = std::filesystem::current_path().c_str();
    // Change to result output directory
    FilesysFun::setCurrentDirectory(simu.getSaveDir());
    double *director(NULL); // TODO: use vector
    set<Simu::SaveFormats> saveFormats = simu.getSaveFormat();

    // determine the iteration counter part of output filename. TODO: use libformat for this
    char numberChar[9];
    sprintf(numberChar, "%08d", simulationState.currentIteration());
    string fileNameNumber(numberChar);
    if (simulationState.state() == RunningState::COMPLETED) { // after completion, output filename with special counter value "final"
        fileNameNumber = "-final";
    }

    if (saveFormats.count(Simu::LCview) || saveFormats.count(Simu::LCviewTXT)) {
        // calculate mesh name with mesh number appended. e.g. mesh.txt -> mesh0.txt
        std::string numberedMeshName(simu.meshName());
        std::stringstream ss;
        ss << simulationState.meshNumber();
        std::string num;
        ss >> num;
        size_t pos = numberedMeshName.find_last_of("."); // position of separator point
        numberedMeshName.insert(pos, num);

        // write the mesh if it has changed since last time
        // TODO: should this be done when the mesh changes instead?
        if (simulationState.meshModified()) {
            ResultIO::writeMesh(geom.getPtrTop(), geom.t, geom.e, geom.getnp(), numberedMeshName);
            simulationState.meshModified(false); // false so we don't need to output mesh until it is modified again
        }

        // Write the actual result file
        if (saveFormats.count(Simu::LCview)) {
            ResultIO::writeLCD_B(geom.getPtrTop(), geom.t, geom.e, &v, &q, currentIteration, currentTime, S0,
                                 numberedMeshName, simu.getSaveDirAbsolutePath());
        }
        if (saveFormats.count(Simu::LCviewTXT)) {
            ResultIO::writeLCD_T(geom.getPtrTop(), geom.t, geom.e, &v, &q, currentIteration, currentTime,
                                 numberedMeshName, simu.getSaveDirAbsolutePath());
        }
    }

    if (saveFormats.count(Simu::RegularVTK)) {
        std::string filename = "regularvtk" + fileNameNumber + ".vtk";
        if (!director) {
            director = tensortovector(q.Values, geom.getnpLC());
        }

        RegularGrid &rGrid = *geom.regularGrid;
        rGrid.writeVTKGrid(filename.c_str(),
                           v.Values,
                           director,
                           geom.getnpLC());
    }

    if (saveFormats.count(Simu::RegularVecMat)) {
        std::string filename = "regularvec" + fileNameNumber + ".m";
        if (!director) {
            director = tensortovector(q.Values, geom.getnpLC());
        }

        RegularGrid &rGrid = *geom.regularGrid;
        rGrid.writeVecMat(filename.c_str(),       // WRITE REGULAR GRID RESULT FILE
                          v.Values,
                          director,
                          geom.getnpLC(),
                          currentTime);
    }

    if (saveFormats.count(Simu::DirStackZ)) {
        std::string filename = "dirstacksz" + fileNameNumber + ".csv";
        if (!director) {
            director = tensortovector(q.Values, geom.getnpLC());
        }
        RegularGrid &rGrid = *geom.regularGrid;
        rGrid.writeDirStackZ(filename.c_str(),
                             director,
                             geom.getnpLC(),
                             currentTime);
    }

    if (saveFormats.count(Simu::CsvUnstructured)) {
        std::string filename = "unstructured.csv." + std::to_string(simulationState.currentIteration());
        ResultIO::writeCsvUnstructured(geom.getPtrTop(), v, q, filename);
    }

    if (saveFormats.count(Simu::VTKUnstructuredAsciiGrid)) {
        std::string fileName = "unstructured" + fileNameNumber + ".vtk";

        if (!director) {
            director = tensortovector(q.Values, geom.getnpLC());
        }

        ResultIO::writeVtkUnstructuredAsciiGrid(
                geom.getPtrTop(), geom.getnp(), geom.getnpLC(), geom.getTetrahedra(), v, q, fileName);
    }

    if (director) { // TODO: use vector
        delete[] director;
    }
    FilesysFun::setCurrentDirectory(currentDirectory); // Go back to execution directory
}


/**
 * THIS IS ONLY CALLED BEFORE SIMULATION STARTS, DOES NOT
 * NEED TO BE AS GENERAL AS handleEvents.
 * TAKES CARE OF:
 *     PRE-REFINEMENT
 *     CALCULATING INITIAL POTENTIALS
 *     OUTPUT result_initial FILE
 */
void handleInitialEvents(SimulationState &simulationState, // non-const since dt may change_
                         EventList &evel,
                         Electrodes &electrodes,
                         Alignment &alignment,
                         Simu &simu,
                         Geometries &geometries,
                         SolutionVectors &solutionvectors,
                         const LC &lc,
                         Settings &settings,
                         SpaMtrix::IRCMatrix &Kpot,
                         SpaMtrix::IRCMatrix &Kq,
                         ResultOutput &resultOutput) {

    int currentIteration = simulationState.currentIteration();
    double currentTime = simulationState.currentTime();
    double timeStep = simulationState.dt();
// IF NEEDS PRE-REFINEMENT. DO IT FIRST
    bool refineMesh = false;

    std::list<Event *> refEvents;    // REFINEMENT EVENTS EXECUTED TOGETHER
    while (evel.eventOccursNow(simulationState)) {
        Event *currentEvent = evel.getCurrentEvent(simulationState); // removes event from queue to be processed
        EventType et = currentEvent->getEventType();
        // REMOVE EVENT FROM LIST AND GET ITS TYPE
        ///EventType et = evel.popCurrentEvent( simu );

        // DEPENDING ON EVENT TYPE, DO STUFF
        switch (et) {
            case (EVENT_SAVE): // INITIAL RESULT IS ALWAYS WRITTEN. SEE BELOW
                delete currentEvent;
                break;
            case (EVENT_SWITCHING):  // SWITCH ELECTRODES
                handleElectrodeSwitching(currentEvent,
                                         electrodes,
                                         *solutionvectors.v,
                                         simu,
                                         simulationState);
                delete currentEvent; // NOT NEEDED ANYMORE
                break;
            case (EVENT_REFINEMENT): // REFINE MESH
                refEvents.push_back(currentEvent);
                refineMesh = true;
                break;
            default:
                throw std::runtime_error(fmt::format("Unknown event type in {}, {}.", __FILE__, __func__));
        }
    }
    if (refineMesh) {
        handlePreRefinement(refEvents,
                            geometries,
                            solutionvectors,
                            simu,
                            simulationState,
                            alignment,
                            electrodes,
                            lc.S0(),
                            Kpot,
                            Kq); // defined in refinementhandler.cpp
    }

    // ALWAYS CALCULATE INITIAL POTENTIAL
    calcpot3d(Kpot,
              solutionvectors.v,
              solutionvectors.q,
              lc,
              *geometries.geom,
              &settings,
              &electrodes);

    // WRITE INITIAL RESULT FILE. ALWAYS!
    /*
    handleResultOutput(simulationState,
                       simu,
                       lc.S0(),
                       *geometries.geom,
                       *solutionvectors.v,
                       *solutionvectors.q);
                       */
    Log::info("Writing initial results");
    resultOutput.writeResults(*geometries.geom, *solutionvectors.v, *solutionvectors.q, simulationState);

    // ADD REOCCURRING EVENTS
    evel.manageReoccurringEvents(currentIteration, currentTime, timeStep);
}

/*!
 * Reduces time step if necessary so that next time step does not skip over next event.
 */
void reduceTimeStep(SimulationState &simulationState, EventList &eventList) {
    double dt = simulationState.dt();
    double currentTime = simulationState.currentTime();
    assert(dt > 0);

    // FIND TIME UNTIL NEXT EVENT
    double tNext = eventList.timeUntilNextEvent(currentTime);
    assert(tNext > 0); // next time step should be in future

    if (tNext < dt) {
        simulationState.dt(tNext);
    }
}

void handleEvents(EventList &evel,
                  Electrodes &electrodes,
                  Alignment &alignment,
                  Simu &simu,
                  SimulationState& simulationState, // non-const, dt may be modified
                  Geometries &geometries,
                  SolutionVectors &solutionvectors,
                  const LC &lc,
                  Settings &settings,
                  SpaMtrix::IRCMatrix &Kpot,
                  SpaMtrix::IRCMatrix &Kq,
                  ResultOutput &resultOutput
) {
    //int currentIteration = simulationState.currentIteration();

// LEAVE IF NO EVENTS LEFT IN QUEUE
    if (!evel.eventsInQueue()) {    // event queue is empty
        evel.manageReoccurringEvents(simulationState.currentIteration(),
                                     simulationState.currentTime(),
                                     simulationState.dt());
        if (simu.simulationMode() == TimeStepping) {
            reduceTimeStep(simulationState, evel);
        }
        return;
    }

    // EVENTS ARE ORDERED BY TIME/ITERATION NUMBER,
    // BUT NOT ACCORDING TO TYPE. HOWEVER, DIFFERENT
    // EVENT TYPES SHOULD ALSO BE EXECUTED IN PARTICULAR ORDER
    // (e.g. UPDATE POTENTIAL VALUES BEFORE SAVING NEW RESULT FILE)
    // USE FOLLOWING FLAGS TO DETERMINE THIS
    bool recalculatePotential = false;
    bool saveResult = false;
    bool refineMesh = false;

    // CHECK WHICH EVENTS ARE OCCURRING *NOW* AND SET CORRESPONDING
    // FLAGS + OTHER PRE-EVENT PROCESSING
    std::list<Event *> refEvents; // STORES REF-EVENTS THAT NEED TO BE EXECUTED

    while (evel.eventOccursNow(simulationState)) {
        // REMOVE EVENT FROM LIST AND GET ITS TYPE
        Event *currentEvent = evel.getCurrentEvent(simulationState); // removes event from queue to be processed
        EventType et = currentEvent->getEventType();

        // DEPENDING ON EVENT TYPE, DO STUFF
        switch (et) {
            case (EVENT_SAVE): // SAVE RESULTS
                saveResult = true;
                delete currentEvent; // NOT NEEDED ANYMORE
                break;
            case (EVENT_SWITCHING):  // SWITCH ELECTRODES
                handleElectrodeSwitching(currentEvent, electrodes, *solutionvectors.v, simu, simulationState);
                delete currentEvent; // NOT NEEDED ANYMORE
                recalculatePotential = true;

                if ((evel.getSaveIter() > 1) || (evel.getSaveTime() > 0)) // OUTPUT RESULT ON SWITCHING ITERATION
                    saveResult = true;

                break;
            case (EVENT_REFINEMENT): // REFINE MESH
                refineMesh = true;
                recalculatePotential = true;
                saveResult = true;
                refEvents.push_back(currentEvent);
                break;
            default:
                throw std::runtime_error(fmt::format("Unknown event type in {}, {}", __FILE__, __func__));
        }
    }

// ADDS REOCCURRING EVENTS TO QUEUE FOR NEXT ITERATION
    evel.manageReoccurringEvents(simulationState.currentIteration(),
                                 simulationState.currentTime(),
                                 simulationState.dt());

// IF MESH REFINEMENT
    if (refineMesh) {
        handleMeshRefinement(refEvents,
                             geometries,
                             solutionvectors,
                             simu,
                             simulationState,
                             alignment,
                             electrodes,
                             lc.S0(),
                             Kpot,
                             Kq); // defined in refinementhandler.cpp
    }

    if (recalculatePotential) {
        calcpot3d(Kpot,
                  solutionvectors.v,
                  solutionvectors.q,
                  lc,
                  *geometries.geom,
                  &settings,
                  &electrodes);
    }

    if (saveResult) {
      resultOutput.writeResults(*geometries.geom, *solutionvectors.v, *solutionvectors.q, simulationState);
    }

    if (simu.simulationMode() == TimeStepping) {
        reduceTimeStep(simulationState, evel);
    }
}//end void HandleEvents


