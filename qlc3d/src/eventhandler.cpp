#include <eventhandler.h>
#include <simulation-state.h>
#include <simulation-time.h>
#include <simulation-adaptive-time-step.h>
#include <filesystem>
#include <util/logging.h>
#include <util/exception.h>
#include <io/result-output.h>
#include <qlc3d.h>
#include <spamtrix_ircmatrix.hpp>
#include <potential/potential-solver.h>

void handleElectrodeSwitching(Event *currentEvent,
                              Electrodes &electr,
                              SolutionVector &v,
                              Simu &simu,
                              SimulationState &simulationState,
                              const Geometry &geom) {
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

    // SET POTENTIAL BOUNDARY CONDITIONS FOR ALL ELECTRODES
    auto potentialsByElectrode = electr.getCurrentPotentials(simulationState.currentTime().getTime());
  v.initialisePotentialBoundaries(geom, potentialsByElectrode);
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
                         SpaMtrix::IRCMatrix &Kq,
                         ResultOutput &resultOutput,
                         PotentialSolver &potentialSolver,
                         SimulationAdaptiveTimeStep &adaptiveTimeStep) {

    int currentIteration = simulationState.currentIteration();
    double timeStep = simulationState.dt();
// IF NEEDS PRE-REFINEMENT. DO IT FIRST
    bool refineMesh = false;

    std::list<Event *> refEvents;    // REFINEMENT EVENTS EXECUTED TOGETHER
    while (evel.eventOccursNow(simulationState)) {
        Event *currentEvent = evel.getCurrentEvent(simulationState); // removes event from queue to be processed
        EventType et = currentEvent->getEventType();

        switch (et) {
            case (EVENT_SAVE):
                delete currentEvent;
                break;
            case (EVENT_SWITCHING):
                handleElectrodeSwitching(currentEvent,
                                         electrodes,
                                         *solutionvectors.v,
                                         simu,
                                         simulationState,
                                         *geometries.geom);
                delete currentEvent;
                break;
            case (EVENT_REFINEMENT):
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
                            Kq); // defined in refinementhandler.cpp
    }

    // ALWAYS CALCULATE INITIAL POTENTIAL
    potentialSolver.solvePotential(*solutionvectors.v, *solutionvectors.q, *geometries.geom);

    Log::info("Writing initial results");
    resultOutput.writeResults(*geometries.geom, *solutionvectors.v, *solutionvectors.q, simulationState);

    // ADD REOCCURRING EVENTS
    evel.manageReoccurringEvents(currentIteration, simulationState.currentTime(), timeStep);
    adaptiveTimeStep.setNextEventTime(evel.nextEventTime());
}


void handleEvents(EventList &evel,
                  Electrodes &electrodes,
                  Alignment &alignment,
                  Simu &simu,
                  SimulationState& simulationState, // non-const, dt may be modified
                  Geometries &geometries,
                  SolutionVectors &solutionvectors,
                  const LC &lc,
                  SpaMtrix::IRCMatrix &Kq,
                  ResultOutput &resultOutput,
                  PotentialSolver &potentialSolver,
                  SimulationAdaptiveTimeStep &adaptiveTimeStep) {

  SimulationTime nextEventTime = evel.nextEventTime();
  if (simulationState.currentTime().greaterThan(nextEventTime)) {
    RUNTIME_ERROR(fmt::format("nextEventTime={} is in the past, currentTime={}", nextEventTime.getTime(), simulationState.currentTime().getTime()));
  }
    // EVENTS ARE ORDERED BY TIME/ITERATION NUMBER,
    // BUT NOT ACCORDING TO TYPE. HOWEVER, DIFFERENT
    // EVENT TYPES SHOULD ALSO BE EXECUTED IN PARTICULAR ORDER
    // (e.g. UPDATE POTENTIAL VALUES BEFORE SAVING NEW RESULT FILE)
    // USE FOLLOWING FLAGS TO DETERMINE THIS
    bool recalculatePotential = false;
    bool saveResult = false;
    bool needsMeshRefinement = false;

    // CHECK WHICH EVENTS ARE OCCURRING *NOW* AND SET CORRESPONDING
    // FLAGS + OTHER PRE-EVENT PROCESSING
    std::list<Event *> refEvents; // STORES REF-EVENTS THAT NEED TO BE EXECUTED

   while (evel.eventOccursNow(simulationState)) {
     Event *currentEvent = evel.getCurrentEvent(simulationState); // removes event from queue to be processed
     EventType et = currentEvent->getEventType();

     switch (et) {
       case (EVENT_SAVE):
         saveResult = true;
         delete currentEvent;
         break;
       case (SAVE_PERIODICALLY_BY_TIME):
         saveResult = true;
         delete currentEvent;
         break;
       case (EVENT_SWITCHING):
         handleElectrodeSwitching(currentEvent, electrodes, *solutionvectors.v, simu, simulationState, *geometries.geom);
         delete currentEvent;
         recalculatePotential = true;

         if ((evel.getSaveIter() > 1) || (evel.getSaveTime() > 0)) { // OUTPUT RESULT ON SWITCHING ITERATION
           saveResult = true;
         }

         break;
       case (EVENT_REFINEMENT):
         needsMeshRefinement = true;
         recalculatePotential = true;
         saveResult = true;
         refEvents.push_back(currentEvent);
         break;
       default:
         RUNTIME_ERROR(fmt::format("Unknown event type {}", et));
     }
   }

// ADDS REOCCURRING EVENTS TO QUEUE FOR NEXT ITERATION
    evel.manageReoccurringEvents(simulationState.currentIteration(),
                                 simulationState.currentTime(),
                                 simulationState.dt());

    if (needsMeshRefinement) {
        bool didRefineMesh = handleMeshRefinement(refEvents,
                             geometries,
                             solutionvectors,
                             simu,
                             simulationState,
                             alignment,
                             electrodes,
                             lc.S0(),
                             Kq); // defined in refinementhandler.cpp
      if (didRefineMesh) {
        potentialSolver.onGeometryChanged();
        simulationState.restrictedTimeStep(true);
      }
    }

    if (recalculatePotential) {
      potentialSolver.solvePotential(*solutionvectors.v, *solutionvectors.q, *geometries.geom);
    }

    if (saveResult) {
      resultOutput.writeResults(*geometries.geom, *solutionvectors.v, *solutionvectors.q, simulationState);
    }

    adaptiveTimeStep.setNextEventTime(evel.nextEventTime());
}


