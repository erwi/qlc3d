#ifndef EVENTHANDLER_H
#define EVENTHANDLER_H
#include <solutionvector.h>
#include <electrodes.h>
#include <geometry.h>
#include <eventlist.h>
#include <meshrefinement.h>
#include <lc.h>
#include <resultio.h>
#include <regulargrid.h>
#include <memory>


namespace SpaMtrix {
class IRCMatrix;
}
class SimulationState;
class ResultOutput;
class PotentialSolver;
class SimulationAdaptiveTimeStep;
// CONVENIENCE STRUCT WITH POINTERS TO THE DIFFERENT GEOMETRY OBJECTS
// NEEDED IN MESH REFINEMENT.
struct Geometries {
    Geometry *geom_orig;    // ORIGINAL, LOADED FROM FILE
    Geometry *geom_prev;    // FROM PREVIOUS REFINEMENT ITERATION
    Geometry *geom;         // CURRENT, WORKING GEOMETRY
    Geometries(): geom_orig(NULL), geom_prev(NULL), geom(NULL) {}
};

struct SolutionVectors {
    SolutionVector *q;      // CURRENT Q-TENSOR
    SolutionVector *v;      // POTENTIAL
    SolutionVectors(): q(NULL), v(NULL) {}
};

void handleInitialEvents(SimulationState &simulationState,
                         EventList &eventList,
                         Electrodes &electr,
                         Alignment &alignment,
                         Simu &simu,
                         Geometries &geometries,
                         SolutionVectors &solutionvectors,
                         const LC &lc,
                         ResultOutput &resultOutput,
                         PotentialSolver &potentialSolver,
                         SimulationAdaptiveTimeStep &simulationAdaptiveTimeStep,
                         std::unique_ptr<RegularGrid> &regGrid);

void handleEvents(EventList &evel,
                  Electrodes &electr,
                  Alignment &alignment,
                  Simu &simu,
                  SimulationState &simulationState,
                  Geometries &geometries,
                  SolutionVectors &solutionvectors,
                  const LC &lc,
                  ResultOutput &resultOutput,
                  PotentialSolver &potentialSolver,
                  SimulationAdaptiveTimeStep &adaptiveTimeStep,
                  std::unique_ptr<RegularGrid> &regGrid);

/** return true/false depending on whether mesh was refined or not */
bool handleMeshRefinement(std::list<Event *> &refEvents,
                          Geometries &geometries,
                          SolutionVectors &solutionvectors,
                          Simu &simu,
                          SimulationState &simulationState,
                          Alignment &alignment,
                          Electrodes &electrodes,
                          double S0,
                          std::unique_ptr<RegularGrid> &regGrid);

void handlePreRefinement(std::list<Event *> &refEvents,
                         Geometries &geometries,
                         SolutionVectors &solutionvectors,
                         Simu &simu,
                         SimulationState &simulationState,
                         Alignment &alignment,
                         Electrodes &electrodes,
                         double S0,
                         std::unique_ptr<RegularGrid> &regGrid);
#endif // EVENTHANDLER_H
