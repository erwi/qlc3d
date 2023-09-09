#ifndef EVENTHANDLER_H
#define EVENTHANDLER_H
#include <solutionvector.h>
#include <electrodes.h>
#include <geometry.h>
#include <eventlist.h>
#include <simu.h>
#include <meshrefinement.h>
#include <lc.h>
#include <solver-settings.h>
#include <resultio.h>

// SPAMTRIX FORWARD DECLARATIONS
namespace SpaMtrix {
class IRCMatrix;
}
class SimulationState;
class ResultOutput;
class PotentialSolver;
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
    SolutionVector *qn;     // PREVIOUS Q-TENSOR
    SolutionVector *v;      // POTENTIAL
    SolutionVectors(): q(NULL), qn(NULL), v(NULL) {}
};

void setElectrodePotentials(EventList &evel,
                            Electrodes &electr,
                            Simu &simu);

void handleInitialEvents(SimulationState &simulationState,
                         EventList &eventList,
                         Electrodes &electr,
                         Alignment &alignment,
                         Simu &simu,
                         Geometries &geometries,
                         SolutionVectors &solutionvectors,
                         const LC &lc,
                         SpaMtrix::IRCMatrix &Kq,
                         ResultOutput &resultOutput,
                         PotentialSolver &potentialSolver);

void handleEvents(EventList &evel,
                  Electrodes &electr,
                  Alignment &alignment,
                  Simu &simu,
                  SimulationState &simulationState,
                  Geometries &geometries,
                  SolutionVectors &solutionvectors,
                  const LC &lc,
                  SpaMtrix::IRCMatrix &Kq,
                  ResultOutput &resultOutput,
                  PotentialSolver &potentialSolver);

/** return true/false depending on whether mesh was refined or not */
bool handleMeshRefinement(std::list<Event *> &refEvents,
                          Geometries &geometries,
                          SolutionVectors &solutionvectors,
                          Simu &simu,
                          SimulationState &simulationState,
                          Alignment &alignment,
                          Electrodes &electrodes,
                          double S0,
                          SpaMtrix::IRCMatrix &Kq);

void handlePreRefinement(std::list<Event *> &refEvents,
                         Geometries &geometries,
                         SolutionVectors &solutionvectors,
                         Simu &simu,
                         SimulationState &simulationState,
                         Alignment &alignment,
                         Electrodes &electrodes,
                         double S0,
                         SpaMtrix::IRCMatrix &Kq);

/*
void handleResultOutput(SimulationState &simulationState,
                        Simu &simu,
                        double S0,
                        Geometry &geom,
                        SolutionVector &v,
                        SolutionVector &q);
*/
#endif // EVENTHANDLER_H
