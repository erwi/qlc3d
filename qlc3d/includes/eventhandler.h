#ifndef EVENTHANDLER_H
#define EVENTHANDLER_H
#include <solutionvector.h>
#include <electrodes.h>
#include <geometry.h>
#include <eventlist.h>
#include <simu.h>
#include <meshrefinement.h>
#include <lc.h>
#include <settings.h>
#include <calcpot3d.h>
#include <resultio.h>

// SPAMTRIX FORWARD DECLARATIONS
namespace SpaMtrix {
class IRCMatrix;
}
class SimulationState;
class ResultOutput;
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
                         EventList &evel,
                         Electrodes &electr,
                         Alignment &alignment,
                         Simu &simu,
                         Geometries &geometries,
                         SolutionVectors &solutionvectors,
                         const LC &lc,
                         Settings &settings,
                         SpaMtrix::IRCMatrix &Kpot,
                         SpaMtrix::IRCMatrix &Kq,
                         ResultOutput &resultOutput);

void handleEvents(EventList &evel,
                  Electrodes &electr,
                  Alignment &alignment,
                  Simu &simu,
                  SimulationState &simulationState,
                  Geometries &geometries,
                  SolutionVectors &solutionvectors,
                  const LC &lc,
                  Settings &settings,
                  SpaMtrix::IRCMatrix &Kpot,
                  SpaMtrix::IRCMatrix &Kq,
                  ResultOutput &resultOutput);

void handleMeshRefinement(std::list<Event *> &refEvents,
                          Geometries &geometries,
                          SolutionVectors &solutionvectors,
                          Simu &simu,
                          SimulationState &simulationState,
                          Alignment &alignment,
                          Electrodes &electrodes,
                          double S0,
                          SpaMtrix::IRCMatrix &Kpot,
                          SpaMtrix::IRCMatrix &Kq);

void handlePreRefinement(std::list<Event *> &refEvents,
                         Geometries &geometries,
                         SolutionVectors &solutionvectors,
                         Simu &simu,
                         SimulationState &simulationState,
                         Alignment &alignment,
                         Electrodes &electrodes,
                         double S0,
                         SpaMtrix::IRCMatrix &Kpot,
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
