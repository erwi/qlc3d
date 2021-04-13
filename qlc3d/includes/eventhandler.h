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
#include <refinfo.h>


// SPAMTRIX FORWARD DECLARATIONS
namespace SpaMtrix {
class IRCMatrix;
}
class SimulationState;
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
                         LC &lc,
                         Settings &settings,
                         SpaMtrix::IRCMatrix &Kpot,
                         SpaMtrix::IRCMatrix &Kq);

void handleEvents(EventList &evel,
                  Electrodes &electr,
                  Alignment &alignment,
                  Simu &simu,
                  SimulationState &simulationState,
                  Geometries &geometries,
                  SolutionVectors &solutionvectors,
                  LC &lc,
                  Settings &settings,
                  SpaMtrix::IRCMatrix &Kpot,
                  SpaMtrix::IRCMatrix &Kq);

void handleMeshRefinement(std::list<Event *> &refEvents,
                          Geometries &geometries,
                          SolutionVectors &solutionvectors,
                          Simu &simu,
                          SimulationState &simulationState,
                          Alignment &alignment,
                          Electrodes &electrodes,
                          LC &lc,
                          SpaMtrix::IRCMatrix &Kpot,
                          SpaMtrix::IRCMatrix &Kq);

void handlePreRefinement(std::list<Event *> &refEvents,
                         Geometries &geometries,
                         SolutionVectors &solutionvectors,
                         Simu &simu,
                         SimulationState &simulationState,
                         Alignment &alignment,
                         Electrodes &electrodes,
                         LC &lc,
                         SpaMtrix::IRCMatrix &Kpot,
                         SpaMtrix::IRCMatrix &Kq);
#endif // EVENTHANDLER_H
