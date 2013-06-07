#ifndef EVENTHANDLER_H
#define EVENTHANDLER_H
#include <solutionvector.h>
#include <electrodes.h>
#include <geometry.h>
#include <eventlist.h>
#include <simu.h>
#include <meshrefinement.h>
//#include <sparsematrix.h>
#include <lc.h>
#include <settings.h>
#include <calcpot3d.h>
#include <resultoutput.h>
#include <refinfo.h>


// SPAMTRIX FORWARD DECLARATIONS
namespace SpaMtrix{
    class IRCMatrix;
}


// CONVENIENCE STRUCT WITH POINTERS TO THE DIFFERENT GEOMETRY OBJECTS
// NEEDED IN MESH REFINEMENT.
struct Geometries
{
    Geometry* geom_orig;    // ORIGINAL, LOADED FROM FILE
    Geometry* geom_prev;    // FROM PREVIOUS REFINEMENT ITERATION
    Geometry* geom;         // CURRENT, WORKING GEOMETRY
    Geometries(): geom_orig(NULL), geom_prev(NULL), geom(NULL){}
};

struct SolutionVectors
{
    SolutionVector* q;      // CURRENT Q-TENSOR
    SolutionVector* qn;     // PREVIOUS Q-TENSOR
    SolutionVector* v;      // POTENTIAL
    SolutionVectors(): q(NULL), qn(NULL), v(NULL){}
};


void setElectrodePotentials(EventList& evel,
                            Electrodes& electr,
                            Simu& simu);


// TAKSES CARE OF EVENTS OCCRRING BEFORE SIMULATION STARTS
// AND PREPARING EVENTLIST FOR SIMULATION
void handleInitialEvents(EventList& evel,      // EVENT LIST
                         Electrodes& electr,   // ELECTRODES WITH POTENTIALS AND TIMING
                         Alignment& alignment, // ANCHORING DATA
                         Simu& simu,
                         Geometries& geometries,   // POINTERS TO CURRENT GEOMETRIES
                         SolutionVectors& solutionvectors, // POINTERS TO SOLUTIONS
                         LC& lc,               // MATERIAL PARAMS.
                         Settings& settings,   // SPARSE SOLVER SETTINGS
                         SpaMtrix::IRCMatrix &Kpot,   // MATRIX FOR POTENTIAL CALCULATION
                         SpaMtrix::IRCMatrix &Kq      // MATRIX FOR Q-TENSOR CALCULATION
                         );


void handleEvents(EventList& evel,      // EVENT LIST
                  Electrodes& electr,   // ELECTRODES WITH POTENTIALS AND TIMING
                  Alignment& alignment, // ANCHORING DATA
                  Simu& simu,
                  Geometries& geometries,   // POINTERS TO CURRENT GEOMETRIES
                  SolutionVectors& solutionvectors, // POINTERS TO SOLUTIONS
                  LC& lc,               // MATERIAL PARAMS.
                  Settings& settings,   // SPARSE SOLVER SETTINGS
                  SpaMtrix::IRCMatrix& Kpot,   // MATRIX FOR POTENTIAL CALCULATION
                  SpaMtrix::IRCMatrix& Kq      // MATRIX FOR Q-TENSOR CALCULATION
                  );

void handleMeshRefinement(std::list<Event*>& refEvents,
                          Geometries& geometries,    // PTRS TO MESHES
                          SolutionVectors& solutionvectors,
                          Simu& simu,
                          Alignment& alignment,
                          Electrodes& electrodes,
                          LC& lc,
                          SpaMtrix::IRCMatrix &Kpot,
                          SpaMtrix::IRCMatrix &Kq
                          );

void handlePreRefinement(std::list<Event*>& refEvents,
                          Geometries& geometries,    // PTRS TO MESHES
                          SolutionVectors& solutionvectors,
                          Simu& simu,
                          Alignment& alignment,
                          Electrodes& electrodes,
                          LC& lc,
                          SpaMtrix::IRCMatrix &Kpot,
                          SpaMtrix::IRCMatrix &Kq
                          );
#endif // EVENTHANDLER_H
