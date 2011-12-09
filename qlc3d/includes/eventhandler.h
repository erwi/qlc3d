#ifndef EVENTHANDLER_H
#define EVENTHANDLER_H
#include <solutionvector.h>
#include <electrodes.h>
#include <geometry.h>
#include <eventlist.h>
#include <simu.h>
#include <meshrefinement.h>
#include <sparsematrix.h>
#include <lc.h>
#include <settings.h>
#include <calcpot3d.h>
#include <resultoutput.h>



void setElectrodePotentials(EventList& evel,
                            Electrodes& electr,
                            Simu& simu);


// TAKSES CARE OF EVENTS OCCRRING BEFORE SIMULATION STARTS
// AND PREPARING EVENTLIST FOR SIMULATION
void handleInitialEvents(EventList& evel,      // EVENT LIST
                        Electrodes& electr,   // ELECTRODES WITH POTENTIALS AND TIMING
                        Simu& simu,
                        SolutionVector& v,    // POTENTIAL SOLUTION
                        Geometry& geom,       // CURRENT MESH
                        MeshRefinement& ref,  // MESH REFINEMENT INFO
                        SparseMatrix* Kpot,   // POTENTIAL CALCULATION MATRIX
                        SolutionVector& q,    // Q-TENSOR
                        LC& lc,               // MATERIAL PARAMS.
                        Settings& settings);  // SPARSE SOLVER SETTINGS


void handleEvents(EventList& evel,      // EVENT LIST
                  Electrodes& electr,   // ELECTRODES WITH POTENTIALS AND TIMING
                  Simu& simu,
                  SolutionVector& v,    // POTENTIAL SOLUTION
                  Geometry& geom,       // CURRENT MESH
                  MeshRefinement& ref,  // MESH REFINEMENT INFO
                  SparseMatrix* Kpot,   // POTENTIAL CALCULATION MATRIX
                  SolutionVector& q,    // Q-TENSOR
                  LC& lc,               // MATERIAL PARAMS.
                  Settings& settings);  // SPARSE SOLVER SETTINGS


#endif // EVENTHANDLER_H
