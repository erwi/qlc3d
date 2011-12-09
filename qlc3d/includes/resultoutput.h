#ifndef RESULTOUTPUT_H
#define RESULTOUTPUT_H


#include <vtkiofun.h>
#include <regulargrid.h>
#include <simu.h>
#include <lc.h>
#include <geometry.h>
#include <solutionvector.h>
#include <meshrefinement.h>

namespace WriteResults
{

// FOLLOWING FUNCTIONS ARE DEFINED IN WriteLCD.cpp

void WriteResult(Simu* simu, 		// Simulation settings
                LC* lc,				// LC material paramters
                Geometry* geom,		// mesh geometry data
                SolutionVector* v,  // potential solution
                SolutionVector* q,  // Q-tensor solution
                MeshRefinement* meshref = NULL); // meshrefinement info. including whether a new mesh has been generated
void CreateSaveDir(Simu& simu);
void ReadLCD_B(Simu* simu, SolutionVector* q);

// REPLACE WRITESETTNGS WITH DIRECT COPY OF
// SETTINGS FILE
//
//void WriteSettings(Simu* simu,
//                   LC* lc,
//                   Boxes* box,
//                   Alignment* alignment ,
//                   Electrodes* electrodes);
}
#endif // RESULTOUTPUT_H
