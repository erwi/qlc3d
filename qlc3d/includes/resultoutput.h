#ifndef RESULTOUTPUT_H
#define RESULTOUTPUT_H


#include <vtkiofun.h>
#include <regulargrid.h>
#include <simu.h>
#include <lc.h>
#include <geometry.h>
#include <solutionvector.h>
#include <meshrefinement.h>


namespace LCviewIO
{

static const char LCVIEW_TEXT_FORMAT_STRING[] = "%i %f %f %f %f %f %f\n";
// FOLLOWING FUNCTIONS ARE DEFINED IN writeLCD_T.cpp

/*
void WriteLCViewResult(Simu* simu, 		// Simulation settings
                LC* lc,				// LC material paramters
                Geometry* geom,		// mesh geometry data
                SolutionVector* v,  // potential solution
                SolutionVector* q); // Q-tensor solution
                //MeshRefinement* meshref = NULL); // meshrefinement info. including whether a new mesh has been generated
*/
/*!
 * Writes the mesh to the specified file name
 * @param p point coordinates
 * @param t tetrahedra
 * @param e triangles
 * @param np numberp of points
 * @param fileName output file name
 */
void writeMesh(double *p, Mesh *t, Mesh *e, idx np, const std::string &fileName);

// WRITES LCVIEW RESULT IN TEXT FORMAT
void writeLCD_T(double *p, Mesh *t, Mesh *e, SolutionVector *v, SolutionVector *q, int currentIteration, double currentTime, const std::string &meshFileName);
// WRITES LCVIEW RESULT IN BINARY FORMAT
void writeLCD_B(double *p,
                Mesh *t, Mesh *e,
                SolutionVector *v, SolutionVector *q,
                int currentIteration,
                double currentTime,
                LC* lc,
                const std::string &meshFileName);
void ReadResult(Simu& simu,         // READS AND LOADS Q-TENSOR VALUES FROM AN EXISTING RESULT FILE
                SolutionVector& q); // TRIES TO FIGURE OUT WHETHER RESULT FILE IS IN TEXT OR BINARY FORMAT




void CreateSaveDir(Simu& simu);
void ReadLCD_B(Simu* simu, SolutionVector* q);
void ReadLCD_T(Simu& simu, SolutionVector& q); // LOADS TEXT FORMAT LCVIEW RESULT FILE
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
