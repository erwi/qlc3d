#ifndef RESULTOUTPUT_H
#define RESULTOUTPUT_H


#include <io/vtkiofun.h>
#include <regulargrid.h>
#include <simu.h>
#include <lc.h>
#include <geometry.h>
#include <solutionvector.h>
#include <meshrefinement.h>


namespace ResultIO {

static const char LCVIEW_TEXT_FORMAT_STRING[] = "%i %f %f %f %f %f %f\n";
/*!
 * Writes the mesh to the specified file name
 * @param p point coordinates
 * @param t tetrahedra
 * @param e triangles
 * @param np numberp of points
 * @param fileName output file name
 */
void writeMesh(double *p, Mesh *t, Mesh *e, idx np, const std::string &fileName);

/**
 * Writes the result in ascii text LcView format.
 * @param p mesh coordinate x,y,z, values interleaved like [x1, y1, z1, x2, y2, z2,... , xn, yn, zn]
 * @param t mesh tetrahedral node indices [4 x <num tets>]
 * @param e mesh triangle node indices [3 x <num tris>]
 * @param v potential solution
 * @param q q-tensor solution
 * @param currentIteration the current iteration number
 * @param currentTime current simulation time [seconds]
 * @param meshFileName the mesh file name that describes the curerent mesh.
 */
void writeLCD_T(double *p,
                Mesh *t,
                Mesh *e,
                SolutionVector *v,
                SolutionVector *q,
                int currentIteration,
                double currentTime,
                const std::string &meshFileName);
// WRITES LCVIEW RESULT IN BINARY FORMAT
void writeLCD_B(double *p,
                Mesh *t, Mesh *e,
                SolutionVector *v, SolutionVector *q,
                int currentIteration,
                double currentTime,
                double S0,
                const std::string &meshFileName);

void ReadResult(Simu& simu,         // READS AND LOADS Q-TENSOR VALUES FROM AN EXISTING RESULT FILE
                SolutionVector& q); // TRIES TO FIGURE OUT WHETHER RESULT FILE IS IN TEXT OR BINARY FORMAT

void CreateSaveDir(Simu& simu);

void ReadLCD_B(Simu* simu, SolutionVector* q);

void ReadLCD_T(Simu& simu, SolutionVector& q); // LOADS TEXT FORMAT LCVIEW RESULT FILE

/**
 * Write a result file containing potential as well as LC director and order parameter on the unstructured mesh points.
 * This file format is compatible, for example, with ParaView. See e.g. https://www.paraview.org/Wiki/ParaView/Data_formats.
 * @param p point coordinates
 * @param v potential solution
 * @param q Q-tensor
 * @param fileName output file name
 */
void writeCsvUnstructured(const double *p, // defined in resultoutput.cpp
                          const SolutionVector &v,
                          const SolutionVector &q,
                          const std::string &fileName);
/**
 *
 * @param p coordinate xyz values
 * @param numPoints number of points in p (= 3 x length of p)
 * @param numLcPoints number of LC points in the geometry (numLcPoints <= numPoints)
 * @param tetMesh
 * @param v
 * @param director the liquid crystal director vector
 * @param fileName
 */
void writeVtkUnstructuredAsciiGrid(
        const double *p,
        size_t numPoints,
        size_t numLcPoints,
        const Mesh &tetMesh,
        const SolutionVector &v,
        const double *director,
        const std::string &fileName);
}
#endif // RESULTOUTPUT_H
