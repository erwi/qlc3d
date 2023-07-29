#ifndef RESULTOUTPUT_H
#define RESULTOUTPUT_H

#include <io/vtkiofun.h>
#include <regulargrid.h>
#include <simu.h>
#include <lc.h>
#include <geometry.h>
#include <solutionvector.h>
#include <meshrefinement.h>

#include <filesystem>

class Coordinates;

namespace ResultIO {

static const char LCVIEW_TEXT_FORMAT_STRING[] = "%i %f %f %f %f %f %f\n";


void ReadResult(Simu& simu,         // READS AND LOADS Q-TENSOR VALUES FROM AN EXISTING RESULT FILE
                SolutionVector& q); // TRIES TO FIGURE OUT WHETHER RESULT FILE IS IN TEXT OR BINARY FORMAT

void ReadLCD_B(Simu* simu, SolutionVector* q);

void ReadLCD_T(Simu& simu, SolutionVector& q); // LOADS TEXT FORMAT LCVIEW RESULT FILE

/**
 * Write a result file containing potential as well as LC director and order parameter on the unstructured mesh points.
 * This file format is compatible, for example, with ParaView. See e.g. https://www.paraview.org/Wiki/ParaView/Data_formats.
 * @param coordinates mesh node coordinates
 * @param v potential solution
 * @param q Q-tensor
 * @param fileName output file name
 */
void writeCsvUnstructured(const Coordinates &coordinates,
                          const SolutionVector &v,
                          const SolutionVector &q,
                          const std::string &fileName);
/**
 *
 * @param mesh node coordinate xyz values
 * @param numPoints number of points in p (= 3 x length of p)
 * @param numLcPoints number of LC points in the geometry (numLcPoints <= numPoints)
 * @param tetMesh
 * @param v potential solution
 * @param q liquid crystal solution
 * @param fileName
 */
void writeVtkUnstructuredAsciiGrid(
        const Coordinates &coordinates,
        size_t numLcPoints,
        const Mesh &tetMesh,
        const SolutionVector &v,
        const SolutionVector &q,
        const std::string &fileName);
}
#endif // RESULTOUTPUT_H
