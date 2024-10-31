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

/**
 * Reads in the data from an LCView result file (either binary or text) and stores it in the q-tensor.
 * It is assumed that the q-tensor size is already known from reading a mesh file and that q is correctly sized.
 */
void ReadResult(const std::string &fileName, SolutionVector& q);

/**
 * Reads in the data from a binary LCView result file and stores it in the q-tensor.
 * It is assumed that the q-tensor size is already known from reading a mesh file and that q is correctly sized.
 */
void readBinaryLcViewResultFile(const std::string &fileName, SolutionVector& q);

/**
 * Loads Q-tensor for an ascii text LCView result file. This assumes that the loaded data size is already known and the
 * q-tensor is correctly sized
 */
void readTextLcViewResultFile(const std::string &fileName, SolutionVector& q); // LOADS TEXT FORMAT LCVIEW RESULT FILE

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
