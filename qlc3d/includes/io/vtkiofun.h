#ifndef VTKIOFUN_H
#define VTKIOFUN_H
//
//  FUNCTIONS FOR INPUT/OUPUT OF VTK FILES
//

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class Mesh;
class SolutionVector;

namespace vtkIOFun {

    extern const char* const ID_STRING;
    enum FileFormat {ASCII , BINARY};

    bool writeID(std::fstream& fid); // writes file identifier

    bool writeHeader( std::fstream& fid,
                      const char* headerString,
                      const FileFormat& format,
                      const unsigned int num_points[3],
                      const double grid_spacing[3],
                      const double origin[3]);


    bool writeScalarData(std::fstream& fid,
                         const unsigned int& np,
                         const char* data_name,
                         const double* data);

    bool writeVectorData(std::fstream& fid,
                         const unsigned int& np,
                         const char* data_name,
                         const double* vec_data1,
                         const double* vec_data2,
                         const double* vec_data3);

    /**
     * Writes VTK unstructured grid compatible file. See https://kitware.github.io/vtk-examples/site/VTKFileFormats/
     * for more info about the format.
     */
    class UnstructuredGridWriter {
        void writePoints(std::ostream &os, size_t numPoints, const double *points) const;
        void writeTetrahedra(std::ostream &os, const Mesh &tetrahedra, size_t numPoints) const;
        void writePotentials(std::ostream &os, size_t numPotentials, const double *potentials) const;
        void writeLiquidCrystal(std::ostream &os, size_t numPoints, size_t numLcPoints, const SolutionVector &q) const;

    public:
        void write(const std::string &fileName,
                   size_t numPoints,
                   size_t numLcPoints,
                   const double *points,
                   const Mesh &tetrahedra,
                   const double *potentials,
                   const SolutionVector &q
                   ) const;
    };
}//end namespace

#endif // VTKIOFUN_H

