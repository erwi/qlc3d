#ifndef VTKIOFUN_H
#define VTKIOFUN_H
//
//  FUNCTIONS FOR INPUT/OUPUT OF VTK FILES
//

#include <fstream>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

class Mesh;
class SolutionVector;
class Coordinates;
namespace qlc3d {
    class Director;
}

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
                         const char* data_name,
                         const std::vector<double>& data);

    bool writeVectorData(std::fstream& fid,
                         const char* data_name,
                         const std::vector<qlc3d::Director>& data);

    /**
     * Writes VTK unstructured grid compatible file. See https://kitware.github.io/vtk-examples/site/VTKFileFormats/
     * for more info about the format.
     *
     * @note For quadratic tetrahedra (TET10), qlc3d stores nodes internally in Gmsh ordering
     * where mid-edge node position [8] = CD midpoint (vertices 2–3) and position [9] = BD midpoint
     * (vertices 1–3). VTK cell type 24 (VTK_QUADRATIC_TETRA) expects the opposite: [8]=BD, [9]=CD.
     * The writer swaps these two positions when emitting quadratic elements so that the output is
     * valid for ParaView / VTK consumers.
     */
    class UnstructuredGridWriter {
        void writePoints(std::ostream &os, const Coordinates &coordinates) const;
        void writeTetrahedra(std::ostream &os, const Mesh &tetrahedra) const;
        void writePotentials(std::ostream &os, size_t numPotentials, const SolutionVector &potentials) const;
        void writeLiquidCrystal(std::ostream &os, size_t numPoints, size_t numLcPoints, const SolutionVector &q) const;

    public:
        void write(const std::filesystem::path &fileName,
                   size_t numLcPoints,
                   const Coordinates &coordinates,
                   const Mesh &tetrahedra,
                   const SolutionVector &potentials,
                   const SolutionVector &q
                   ) const;

        static constexpr unsigned int CELL_TYPE_LINEAR_TETRAHEDRON = 10;
        static constexpr unsigned int CELL_TYPE_QUADRATIC_TETRAHEDRON = 24;
    };
}//end namespace

#endif // VTKIOFUN_H

