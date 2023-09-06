
#include <io/vtkiofun.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <mesh.h>
#include <assert.h>
#include <solutionvector.h>
#include <lc-representation.h>

namespace vtkIOFun {
    const char* const ID_STRING = "# vtk DataFile Version 3.0";

    bool fidOK(std::fstream& fid) {
        return fid.is_open();
    }

    bool writeID(std::fstream &fid) {
        if ( !fidOK(fid) )
            return false;

        fid << ID_STRING << std::endl;
        return true;
    }

    bool writeHeader(std::fstream &fid,
                     const char* headerString,
                     const FileFormat& format,
                     const unsigned int num_points[3],
                     const double grid_spacing[3],
                     const double origin[3]) {
        if ( !fidOK(fid) )
            return false;

        fid << headerString << "\n";

// CHOOSE ASCII / BINARY
        switch( format ) {
            case ASCII:
                fid << "ASCII" << std::endl;
                break;
            case BINARY:
                fid << "BINARY" << std::endl;
                break;
            default:
                return false;

        }
        int np = num_points[0] * num_points[1] * num_points[2];
        fid <<"DATASET STRUCTURED_POINTS"<<std::endl;
        fid <<"DIMENSIONS "<<num_points[0]<<" "<<num_points[1]<<" "<<num_points[2]<<std::endl;
        fid <<"ORIGIN "<<origin[0]<<" "<<origin[1]<<" "<<origin[2]<<std::endl;
        fid <<"SPACING "<<grid_spacing[0]<<" "<<grid_spacing[1]<<" "<<grid_spacing[2]<<std::endl;
        fid <<"POINT_DATA "<< np << std::endl;
        return true;
    }

    bool writeScalarData(std::fstream &fid,
                         const char *data_name,
                         const std::vector<double> &data) {
// APPENDS SCALAR DATA TO END OF FILE
// NON-LC VALUES ARE SET TO 0 (VTK DOES NOT LIKE NaN)

        if (!fidOK(fid))
            return false;

        fid <<"SCALARS "<< data_name <<" double 1" << std::endl;
        fid <<"LOOKUP_TABLE default"<<std::endl;

        for (unsigned int i = 0 ; i < data.size(); i++) {
            if (data[i]!= data[i]) {// IF NaN
              fid << 0 << " ";
            } else {
              fid << data[i] << " ";
            }
        }

        fid<<std::endl;
        return true;
    }

    bool writeVectorData(std::fstream &fid, const char *data_name, const std::vector<qlc3d::Director> &data) {
        if (!fidOK(fid)) {
          return false;
        }

        fid << "VECTORS "<< data_name <<" double"<<std::endl;

        for (size_t i = 0; i < data.size(); i++)
          if (std::isnan(data[i].S())) {
            fid << 0 << " " << 0 << " " << 0 << std::endl;
          } else {
            fid << data[i].nx() << " " << data[i].ny() << " " << data[i].nz() << std::endl;
          }

        return true;
    }

    void UnstructuredGridWriter::write(const std::filesystem::path &fileName,
                                       size_t numLcPoints,
                                       const Coordinates &coordinates,
                                       const Mesh &tetrahedra,
                                       const SolutionVector &potentials,
                                       const SolutionVector &q) const {
        assert(numLcPoints <= coordinates.size());

        using namespace std;

        ofstream os(fileName, fstream::out);

        os << "# vtk DataFile Version 2.0\n";
        os << "qlc3d result\n";
        os << "ASCII\n";
        os << "DATASET UNSTRUCTURED_GRID\n";

        writePoints(os, coordinates);
        writeTetrahedra(os, tetrahedra, coordinates.size());

        writePotentials(os, coordinates.size(), potentials);
        writeLiquidCrystal(os, coordinates.size(), numLcPoints, q);

        os.close();
    }

    void UnstructuredGridWriter::writePoints(std::ostream &os, const Coordinates &coordinates) const {
        os << "\n";
        os << "POINTS " << coordinates.size() << " float\n";

        for (unsigned int i = 0; i < coordinates.size(); i++) {
          auto &p = coordinates.getPoint(i);
            os << p.x() << " " << p.y() << " " << p.z() << "\n";
        }
    }

    void UnstructuredGridWriter::writeTetrahedra(std::ostream &os, const Mesh &tetrahedra, size_t numPoints) const {
        size_t numTetrahedra = tetrahedra.getnElements();
        size_t numNodes = tetrahedra.getnNodes();
        size_t arrayLength = numTetrahedra * (numNodes + 1); // length of array required to store cell data

        assert(numNodes == 4);

        os << "\n";
        os << "CELLS " << numTetrahedra << " " << arrayLength << "\n";
        for (unsigned int i = 0; i < numTetrahedra; i++) {
            size_t n1 = tetrahedra.getNode(i, 0);
            size_t n2 = tetrahedra.getNode(i, 1);
            size_t n3 = tetrahedra.getNode(i, 2);
            size_t n4 = tetrahedra.getNode(i, 3);

            assert(n1 < numPoints && n2 < numPoints && n3 < numPoints && n4 < numPoints);

            os << numNodes << " " << n1 << " " << n2 << " " << n3 << " " << n4 << "\n";
        }
        os << "\n";
        os << "CELL_TYPES " << numTetrahedra << "\n";
        const int CELL_TYPE_TETRAHEDRON = 10;
        for (unsigned int i = 0; i < numTetrahedra; i++) {
            os << CELL_TYPE_TETRAHEDRON << "\n";
        }
    }

    void UnstructuredGridWriter::writePotentials(std::ostream &os, size_t numPotentials, const SolutionVector &potentials) const {
        os << "\n";
        os << "POINT_DATA " << numPotentials << "\n";
        os << "SCALARS potential float 1\n";
        os << "LOOKUP_TABLE default\n";

        for (unsigned int i = 0; i < numPotentials; i++) {
            os << potentials.getValue(i) << "\n";
        }
    }

    void UnstructuredGridWriter::writeLiquidCrystal(std::ostream &os, size_t numPoints, size_t numLcPoints, const SolutionVector &q) const {
        os << "\n";
        os << "VECTORS director float" << "\n";
        for (unsigned int i = 0; i < numLcPoints; i++) {
            qlc3d::Director n = q.getDirector(i);
            os << n.nx() << " " << n.ny() << " " << n.nz() << "\n";
        }

        // for dielectric regions, write director with zero length.
        for (unsigned int i = numLcPoints; i < numPoints; i++) {
            os << "0 0 0\n";
        }

        // write order parameter
        os << "\n";
        os << "SCALARS S float 1\n";
        os << "LOOKUP_TABLE default\n";
        for (unsigned int i = 0; i < numLcPoints; i++) {
            qlc3d::Director n = q.getDirector(i);
            os << n.S() << "\n";
        }

        // for dielectric regions, write S = 1. This makes it easier to find low order regions than using 0
        for (unsigned int i = numLcPoints; i < numPoints; i++) {
            os << "1\n";
        }
    }
}// end namespace vtkIOFun


