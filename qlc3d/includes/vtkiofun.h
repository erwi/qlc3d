<<<<<<< HEAD
#ifndef VTKIOFUN_H
#define VTKIOFUN_H
//
//  FUNCTIONS FOR INPUT/OUPUT OF VTK FILES
//

#include <fstream>
#include <iostream>
#include <string>
namespace vtkIOFun
{

    extern const char* const ID_STRNIG;
    enum FileFormat {ASCII , BINARY};




    bool writeID(std::fstream& fid); // writes file identifier

    bool writeHeader( std::fstream& fid,
                      const char* headerString);

    bool writeFileFormat(  std::fstream& fid,
                           const FileFormat& format);

    bool writeDatasetFormat(std::fstream& fid,
                      const int& nx, const int& ny, const int& nz,
                      const double& ox, const double& oy, const double& oz,     // origin coords
                      const double& sx, const double& sy, const double& sz );   // grid spacing

    bool writeScalarData(std::fstream& fid,
                         const unsigned int& np,
                         const char* data_name,
                         const double* data);
}//end namespace

#endif // VTKIOFUN_H
=======
#ifndef VTKIOFUN_H
#define VTKIOFUN_H
//
//  FUNCTIONS FOR INPUT/OUPUT OF VTK FILES
//

#include <fstream>
#include <iostream>
#include <string>
namespace vtkIOFun
{

    extern const char* const ID_STRNIG;
    enum FileFormat {ASCII , BINARY};




    bool writeID(std::fstream& fid); // writes file identifier

    bool writeHeader( std::fstream& fid,
                      const char* headerString,
                      const FileFormat& format,
                      const int num_points[3],
                      const double grid_spacing[3]);


    bool writeScalarData(std::fstream& fid,
                         const unsigned int& np,
                         const char* data_name,
                         const double* data);

    bool writeVectorData(std::fstream& fid,
                         const unsigned int& np,
                         const char* data_name,
                         const double* vec_data);

}//end namespace

#endif // VTKIOFUN_H
>>>>>>> 5385d0f5615cb5b4d1e7b2fb7506a9397ef2a0db
