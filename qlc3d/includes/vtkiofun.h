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
