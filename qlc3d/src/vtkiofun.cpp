#include <vtkiofun.h>

namespace vtkIOFun
{


const char* const ID_STRNIG = "# vtk DataFile Version 3.0";

bool fidOK(std::fstream& fid)
{
    return fid.is_open();
}


bool writeID(std::fstream &fid)
{
    if ( !fidOK(fid) )
        return false;


    fid << ID_STRNIG << std::endl;
    return true;
}

bool writeHeader(std::fstream &fid,
                 const char* headerString,
                 const FileFormat& format,
                 const int num_points[3],
                 const double grid_spacing[3]
                 )
{
    if ( !fidOK(fid) )
        return false;

    fid << headerString << "\n";

// CHOOSE ASCII / BINARY
    switch( format )
    {
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
    fid <<"ORIGIN "<<0<<" "<<0<<" "<<0<<std::endl;
    fid <<"SPACING "<<grid_spacing[0]<<" "<<grid_spacing[1]<<" "<<grid_spacing[2]<<std::endl;
    fid <<"POINT_DATA "<< np << std::endl;
    return true;
}


bool writeScalarData(std::fstream &fid,
                     const unsigned int &np,    // NUMBER OF REGULAR GRID POINTS
                     const char *data_name,     //
                     const double *data)        //
{
// APPENDS SCALAR DATA TO END OF FILE


    if (!fidOK(fid))
        return false;

    fid <<"SCALARS "<< data_name <<" double 1" << std::endl;
    fid <<"LOOKUP_TABLE default"<<std::endl;

    for (unsigned int i = 0 ; i < np ; i++ )
        fid << data[i] <<" ";

    fid<<std::endl;
    return true;
}

bool writeVectorData(std::fstream &fid,
                     const unsigned int &np, //  NUMBER OF REGULAR GRID POINTS
                     const char *data_name,
                     const double *vec_data)
{
// APPENDS VECTOR DATA TO END OF FILE

    if ( !fidOK(fid) )
        return false;

    fid << "VECTORS "<< data_name <<" double"<<std::endl;

    for (size_t i = 0 ; i < np ; i++)
        fid << vec_data[i] <<" "<< vec_data[i+np] <<" "<< vec_data[i+2*np] << std::endl;

    return true;
}






}// end namespace vtkIOFun


