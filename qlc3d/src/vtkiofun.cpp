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
                 const char* headerString )
{
    if ( !fidOK(fid) )
        return false;

    fid << headerString << "\n";

    return true;
}

bool writeFileFormat(std::fstream &fid,
                     const FileFormat &format)
{
    if ( !fidOK(fid) )
        return false;

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

    return true;

}

bool writeDatasetFormat(std::fstream &fid,
                  const int &nx, const int &ny, const int &nz,
                  const double &ox, const double &oy, const double &oz,
                  const double &sx, const double &sy, const double &sz)
{

    if (!fidOK(fid))
        return false;

    fid <<"DATASET STRUCTURED_POINTS"<<std::endl;
    fid <<"DIMENSIONS "<<nx<<" "<<ny<<" "<<nz<<std::endl;
    fid <<"ORIGIN "<<ox<<" "<<oy<<" "<<oz<<std::endl;
    fid <<"SPACING "<<sx<<" "<<sy<<" "<<sz<<std::endl;

    return true;

}

bool writeScalarData(std::fstream &fid,
                     const unsigned int &np,
                     const char *data_name,
                     const double *data)
{
    if (!fidOK(fid))
        return false;

    fid <<"POINT_DATA "<< np << std::endl;
    fid <<"SCALARS "<< data_name <<" double 1" << std::endl;
    fid <<"LOOKUP_TABLE default"<<std::endl;

    for (unsigned int i = 0 ; i < np ; i++ )
        fid << data[i] <<" ";

    fid<<std::endl;


    return true;
}





}// end namespace vtkIOFun


