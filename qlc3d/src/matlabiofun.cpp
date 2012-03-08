#include <matlabiofun.h>


bool MatlabIOFun::writeNumberArray(std::ofstream &fid,
                              const char *varName,
                              const double *values,
                              const idx n)
{
    // WRITES A SINGLE ARRAY OF NUMBERS

    if ( !fid.good() )
        return false;

    fid << varName << "=[";

    for (idx i = 0; i < n-1 ; i++ )
    {
        fid <<values[i]<<",";
    }
    fid <<values[n-1] <<"];"<<std::endl;

    return true;
}




bool MatlabIOFun::writeNumberColumns(std::ofstream &fid,
                                     const char *varName,
                                     const double *values,
                                     const idx nx,
                                     const idx ny,
                                     const idx nz)
{

    if ( !fid.good() )
        return false;

    fid << varName << " = [ "<<std::endl;

    for (idx ix = 0 ; ix < nx ; ix++)
    {
        for (idx iy = 0 ; iy < ny ; iy ++)
        {
            fid <<"[";
            idx iz, i;
            for (iz = 0 ; iz < nz-1 ; iz++)
            {
                i = iz*nx*ny + iy*nx +ix;   // VALUES ARE ORDERED DIFFERENTLY
                double val = values[i];
                fid << val << ",";
            }

            i = iz*nx*ny + iy*nx +ix;
            double val = values[i];
            fid << val <<"];" << std::endl;
        }// end for iy
    }// end for ix
    fid << "];" << std::endl;
return true;
}
