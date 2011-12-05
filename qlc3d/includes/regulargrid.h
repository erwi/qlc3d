
#ifndef REGULARGRID_H
#define REGULARGRID_H
#include <geometry.h>
#include <vector>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vtkiofun.h>
class RegularGrid {
public:
    enum LookupType{ OK, NOT_LC, NOT_FOUND};
private:


    // A UNSTRCTURED NODES -> REGULAR NODE INTERPOLATION LOOKUP TABLE
    struct lookup{
        LookupType  type;    // descriptor of tet element containing sough node
        unsigned int ind[4]; // index to nodes from which to interpolate ( = tet corner nodes)
        double weight[4];    // weights applied to each value ( = tet local coords)

    };


    unsigned int nx_;    // number of points in each direction
    unsigned int ny_;
    unsigned int nz_;
    unsigned int npr_;   // total number of regular points

    double dx_;          // grid spacing in each direction
    double dy_;
    double dz_;

    double xLimits_[2];  // min and max values
    double yLimits_[2];
    double zLimits_[2];

    // Return position of i'th regular grid coordinate
    double getGridX(const unsigned int& xi) const {return xLimits_[0] + xi*dx_;}
    double getGridY(const unsigned int& yi) const {return yLimits_[0] + yi*dy_;}
    double getGridZ(const unsigned int& zi) const {return zLimits_[0] + zi*dz_;}

    std::vector <lookup> lookupList; // pre-calculated index-weight values for each node

    bool generateLookupList(Geometry& geom);

    double interpolateNode(const double* valueIn,
                            const lookup& L) const;

public:

    static const unsigned int NOT_AN_INDEX; // A MAGIC NUMBER INDICATING A NODE NOT IN THE GRID

    RegularGrid();

    // CREATES A TET MESH -> REGULAR GRID LOOKUP
    bool createFromTetMesh(const int& nx,   //number of points in each direction
                           const int& ny,
                           const int& nz,
                           Geometry& geom); // underlying tet mesh

    // INTERPOLATES A SCALAR VALUE
    void interpolateToRegular( const double*& sclrIn,    // input values
                               double*& sclrOut);        // output values

    // INTERPOLATES A VECTOR VALUE
    void interpolateToRegular( const double*& vecIn,    // input irregular
                               double*& vecOut,         // output regular
                               const size_t& np);       // number of nodes in irregular input

    // ==============================================
    //
    // FILE OUTPUT FUNCTIONS
    //
    // ==============================================
    bool writeVTKGrid(const char* filename,
                      const double* pot,        // POTENTIAL DATA
                      const double* n,          // DIRECTOR DATA (INCLUDING S)
                      const size_t& npLC );     // NUMBER OF LC NODES


};

#endif // REGULARGRID_H

