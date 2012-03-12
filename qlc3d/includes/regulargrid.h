
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
#include <matlabiofun.h>
#include <globals.h>

class Geometry; // NEED TO DECLARE EXISTENCE OF Geometry CLASS HERE TO AVOID CIRCULAR #includes

class RegularGrid {
public:
    enum LookupType{ OK, NOT_LC, NOT_FOUND}; // TYPE OF REGULAR FRID NODE
private:


    // A UNSTRCTURED NODES -> REGULAR NODE INTERPOLATION LOOKUP TABLE
    struct lookup{
        LookupType  type;    // descriptor of tet element containing sough node
        idx ind[4]; // index to nodes from which to interpolate ( = tet corner nodes)
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
    std::vector <lookup> lookupList; // pre-calculated index-weight values for each node
// ---------------------------------------------------------------------------
// MEMBER FUNCTIONS
    // Return position of i'th regular grid coordinate
    double getGridX(const idx& xi)const;// const {return xLimits_[0] + xi*dx_;}
    double getGridY(const idx& yi)const;// const {return yLimits_[0] + yi*dy_;}
    double getGridZ(const idx& zi)const;// const {return zLimits_[0] + zi*dz_;}

    bool generateLookupList(Geometry* geom);

    double interpolateNode(const double* valueIn,
                            const lookup& L) const;

    // INTEPOLATES DIRECTOR, CHECKING HEAD-TAIL
    void interpolateDirNode(const double* dirIn,
                           double dirOut[3],
                           const lookup& L,
                           const idx npLC) const;


    void printLookup(const lookup& lu) const;

public:

    static const idx NOT_AN_INDEX; // A MAGIC NUMBER INDICATING A NODE NOT IN THE GRID
    static const idx MAX_SIZE_T;
    RegularGrid();
    RegularGrid(const RegularGrid& rg);

    // CREATES A TET MESH -> REGULAR GRID LOOKUP
    bool createFromTetMesh(const unsigned int& nx,   //number of points in each direction
                           const unsigned int& ny,
                           const unsigned int& nz,
                           Geometry* geom); // underlying tet mesh


    // INTERPOLATES A VALUE
    void interpolateToRegular( const double* valIn,    // input irregular
                               double*& valOut,         // output regular
                               const idx maxnp = MAX_SIZE_T); // number of nodes in irregular input (npLC or np)
    void interpolateDirToRegular( const double* vecIn,
                                  double*& vecOut,
                                  const idx npLC );
    // ==============================================
    //
    // FILE OUTPUT FUNCTIONS
    //
    // ==============================================
    bool writeVTKGrid(const char* filename,
                      const double* pot,        // POTENTIAL DATA
                      const double* n,          // DIRECTOR DATA (INCLUDING S)
                      const idx npLC );     // NUMBER OF LC NODES

    bool writeVecMat(const char* filename,
                     const double* pot,
                     const double* n,
                     const idx npLC,
                     const double time = 0 );
};

#endif // REGULARGRID_H

