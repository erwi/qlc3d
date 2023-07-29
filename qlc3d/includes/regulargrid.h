
#ifndef REGULARGRID_H
#define REGULARGRID_H
#include <geometry.h>
#include <vector>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <io/vtkiofun.h>
#include <io/matlabiofun.h>
#include <globals.h>
#include <assert.h>

class Geometry; // FORWARD DECLARATION
class SolutionVector;
namespace qlc3d {
  class Director;
};

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
    unsigned int numRegularPoints_;   // total number of regular points

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

    // CONVERSION FROM LINEAR INDEXING WHICH IS USED TO STORE ALL LOOKUP STRUCTS
    // IN A VECTOR TO POSITIONAL INDEXING GINVING THE GRID POINT POSITION IN
    // X,Y,Z DIMENSIONS
    void linearToGridIndex(const idx li, idx &xi, idx &yi, idx &zi);
    // CALCULATE ARRAY POSITION FROM GRID X,Y AND Z INDEXES
    idx gridToLinearIndex(const idx xi, const idx yi, const idx zi);

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

    [[nodiscard]] std::vector<double> interpolateToRegular(const SolutionVector &pot) const;
    [[nodiscard]] std::vector<double> interpolateToRegularS(const std::vector<qlc3d::Director> &dir) const;
    [[nodiscard]] std::vector<qlc3d::Director> interpolateToRegularDirector(const std::vector<qlc3d::Director> &dir) const;
    // ==============================================
    //
    // FILE OUTPUT FUNCTIONS
    //
    // ==============================================
    bool writeVTKGrid(const std::filesystem::path &fileName,
                      const SolutionVector &pot,
                      const std::vector<qlc3d::Director> &dir);

    bool writeVecMat(const std::filesystem::path &filename,
                     const SolutionVector &pot,
                     const std::vector<qlc3d::Director> &dir,
                     double time = 0 );

  /**
  * Writes output in a comma separated values text file. Only the Director component values will be written.
  * Each row in file corresponds to a column along the z-axis, where the director components are interleaved in order nx,ny,nz, nx,ny,nz ...
  * The positions of the columns within the structure increase in rows along the x-axis, then incrementing the y-position at the end of each row.
  * Additionally, on the first row, the number of points in x,y,and z directions and current simulation time
  * are printed, with director data starting on second row.
  */
  bool writeDirStackZ(const std::filesystem::path &filename,
                      const std::vector<qlc3d::Director> &dir,
                      double time = 0);
};

#endif // REGULARGRID_H

