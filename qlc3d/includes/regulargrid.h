#ifndef REGULARGRID_H
#define REGULARGRID_H
#include <geometry.h>
#include <vector>
#include <limits>
class RegularGrid {

    // A UNSTRCTURED NODES -> REGULAR NODE INTERPOLATION LOOKUP TABLE
    struct lookup{
        unsigned int ind[4]; // index to nodes from which to interpolate ( = tet corner nodes)
        double weight[4];       // weights applied to each value ( = tet local coords)
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

public:

    static const unsigned int NOT_AN_INDEX; // A MAGIC NUMBER INDICATING A NODE NOT IN THE GRID

    RegularGrid();

    bool createFromTetMesh(const int& nx, const int& ny, const int& nz, //number of points in each direction
                           Geometry& geom); // underlying tet mesh

};

#endif // REGULARGRID_H
