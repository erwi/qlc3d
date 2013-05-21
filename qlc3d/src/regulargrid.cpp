
#include <regulargrid.h>

const idx RegularGrid::NOT_AN_INDEX = std::numeric_limits<idx>::max();
const idx RegularGrid::MAX_SIZE_T = std::numeric_limits<idx>::max();

RegularGrid::RegularGrid():
    nx_(0), ny_(0), nz_(0),
    npr_(0),
    dx_(0), dy_(0), dz_(0)

{
    xLimits_[0] = 0; xLimits_[1] = 0;
    yLimits_[0] = 0; yLimits_[1] = 0;
    zLimits_[0] = 0; zLimits_[1] = 0;

}

RegularGrid::RegularGrid(const RegularGrid &rg):
    nx_(rg.nx_), ny_(rg.ny_), nz_(rg.nz_),
    npr_(rg.npr_),
    dx_(rg.dx_), dy_(rg.dy_), dz_(rg.dz_)
{
    xLimits_[0] = rg.xLimits_[0];
    xLimits_[1] = rg.xLimits_[1];
    yLimits_[0] = rg.yLimits_[0];
    yLimits_[1] = rg.yLimits_[1];
    zLimits_[0] = rg.zLimits_[0];
    zLimits_[1] = rg.zLimits_[1];

    lookupList.insert(lookupList.begin(), rg.lookupList.begin() , rg.lookupList.end() );
}



double RegularGrid::getGridX(const unsigned int &xi) const
{
    if (nx_ == 1)
        return 0.5*(xLimits_[0] + xLimits_[1]);
    else
        return xLimits_[0] + xi*dx_;
}

double RegularGrid::getGridY(const unsigned int &yi) const
{
    if (ny_ == 1)
        return 0.5 * (yLimits_[0] + yLimits_[1] );
    else
        return yLimits_[0] + yi*dy_;
}

double RegularGrid::getGridZ(const unsigned int &zi) const
{
    if (nz_ == 1)
        return 0.5 * (zLimits_[0] + zLimits_[1] );
    else
        return zLimits_[0] + zi*dz_;
}

// CONVERSION FROM LINEAR INDEXING WHICH IS USED TO STORE ALL LOOKUP STRUCTS
// IN A VECTOR TO POSITIONAL INDEXING GINVING THE GRID POINT POSITION IN
// X,Y,Z DIMENSIONS
void RegularGrid::linearToGridIndex(const idx li, idx &xi, idx &yi, idx &zi)
{
#ifdef DEBUG
    assert(li < this->lookupList.size() ); // MAKE SURE POINT IS IN LIST
#endif

    // REGULAR GRID POINTS ARE ARRANGED STARTING FORM xmin,ymin,zmin
    // THE INDEX INCREASES FIRST IN X-DIRECTION, THEN Y-DIRECTION AND
    // FINALLY Z-DIRECTION


    idx ppxy = nx_*ny_; // POINTS PER X-Y PLANE

    zi = li / ppxy;     // GET Z-LEVEL

    idx  c = li % ppxy;   // NUMBER OF POINT IN INCOMPLETE LAYER

    yi = c / nx_;
    xi = c % nx_;

}

//void RegularGrid::gridToLinearIndex(const idx xi, const idx yi, const idx zi)
idx RegularGrid::gridToLinearIndex(const idx xi, const idx yi, const idx zi)
{
    // CALCULATE ARRAY POSITION FROM GRID X,Y AND Z INDEXES xi, yi, zi
#ifdef DEBUG
    assert(xi < nx_);
    assert(yi < ny_);
    assert(zi < nz_);
#endif
    const idx nxyp = nx_*ny_; // NUMBER OF NODES IN X-Y PLANE

    return xi + yi*nx_ + zi*nxyp;
}


bool RegularGrid::createFromTetMesh(const unsigned int &nx,
                                    const unsigned int &ny,
                                    const unsigned int &nz,
                                    Geometry *geom)
{
    // CREATES INTEPOLATION TABLE FROM A TETRAHEDRAL MESH DESCRIBED BY geom
    // SO THAT FAST INTEPOLATION CAN BE PERFORMED LATER ON

    xLimits_[0] = geom->getXmin();   xLimits_[1] = geom->getXmax();
    yLimits_[0] = geom->getYmin();   yLimits_[1] = geom->getYmax();
    zLimits_[0] = geom->getZmin();   zLimits_[1] = geom->getZmax();

    // LIMIT MIN NUMBER OF NODES TO 1 PER DIMENSION
    nx_ = nx == 0 ? 1 : nx;
    ny_ = ny == 0 ? 1 : ny;
    nz_ = nz == 0 ? 1 : nz;

    npr_ = nx_*ny_*nz_;


    dx_ = ( xLimits_[1] - xLimits_[0] ) / ( nx_ - 1 );
    dy_ = ( yLimits_[1] - yLimits_[0] ) / ( ny_ - 1 );
    dz_ = ( zLimits_[1] - zLimits_[0] ) / ( nz_ - 1 );

    // SPECIAL CASE, WHEN ONLY A SINGLE NODE IN A DIRECTION IS REQUIRED, MAKE dx WHOLE WIDTH OF STRUCTURE
    if ( nx_ == 1 ) dx_ = xLimits_[1] - xLimits_[0];
    if ( ny_ == 1 ) dy_ = yLimits_[1] - yLimits_[0];
    if ( nz_ == 1 ) dz_ = zLimits_[1] - zLimits_[0];

    generateLookupList(geom);
    //validateLookupVector();
    return true;
}


bool RegularGrid::generateLookupList(Geometry *geom)
{
    // FILLS IN VALUES FOR THE INTEPOLATION LOOKUP TABLES

    double* coords = new double[npr_*3]; // allocate temporary memory for regular grid coordinates

    size_t cc = 0;    // coordinate counter
    for (unsigned int k = 0 ; k < nz_ ; k++ )// loop over z
    {
        double z = getGridZ(k);

        for (unsigned int j = 0 ; j < ny_ ; j++ ) // loop over y
        {
            double y = getGridY( j );

            for ( unsigned int i = 0 ; i < nx_ ; i++, cc++) // loop over x
            {
                double x = getGridX( i );

                size_t ind = 3*cc;
                coords[ ind + 0 ] = x;
                coords[ ind + 1 ] = y;
                coords[ ind + 2 ] = z;
            }
        }// end loop over y
    }// end loop over z

    // GENERATE INDEX TO TETS THAT CONTAIN EARCH REGULAR GRID POINT
    // INDEX WILL HAVE SPECIAL VALUE Geom::NOT_AN_INDEX, IF POITN WAS NOT FOUND
    // THIS MAY HAPPEN WHEN THE UNDERLYING TET MESH IS NOT A CUBE
    std::vector< unsigned int > indT; // index to tet containing a regular coordinate
    geom->genIndToTetsByCoords(indT,
                               coords,
                               npr_,
                               false, // do NOT terminate app if a coord is not found
                               false );//do NOT require LC element (although it should be preferred, add this option later)

    // NOW CALCULATE WEIGHTS AND NODE INDEXES FOR EACH REGULAR GRID POINT
    lookupList.clear();
    lookupList.reserve( npr_ );
    Mesh* t = geom->t;   // TETRAHEDRAL MESH. OBVIOUSLY THIS IS EVIL...
    double* p = geom->getPtrTop();
    for (idx i = 0; i < npr_ ; i++)
    {

        lookup lu;                  // NEW LOOKUP TABLE ENTRY
        lu.type = RegularGrid::OK;  // INITIALISE TO GOOD

        if ( indT[i] != Geometry::NOT_AN_INDEX ) // IF CONTAINING TET ELEMENT WAS FOUND
        {
            double* coord = &coords[3*i]; // pointer to coordinates of i'th regular point
            // SET INDEXES TO NEIGHBOURING VERTEXES
            lu.ind[0] = t->getNode( indT[i], 0 );
            lu.ind[1] = t->getNode( indT[i], 1 );
            lu.ind[2] = t->getNode( indT[i], 2 );
            lu.ind[3] = t->getNode( indT[i], 3 );

            // CALCULATE NEIGHBOUR NODE WEIGHTS - THESE ARE THE LOCAL COORDINATES
            // OF THE VERTEXES OF THE TET CONTAINING THE REGULAR POINT
            t->CalcLocCoords( indT[i], p, coord, &lu.weight[0] );

            if (t->getMaterialNumber( indT[i]) > MAT_DOMAIN7 )
                lu.type = RegularGrid::NOT_LC;

        }
        else // CONTAINING TET ELEMENT WAS NOT FOUND
        {
            lu.type = RegularGrid::NOT_FOUND;   // containing element not found

            lu.ind[0] = NOT_AN_INDEX; lu.ind[1] = NOT_AN_INDEX;
            lu.ind[2] = NOT_AN_INDEX; lu.ind[3] = NOT_AN_INDEX;
            lu.weight[0] = 0; lu.weight[1] = 0;
            lu.weight[2] = 0; lu.weight[3] = 0;
        }
        lookupList.push_back( lu );
    }


    delete [] coords;
    return true;


}
double RegularGrid::interpolateNode(const double* valuesIn,
                                    const RegularGrid::lookup& L)const
{
    // USES PRE-CALCULATED WEIGHTS TO INTERPOLATE A SINGLE VALUE
    // WITHIN A SINGLE TET-ELEMENT

    //
    //  SOMETIMES REGULAR GRID NODES ARE NOT FOUND
    //  (NUMERICAL NOISE, HOLE IN MESH, NON CUBOIDAL MESH etc. )
    //  DEAL WITH IT

    if ( L.type == RegularGrid::OK )
    {
        return valuesIn[ L.ind[0] ]*L.weight[0] +
                valuesIn[L.ind[1] ]*L.weight[1] +
                valuesIn[L.ind[2] ]*L.weight[2] +
                valuesIn[L.ind[3] ]*L.weight[3];

    }
    else
    {
        return std::numeric_limits<double>::quiet_NaN(); // OUTPUTS NaN
    }




}

void RegularGrid::interpolateDirNode(const double* vecin,
                                     double* dirout,
                                     const RegularGrid::lookup& L,
                                     const idx npLC)const
{
    // INTERPOLATE DIRECTOR TO A NODE TAKING INTO ACCOUNT HEAD-TAIL SYMMETRY.
    // USE DIRECTOR OF FIRST NODE AS REFERENCE AND MAKE SURE DOT-PRODUCTS WITH
    // OTHER 3 NODES IS POSITIVE

    // LOCAL COPIES OF DIRECTOR AT 4 ELEMENT CORNER NODES
    double n1[3] = { vecin[L.ind[0]] , vecin[ L.ind[0]+npLC ] , vecin[L.ind[0] + 2*npLC]};
    double n2[3] = { vecin[L.ind[1]] , vecin[ L.ind[1]+npLC ] , vecin[L.ind[1] + 2*npLC]};
    double n3[3] = { vecin[L.ind[2]] , vecin[ L.ind[2]+npLC ] , vecin[L.ind[2] + 2*npLC]};
    double n4[3] = { vecin[L.ind[3]] , vecin[ L.ind[3]+npLC ] , vecin[L.ind[3] + 2*npLC]};

    // 3 DOT PRODUCTS WITH REFERENCE DIRECTOR
    double dots[3] = {
        n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2],
        n1[0]*n3[0] + n1[1]*n3[1] + n1[2]*n3[2],
        n1[0]*n4[0] + n1[1]*n4[1] + n1[2]*n4[2] };

    // REDUCE TO SIGN ONLY
    dots[0] = dots[0] >= 0 ? 1.0 : -1.0;
    dots[1] = dots[1] >= 0 ? 1.0 : -1.0;
    dots[2] = dots[2] >= 0 ? 1.0 : -1.0;


    // MULTIPLY EACH DIRECTOR BY SIGN
    n2[0]*=dots[0]; n2[1]*=dots[0]; n2[2]*= dots[0];
    n3[0]*=dots[1]; n3[1]*=dots[1]; n3[2]*= dots[1];
    n4[0]*=dots[2]; n4[1]*=dots[2]; n4[2]*= dots[2];

    // INTERPOLATE
    dots[0] = n1[0]*L.weight[0] + n2[0]*L.weight[1] + n3[0]*L.weight[2] + n4[0]*L.weight[3]; //temp
    dots[1] = n1[1]*L.weight[0] + n2[1]*L.weight[1] + n3[1]*L.weight[2] + n4[1]*L.weight[3];
    dots[2] = n1[2]*L.weight[0] + n2[2]*L.weight[1] + n3[2]*L.weight[2] + n4[2]*L.weight[3];

    // MAINTAIN UNIT LENGTH OF DIRECTOR - NORMALISE IT
    double len = dots[0]*dots[0] + dots[1]*dots[1] + dots[2]*dots[2];
    len = sqrt(len);

    dirout[0] = dots[0] / len;
    dirout[1] = dots[1] / len;
    dirout[2] = dots[2] / len;

}




void RegularGrid::printLookup(const lookup &lu) const
{
    std::cout << "Type: "<< lu.type << std::endl;
    std::cout << "ind: "<<lu.ind[0]<<", "<<lu.ind[1]<<", "<<lu.ind[2]<<", "<<lu.ind[3]<<std::endl;
    std::cout << "weight: "<<lu.weight[0]<<", "<<lu.weight[1]<<", "<<lu.weight[2]<<", "<<lu.weight[3]<<std::endl;
}

void RegularGrid::interpolateToRegular(const double *valIn,
                                       double *&valOut,
                                       const idx np)

{
    // INTERPOLATES A VARIABLE TO REGULAR GRID
    if (!npr_)
    {
        printf("error in %s, Regular grid doesn't seem to be initialised.\n", __func__);
        printf("grid count in x,y,z is %u, %u, %u - bye!\n", this->nx_, this->ny_, this->nz_);
        exit(1);
    }

    for ( idx i = 0 ; i < lookupList.size() ; i++ )
    {
        lookup L = lookupList[i];

        // IF TRYING TO INTEPOLATE LC PARAM TO A NON-LC GRID NODE
        if ( (L.type == RegularGrid::NOT_LC ) &&
             (np < MAX_SIZE_T ) )
        {
            valOut[i] = std::numeric_limits<double>::quiet_NaN(); // OUTPUTS NaN
        }
        else
        {
            valOut[i] = interpolateNode( valIn, L);
        }

    }
}






void RegularGrid::interpolateDirToRegular(const double *vecIn,
                                          double *&vecOut,
                                          const idx npLC)
{
    // INTERPOLATES DIRECTOR TO REGULAR GRID.
    // DOES DIRECTOR SWAPPING WITHIN ELEMENT TO MAKE SURE THAT
    // ALL ELEMENT ARE ORIENTED IN SAME(ISH) DIRECTION.
    // THIS IS NECESSARY TO MAINTAIN UNIT LENGTH OF DIRECTOR
    // DIRECTOR COMPONENTS ARE ORDERED AS nx,nx,nx..., ny,ny,ny... nz,nz,nz...
    if (!npr_)
    {
        printf("error in %s, Regular grid doesn't seem to be initialised.\n", __func__);
        printf("grid count in x,y,z is %u, %u, %u - bye!\n", this->nx_, this->ny_, this->nz_);
        exit(1);
    }
    for ( idx i = 0 ; i < npr_ ; i++ )
    {

        lookup L = lookupList[i];

        // IF TRYING TO INTEPOLATE LC PARAM TO A NON-LC GRID NODE
        if (L.type == RegularGrid::OK) // LC NODE
        {
            double dir[3];
            interpolateDirNode( vecIn,dir,L, npLC);
            vecOut[i + 0*npr_] = dir[0];
            vecOut[i + 1*npr_] = dir[1];
            vecOut[i + 2*npr_] = dir[2];
        }
        else
            //if ( (L.type == RegularGrid::NOT_LC ) )
        {
            vecOut[i +0*npr_] = std::numeric_limits<double>::quiet_NaN();
            vecOut[i +1*npr_] = std::numeric_limits<double>::quiet_NaN();
            vecOut[i +2*npr_] = std::numeric_limits<double>::quiet_NaN();
        }

        //else
        //{
        //    double dir[3];
        //    interpolateDirNode( vecIn,dir,L, npLC);
        //    vecOut[i + 0*npr_] = dir[0];
        //    vecOut[i + 1*npr_] = dir[1];
        //    vecOut[i + 2*npr_] = dir[2];
        //}

    }
}

void RegularGrid::validateLookupVector()
{
    // VALIDATES LOOKUP VECTOR TO MAKE IT IS OK AND WON'T CAUSE
    // TROUBLE LATER ON.
    // TERMINATES PROGRAM AND TRIES TO PRINT INFO ON WHAT WENT WRONG
    // IF LOOKUP VECTOR CONTAINS INVALID ENTRIES
    using std::cout;
    using std::endl;
    for (idx i = 0 ; i < this->lookupList.size() ; ++i)
    {
        const lookup &L = lookupList[i];

        if (L.type != OK)
        {
            cout << "\nerror in regular grid lookup table. \nEntry " << i << " not found\n"<<endl;

            cout << "total number of regular grid points : " << npr_ << endl;
            cout << "number of points in each direction : " << nx_ << ", " <<
                    ny_ << ", " << nz_ << endl;

            cout << "Grid bounds (min,max) :" << endl;
            cout << "\t x: "<< xLimits_[0] << ", " << xLimits_[1] << endl;
            cout << "\t y: "<< yLimits_[0] << ", " << yLimits_[1] << endl;
            cout << "\t z: "<< zLimits_[0] << ", " << zLimits_[1] << endl;

            idx ix, iy, iz;
            this->linearToGridIndex(i, ix, iy, iz);
            cout << "Point grid indexes (nx,ny,nz) :"<< ix <<", "
                 << iy <<", " << iz << endl;
            cout << "Point coordinates (x,y,z) : "<< getGridX(ix) << ", "<<
                    getGridY(iy) << ", " << getGridZ(iz) << endl;

            this->printLookup( L );


            cout << "If values above look OK, try reducing the coordinate comparison accuracy"<< endl;
            cout << "Current coordinate comparison accuracy is " << qlc3d_GLOBALS::GLOBAL_COORD_EPS << endl;


            cout << "Otherwise, this may be due to a problem in the mesh - bye!" << endl;
            exit(1);
        }
    }


}








bool RegularGrid::writeVTKGrid(const char *filename,
                               const double *pot,
                               const double *n,
                               const idx npLC)
{

    // WRITES POTENTIAL, ORDER PARAMETER AND DIRECTOR ONTO VTK REGULAR GRID FILE

    if (npr_ == 0 ){
        printf("error in %s, Regular grid doesn't seem to be initialised.\n", __func__);
        printf("grid count in x,y,z is %u, %u, %u - bye!\n", this->nx_, this->ny_, this->nz_);
        exit(1);
    }

    std::fstream fid;
    fid.open( filename, std::fstream::out );

    if (!fid.is_open() ){    // EXIT IF COULDN'T OPEN FILE
        return false;
    }


    double* regU = new double[ npr_ ];  // TEMPORARY STORAGE FOR INTERPOLATED VALUES
    double* regN = new double[ 3*npr_];
    double* regS = new double[ npr_];
    const double* S = n+3*npLC;   // POINTS TO START OF IRREGULAR S

    interpolateToRegular( pot , regU );     // IRREGULAR TO REGULAR CONVERSION
    interpolateDirToRegular( n , regN, npLC );
    interpolateToRegular( S , regS, npLC );


    idx num_points[3] = {nx_, ny_, nz_};
    double grid_spacing[3] = {dx_, dy_, dz_};
    double origin[3] = { getGridX(0), getGridY(0), getGridZ(0)};
    vtkIOFun::writeID( fid );

    vtkIOFun::writeHeader( fid,
                           "header string",
                           vtkIOFun::ASCII,
                           num_points,
                           grid_spacing,
                           origin);

    vtkIOFun::writeScalarData( fid , npr_, "potential", regU);

    vtkIOFun::writeScalarData( fid , npr_, "S", regS);

    vtkIOFun::writeVectorData( fid, npr_, "director", regN, regN+npr_, regN+2*npr_);
    printf("%s VTK\n", __func__); fflush(stdout);
    delete [] regU;
    delete [] regN;
    delete [] regS;

    return true;

}



bool RegularGrid::writeVecMat(const char *filename,
                              const double *pot,
                              const double *n,
                              const idx npLC,
                              const double time)
{
    // WRITES OUTPUT IN A MATLAB FILE.
    // VALUES ARE WRITTEN IN 2D MATRIXES, WHERE EACH ROW CORRESPONDS TO A
    // COLUMN OF VALUES Zmin->Zmax  IN THE MODELLED STRUCTURE

    std::ofstream fid(filename);
    if ( !fid.good() )
        return false;

    fid << "grid_size = ["<<nx_<<","<<ny_<<","<<nz_<<"];"<<std::endl;
    fid << "current_time = " << time << ";" << std::endl;

    double* regU = new double[ npr_ ];  // TEMPORARY STORAGE FOR INTERPOLATED VALUES
    double* regN = new double[ 3*npr_]; // INTERPOLATED DIRECTOR
    double* regS = new double[ npr_];   // ORDER PARAMETER

    const double* S = n+3*npLC;   // START OF IRREGULAR S

    interpolateToRegular( pot , regU );     // IRREGULAR TO REGULAR CONVERSION
    interpolateToRegular( S , regS, npLC );
    interpolateDirToRegular( n , regN, npLC );

    //*
    for (idx i = 0 ; i < npr_ ; i++)
    {
        if (regN[i]==regN[i])
        {
            double nx = regN[i];
            double ny = regN[i+npr_];
            double nz = regN[i+2*npr_];
            double len = nx*nx + ny*ny + nz*nz;

            if ( fabs( len-1.0 ) > 0.05)
            {
                printf("len[%i]  = %e\n", i, len);
                exit(1);
            }
        }
    }
    //*/

    MatlabIOFun::writeNumberColumns(fid,
                                    "V",
                                    regU,
                                    nx_, ny_, nz_);
    MatlabIOFun::writeNumberColumns(fid,
                                    "nx",
                                    regN,
                                    nx_, ny_, nz_ );
    MatlabIOFun::writeNumberColumns(fid,
                                    "ny",
                                    regN+npr_,
                                    nx_, ny_, nz_ );
    MatlabIOFun::writeNumberColumns(fid,
                                    "nz",
                                    regN+2*npr_,
                                    nx_, ny_, nz_ );
    MatlabIOFun::writeNumberColumns(fid,
                                    "S",
                                    regS,
                                    nx_, ny_, nz_ );


    delete [] regU;
    delete [] regN;
    delete [] regS;
    fid.close();

    return true;

}

bool RegularGrid::writeDirStackZ(const char *filename,
                                 const double *n,
                                 const idx npLC,
                                 const double time)
{
    /*!
  Writes output in a comma separated values text file. Only the Director
  component values will be written.
  Each row in file corresponds to a column along the z-axis, where the
  director components are interleaved in order nx,ny,nz, nx,ny,nz ...
  The positions of the columns within the structure increase in rows along
  the x-axis, then incrementing the y-position at the end of each row.
  Additionally, on the first row, the number of points in x,y,and z
  directions and current simulation time
  are printed, with director data starting on second row.
  */

    std::ofstream fid(filename);
    if( !fid.good() )
        return false;

    // FIRST WRITE GRID SIZE AND CURRENT TIME
    fid << nx_ <<','<< ny_<<',' << nz_ <<','<< time << std::endl;

    // INTEPOLATE DIRECTOR TO REGULAR GRID
    double *regN = new double[3*npr_];
    interpolateDirToRegular(n, regN, npLC);

    for (idx y = 0 ; y < ny_ ; y++)
    {
        for (idx x = 0 ; x < nx_ ; x++)
        {
            //double os = x*nz_ + n*nx_*nz_; // OFFSET TO START OF THIS COLUMN
            for (idx z = 0 ; z < nz_ ; z++)
            {
                idx i = gridToLinearIndex(x, y, z);
                double nx = regN[i];
                double ny = regN[i+npr_];
                double nz = regN[i+2*npr_];

                if (z > 0) fid <<",";

                // DIRECTO IN DIELECTRICS IS NaN
                if (nx==nx)
                   fid << nx<<","<<ny<<","<<nz;
                else
                    fid <<"NaN,NaN,NaN";
            }// end for zcolumn
            fid <<std::endl;
        }// end for x
    }// end for y

    fid.close();
    delete [] regN;
    return true;
}
