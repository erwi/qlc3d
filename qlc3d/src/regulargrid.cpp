
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
// THIS MAY HAPPEN WHEN THE UNDERLYING TET MESH IS NOT A QUBE
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

    return valuesIn[ L.ind[0] ]*L.weight[0] +
            valuesIn[L.ind[1] ]*L.weight[1] +
            valuesIn[L.ind[2] ]*L.weight[2] +
            valuesIn[L.ind[3] ]*L.weight[3];


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


bool RegularGrid::writeVTKGrid(const char *filename,
                               const double *pot,
                               const double *n,
                               const idx npLC)
{
// WRITES POTENTIAL, ORDER PARAMETER AND DIRECTOR ONTO VTK REGULAR GRID FILE

    if (npr_ == 0 )
    {
        printf("error in %s, Regular grid doesn't seem to be initialised.\n", __func__);
        printf("grid count in x,y,z is %u, %u, %u - bye!\n", this->nx_, this->ny_, this->nz_);
        exit(1);

    }

    std::fstream fid;
    fid.open( filename, std::fstream::out );

    if (!fid.is_open() )    // EXIT IF COULDN'T OPEN FILE
        return false;


    double* regU = new double[ npr_ ];  // TEMPORARY STORAGE FOR INTERPOLATED VALUES
    double* regNx = new double[ npr_];
    double* regNy = new double[ npr_];
    double* regNz = new double[ npr_];
    double* regS = new double[ npr_];

    const double* ny = n+1*npLC;
    const double* nz = n+2*npLC;
    const double* S = n+3*npLC;   // START OF IRREGULAR S


    interpolateToRegular( pot , regU );     // IRREGULAR TO REGULAR CONVERSION
    interpolateToRegular( n , regNx, npLC );
    interpolateToRegular( ny , regNy, npLC );
    interpolateToRegular( nz , regNz, npLC );
    interpolateToRegular( S , regS, npLC );



    int num_points[3] = {nx_, ny_, nz_};
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

    vtkIOFun::writeVectorData( fid, npr_, "director", regNx, regNy, regNz);
    printf("%s VTK\n", __func__); fflush(stdout);
    delete [] regU;
    delete [] regNx;
    delete [] regNy;
    delete [] regNz;
    delete [] regS;

    return true;

}


