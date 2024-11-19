
#include <regulargrid.h>
#include <util/logging.h>
#include <util/exception.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <solutionvector.h>
#include <lc-representation.h>
#include <filesystem>


const idx RegularGrid::NOT_AN_INDEX = std::numeric_limits<idx>::max();
const idx RegularGrid::MAX_SIZE_T = std::numeric_limits<idx>::max();

RegularGrid::RegularGrid():
        nx_(0), ny_(0), nz_(0),
        numRegularPoints_(0),
        dx_(0), dy_(0), dz_(0)

{
    xLimits_[0] = 0; xLimits_[1] = 0;
    yLimits_[0] = 0; yLimits_[1] = 0;
    zLimits_[0] = 0; zLimits_[1] = 0;

}

RegularGrid::RegularGrid(const RegularGrid &rg):
        nx_(rg.nx_), ny_(rg.ny_), nz_(rg.nz_),
        numRegularPoints_(rg.numRegularPoints_),
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
    auto &bounds = geom->getBoundingBox();
    xLimits_[0] = bounds.getXMin(); xLimits_[1] = bounds.getXMax();
    yLimits_[0] = bounds.getYMin(); yLimits_[1] = bounds.getYMax();
    zLimits_[0] = bounds.getZMin(); zLimits_[1] = bounds.getZMax();

    // LIMIT MIN NUMBER OF NODES TO 1 PER DIMENSION
    nx_ = nx == 0 ? 1 : nx;
    ny_ = ny == 0 ? 1 : ny;
    nz_ = nz == 0 ? 1 : nz;

  numRegularPoints_ = nx_ * ny_ * nz_;


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


bool RegularGrid::generateLookupList(Geometry *geom) {
  size_t cc = 0;    // coordinate counter
  std::vector<Vec3> gridPoints;
  for (unsigned int k = 0; k < nz_; k++) {// loop over z
    double z = getGridZ(k);

    for (unsigned int j = 0; j < ny_; j++) { // loop over y
      double y = getGridY(j);
      for (unsigned int i = 0; i < nx_; i++, cc++) { // loop over x
        double x = getGridX(i);
        gridPoints.emplace_back(x, y, z);
      }
    }// end loop over y
  }// end loop over z
  Coordinates gridCoordinates(std::move(gridPoints));

  // GENERATE INDEX TO TETS THAT CONTAIN EARCH REGULAR GRID POINT
  // INDEX WILL HAVE SPECIAL VALUE Geom::NOT_AN_INDEX, IF POITN WAS NOT FOUND
  // THIS MAY HAPPEN WHEN THE UNDERLYING TET MESH IS NOT A CUBE
  std::vector<unsigned int> indT; // index to tet containing a regular coordinate
  geom->genIndToTetsByCoords(indT,
                             gridCoordinates,
                             false, // do NOT terminate app if a coord is not found
                             false);//do NOT require LC element (although it should be preferred, add this option later)

  // NOW CALCULATE WEIGHTS AND NODE INDEXES FOR EACH REGULAR GRID POINT
  lookupList.clear();
  lookupList.reserve(numRegularPoints_);
  const Mesh &t = geom->getTetrahedra();
  //double* p = geom->getPtrTop();
  for (idx i = 0; i < numRegularPoints_; i++) {

    lookup lu;                  // NEW LOOKUP TABLE ENTRY
    lu.type = RegularGrid::OK;  // INITIALISE TO GOOD

    if (indT[i] != Geometry::NOT_AN_INDEX) { // IF CONTAINING TET ELEMENT WAS FOUND
      Vec3 targetPoint = gridCoordinates.getPoint(i);

      // SET INDEXES TO NEIGHBOURING VERTEXES
      lu.ind[0] = t.getNode(indT[i], 0);
      lu.ind[1] = t.getNode(indT[i], 1);
      lu.ind[2] = t.getNode(indT[i], 2);
      lu.ind[3] = t.getNode(indT[i], 3);

      // CALCULATE NEIGHBOUR NODE WEIGHTS - THESE ARE THE LOCAL COORDINATES
      // OF THE VERTEXES OF THE TET CONTAINING THE REGULAR POINT
      t.calcLocCoords(indT[i], geom->getCoordinates(), targetPoint, &lu.weight[0]);

      if (t.getMaterialNumber(indT[i]) > MAT_DOMAIN7) {
        lu.type = RegularGrid::NOT_LC;
      }

    } else {// CONTAINING TET ELEMENT WAS NOT FOUND
      lu.type = RegularGrid::NOT_FOUND;   // containing element not found

      lu.ind[0] = NOT_AN_INDEX;
      lu.ind[1] = NOT_AN_INDEX;
      lu.ind[2] = NOT_AN_INDEX;
      lu.ind[3] = NOT_AN_INDEX;
      lu.weight[0] = 0;
      lu.weight[1] = 0;
      lu.weight[2] = 0;
      lu.weight[3] = 0;
    }
    lookupList.push_back(lu);
  }
  return true;
}
double RegularGrid::interpolateNode(const double* valuesIn,
                                    const RegularGrid::lookup& L) const
{
    // USES PRE-CALCULATED WEIGHTS TO INTERPOLATE A SINGLE VALUE
    // WITHIN A SINGLE TET-ELEMENT

    //
    //  SOMETIMES REGULAR GRID NODES ARE NOT FOUND
    //  (NUMERICAL NOISE, HOLE IN MESH, NON CUBOIDAL MESH etc. )
    //  DEAL WITH IT

    if ( L.type == RegularGrid::OK ) {
        return valuesIn[ L.ind[0] ]*L.weight[0] +
                valuesIn[L.ind[1] ]*L.weight[1] +
                valuesIn[L.ind[2] ]*L.weight[2] +
                valuesIn[L.ind[3] ]*L.weight[3];
    }
    else {
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

void RegularGrid::interpolateToRegular(const double *valIn,
                                       double *&valOut,
                                       const idx np) {
    // INTERPOLATES A VARIABLE TO REGULAR GRID
    if (!numRegularPoints_) {
        RUNTIME_ERROR("Regular grid is not initialised.");
    }

    for (idx i = 0; i < lookupList.size(); i++ ) {
        lookup L = lookupList[i];

        // If trying to interpolate to a non-LC region
        if ((L.type == RegularGrid::NOT_LC) && (np < MAX_SIZE_T)) {
            valOut[i] = std::numeric_limits<double>::quiet_NaN(); // OUTPUTS NaN
        }
        else {
            valOut[i] = interpolateNode( valIn, L);
        }
    }
}

std::vector<double> RegularGrid::interpolateToRegular(const SolutionVector &pot) const {
  std::vector<double> regPot;
  for (auto &l : lookupList) {
    if (l.type == RegularGrid::OK) {
      regPot.push_back(
              pot.getValue(l.ind[0]) * l.weight[0] +
              pot.getValue(l.ind[1]) * l.weight[1] +
              pot.getValue(l.ind[2]) * l.weight[2] +
              pot.getValue(l.ind[3]) * l.weight[3]);
    }
    else {
      regPot.push_back(std::numeric_limits<double>::quiet_NaN());
    }
  }
  return regPot;
}

std::vector<double> RegularGrid::interpolateToRegularS(const std::vector<qlc3d::Director> &dir) const {
  std::vector<double> regS;
  for (auto &l : lookupList) {
    if (l.type == RegularGrid::NOT_LC || l.type == RegularGrid::NOT_FOUND) {
      regS.push_back(std::numeric_limits<double>::quiet_NaN());
    } else {
      regS.push_back(
              dir[l.ind[0]].S() * l.weight[0] +
              dir[l.ind[1]].S() * l.weight[1] +
              dir[l.ind[2]].S() * l.weight[2] +
              dir[l.ind[3]].S() * l.weight[3]);
    }
  }
  return regS;
}

std::vector<qlc3d::Director> RegularGrid::interpolateToRegularDirector(const std::vector<qlc3d::Director> &dir) const {
  std::vector<qlc3d::Director> regDir;
  for (auto &l : lookupList) {
    if (l.type == RegularGrid::NOT_LC || l.type == RegularGrid::NOT_FOUND) {
      regDir.push_back(qlc3d::Director({1, 0, 0}, std::numeric_limits<double>::quiet_NaN())); // S = NaN in non-LC regions
    } else {

      double nx = dir[l.ind[0]].nx() * l.weight[0] +
                  dir[l.ind[1]].nx() * l.weight[1] +
                  dir[l.ind[2]].nx() * l.weight[2] +
                  dir[l.ind[3]].nx() * l.weight[3];

      double ny = dir[l.ind[0]].ny() * l.weight[0] +
                  dir[l.ind[1]].ny() * l.weight[1] +
                  dir[l.ind[2]].ny() * l.weight[2] +
                  dir[l.ind[3]].ny() * l.weight[3];

      double nz = dir[l.ind[0]].nz() * l.weight[0] +
                  dir[l.ind[1]].nz() * l.weight[1] +
                  dir[l.ind[2]].nz() * l.weight[2] +
                  dir[l.ind[3]].nz() * l.weight[3];

      double S = dir[l.ind[0]].S() * l.weight[0] +
                 dir[l.ind[1]].S() * l.weight[1] +
                 dir[l.ind[2]].S() * l.weight[2] +
                 dir[l.ind[3]].S() * l.weight[3];

      Vec3 n = {nx, ny, nz};
      n.normalize();
      qlc3d::Director d(n, S);
      regDir.push_back(d);
    }
  }
  return regDir;
}


void RegularGrid::interpolateDirToRegular(const double *vecIn,
                                          double *&vecOut,
                                          const idx npLC) {
    // INTERPOLATES DIRECTOR TO REGULAR GRID.
    // DOES DIRECTOR SWAPPING WITHIN ELEMENT TO MAKE SURE THAT
    // ALL ELEMENT ARE ORIENTED IN SAME(ISH) DIRECTION.
    // THIS IS NECESSARY TO MAINTAIN UNIT LENGTH OF DIRECTOR
    // DIRECTOR COMPONENTS ARE ORDERED AS nx,nx,nx..., ny,ny,ny... nz,nz,nz...
    if (!numRegularPoints_) {
        throw std::runtime_error(fmt::format("Regular grid is not intialised in {}, {}", __FILE__, __func__ ));
    }
    for (idx i = 0; i < numRegularPoints_; i++) {
        lookup L = lookupList[i];

        // If LC node
        if (L.type == RegularGrid::OK) {
            double dir[3];
            interpolateDirNode( vecIn,dir,L, npLC);
            vecOut[i + 0 * numRegularPoints_] = dir[0];
            vecOut[i + 1 * numRegularPoints_] = dir[1];
            vecOut[i + 2 * numRegularPoints_] = dir[2];
        }
        else {
            vecOut[i + 0 * numRegularPoints_] = std::numeric_limits<double>::quiet_NaN();
            vecOut[i + 1 * numRegularPoints_] = std::numeric_limits<double>::quiet_NaN();
            vecOut[i + 2 * numRegularPoints_] = std::numeric_limits<double>::quiet_NaN();
        }
    }
}

/**
 * Writes potential, director and order parameted to a VTK grid. TODO: this should be in part of a separate resul IO
 * file/class
*/
bool RegularGrid::writeVTKGrid(const std::filesystem::path &fileName,
                               const SolutionVector &pot,
                               const std::vector<qlc3d::Director> &dir) {

    if (numRegularPoints_ == 0) {
      RUNTIME_ERROR("Regular grid is not initialised.");
    }
    std::fstream fid;
    fid.open( fileName, std::fstream::out );

    if (!fid.is_open() ) {    // EXIT IF COULDN'T OPEN FILE // TODO: throw exception instead?
        return false;
    }
    std::vector<double> regPot = interpolateToRegular(pot);
    std::vector<double> regS = interpolateToRegularS(dir);
    std::vector<qlc3d::Director> regN = interpolateToRegularDirector(dir);

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

    vtkIOFun::writeScalarData(fid, "potential", regPot);
    vtkIOFun::writeScalarData(fid,  "S", regS);
    vtkIOFun::writeVectorData(fid, "director", regN);
    fid.close();
    return true;
}



bool RegularGrid::writeVecMat(const std::filesystem::path &fileName,
                              const SolutionVector &pot,
                              const std::vector<qlc3d::Director> &dir,
                              const double time) {
  // WRITES OUTPUT IN A MATLAB FILE.
  // VALUES ARE WRITTEN IN 2D MATRIXES, WHERE EACH ROW CORRESPONDS TO A
  // COLUMN OF VALUES Zmin->Zmax  IN THE MODELLED STRUCTURE

  std::ofstream fid(fileName);
  if (!fid.good()) {
    return false;
  }

  fid << "grid_size = ["<<nx_<<","<<ny_<<","<<nz_<<"];"<<std::endl;
  fid << "current_time = " << time << ";" << std::endl;

  std::vector<double> regPot = interpolateToRegular(pot);
  std::vector<double> regS = interpolateToRegularS(dir);
  std::vector<qlc3d::Director> regN = interpolateToRegularDirector(dir);

  std::vector<double> nx, ny, nz;
  nx.resize(regS.size(), 0);
  ny.resize(regS.size(), 0);
  nz.resize(regS.size(), 0);

  for (auto &d : regN) {
    nx.push_back(d.nx());
    ny.push_back(d.ny());
    nz.push_back(d.nz());
  }

  MatlabIOFun::writeNumberColumns(fid, "V", regPot, nx_, ny_, nz_);
  MatlabIOFun::writeNumberColumns(fid, "nx", nx,nx_, ny_, nz_ );
  MatlabIOFun::writeNumberColumns(fid,"ny", ny,nx_, ny_, nz_ );
  MatlabIOFun::writeNumberColumns(fid,"nz", nz,nx_, ny_, nz_ );
  MatlabIOFun::writeNumberColumns(fid,"S",regS,nx_, ny_, nz_ );

  fid.close();

  return true;
}

bool RegularGrid::writeDirStackZ(const std::filesystem::path &fileName,
                                 const std::vector<qlc3d::Director> &dir,
                                 double time) {
  std::ofstream fid(fileName);
  if (!fid.good()) {
    return false;
  }

  // First line is grid size
  fid << nx_ <<','<< ny_<<',' << nz_ <<','<< time << std::endl;

  std::vector<qlc3d::Director> regN = interpolateToRegularDirector(dir);

  for (idx y = 0; y < ny_; y++) {
    for (idx x = 0; x < nx_; x++) {
      for (idx z = 0; z < nz_; z++) {
        idx i = gridToLinearIndex(x, y, z);

        if (z > 0) {
          fid <<",";
        }

        if (std::isnan(regN[i].S())) {
          fid << "NaN,NaN,NaN";
        } else {
          fid << regN[i].nx() << "," << regN[i].ny() << "," << regN[i].nz();
        }
      }
      fid << std::endl;
    }
  }
  fid.close();
  return true;
}
