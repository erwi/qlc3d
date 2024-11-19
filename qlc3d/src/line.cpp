#include "../includes/line.h"
#include <geom/coordinates.h>
#include <geom/vec3.h>

bool Line::operator<(const Line& other) const {
    // L[0] is always smaller than L[1] (see constructor)
    // a line is smaller if its node number ( L[0] or L[1] )
    // is smaller than the other lines node number
    if ( this->L[0] < other.L[0] ) return true;
    if ( ( this->L[0] == other.L[0]) && (this->L[1] < other.L[1] ) ) return true;

    return false; // default
}

bool Line::operator ==( const Line& other) const {
    // same line if same objects or if same node numbersR
    if ( this ==  &other) 
		return true;  // same object
    if ( ( this->L[0] == other.L[0]) && (this->L[1] == other.L[1] ) ) 
		return true; // same node numbers

    return false; // default
}

Line::Line() { }

Line::Line(const int& a,const int& b){
    if (a<b) {
	L[0] = a;
	L[1] = b;
    }
    else {
	L[0] = b;
	L[1] = a;
    }
}

bool Line::isOnFrontSurface(Geometry* geom) {
  auto &bounds = geom->getBoundingBox();
  double y1 = geom->getCoordinates().getPoint(L[0]).y();
  double y2 = geom->getCoordinates().getPoint(L[1]).y();
  return (y1 == bounds.getYMin() && y2 == bounds.getYMin());
}
		
bool Line::isOnBackSurface(Geometry* geom) {
  auto &bounds = geom->getBoundingBox();
  double y1 = geom->getCoordinates().getPoint(L[0]).y();
  double y2 = geom->getCoordinates().getPoint(L[1]).y();
  return (y1 == bounds.getYMax() && y2 == bounds.getYMax());
}

bool Line::isOnRightSurface(Geometry* geom) {
  auto &bounds = geom->getBoundingBox();
  double x1 = geom->getCoordinates().getPoint(L[0]).x();
  double x2 = geom->getCoordinates().getPoint(L[1]).x();
  return (x1 == bounds.getXMax() && x2 == bounds.getXMax());
}

bool Line::isOnLeftSurface(Geometry* geom) {
  auto &bounds = geom->getBoundingBox();
  double x1 = geom->getCoordinates().getPoint(L[0]).x();
  double x2 = geom->getCoordinates().getPoint(L[1]).x();
  return (x1 == bounds.getXMin() && x2 == bounds.getXMin());
}

bool Line::isOnTopSurface(Geometry* geom) {
  auto &bounds = geom->getBoundingBox();
  double z1 = geom->getCoordinates().getPoint(L[0]).z();
  double z2 = geom->getCoordinates().getPoint(L[1]).z();
  return (z1 == bounds.getZMax() && z2 == bounds.getZMax());
}

bool Line::isOnBottomSurface(Geometry* geom) {
  auto &bounds = geom->getBoundingBox();
  double z1 = geom->getCoordinates().getPoint(L[0]).z();
  double z2 = geom->getCoordinates().getPoint(L[1]).z();
  return (z1 == bounds.getZMin() && z2 == bounds.getZMin());
}
	
bool Line::isTopBottomCornerLine(Geometry* geom){
    // is corner line if on two connected vertical surfaces simultaneously
    if(  ( isOnFrontSurface(geom) &&  ( isOnLeftSurface(geom) || isOnRightSurface(geom) ) ) || //front-left or front-right corner ...
	 ( isOnBackSurface(geom)  &&  ( isOnLeftSurface(geom) || isOnRightSurface(geom) ) ) )  // ... or back-left or back-right corner
	return true;
    else
	return false;
}
	
bool Line::isFrontBackCornerLine(Geometry* geom){
	// is corner line if on a side and top/bottom surface simultaneously
	if ( ( isOnLeftSurface(geom) && ( isOnTopSurface(geom) || isOnBottomSurface(geom) ) ) || // left-top or left-bottom edge
	     ( isOnRightSurface(geom) &&( isOnTopSurface(geom) || isOnBottomSurface(geom) ) ) )  // ... or on right-top or right-bottom
	    return true;
	else
	return false;
}
	
bool Line::isLeftRightCornerLine(Geometry* geom){
	// if is on front/back and top/bottom simultaneously
	if ( ( isOnFrontSurface(geom) && ( isOnTopSurface(geom) || isOnBottomSurface(geom) ) ) || // front-top or front-bottom ...
		 ( isOnBackSurface(geom)  && ( isOnTopSurface(geom) || isOnBottomSurface(geom) ) ) )  // ... or on back-top or back-bottom
	    return true;
	else
	    return false;
}

bool Line::isTranslationOf( Line& L2, Geometry* geom, double* dir)
{
	/*! Checks whether this line is a translation of line L2. 
	 * dir determines which dimensions are ignored from comparison. e.g.
	 * dir = [0,0,1] only compares for z-components of the lines
	 * dir = [1,0,1] makes sure x and z components match (i.e. shift along y-axis)
	 */
	
	
	if (L2 == *this) return false; // avoid selfs
	
	// The order of nodes is not guaranteed. Check for parallel lines by 
	// taking the dot product of the two vectors
  auto& coords = geom->getCoordinates();
	double eps = 1e-5;

  Vec3 p11 = coords.getPoint(L[0]);
  Vec3 p12 = coords.getPoint(L[1]);

  Vec3 p21 = coords.getPoint(L2.L[0]);
  Vec3 p22 = coords.getPoint(L2.L[1]);

  // calculate vectors along the lines
  Vec3 v1 = p12 - p11;
  Vec3 v2 = p22 - p21;

  double l1 = v1.norm();
  double l2 = v2.norm();
  if (std::abs(l1- l2) > eps) {
    return false; // not same length
  }

  // check if parallel
  double dotNormalised = std::abs(v1.dot(v2) / (l1*l2)); // should be 1 for parallel or opposite vectors
  if (std::abs(1 - dotNormalised) > eps) {
    return false; // not parallel
  }

  // lines are parallel and of same length. Check if they are translations of each other (I don't understand this part anymore)
  // calculate differences between line end points
  Vec3 diff1 = p11;
  Vec3 diff2 = p12;

  // calculate shifts. differences depend on line orientations
  if (v1.dot(v2) > 0) {// same direction
    diff1-= p21; // difference between nodes 1
    diff2-= p22; // difference between nodes 2
  }
  else { // opposite directions
    diff1-= p22; // difference between nodes 1
    diff2-= p21; // difference between nodes 2
  }

  // set x,y,z components of shift to zero, as defined in dir[0]->dir[2]
  double x = (std::abs(diff1.x()) + std::abs(diff2.x())) * dir[0];
  double y = (std::abs(diff1.y()) + std::abs(diff2.y())) * dir[1];
  double z = (std::abs(diff1.z()) + std::abs(diff2.z())) * dir[2];

  //printf("%f,%f,%f\n", diff1[0], diff1[1], diff1[2]);
  return (x + y + z) <= eps;
}
	
// corners along z
bool Line::isCorn0(Geometry* geom) // xmin, ymin
{
	return (isOnLeftSurface(geom) && isOnFrontSurface(geom) );
}

bool Line::isCorn1(Geometry* geom) // xmax, ymin
{
	return (isOnFrontSurface(geom) && ( isOnRightSurface(geom) ) );
}
bool Line::isCorn2(Geometry* geom) // xmax, ymax
{
	return ( isOnRightSurface(geom) && isOnBackSurface(geom) );
}
bool Line::isCorn3(Geometry* geom) // xmin, ymax
{
	return ( isOnLeftSurface(geom) && isOnBackSurface(geom) );
}

	// corners along x
bool Line::isCorna(Geometry* geom) // ymin, zmin
{
	return ( isOnFrontSurface(geom) && isOnBottomSurface(geom) );
}
bool Line::isCornb(Geometry* geom) // ymax, zmin
{
	return ( isOnBottomSurface(geom) && isOnBackSurface(geom) );
}
bool Line::isCornc(Geometry* geom) // ymax, zmax
{
	return (isOnBackSurface(geom) && isOnTopSurface(geom) );
}
bool Line::isCornd(Geometry* geom) // ymin, zmax
{
	return (isOnTopSurface(geom) && isOnFrontSurface(geom) );
}

// corners along y
bool Line::isCornA(Geometry* geom) // xmin, zmin
{
	return (isOnLeftSurface(geom) && isOnBottomSurface(geom) );
}
bool Line::isCornB(Geometry* geom) // xmax, zmin
{
	return (isOnBottomSurface(geom) && isOnRightSurface(geom) );
}
bool Line::isCornC(Geometry* geom) // xmax, zmax
{
	return (isOnRightSurface(geom) && isOnTopSurface(geom) );
}

bool Line::isCornD(Geometry* geom) // xmin, zmax
{
	return (isOnTopSurface(geom) && isOnLeftSurface(geom) );
}

