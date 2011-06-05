#include "../includes/line.h"

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

Line::Line(){
    // cout << "error- Line::Line() should never be called. Use Line::Line( a, b ) instead - bye!" << endl;
    // exit(1);
}

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


void Line::PrintLine(){
    printf("[%i,%i]\n", L[0], L[1]);
}
void Line::PrintLine(Geometry* geom)
{
	printf("%i,%i=[%f,%f,%f], [%f,%f,%f]\n", L[0], L[1], 
		geom->getpX( L[0] ),	geom->getpY( L[0] ),	geom->getpZ( L[0] ),
		geom->getpX( L[1] ),	geom->getpY( L[1] ),	geom->getpZ( L[1] ) );
}

bool Line::isBoundingBoxLine(Geometry* geom){// returns true if this line is on a surface of the external boundix box of the modelling window
    double x1 = geom->getpX(L[0]);
    double y1 = geom->getpY(L[0]);
    double z1 = geom->getpZ(L[0]);
			
    double x2 = geom->getpX(L[1]);
    double y2 = geom->getpY(L[1]);
    double z2 = geom->getpZ(L[1]);
			
    // if left side
    if ( (x1 == geom->getXmin() ) && (x2 == geom->getXmin() ) )
	return true;
    // if right side
    else
    if ( (x1 == geom->getXmax() ) && (x2 == geom->getXmax() ) )
	return true;
    // if front
    else
    if ( (y1 == geom->getYmin() ) && (y2 == geom->getYmin() ) )
	return true;
    // if back
    else
    if ( (y1 == geom->getYmax() ) && (y2 == geom->getYmax() ) )
	return true;
    //if bottom
    else
    if ( (z1 == geom->getZmin() ) && (z2 == geom->getZmin() ) )
	return true;
    //if top
    else
    if ( (z1 == geom->getZmax() ) && (z2 == geom->getZmax() ) )
	return true;
    // otherwise NO
    else
	return false;
    }// end isBoundingBoxLine()

bool Line::isOnFrontSurface(Geometry* geom){
				
    double y1 = geom->getpY(L[0]);
    double y2 = geom->getpY(L[1]);

    if ( ( y1==geom->getYmin() ) && (y2 == geom->getYmin() ) )
	return true;
    else
	return false;
}
		
bool Line::isOnBackSurface(Geometry* geom){
    double y1 = geom->getpY(L[0]);
    double y2 = geom->getpY(L[1]);
    if ( ( y1==geom->getYmax() ) && (y2 == geom->getYmax() ) )
	return true;
    else
	return false;
}

bool Line::isOnRightSurface(Geometry* geom){
    double x1 = geom->getpX(L[0]);
    double x2 = geom->getpX(L[1]);
    if ( (x1 == geom->getXmax() ) && (x2 == geom->getXmax() ))
	return true;
    else
    	return false;
		
}
bool Line::isOnLeftSurface(Geometry* geom){
    double x1 = geom->getpX(L[0]);
    double x2 = geom->getpX(L[1]);
    if ( (x1 == geom->getXmin() ) && (x2 == geom->getXmin() ))
	return true;
    else
	return false;
}
bool Line::isOnTopSurface(Geometry* geom){
    double z1 = geom->getpZ(L[0]);
    double z2 = geom->getpZ(L[1]);
    if ( (z1 == geom->getZmax() ) && (z2 == geom->getZmax() ) )
	return true;
    else
	return false;
}
bool Line::isOnBottomSurface(Geometry* geom){
    double z1 = geom->getpZ(L[0]);
    double z2 = geom->getpZ(L[1]);
    if ( (z1 == geom->getZmin() ) && (z2 == geom->getZmin() ) )
	return true;
    else
	return false;
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
	
	
bool Line::isPeriodicTo(Geometry* geom, Line* L2)
{// compares if this line is periodic (on opposite face) to line L2
	// if comparison with self, return false
	if ( (L[0] == L2->L[0]) && (L[1] == L2->L[1]) )
	    return false;

	double eps = 1e-6; // accuracy for cordinate comparison (sometimes numerical noise from GiD)
	// get coordinates for this line
	double x1 = geom->getpX(L[0]);
	double y1 = geom->getpY(L[0]);
	double z1 = geom->getpZ(L[0]);
	double x2 = geom->getpX(L[1]);
	double y2 = geom->getpY(L[1]);
	double z2 = geom->getpZ(L[1]);
			
	// get coordintes for L2
	double X1 = geom->getpX(L2->L[0]);
	double Y1 = geom->getpY(L2->L[0]);
	double Z1 = geom->getpZ(L2->L[0]);
	double X2 = geom->getpX(L2->L[1]);
	double Y2 = geom->getpY(L2->L[1]);
	double Z2 = geom->getpZ(L2->L[1]);


// CASE1: ONLY FRONT BACK PERIODIC
	if ( ( geom->getfront_back_is_periodic()) && (!geom->getleft_right_is_periodic()) && (!geom->gettop_bottom_is_periodic()) ){
	    if (  isOnFrontSurface(geom) ||  isOnBackSurface(geom)  ){
	    // compare x and z coordinates
		if ( ( fabs(x1 - X1) < eps ) && ( fabs(z1 - Z1) < eps) && (fabs(x2 - X2) < eps) && (fabs(z2 - Z2) < eps) )
		    return true;
	    }
	}// end CASE1
			
// CASE2: FRONT/BACK and LEFT RIGHT ARE PERIODIC
	else
	if ( (geom->getfront_back_is_periodic() ) && (geom->getleft_right_is_periodic() ) && (!geom->gettop_bottom_is_periodic() ) ){
	    // if both lines are corner lines
	    if (isTopBottomCornerLine(geom) ){
		if (L2->isTopBottomCornerLine(geom) ){
		// both lines are corner linews -> only need to compare z-coordinates
		if ( ((fabs(z1-Z1)< eps)|| (fabs(z1-Z2) < eps) ) && ( ( fabs(z2 - Z2) < eps ) || ( fabs(z2 - Z1 ) < eps ) ) )
		    return true;
		}
	    }
	else
	if ( isOnFrontSurface(geom) || isOnBackSurface(geom) ){ // front-back
	    // compare x and z coordinates (same as CASE1)
	    double dx1 = mymin( fabs(x1-X1) , fabs(x1-X2) );
	    double dx2 = mymin( fabs(x2-X1) , fabs(x2-X2) );
	    double dz1 = mymin( fabs(z1-Z1) , fabs(z1-Z2) );
	    double dz2 = mymin( fabs(z2-Z1) , fabs(z2-Z2) );
		
	    //if ( (L[0] == 4) && (L[1] == 5) ){
	    //	if (L2->L[0] == 2)
	    //    L2->PrintLine();
	    //	}
	    if ( (dx1<eps) && (dx2<eps) && (dz1<eps) && (dz2<eps) )
		return true;
	    }
	    else
	    if (isOnLeftSurface(geom) || isOnRightSurface(geom) ){ // left-right
		// compare y and z coordinates
		double dy11 = fabs(y1-Y1); // consider both orientations of lines
		double dy21 = fabs(y2-Y2);
		double dy12 = fabs(y1-Y2);
		double dy22 = fabs(y2-Y1);
				
		double dz11 = fabs(z1-Z1);
		double dz21 = fabs(z2-Z2);
		double dz12 = fabs(z1-Z2);
		double dz22 = fabs(z2-Z1);
				
	    if ( (dy11<eps || dy12<eps) && (dy21<eps||dy22<eps) && (dz11<eps||dz12<eps) && (dz21<eps||dz22<eps) )
		return true;
			
	    }
	    else{
		// line is not on front/back/left/right
		printf("error - Line::IsPEriodicTo - CASE2 - bye\n");
		exit(1);
	    }
			
	}// end CASE2
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
	double eps = 1e-5;
	double* p11 = geom->getPtrTop(this->L[0]);
	double* p12 = geom->getPtrTop(this->L[1]);
	
	double* p21 = geom->getPtrTop(L2.L[0]);
	double* p22 = geom->getPtrTop(L2.L[1]);
	
	double v1[3] = {p11[0] - p12[0] ,
					p11[1] - p12[1] ,
					p11[2] - p12[2]  };

	double v2[3] = {p21[0] - p22[0],
					p21[1] - p22[1],
					p21[2] - p22[2] };
	// line (vector) lengths
	double l1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2] ); 
	double l2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2] );
				
	double dot = ( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] )/l1/l2;				
	
	if ( (fabs( fabs(dot) -1.0 ) <= eps ) &&// if v1 and v2 are parallel AND
	     (fabs(l1-l2) <= eps ) ) //if v1 and v2 are same length
	{ 
		// calculate differences between line end points
		double diff1[3] = { p11[0], p11[1], p11[2] };
		double diff2[3] = { p12[0], p12[1], p12[2] };
		
		// calculate shifts. differences depend on line orientations
		if (dot > 0) // same direction
		{
			diff1[0] -= p21[0]; // difference between nodes 1
			diff1[1] -= p21[1];
			diff1[2] -= p21[2];
			
			diff2[0] -= p22[0]; // difference between nodes 2
			diff2[1] -= p22[1];
			diff2[2] -= p22[2];
		}
		else // opposite directions
		{
			diff1[0] -= p22[0];
			diff1[1] -= p22[1];
			diff1[2] -= p22[2];
			
			diff2[0] -= p21[0];
			diff2[1] -= p21[1];
			diff2[2] -= p21[2];
		}
		
		// set x,y,z components of shift to zero, as defined in dir[0]->dir[2]		
		diff1[0] = ( fabs(diff1[0]) + fabs(diff2[0] ) )*dir[0]; 
		diff1[1] = ( fabs(diff1[1]) + fabs(diff2[1] ) )*dir[1];
		diff1[2] = ( fabs(diff1[2]) + fabs(diff2[2] ) )*dir[2];
		
		//printf("%f,%f,%f\n", diff1[0], diff1[1], diff1[2]);
		if ( ( diff1[0]+diff1[1]+diff1[2] ) <= eps  )
		{
			return true;
		}
	}// end if same length and parallel
	return false;
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

