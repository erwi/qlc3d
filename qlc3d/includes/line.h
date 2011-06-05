#ifndef LINE_H
#define LINE_H


#include <geometry.h>

#include <vector>
#include <algorithm>



// Line is a simple line-element class joining two nodes (coordinates).
// It contains indexes to the two coordinates, their position in the 
// array of coordinates, p. 
class Line{
    public:
	int L[2];
	Line(); // dont use this, declared here only to keep STL happy
	Line(const int& a, const int& b);
	void PrintLine(); // prints node numbers only
	void PrintLine(Geometry* geom);  // prints node numbers and coordinates
	
	bool isOnFrontSurface(Geometry* geom);
	bool isOnBackSurface(Geometry* geom);
	bool isOnRightSurface(Geometry* geom);
	bool isOnLeftSurface(Geometry* geom);
	bool isOnTopSurface(Geometry* geom);
	bool isOnBottomSurface(Geometry* geom);
	bool isTopBottomCornerLine(Geometry* geom);
	bool isFrontBackCornerLine(Geometry* geom);
	bool isLeftRightCornerLine(Geometry* geom);
	bool isPeriodicTo(Geometry* geom, Line* L2);

	bool isTranslationOf( Line& L2, Geometry* geom, double* dir);// checks whether this line is a 
																//translation of L2 along drection 
																// dir = [1,0,0] ->compare x only, ignore y,z
																// [1,1,0] -> compare x and y only, ignore z
	// corners along z
	bool isCorn0(Geometry* geom); // xmin, ymin
	bool isCorn1(Geometry* geom); // xmax, ymin
	bool isCorn2(Geometry* geom); // xmax, ymax
	bool isCorn3(Geometry* geom); // xmin, ymax
	// corners along x
	bool isCorna(Geometry* geom); // ymin, zmin
	bool isCornb(Geometry* geom); // ymax, zmin
	bool isCornc(Geometry* geom); // ymax, zmax
	bool isCornd(Geometry* geom); // ymin, zmax
	// corners along y
	bool isCornA(Geometry* geom); // xmin, zmin
	bool isCornB(Geometry* geom); // xmax, zmin
	bool isCornC(Geometry* geom); // xmax, zmax
	bool isCornD(Geometry* geom); // xmin, zmax
	
	
			
	bool isBoundingBoxLine(Geometry* geom);

	bool operator<(const Line& other) const;
	bool operator==(const Line& other) const;
	double mymin(double x, double y){
	    return x <= y ? x:y;
	    
	}
}; // end class line


inline void uniquefy_line_vector( std::vector<Line>& v)
{/*! removes duplicate entries of vectors of lines
	SORT, REORDER, RESIZE
*/
	sort( v.begin() , v.end() ); // SORT lines
	std::vector<Line> :: iterator u;
	u = unique( v.begin() , v.end() ); // REORDER repeated lines after u
	v.erase(u, v.end() ); // RESIZE to only include unique lines
}


#endif
