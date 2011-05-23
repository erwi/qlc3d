#ifndef LINE_H
#define LINE_H


#include <geometry.h>

class Line{
    public:
	int L[2];
	Line(); // dont use this, declared here only to keep STL happy
	Line(const int& a, const int& b);
	void PrintLine();
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
			
	bool isBoundingBoxLine(Geometry* geom);

	bool operator<(const Line& other) const;
	bool operator==(const Line& other) const;
	double mymin(double x, double y){
	    return x <= y ? x:y;
	    //if (x<=y)
	//	return x;
	  //  else
		//return y;
	}
		
	
		
	}; // end class line

	


    //static bool operator<(const Line a, const Line b)
//	{
//		if ((a.L[0] < b.L[0]) && (a.L[0] < a.L[1]))
//			return true;
//		else
//		if ((a.L[0] == b.L[0]) && (a.L[1] < b.L[1]))
//			return true;
//
//		return false;
//	}
//
//	static bool operator==(const Line a, const Line b)
//	{
//		if ( (a.L[0] == b.L[0]) && ( a.L[1] == b.L[1]) )
//			return true;
//
//		return false;
//	}


#endif
