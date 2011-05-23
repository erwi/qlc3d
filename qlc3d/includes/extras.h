#ifndef EXTRAS_H
#define EXTRAS_H

//#include <geometry.h>
//#include <solutionvector.h>
#include <qlc3d.h>
#include <stdio.h>

int FindElemByCoord(Geometry* geom , double px, double py, double pz);
void extras(Geometry* geom, SolutionVector* sol);

#endif
