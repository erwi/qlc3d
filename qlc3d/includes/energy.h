#ifndef ENERGY_H
#define ENERGY_H
#include <qlc3d.h>
//#include <simu.h>

class Simu;
class SolutionVector;
void CalculateFreeEnergy(FILE* fid, Simu* simu, LC* lc, Geometry* geom, SolutionVector* v, SolutionVector* q );

void closeEnergyFile(FILE* fid, Simu& simu);
#endif
