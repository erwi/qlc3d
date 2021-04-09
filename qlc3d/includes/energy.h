#ifndef ENERGY_H
#define ENERGY_H
#include <qlc3d.h>
class Simu;
class SolutionVector;
void CalculateFreeEnergy(FILE* fid,
                         int currentIteration,
                         double currentTime,
                         LC* lc,
                         Geometry* geom,
                         SolutionVector* v,
                         SolutionVector* q);
void closeEnergyFile(FILE* fid, Simu& simu);
#endif
