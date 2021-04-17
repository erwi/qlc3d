#ifndef CALCPOT3D_H
#define CALCPOT3D_H
#include <solutionvector.h>
#include <lc.h>
#include <geometry.h>
#include <settings.h>
#include <electrodes.h>
namespace SpaMtrix {
class IRCMatrix;
}
void calcpot3d(SpaMtrix::IRCMatrix &Kpot,
               SolutionVector *v,
               SolutionVector *q,
               const LC &lc,
               Geometry &geom,
               Settings *settings,
               Electrodes *electrodes);
#endif // CALCPOT3D_H
