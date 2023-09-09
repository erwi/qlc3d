#ifndef CALCPOT3D_H
#define CALCPOT3D_H
#include <solutionvector.h>
#include <lc.h>
#include <geometry.h>
#include <solver-settings.h>
#include <electrodes.h>
namespace SpaMtrix {
class IRCMatrix;
}
void calcpot3d(SpaMtrix::IRCMatrix &Kpot,
               SolutionVector &v,
               const SolutionVector &q,
               const LC &lc,
               const Geometry &geom,
               const SolverSettings &settings,
               const Electrodes &electrodes);
#endif // CALCPOT3D_H
