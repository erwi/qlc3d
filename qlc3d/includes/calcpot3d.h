#ifndef CALCPOT3D_H
#define CALCPOT3D_H

#include <sparsematrix.h>
#include <solutionvector.h>
#include <lc.h>
#include <geometry.h>
#include <settings.h>
#include <electrodes.h>

// SPAMTRIX INCLUDES
//#include <ircmatrix.h>

namespace SpaMtrix{
    class IRCMatrix;
}

void calcpot3d(SpaMtrix::IRCMatrix &Kpot,
               SolutionVector *v,
               SolutionVector *q,
               LC* lc,
               Geometry& geom,
               Settings* settings,
               Electrodes* electrodes);


#endif // CALCPOT3D_H
