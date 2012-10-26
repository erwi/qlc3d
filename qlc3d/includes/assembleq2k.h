#ifndef ASSEMBLEQ2K_H
#define ASSEMBLEQ2K_H
#include <ircmatrix.h>
#include <solutionvector.h>
#include <geometry.h>
#include <lc.h>
#include <simu.h>
#include <alignment.h>


void assemble_volumes2K(IRCMatrix &K,
                        Vector &L,
                        const SolutionVector& q,
                        const SolutionVector& v,
                        const LC& mat_par,
                        const Simu& simu,
                        const Mesh& t,
                        double* p
                        );

void assemble_prev_rhs_K2(Vector &RHS,
                         SolutionVector& qn,
                         SolutionVector& v,
                         LC& mat_par,
                         Simu& simu,
                         Geometry& geom );

#endif // ASSEMBLEQ2K_H

