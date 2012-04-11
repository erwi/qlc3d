#ifndef ASSEMBLEQ2K_H
#define ASSEMBLEQ2K_H
#include <sparsematrix.h>
#include <solutionvector.h>
#include <geometry.h>
#include <lc.h>
#include <simu.h>
#include <alignment.h>


void assemble_volumes2K(SparseMatrix& K,
                        double* L,
                        SolutionVector& q,
                        SolutionVector& v,
                        LC& mat_par,
                        Simu& simu,
                        Mesh& t,
                        double* p
                        );


#endif // ASSEMBLEQ2K_H

