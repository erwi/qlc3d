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
                        const SolutionVector& q,
                        const SolutionVector& v,
                        const LC& mat_par,
                        const Simu& simu,
                        const Mesh& t,
                        double* p
                        );



#endif // ASSEMBLEQ2K_H

