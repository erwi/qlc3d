#include <settings.h>
#include <simu.h>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_iterativesolvers.hpp>
#include <spamtrix_luincpreconditioner.hpp>
#include <spamtrix_diagpreconditioner.hpp>

/*
void solve_QTensor(SpaMtrix::IRCMatrix &K,
                   SpaMtrix::Vector &b,
                   SpaMtrix::Vector &x,
                   const Simu &simu,
                   const Settings &settings)
{
    // GMRES settings...
    idx maxiter 	= settings.getQ_GMRES_Maxiter();
    idx restart 	= settings.getQ_GMRES_Restart();
    real toler      = settings.getQ_GMRES_Toler();
    SpaMtrix::IterativeSolvers solver(maxiter, restart, toler);

    if (simu.getQMatrixSolver() == Simu::PCG ){
        printf("PCG...");
        SpaMtrix::DiagPreconditioner M(K);
        if (!solver.pcg(K,x,b,M) ){
            printf("PCG did not converge in %i iterations\nTolerance achieved is %f\n", solver.maxIter, solver.toler);
        }
    }
    else if (simu.getQMatrixSolver() == Simu::GMRES){
        printf("GMRES...");
        SpaMtrix::DiagPreconditioner M(K);
        if (!solver.gmres(K,x,b,M)){
            printf("GMRES did not converge in %i iterations \nTolerance achieved is %f\n",solver.maxIter,solver.toler);
        }
    }
}
*/
