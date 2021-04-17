
#include <math.h>
#include <time.h>
#include <omp.h>
#include <qlc3d.h>
#include <cstdio>
#include <simulation-state.h>


// SpaMtrix INCLUDES
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_tickcounter.hpp>
#include <spamtrix_iterativesolvers.hpp>
#include <spamtrix_diagpreconditioner.hpp>

Simu::QMatrixSolvers selectQMatrixSolver(const Simu &simu, const LC &lc) {
    // SELECTS WHICH MATRIX SOLVER TO USE FOR Q-TENSOR
    // EQUATIONS. IF MATRIX IS SYMMETRIC USE PCG, ELSE
    // GMRES

    // IF SOLVER HAS ALREADY BEEN CHOSEN IN SETTINGS
    // FILE, DON'T DO ANYTHING
    if (simu.getQMatrixSolver()!=Simu::Auto) {
        return simu.getQMatrixSolver();
    }
    // SINGLE ELASTIC COEFF. EQUATIONS -> SYMMETRIC
    // CHIRALITY -> NON-SYMMETRIC
    bool isSymmetric = true;
    if (lc.K11() != lc.K22() ) { isSymmetric = false; }
    if (lc.K11() != lc.K33() ) { isSymmetric = false; }
    if (lc.p0() != 0.0) { isSymmetric = false; }

    return isSymmetric ? Simu::PCG : Simu::GMRES;
}

void solve_QTensor(SpaMtrix::IRCMatrix &K,
                   SpaMtrix::Vector &b,
                   SpaMtrix::Vector &x,
                   const Simu &simu,
                   const Settings &settings,
                   const LC &lc)
{
    SpaMtrix::DiagPreconditioner M(K);
    // Iterative solvers' settings...
    idx maxiter 	= settings.getQ_GMRES_Maxiter();
    idx restart 	= settings.getQ_GMRES_Restart();
    real toler      = settings.getQ_GMRES_Toler();
    SpaMtrix::IterativeSolvers solver(maxiter, restart, toler);
    Simu::QMatrixSolvers solverType = selectQMatrixSolver(simu, lc);

    if (solverType == Simu::PCG ){
        printf("PCG...");
        if (!solver.pcg(K,x,b,M) ){
            printf("PCG did not converge in %i iterations\nTolerance achieved is %f\n", solver.maxIter, solver.toler);
        }
    }
    else if (solverType == Simu::GMRES){
        printf("GMRES...");
        if (!solver.gmres(K,x,b,M)){
            printf("GMRES did not converge in %i iterations \nTolerance achieved is %f\n",solver.maxIter,solver.toler);
        }
    }
}


void setThreadCount(unsigned int nt)
{
#ifndef DEBUG
    omp_set_num_threads(nt);
#endif
}

void updateSolutionVector(SolutionVector &q,
                          const SpaMtrix::Vector &dq,
                          double &maxdq,
                          double &damping,
                          const Simu &simu )
{
    // UPDATES SOLUTION VECTOR q = q + dq
    // USES getEquNode REDIRECTION
    // FOR PERIODIC AND FIXED NODES

    // CALCULATE DAMPING COEFFICIENT IF GOING FOR STEADY STATE WITH NEWTON METHOD
    maxdq = fabs(dq[0]);
    if (simu.simulationMode() == SteadyState) {
        // FIND MAXIMUM CHANGE
        for (idx i = 0 ; i < dq.getLength() ; i++){
            maxdq = maxdq > fabs(dq[i]) ? maxdq:fabs(dq[i]);
        }

        if ( maxdq>simu.getTargetdQ() ){ // IF DAMPING IS NEEDED
            damping = simu.getTargetdQ() / maxdq;
        }
    }

    const idx npLC = q.getnDoF();
    for (unsigned int i = 0 ; i < 5 ; i++){   // LOOP OVER DIMENSIONS
        for (idx j = 0; j < npLC ; j++){ // LOOP OVER EACH NODE IN DIMENSION i
            const idx n = j + i*npLC;
            const idx effDoF = q.getEquNode(n);

            // EQUIVALENT DOF OF FIXED NODES ARE LABELLED AS "NOT_AN_INDEX"
            if (effDoF < NOT_AN_INDEX ){
                const double dqj = dq[ effDoF ];
                q.Values[n] -= dqj ;
                // KEEP TRACK OF LARGEST CHANGE IN Q-TENSOR
                //maxdq = fabs(dqj) > fabs(maxdq) ? dqj:maxdq;
                maxdq = fabs(dqj) > maxdq? fabs(dqj):maxdq;
            }
          //  cout << maxdq << "," << damping << endl;
        }// end for j
    }// end for i
}

double calcQ3d(SolutionVector *q,   // current Q-tensor
               SolutionVector *qn,  // previous Q-tensor
               SolutionVector *v,   // potential
               Geometry& geom,
               LC *mat_par,
               Simu* simu,
               SimulationState &simulationState,
               SpaMtrix::IRCMatrix &K,
               Settings* settings,
               Alignment* alignment)
{
    const idx numCols = K.getNumCols();
    double maxdq = 10;
    double maxdq_initial = 0;

    const idx isTimeStepping = simu->simulationMode() == TimeStepping ? 1:0;
    double timeStep = simulationState.dt();
    SpaMtrix::Vector RHS(isTimeStepping * numCols); // THIS IS ZERO SIZED IF NOT TIME STEPPING
    SpaMtrix::Vector dq(numCols);
    SpaMtrix::Vector L(numCols);

    int newton_iter = 0;    // COUNTER FOR NEWTON ITERATIONS

    // SAVE Q FROM PREVIOUS TIME STEP
    if (simu->simulationMode() == TimeStepping) {
        qn->setValuesTo(*q);
    }

    // CREATE MILLISECOND ACCURACY TIMER
    TickCounter <std::chrono::milliseconds> timer;
    timer.start();

    bool LOOP = true;
    while ( LOOP ){ // Q-TENSOR SOLUTION LOOP
        newton_iter++;
        timer.reset();
        // SET NUMBER OF THREADS USED IN ASSEMBLY
        setThreadCount(simu->getAssemblyThreadCount() );

        // ASSEMBLE RHS CONTRIBUTION FROM PREVIOUS TIME-STEP, IF NEEDED
        if ( (isTimeStepping) && (newton_iter == 1) ){
            RHS = 0;
            assemble_prev_rhs(RHS, *qn, *v, *mat_par, timeStep, geom );
        }
        // CLEAR MATRIX, L AND dq
        K = 0;
        L = 0;
        dq = 0;

        //======================================
        //  MATRIX ASSEMBLY
        //======================================

        std::cout << " " <<newton_iter << " Assembly..."; fflush(stdout);
        // ASSEMBLE MATRIX AND RHS (L)
        assembleQ(K, L, q, v, geom.t, geom.e, geom.getPtrTop(), mat_par, timeStep, alignment, geom.getPtrToNodeNormals());

        if (isTimeStepping) { // make Non-linear Crank-Nicholson RHS
            for (idx i = 0; i < numCols; i++) {
                L[i] += RHS[i];
            }
        }

        printf("OK %1.3es.", (float) timer.getElapsed() / 1000.0 );
        timer.reset();

        //======================================
        //  MATRIX SOLUTION
        //======================================
        // SET SOLVER THREAD COUNT
        setThreadCount(simu->getMatrixSolverThreadCount());
        // SOLVES Ax = b MATRIX PROBLEM
        solve_QTensor(K, L, dq, *simu, *settings, *mat_par);
        double damping = 1.0; // steady state damping coeff. calculated in updateSolutionVector
        updateSolutionVector(*q, dq, maxdq,damping,*simu); // q += damping*dq , taking into account periodic and fixed nodes

        if (newton_iter == 1) {
            maxdq_initial = maxdq; // maxdq_initial is needed elsewhere to adjust time-step size
        }

        // PRINT SOLUTION TIME
        printf("OK %1.3es.\tdQ = %1.3e\n", (float) timer.getElapsed() / 1000.0, maxdq);
        if (damping < 1.0){ // if damped, display by how much
            printf("damping=%1.3e\n", damping);
        }
        fflush( stdout );

        // PANIC!! if looks like no convergence
        if (newton_iter > settings->getQ_Newton_Panic_Iter() ){
            printf("Newton in distress!! - reducing time step by a factor of %f\n", settings->getQ_Newton_Panic_Coeff() );
            double newTimeStep = settings->getQ_Newton_Panic_Coeff() * simulationState.dt();
            simulationState.dt(newTimeStep);
            printf("new time step is : %e s.\n", simulationState.dt());
            newton_iter = 0;
            q->setValuesTo(*qn);
        }// end if PANIC!!!

        // DETERMINE WHETHER NEWTON LOOP IS DONE
        if (!isTimeStepping) {// IF dt == 0, GOING FOR STEADY STATE AND ONLY DOING ONE ITERATION -> NO LOOPS
            LOOP = false;
        } else if ( fabs(maxdq) < simu->getMaxError() ) { // EXIT IF ACCURATE ENOUGH
            LOOP = false;
        }
    } //end while LOOP
    return maxdq_initial;
}

