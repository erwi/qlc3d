
#include <math.h>
#include <time.h>
#include <omp.h>
#include <qlc3d.h>
#include <cstdio>
#include <assembleq2k.h>


// SpaMtrix INCLUDES
#include <ircmatrix.h>
#include <vector.h>
#include <tickcounter.h>

// THESE TWO ARE DEFINED IN solve_pcg.cpp
void solve_pcg(SpaMtrix::IRCMatrix &K, SpaMtrix::Vector &b, SpaMtrix::Vector &x ,Settings* settings);
void solve_gmres(SpaMtrix::IRCMatrix &K, SpaMtrix::Vector &b, SpaMtrix::Vector &x ,Settings* settings);

double calcQ3d(SolutionVector *q,   // current Q-tensor
               SolutionVector *qn,  // previous Q-tensor
               SolutionVector *v,   // potential
               Geometry& geom,
               LC *mat_par,
               Simu* simu,
               SpaMtrix::IRCMatrix &K,
               Settings* settings,
               Alignment* alignment)
//double* NodeNormals)
{
    const idx numCols = K.getNumCols();
    double maxdq = 10;
    double maxdq_previous = 10;
    double maxdq_initial = 0;

    const idx npLC = q->getnDoF();

    const idx isTimeStepping = simu->getdt() > 0 ? 1:0;

    SpaMtrix::Vector RHS(isTimeStepping*numCols); // THIS IS ZERO SIZED IF NOT TIME STEPPING
    SpaMtrix::Vector dq(numCols);
    SpaMtrix::Vector L(numCols);

    int newton_iter = 0;    // COUNTER FOR NEWTON ITERATIONS

    // SAVE Q FROM PREVIOUS TIME STEP
    if (simu->getdt() >0){
        qn->setValuesTo(*q);
    }

    // CREATE MILLISECOND ACCURACY TIMER
    TickCounter <std::chrono::milliseconds> timer;
    timer.start();

    bool LOOP = true;
    while ( LOOP ) // Q-TENSOR SOLUTION LOOP
    {
        newton_iter++;
        timer.reset();
        // ASSEMBLE RHS CONTRIBUTION FROM PREVIOUS TIME-STEP, IF NEEDED
        if ( (isTimeStepping) && (newton_iter == 1) )
        {
            RHS = 0;

            // CHOOSE WHICH VERSION TO CALL DEPENDING ON FORMULATION USED
            if (mat_par->PhysicsFormulation == LC::Nematic )
            {
                assemble_prev_rhs(RHS, *qn, *v, *mat_par, *simu, geom );
            }
            else if (mat_par->PhysicsFormulation == LC::BluePhase )
            {
                //printf(" BLUE PAHSE FORMULATION DISABLED IN %s\n",__func__);
                //exit(1);
                assemble_prev_rhs_K2(RHS, *qn, *v, *mat_par, *simu, geom);
            }
        }


        // CLEAR MATRIX, L AND dq
        K = 0; //->setAllValuesTo(0);
        L = 0;
        dq = 0;

        std::cout << " " <<newton_iter << " Assembly..."; fflush(stdout);
        // ASSEMBLE MATRIX AND RHS (L)
        assembleQ(K, L, q, v, geom.t, geom.e, geom.getPtrTop(), mat_par, simu, settings, alignment, geom.getPtrToNodeNormals());

        printf("OK %1.3es.", (float) timer.getElapsed() / 1000.0 );
        timer.reset();

        if (isTimeStepping) // make Non-linear Crank-Nicholson RHS
        {
#ifndef DEBUG
#pragma omp parallel for
#endif
            for (size_t i = 0 ; i < numCols ; i++)
                L[i] += RHS[i];
        }

        // SOLUTION
        if (settings->getQ_Solver() == Q_SOLVER_PCG)// use PCG for Q solution
        {
            printf("PCG..."); fflush( stdout );
            solve_pcg(K,L,dq,settings);
        }
        else if (settings->getQ_Solver() == Q_SOLVER_GMRES)// use GMRES for Q solution
        {
            printf("GMRES..."); fflush( stdout );
            solve_gmres(K,L,dq,settings);
        }

        maxdq = 0;

        // UPDATE SOLUTION VECTOR - USES getEquNode REDIRECTION FOR PERIODIC NODES
        maxdq_previous = maxdq;
        for (unsigned int i = 0 ; i < 5 ; i++)   // LOOP OVER DIMENSIONS
        {
            for (idx j = 0; j < npLC ; j++) // LOOP OVER EACH NODE IN DIMENSION i
            {
                const idx n = j + i*npLC;
                const idx effDoF = q->getEquNode(n);

                // EQUIVALENT DOF OF FIXED NODES ARE LABELLED AS "NOT_AN_INDEX"
                if (effDoF < NOT_AN_INDEX )
                {
                    const double dqj = dq[ effDoF ];
                    q->Values[n] -= dqj ;

                    // KEEP TRACK OF LARGEST CHANGE IN Q-TENSOR
                    maxdq = fabs(dqj) > fabs(maxdq) ? dqj:maxdq;
                }
            }// end for j
        }// end for i

        if (newton_iter==1) maxdq_initial = maxdq; // maxdq_initial is needed elsewhere to adjust time-step size

        // PRINT SOLUTION TIME
        printf("OK %1.3es.\tdQ = %1.3e\n", (float) timer.getElapsed() / 1000.0, maxdq);
        fflush( stdout );

        // PANIC!! if looks like no convergence
        if (newton_iter > settings->getQ_Newton_Panic_Iter() )
        {
            printf("Newton in distress!! - reducing time step by a factor of %f\n", settings->getQ_Newton_Panic_Coeff() );
            simu->setdt(settings->getQ_Newton_Panic_Coeff() * simu->getdt() );
            printf("new time step is : %e s.\n",simu->getdt());
            newton_iter = 0;
            q->setValuesTo(*qn);
        }// end if PANIC!!!

        // DETERMINE WHETHER NEWTON LOOP IS DONE
        if (!isTimeStepping) // IF dt == 0, GOING FOR STEADY STATE AND ONLY DOING ONE ITERATION -> NO LOOPS
            LOOP = false;
        else if ( fabs(maxdq) < simu->getMaxError() ) // EXIT IF ACCURATE ENOUGH
            LOOP = false;

    }//end while LOOP
    return maxdq_initial;
}

