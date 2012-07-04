
#include <math.h>
#include <time.h>
#include <omp.h>
#include <qlc3d.h>
#include <cstdio>
#include <assembleq2k.h>
double calcQ3d(SolutionVector *q,   // current Q-tensor
               SolutionVector *qn,  // previous Q-tensor
               SolutionVector *v,   // potential
               Geometry& geom,
               LC *mat_par,
               Simu* simu,
               SparseMatrix* K,
               Settings* settings,
               Alignment* alignment)
//double* NodeNormals)
{

    double maxdq = 10;
    double maxdq_previous = 10;
    double maxdq_initial = 0;

    idx npLC = q->getnDoF();
    double* RHS(NULL);		// right hand side vector for Crank-Nicholson.

    idx numCols = (idx) K->cols;
    double* dq = (double*) malloc( numCols * sizeof(double) );
    double*  L = (double*) malloc( numCols * sizeof(double) );
    if ( simu->getdt() > 0 )
        RHS =  (double*) malloc( numCols * sizeof(double) ); // ONLY NEEDED IF TIME-STEPPING

    int newton_iter = 0;    // COUNTER FOR NEWTON ITERATIONS

    // SAVE Q FROM PREVIOUS TIME STEP
    if (simu->getdt() >0){
        qn->setValuesTo(*q);
    }

    bool LOOP = true;

    while ( LOOP ) // Q-TENSOR SOLUTION LOOP
    {
        newton_iter++;
        clock_t time1 = clock();

        // ASSEMBLE RHS CONTRIBUTION FROM PREVIOUS TIME-STEP, IF NEEDED
        if ( (newton_iter == 1) && (simu->getdt() > 0 ) )
        {
            memset(RHS, 0, numCols*sizeof( double ) );

            // CHOOSE WHICH VERSION TO CALL DEPENDING ON FORMULATION USED
            if (mat_par->PhysicsFormulation == LC::Nematic )
            {
                assemble_prev_rhs(RHS, *qn, *v, *mat_par, *simu, geom );
            }
            else
            if (mat_par->PhysicsFormulation == LC::BluePhase )
            {
                printf(" BLE PAHSE FORMULATION DISABLED IN %s\n",__func__);
                exit(1);
                assemble_prev_rhs_K2(RHS, *qn, *v, *mat_par, *simu, geom);
            }
        }


        // CLEAR MATRIX, RHS AND dq
        K->setAllValuesTo(0);
        memset(L  , 0 , numCols * sizeof(double) );
        memset(dq , 0 , numCols * sizeof(double) );

        std::cout << " " <<newton_iter << " Assembly..."; fflush(stdout);
        // ASSEMBLE MATRIX AND RHS
        assembleQ(K, L, q, v, geom.t, geom.e, geom.getPtrTop(), mat_par, simu, settings, alignment, geom.getPtrToNodeNormals());

#ifdef DEBUG
        K->DetectZeroDiagonals();
#endif


        float elapsed  = 0;
        elapsed = ( (float) clock() - (float) time1 ) / (float) CLOCKS_PER_SEC; // get assembly time
        time1 = clock(); // used for solver timing next.

        printf("OK %1.3es. ", elapsed ); fflush(stdout);

        if (simu->getdt() > 0) // make Non-linear Crank-Nicholson RHS
        {
            //#ifndef DEBUG
            //#pragma omp parallel for
            //#endif
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
                idx n = j + i*npLC;
                idx effDoF = q->getEquNode(n);

                // EQUIVALENT DOF OF FIXED NODES ARE LABELLED AS "NOT_AN_INDEX"
                if (effDoF < NOT_AN_INDEX )
                {
                    double dqj = dq[ effDoF ];
                    q->Values[n] += dqj ;
                    if (fabs( dqj ) > fabs(maxdq) ) // KEEP TRACK OF LARGEST CHANGE IN Q-TENSOR
                        maxdq = dqj;
                }
            }// end for j
        }// end for i



        if (newton_iter==1) maxdq_initial = maxdq; // maxdq_initial is needed elsewhere to adjust time-step size

        // PRINT SOLUTION TIME
        elapsed = ( (float) clock() - (float) time1 ) / (float) CLOCKS_PER_SEC ;


        //elapsed = 0; maxdq = 0; //valgrind
        printf("OK %1.3es.\tdQ = %1.3e\n", elapsed, maxdq);
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
        if (simu->getdt() == 0) // IF dt == 0, GOING FOR STEADY STATE AND ONLY DOING ONE ITERATION -> NO LOOPS
        {
            LOOP = false;
        }
        else
            if ( fabs(maxdq) < simu->getMaxError() ) // EXIT IF ACCURATE ENOUGH
            {
                LOOP = false;
            }// end if
    }//end while dq > MaxError  ---- ENDS NEWTON LOOP


    free(dq);
    free(L);
    if (RHS) free(RHS); // free if using time stepping
    return maxdq_initial;
}

