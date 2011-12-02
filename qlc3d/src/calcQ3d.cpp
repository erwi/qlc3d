
#include <math.h>
#include <time.h>
#include <omp.h>
#include <qlc3d.h>
#include <cstdio>
double calcQ3d(SolutionVector *q,   // current Q-tensor
	       SolutionVector* qn,  // previous Q-tensor
	       SolutionVector *v,   // potential
	       Mesh *t,		    // tet-mesh
	       Mesh *e,		    // tri-mesh
	       double *p,	    // node coordinates
	       LC *mat_par,
	       Simu* simu,
	       SparseMatrix* K,
               Settings* settings,
	       Alignment* alignment,
	       double* NodeNormals)
{

    double maxdq = 10;
    double maxdq_previous = 10;
    double maxdq_initial = 0;

    int npLC = q->getnDoF();
    double* RHS(NULL);		// right hand side vector for Crank-Nicholson.

    size_t numCols = (size_t) K->cols;
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

    while ( LOOP){// Q-TENSOR SOLUTION LOOP
        newton_iter++;
        clock_t time1 = clock();

        // ASSEMBLE RHS CONTRIBUTION FROM PREVIOUS TIME-STEP
        if ( (newton_iter == 1) && (simu->getdt() > 0 ) )
        {
            memset(RHS, 0, numCols*sizeof( double ) );
            assemble_prev_rhs(RHS, *qn, *v, *t, *e, p , *mat_par, *simu);
        }

        if ( simu->IsAssembleMatrix() )
            // CLEAR MATRIX, RHS AND dq
            K->setAllValuesTo(0);
            memset(L  , 0 , numCols * sizeof(double) );
            memset(dq , 0 , numCols * sizeof(double) );

            std::cout << " " <<newton_iter << " Assembly...";
            // ASSEMBLE MATRIX AND RHS
            assembleQ(K, L, q, v, t, e, p, mat_par, simu, settings, alignment, NodeNormals);

            #ifdef DEBUG
                K->DetectZeroDiagonals();
            #endif

/*
    // DEBUG BLOCK
     {
     printf("Num fixed nodes = %i\n", q->getnFixed() );
     for (int i = 0 ; i < q->getnDoF() ; i++ )
      {
          if
              printf("L[%i] = %e, K[%i,%i] = %e\n",i, L[i], i,i,K->sparse_get(i,i));

       }
     }//*/
    // END DEBUG BLOCK



            float elapsed  = 0;
            elapsed = ( (float) clock() - (float) time1 ) / (float) CLOCKS_PER_SEC; // get assembly time
            time1 = clock(); // used for solver timing next.

            printf("OK %1.3es. ", elapsed );
            fflush(stdout);

            if (simu->getdt() > 0) // make Non-linear Crank-Nicholson RHS
            {
                #pragma omp parallel for
                for (unsigned int i = 0 ; i < (unsigned int) K->cols ; i++)
                    L[i] += RHS[i];
            }

            // SOLUTION
            if (settings->getQ_Solver() == Q_SOLVER_PCG)// use PCG for Q solution
            {
                printf("PCG...");
                fflush( stdout );
                solve_pcg(K,L,dq,settings);
            }
            else if (settings->getQ_Solver() == Q_SOLVER_GMRES)// use GMRES for Q solution
            {
                printf("GMRES...");
                fflush( stdout );
                solve_gmres(K,L,dq,settings);
            }
 /*
            {// DEBUG

                printf("\n");
                for (int i = 0 ; i < 10 ; i++)
                {
                    printf("dq[%i] = %e\n", i, dq[i]);
                }

            }// END DEBUG
  */




            maxdq = 0;

            // UPDATE SOLUTION VECTOR - USES getEquNode REDIRECTION FOR PERIODIC NODES
            maxdq_previous = maxdq;
            //int indMax = 0;
            for (int i = 0 ; i < 5 ; i++)   // LOOP OVER DIMENSIONS
            {
                for (int j = 0; j < npLC ; j++) // LOOP OVER EACH NODE IN DIMENSION i
                {
                    int n = j + i*npLC;
                    double dqj = dq[ q->getEquNode(n) ];

                    //#ifdef DEBUG
                    if ( fabs( dqj ) > fabs(maxdq) )// find max dQ
                    {
                        maxdq = dqj;
                        //indMax = j;
                        //printf("\nindMax = %i, dim %i\n", indMax, i);
                    }
                    //#endif

                    q->Values[n] += dqj ;
                }// end for j
            }// end for i

                /// DEBUG BLOCK
                /*{
                printf("\n----------------------------------------------\n");
                printf("indMax = %i, maxdq = %e\n", indMax, maxdq);
                printf("a"); fflush(stdout );
                printf("p[indMax] = [%e,%e,%e]\n", p[indMax*3], p[indMax*3+1], p[indMax*3+2] );

                int equMax = q->getEquNode(indMax);
                printf("equMax = %i, at p = [%e,%e,%e]\n", equMax, p[equMax*3], p[equMax*3+1], p[equMax*3+2]);
                printf("\n----------------------------------------------\n");
                printf("a"); fflush(stdout );
                }*/
                /// END DEBUG BLOCK

                if (newton_iter==1) maxdq_initial = maxdq; // maxdq_initial is needed elsewhere to adjust time-step size

		// PRINT SOLUTION TIME
		elapsed = ( (float) clock() - (float) time1 ) / (float) CLOCKS_PER_SEC ;


		//elapsed = 0; maxdq = 0; //valgrind
		printf("OK %1.3es.\tdQ = %1.3e\n", elapsed, maxdq);
		fflush( stdout );

		// PANIC!! if looks like no convergence
		if (newton_iter > settings->getQ_Newton_Panic_Iter() ){
			printf("Newton in distress!! - reducing time step by a factor of %f\n", settings->getQ_Newton_Panic_Coeff() );
			simu->setdt(settings->getQ_Newton_Panic_Coeff() * simu->getdt() );
			printf("new time step is : %f ms.\n",simu->getdt() * 1e3);
			newton_iter = 0;
			q->setValuesTo(*qn);
		}// end if PANIC!!!


// DETERMINE WHETHER NEWTON LOOP IS DONE
		if (simu->getdt() == 0) // IF dt == 0, GOING FOR STEADY STATE AND ONLY DOING ONE ITERATION -> NO LOOPS
			LOOP = false;
		else
        if ( fabs(maxdq) < simu->getMaxError() ){ // EXIT IF ACCURATE ENOUGH
            //printf("maxdq = %f, maxError = %f\n", fabs(maxdq), simu->getMaxError() );
            LOOP = false;
        }// end if
	}//end while dq > MaxError  ---- ENDS NEWTON LOOP

	if ( simu->getdt() > 0 ){
	    simu->IncrementCurrentTime();
    }// end if dt>0

	simu->IncrementCurrentIteration();

	free(dq);
	free(L);
	if (RHS) free(RHS); // free if using time stepping
	return maxdq_initial;
}

