
#include <math.h>
#include <time.h>
#include <omp.h>
#include <qlc3d.h>
#include <cstdio>
#ifndef NO_QT
    #include <QTime>
#endif
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

    double* dq = (double*) malloc(q->getnFreeNodes() * q->getnDimensions() * sizeof(double) );
    double*  L = (double*) malloc(q->getnFreeNodes() * q->getnDimensions() * sizeof(double) );

    double* RHS = NULL;		// right hand side vector for Crank-Nicholson.
    if ( simu->getdt() > 0 )	// only allocate if needed
        RHS =  (double*) malloc(K->cols * sizeof(double) );



    int newton_iter = 0;

// SAVE Q FROM PREVIOUS TIME STEP
    // qn = q_old
    if (simu->getdt() >0){
        qn->setValuesTo(*q);
    }

    bool LOOP = true;

    while ( LOOP){// Q-TENSOR SOLUTION LOOP
	newton_iter++;
	// ASSEMBLY
	clock_t time1 = clock();
	#ifndef NO_QT // if using Qt
	    QTime time;
	    time.start();
	#endif
	// ASSEMBLE RHS CONTRIBUTION FROM PREVIOUS TIME-STEP
	if ( (newton_iter == 1) && (simu->getdt() > 0 ) ){
	    memset(RHS, 0, K->cols * sizeof(double) );
	    assemble_prev_rhs(RHS, *qn, *v, *t, *e, p , *mat_par, *simu);
	}



	    if ( simu->IsAssembleMatrix() )
		// CLEAR MATRIX, RHS AND dq
		K->setAllValuesTo(0);
		memset(L  , 0 , K->cols * sizeof(double) );
		memset(dq , 0 , K->cols * sizeof(double) );

		std::cout << " " <<newton_iter << " Assembly...";
		// ASSEMBLE MATRIX AND RHS
        assembleQ(K, L, q, v, t, e, p, mat_par, simu, settings, alignment, NodeNormals);

        //for (int i = 0 ; i < K->cols ; i++)
        //    printf("L[%i] = %e\n",i,L[i]);



		#ifdef DEBUG
		    K->DetectZeroDiagonals();
		#endif

		float elapsed  = 0;
		elapsed = ( (float) clock() - (float) time1 ) / (float) CLOCKS_PER_SEC; // get assembly time
		time1 = clock(); // used for solver timing next.
		#ifndef NO_QT // if using Qt
		    elapsed = ( (float) time.elapsed() ) / 1000.0;
		    time.restart();
		#endif
		printf("OK %1.3es. ", elapsed );
		fflush(stdout);

		if (simu->getdt() > 0){ // make Non-linear Crank-Nicholson RHS
		    #pragma omp parallel for
		    for (unsigned int i = 0 ; i < (unsigned int) K->cols ; i++){
			L[i] += RHS[i];
		    }

		}



		// SOLUTION
		if (settings->getQ_Solver() == Q_SOLVER_PCG){			// use PCG for Q solution
			printf("PCG...");
			fflush( stdout );
			solve_pcg(K,L,dq,settings);
		}
		else if (settings->getQ_Solver() == Q_SOLVER_GMRES){		// use GMRES for Q solution
			printf("GMRES...");
			fflush( stdout );
			solve_gmres(K,L,dq,settings);
		}

		maxdq = 0;

		// UPDATE SOLUTION VECTOR - USES getEquNode REDIRECTION FOR PERIODIC NODES
		maxdq_previous = maxdq;
		for (int i = 0 ; i < 5*npLC ; i++){
            //printf("dq[%i] = %e\n", i , dq[i]);
            if ( fabs(dq[q->getEquNode(i)]) > fabs(maxdq) )
			    maxdq = dq[q->getEquNode(i)]; // find max dQ

			q->Values[i] +=dq[q->getEquNode(i)];
		}

		if (newton_iter==1) maxdq_initial = maxdq; // maxdq_initial is needed elsewhere to adjust time-step size

		// PRINT SOLUTION TIME
		elapsed = ( (float) clock() - (float) time1 ) / (float) CLOCKS_PER_SEC ;
		#ifndef NO_QT // is using Qt
		    elapsed = ( (float) time.elapsed() ) / 1000.0 ;
		#endif

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

