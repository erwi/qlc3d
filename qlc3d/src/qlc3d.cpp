	
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#include <qlc3d.h>
#include <simu.h>
#include <electrodes.h>
#include <lc.h>
#include <box.h>
#include <alignment.h>
#include <material_numbers.h>
#include <refinement.h>

// make use of QT optional
#ifndef NO_QT
    #include <QtCore/QCoreApplication>
    #include <QTime>
    #include <QString>
#endif



int minnode(int *t, int dim1, int dim2){ // what does this do?
	int min = (int)1e9;
	for (int i =0 ; i < dim1 * dim2 ; i++ )
		if (t[i]<min) min = t[i];
	return min;
}
void AdjustTimeStep(Simu& simu, const double& maxdq){
	if (simu.getdt()>0){

		// ADAPT TIME STEP ACCORDING TO CONVERGENCE

		double R = simu.getTargetdQ() / fabs(maxdq);

		double Rf[4];
		simu.getdtFunction( &Rf[0] );

		double Rmin  = Rf[0]; // S = 0
		double RLmin = Rf[1]; // S = R (min)
		double RLmax = Rf[2]; // S = R (max)
		double Rmax  = Rf[3]; // S = 2

		double Smax = 2;
		double S = 1;

		// LINEAR REGION
		if ( (R < RLmax) && (R > RLmin) )
		{
		    //cout << "linear" << endl;
		    S = R;
		}
		else
		// BELOW LINEAR
		    if( R <= RLmin){
			double  k = RLmin / (RLmin - Rmin);
			double  dR = (RLmin - R);
			S = RLmin - k * dR;
		    }
		else
		// above LINEAR
		    if ( R >= RLmax){
			double k = (Smax - RLmax) / (Rmax - RLmax);
			double dr = R - RLmax;
			S = RLmax + k * dr;
			if (S>Smax) S = Smax;
		}

		cout << "S: " << S << endl;
		double newdt = simu.getdt()* S;
		if ( newdt < simu.getMindt() ) newdt = simu.getMindt();
		//cout << "newdt: " << newdt << endl;
		simu.setdt(newdt); // min/max limits taken care of here
	}

	simu.setCurrentChange(maxdq); // updates maximum change value for this iteration
}
void copyTo(double* arr, const SolutionVector& sv){ // copies all values from solution vector to array
	unsigned int n = sv.getnDoF() * sv.getnDimensions();
	#pragma omp parallel for
	for (unsigned int i = 0 ; i < n ; i++)
		arr[i] = sv.getValue(i);
}
double maxDiff(const double* arr, const SolutionVector& sv){// finds the maximum difference between array and solutionvector 
	unsigned int n = sv.getnDoF() * sv.getnDimensions();
	double md = 0;
	for (unsigned int i = 0 ; i < n ; i++){
		double diff = fabs( arr[i] - sv.getValue(i) );
		if (diff>md)
			md = diff;
	}

	return md;
}
double ConsistencyLoop(SolutionVector& v, SolutionVector& q , SolutionVector& qn,
					   double* v_cons, double* q_cons, double* qn_cons,
					   Geometry& geom1,
					   Simu& simu, LC& lc, Settings& settings,
					   Alignment& alignment, Electrodes& electrodes,
                       SparseMatrix* Kq, SparseMatrix* Kpot ){
	double consistency = 1e99; // arb. bignumber!
	double maxdq = 0;

	// in some cases Pot Cons loops are not needed (if gong for steady state or PotCons == Off)
	bool isPotCons = ( simu.getdt() > 0 ) && ( simu.getPotCons() != Off );
    isPotCons = false; // turn of consisteny loop
	// DO THIS LOOP UNTIL POTENTIAL CONSISTENCY IS ACHIEVED
	while( consistency > simu.getTargetPotCons() ){
		maxdq = 0;
		// make copies for potential consistency
		if (isPotCons){
			copyTo(v_cons , v);
			copyTo(q_cons , q);
			copyTo(qn_cons, qn);
		}

		// Q-Tensor
        maxdq = calcQ3d(&q,&qn,&v,geom1.t,geom1.e,geom1.getPtrTop(),&lc, &simu, Kq, &settings, &alignment, geom1.getPtrToNodeNormals());
		// Potential for consistency check

		calcpot3d(Kpot,&v, &q, &lc, geom1.t,geom1.e, geom1.getPtrTop(), &settings, &electrodes);


		if (!isPotCons) break; // can leave if no PotCons is needed


		// get max potential difference between before and after Q-tensor calculation
		consistency = maxDiff(v_cons, v);

		if ( consistency > simu.getTargetPotCons() ){
			cout << "\tpotential consistency :" << consistency <<" repeating iteration"<< endl;
			q.setValuesTo(q_cons);
			qn.setValuesTo(qn_cons);
		}
	}// end while

	return maxdq;

}



void HandleEvents(EventList* eventlist, Simu* simu, SolutionVector* v, Electrodes* electrodes, Mesh* ee){
	if (simu->getdt() != 0){// if time stepping

	double te =  eventlist->getNextEventTime() - simu->getCurrentTime()  ;
//	double next = eventlist->getNextEventTime();
//	double now = simu->getCurrentTime();
	//cout << "next, now , te: " <<next<<","<<now<<"," << te<<endl;

	while (( fabs(te) < 1e-20) ){ // if time to next event is ridiculously small
		//cout << "EVENT"<< endl;
	    if  ( eventlist->getNextEventNumber() == EVENT_SWITCHING ){ // IF SWITCHING EVENT
			//cout << "SWITCHING" << endl;
			if (simu->getCurrentIteration() > 0){
				//cout << "reducing time step to: " << simu->getMindt() << endl;
				simu->setdt( simu->getMindt() );	// reduce time step
				}
			v->setFixedNodesPot(electrodes,ee,simu->getCurrentTime() );  // set boundary conditions
			eventlist->RemoveFirstEvent();	// remove event from list

			if (eventlist->getLength() == 0) break; // leave if no more events to check

			if ( electrodes->getnElectrodes() >= 2 ) {// check if potential calculation needed
				double pot1 = electrodes->E[0]->Potential[0];
				bool calcpot = false;

				for (int i = 0 ; i < electrodes->getnElectrodes() ; i ++){// compare all electrode potentials
					if (electrodes->E[i]->Potential[0] != pot1){ // iff all same, no calculation needed
						calcpot = true;
						break;
					}
				}
				if (!calcpot){// if not needed, set all values to same
					electrodes->setCalcPot(calcpot);
					v->setValuesTo(pot1);
					}
			}// end if potential calc needed

	    }// end if EVENT_SWITCHING

	te = fabs( eventlist->getNextEventTime() - simu->getCurrentTime() ) ;
	}
	if (( te > 0  ) && (te < simu->getdt() )){ // adjust time step if event will occur during next step
		simu->setdt(te);
		printf("adjusting time step to %f ms.\n", simu->getdt() * 1e3);
	}
	}// end if time stepping
}//end void HandleEvents

int main(int argc, char* argv[]){

#ifndef NO_QT
	int nin = 1;
	char* chin[1];
	QCoreApplication a(nin,chin);//argc, argv);
#endif

    printf("\n\n\n");
    printf("=============================================================\n");
    printf("=                                                           =\n");
    printf("=                           QLC3D                           =\n");
    printf("=                        %s                        =\n",__DATE__);
    printf("= E. J. Willman, R. L. James, S. E. Day, F. A. Fernandez    =\n");
    printf("=============================================================\n");
    printf("\n\n\n");


//*
// Simulation settings (material parameters etc.)
	Simu		simu;
	Electrodes	electrodes;
	LC			lc;
	Boxes		boxes;
	Alignment	alignment;
	MeshRefinement  meshrefinement;
	string settings_filename = "./meshes/test.txt"; 	// default settings file, loaded when no i/p args.
	if ( argc >= 2){
		settings_filename.clear();
		settings_filename = argv[1];  // change if set by command line argument
	}
        ReadSettings(settings_filename,&simu,&lc,&boxes,&alignment,&electrodes, &meshrefinement);

	FILE* Energy_fid = NULL; // file for outputting free energy

	// Solver settings	(choose solver, preconditioner etc.)
	Settings settings = Settings();

	ReadSolverSettings("solver.qfg", &settings);


	// ================================================================
	//	CREATE GEOMETRY
        //	NEED 3 GEOMETRY OBJECTS WHEN USING MESH REFINEMENT
	// ================================================================
    Geometry geom1 = Geometry();	    // working geometry
    Geometry geom_orig = Geometry();    // original, loaded from file
    Geometry geom_prev = Geometry();    // geometry from previous ref. iteration
    prepareGeometry(geom_orig, simu);   // mesh file is read and geometry is loaded in this function (in inits.cpp)
    geom_prev.setTo( &geom_orig);	    // for first iteration, geom_prev = geom_orig
    geom1.setTo( &geom_orig);
    Refine(geom_orig, geom_prev, geom1 , &meshrefinement);

    geom_orig.setTo( &geom1);
    geom_prev.setTo( &geom_orig);	    // for first iteration, geom_prev = geom_orig
    geom1.setTo( &geom_orig);

    //geom1.getTotalSize();
//*/

// ==============================================
//
//	POTENTIAL SOLUTION DATA
//
//================================================

    EventList eventlist = EventList();
    eventlist.setElectrodeEvents(&electrodes);

    cout << "Creating V...";
    SolutionVector v( geom1.getnp() );
    v.setFixedNodesPot(&electrodes , geom1.e , simu.getCurrentTime());
    v.setPeriodicEquNodes(&geom1 ); // periodic nodes
    v.EnforceEquNodes(); // makes sure values at periodic boundaries match
	cout << "OK"<<endl;

// =============================================================
//
// 	SET LC INITIAL CONDITIONS
//
//==============================================================
    // create vector for 5 * npLC Q-tensor components
    cout << "Creating Q..." ;
    SolutionVector q(geom1.getnpLC(),5);    //  Q-tensor for current time step
    SolutionVector qn(geom1.getnpLC(),5);   //  Q-tensor from previous time step

    SetVolumeQ(&q, &lc, &boxes, geom1.getPtrTop());
	setSurfacesQ(&q, &alignment, &lc, &geom1);
		

//  LOAD Q FROM RESULT FILE
    if (!simu.getLoadQ().empty() ){
        ReadLCD_B(&simu,&q);
        setStrongSurfacesQ(&q, &alignment, &lc, &geom1); // over writes surfaces with user specified values
	}
    q.setFixedNodesQ(&alignment, geom1.e);  // set fixed surface anchoring
    q.setPeriodicEquNodes(&geom1);          // periodic nodes
    
	//q.PrintEquNodes();
	//exit(1);
	
	
	q.EnforceEquNodes();					// makes sure values at periodic boundaies match
	qn=q;                                   // q-previous = q-current in first iteration
    cout << "OK" << endl;                   // Q-TENSOR CREATED OK
 

// autoref( geom_orig, geom_prev, geom1, q,qn, v, meshrefinement, simu, alignment , electrodes, lc);

// Make matrices for potential and Q-tensor..
    cout << "Creating sparse matrix for potential...";
    SparseMatrix* Kpot = createSparseMatrix(geom1 , v);
    cout << "potential matrix OK"<< endl;

    cout << "Creating matrix for Q-tensor..." << endl;
    SparseMatrix* Kq = createSparseMatrix(geom1, q, MAT_DOMAIN1);
    cout << "Q-tensor matrix OK" << endl;

//********************************************************************
//*
//*		Save Initial configuration and potential
//*
//********************************************************************

	calcpot3d(Kpot,&v,&q,&lc,geom1.t,geom1.e, geom1.getPtrTop(), &settings, &electrodes);

// writes a copy of settings file in results directory
	printf("Saving a copy of the settings file...");
	CreateSaveDir(&simu);
	simu.setCurrentTime(0);
	WriteSettings(&simu, &lc, &boxes, &alignment, &electrodes);

	printf("\nSaving starting configuration (iteration -1)...\n");
	WriteResult(&simu, &lc , &geom1, &v, &q);
	printf("OK\n");

    //return 0;
	Energy_fid = createOutputEnergyFile(simu); // done in inits


	// TEMPORARY ARRAYS, USED IN POTENTIAL CONSISTENCY CALCULATIONS
    simu.setPotCons( Off );
    double* v_cons = NULL ;
	double* q_cons = NULL ;
	double* qn_cons= NULL ;
	if ( (simu.getPotCons() != Off) && (simu.getdt() > 0 ) ){
            v_cons  = (double*) malloc(v.getnDoF() * sizeof(double) );
            q_cons  = (double*) malloc(q.getnDoF() * q.getnDimensions() * sizeof(double) );
            qn_cons = (double*) malloc(q.getnDoF() * q.getnDimensions() * sizeof(double) );
	}

//===============================================
//	 Start simulation
//
//===============================================

	simu.setCurrentTime(0);
	time_t t1, t2;
	time(&t1);
	time(&t2);
	double maxdq = 0;
	// MAIN LOOP - while simulation is running
	#ifndef NO_QT
	QTime timer;
	timer.start();
	#endif
	do{
	#ifndef NO_QT
		timer.restart();
	#endif

            time(&t2);
            printf("=======================================================================\n");
            printf("Iteration %i, Time = %1.3es. Real time = %1.2es. dt = %1.3es. \n\b\b",simu.getCurrentIteration(),simu.getCurrentTime(),(float) t2-t1,simu.getdt());
            fflush(stdout);

            CalculateFreeEnergy(Energy_fid, &simu, &lc, &geom1, &v, &q);

            //EVENTS (ELECTRODES SWITCHING ETC.)
            HandleEvents(&eventlist, &simu, &v, &electrodes, geom1.e);

            // CALCULATES Q-TENSOR AND POTENTIAL
            maxdq = ConsistencyLoop(v,q,qn,v_cons,q_cons,qn_cons,
                                        geom1, simu, lc, settings,
                                        alignment, electrodes,
                                        Kq, Kpot);

            AdjustTimeStep(simu, maxdq);

            // SAVE
            WriteResult(&simu, &lc , &geom1, &v, &q);
            printf("\n");
/* AUTOREFINEMENT
        if (( (meshrefinement.isRefinementIteration( simu.getCurrentIteration() ) ) ))  // if this is a refinement iteration
             //|| (simu.getCurrentIteration() == 1) ) &&                               // or if first iteration this should be done before start of simulation
            //(needsRefinement(geom1,q,meshrefinement) ) ) {
            {
            delete Kpot;
            delete Kq;
            autoref(geom_orig, geom_prev, geom1, q, qn, v,  meshrefinement, simu,alignment, electrodes, lc);
            Kpot = createSparseMatrix(geom1, v );
            Kq   = createSparseMatrix(geom1, q, MAT_DOMAIN1);
            calcpot3d(Kpot,&v,&q,&lc,geom1.t,geom1.e, geom1.getPtrTop(), &settings, &electrodes);

            if (simu.getdt() > 0 )// if not steady state
            {
                simu.setdt( simu.getMindt() );
            }

            WriteResult(&simu, &lc, &geom1, &v, &q, &meshrefinement);
        }
//*/

	// LAZY....
	#ifndef NO_QT
		cout << "Total iteration duration: " << (float) timer.elapsed()/1000<<"s." << endl;
	#endif

	// if need end-refinement
    if ( (!simu.IsRunning() ) && (needsEndRefinement( geom1, q , meshrefinement ) ) )
	{
        printf("End refinement\n");

        delete Kpot;
        delete Kq;
                    
	  	autoref(geom_orig, geom_prev, geom1, q, qn, v,  meshrefinement, simu,alignment, electrodes, lc);
        Kpot = createSparseMatrix(geom1, v );
        Kq   = createSparseMatrix(geom1, q, MAT_DOMAIN1);
        calcpot3d(Kpot,&v,&q,&lc,geom1.t,geom1.e, geom1.getPtrTop(), &settings, &electrodes);

        WriteResult(&simu, &lc, &geom1, &v, &q, &meshrefinement);

        if ( simu.getdt() > 0.0 )
        {
            simu.setdt( simu.getMindt() );
        }
    }// end if need end-refinement

	}while ( simu.IsRunning() ); // end MAIN LOOP - while simulation is runnning








    printf("\nSaving final result file...\n");
    simu.setCurrentIteration( SIMU_END_SIMULATION );
    WriteResult(&simu, &lc , &geom1, &v, &q);
    printf("OK\n");

    closeEnergyFile(Energy_fid , simu);

    delete Kpot;
    delete Kq;

#ifndef NO_QT
	a.quit();
#endif



    if (v_cons) free(v_cons);
    if (q_cons) free(q_cons);
    if (qn_cons) free(qn_cons);
    return (0);
}
// end main
