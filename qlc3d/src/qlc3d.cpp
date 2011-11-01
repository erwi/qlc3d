	
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#include <filesysfun.h> // functions for handling directories

#include <qlc3d.h>
#include <simu.h>
#include <electrodes.h>
#include <lc.h>
#include <box.h>
#include <alignment.h>
#include <material_numbers.h>
#include <refinement.h>



int minnode(int *t, int dim1, int dim2){ // what does this do?
	int min = (int)1e9;
	for (int i =0 ; i < dim1 * dim2 ; i++ )
		if (t[i]<min) min = t[i];
	return min;
}
void AdjustTimeStep(Simu& simu, const double& maxdq){

    // if electrode switching etc. has just happened, dont adapt time step
    if ( simu.restrictedTimeStep ) return;

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

        cout << " scaling dt by: " << S << endl;
        double newdt = simu.getdt()* S;
        if ( newdt < simu.getMindt() ) newdt = simu.getMindt();
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
        //double consistency = 1e99; // arb. bignumber!
	double maxdq = 0;

	// in some cases Pot Cons loops are not needed (if gong for steady state or PotCons == Off)
	bool isPotCons = ( simu.getdt() > 0 ) && ( simu.getPotCons() != Off );
        isPotCons = false; // turn of consisteny loop
        calcpot3d(Kpot,&v, &q, &lc, geom1.t,geom1.e, geom1.getPtrTop(), &settings, &electrodes);
        maxdq = calcQ3d(&q,&qn,&v,geom1.t,geom1.e,geom1.getPtrTop(),&lc, &simu, Kq, &settings, &alignment, geom1.getPtrToNodeNormals());

	
	// DO THIS LOOP UNTIL POTENTIAL CONSISTENCY IS ACHIEVED
	/*
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
*/
	return maxdq;

}



void HandleEvents(EventList* eventlist, Simu* simu, SolutionVector* v, Electrodes* electrodes, Mesh* ee){
    if (simu->getdt() != 0){// if time stepping

    // time until next event
    double te =  eventlist->getNextEventTime() - simu->getCurrentTime()  ;
    simu->restrictedTimeStep = false; // no restriction needed on time step adaptation

    while ( te < 1e-20 ){ // if time to next event is ridiculously small (or has just passed, te < 0 )

        if  ( eventlist->getNextEventNumber() == EVENT_SWITCHING ){ // IF SWITCHING EVENT

            simu->restrictedTimeStep = true; // time step shoul be set to minimum after switchign
            simu->setdt( simu->getMindt() );

            v->setFixedNodesPot(electrodes,ee,simu->getCurrentTime() );  // set potential boundary conditions for all electrodes
            eventlist->RemoveFirstEvent();	// remove event from list

            // check if potential calculation needed
            bool needsCalcpot = false;
            double currentPot = 0;
            if ( electrodes->getnElectrodes() >= 2 ) {
                currentPot = electrodes->E[0]->getCurrentPotential();

                for (int i = 0 ; i < electrodes->getnElectrodes() ; i ++)// for each electrode
                {
                    // if differing potentials are applied to any electrode, potential must be calculated
                    if (electrodes->E[i]->getCurrentPotential() != currentPot)
                    {
                        needsCalcpot = true;
                        break;
                    }
                }
            }// end if potential calc needed
			
            electrodes->setCalcPot( needsCalcpot ); // sets flag whether to calculate potential
            if (! needsCalcpot) v->setValuesTo( currentPot ); // if no calculation needed, all values will be same
            if (eventlist->getLength() == 0) break; // leave if no more events to check
        }// end if EVENT_SWITCHING

        te = fabs( eventlist->getNextEventTime() - simu->getCurrentTime() ) ;

        }// end while next event

	// adjust time step if an event will occur during next time-step
	if (( te > 0  ) && (te < simu->getdt() )){ 
		simu->setdt(te);
		printf("adjusting time step to %e s. due to event\n", simu->getdt() );
	}
	}// end if time stepping



}//end void HandleEvents

int main(int argc, char* argv[]){


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
    simu.setCurrentDir( getCurrentDirectory() );
    if ( argc >= 2)
    {
        settings_filename.clear();
        settings_filename = argv[1];  // change if set by command line argument
    }
    if (argc >= 3 )// if working directory is defined
    {
        simu.setCurrentDir( argv[2] );
    }

    if ( !setCurrentDirectory( simu.getCurrentDir() ) )
    {
        printf("error - could not set working directory to:\n%s\nbye!", simu.getCurrentDir().c_str() );
        exit(1);
    }

    printf("current working directory:\n%s\n", getCurrentDirectory().c_str() );
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

// ==============================================
//
//	POTENTIAL SOLUTION DATA
//
//================================================

    EventList eventlist = EventList();
    eventlist.setElectrodeEvents(&electrodes);
	eventlist.printEventList();
	
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
    if (!simu.getLoadQ().empty() )
    {
        ReadLCD_B(&simu,&q);
        setStrongSurfacesQ(&q, &alignment, &lc, &geom1); // over writes surfaces with user specified values
    }
    q.setFixedNodesQ(&alignment, geom1.e);  // set fixed surface anchoring
    q.setPeriodicEquNodes(&geom1);          // periodic nodes
	
    q.EnforceEquNodes();		    // makes sure values at periodic boundaies match
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

	do{

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



            // SAVE
            WriteResult(&simu, &lc , &geom1, &v, &q);

            AdjustTimeStep(simu, maxdq); // also increments current time

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
	//#ifndef NO_QT
	//	cout << "Total iteration duration: " << (float) timer.elapsed()/1000<<"s." << endl;
	//#endif

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

//#ifndef NO_QT
//	a.quit();
//#endif



    if (v_cons) free(v_cons);
    if (q_cons) free(q_cons);
    if (qn_cons) free(qn_cons);
    return (0);
}
// end main
