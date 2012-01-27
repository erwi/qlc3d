	
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
#include <regulargrid.h>
#include <string>
#include <eventhandler.h>
#include <calcpot3d.h>
#include <resultoutput.h>
int minnode(int *t, int dim1, int dim2){ // what does this do?
	int min = (int)1e9;
	for (int i =0 ; i < dim1 * dim2 ; i++ )
		if (t[i]<min) min = t[i];
	return min;
}
void adjustTimeStepSize(Simu& simu, const double& maxdq){

      //simu.IncrementCurrentIteration();

    // if electrode switching etc. has just happened, dont adapt time step this time
    // but set switch to false to allow adjutment starting next step
    if ( simu.restrictedTimeStep )
    {
        simu.restrictedTimeStep = false;
        return;
    }

    if (simu.getdt()>0){

        // INCREMENT CURRENT TIME BEFORE CHANGING TIME STEP SIZE
        // dt SHOULD BE CORRECT HERE TO MAKE SURE THAT A TIME STEP
        // COINCIDES WITH AN EVENT, IF ANY.

        /// simu.IncrementCurrentTime();

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

//void incrementIteration(Simu& simu)
//{
//    simu.IncrementCurrentIteration();
//    simu.setdt( simu.getdt() + simu.getCurrentTime() );
//}


// THESE SHOULD NOT BE HERE!!
void copyTo(double* arr, const SolutionVector& sv){ // copies all values from solution vector to array
	unsigned int n = sv.getnDoF() * sv.getnDimensions();
	#pragma omp parallel for
	for (unsigned int i = 0 ; i < n ; i++)
		arr[i] = sv.getValue(i);
}

// THESE SHOULD NOT BE HERE!!
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
double updateSolutions(SolutionVector& v, SolutionVector& q , SolutionVector& qn,
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
        calcpot3d(Kpot,&v, &q, &lc, geom1, &settings, &electrodes);
        maxdq = calcQ3d(&q,&qn,&v,geom1,&lc, &simu, Kq, &settings, &alignment );

	return maxdq;

}



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
    Simu	simu;
    Electrodes	electrodes;
    LC		lc;
    Boxes	boxes;
    Alignment	alignment;
    MeshRefinement  meshrefinement;
    RegularGrid regGrid;
    EventList eventlist;
    string settings_filename = "./meshes/test.txt"; 	// default settings file, loaded when no i/p args.
    simu.setCurrentDir( FilesysFun::getCurrentDirectory() );




    if ( argc >= 2)
    {
        settings_filename.clear();
        settings_filename = argv[1];  // change if set by command line argument
    }
    if (argc >= 3 )// if working directory is defined
    {
        simu.setCurrentDir( argv[2] );
    }

    if ( !FilesysFun::setCurrentDirectory( simu.getCurrentDir() ) )
    {
        printf("error - could not set working directory to:\n%s\nbye!", simu.getCurrentDir().c_str() );
        exit(1);
    }

    printf("current working directory:\n%s\n", FilesysFun::getCurrentDirectory().c_str() );

    ReadSettings(settings_filename, // CHANGE POINTERS TO REFS.
                 &simu,
                 lc,
                 &boxes,
                 &alignment,
                 &electrodes,
                 &meshrefinement,
                 eventlist);


    FILE* Energy_fid = NULL; // file for outputting free energy

    // Solver settings	(choose solver, preconditioner etc.)
    Settings settings = Settings();
    ReadSolverSettings("solver.qfg", &settings);


    // ================================================================
    //	CREATE GEOMETRY
    //	NEED 3 GEOMETRY OBJECTS WHEN USING MESH REFINEMENT
    // ================================================================
    Geometry geom1 = Geometry();	// working geometry
    Geometry geom_orig = Geometry();    // original, loaded from file
    Geometry geom_prev = Geometry();    // geometry from previous ref. iteration
    prepareGeometry(geom_orig, simu, alignment);   // mesh file is read and geometry is loaded in this function (in inits.cpp)
    geom_prev.setTo( &geom_orig);	    // for first iteration, geom_prev = geom_orig
    geom1.setTo( &geom_orig);
    Refine(geom_orig, geom_prev, geom1 , &meshrefinement);

    geom_orig.setTo( &geom1);
    geom_prev.setTo( &geom_orig);	    // for first iteration, geom_prev = geom_orig
    geom1.setTo( &geom_orig);

    // SET CONVENIENCE STRUCT OF POINTERS
    Geometries geometries;
    geometries.geom = &geom1;
    geometries.geom_orig = &geom_orig;
    geometries.geom_prev = &geom_prev;

// ==============================================
//
//	POTENTIAL SOLUTION DATA
//
//================================================

    cout << "Creating V...";
    SolutionVector v( geom1.getnp() );

    v.allocateFixedNodesArrays( geom1 );

    v.setPeriodicEquNodes( &geom1 ); // periodic nodes

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
        WriteResults::ReadLCD_B(&simu,&q);
        setStrongSurfacesQ(&q, &alignment, &lc, &geom1); // over writes surfaces with user specified values
    }
    q.setFixedNodesQ(&alignment, geom1.e);  // set fixed surface anchoring
    q.setPeriodicEquNodes(&geom1);          // periodic nodes
    q.EnforceEquNodes(geom1);		    // makes sure values at periodic boundaies match
    qn=q;                                   // q-previous = q-current in first iteration
    cout << "OK" << endl;                   // Q-TENSOR CREATED OK

    // SET CONVENIENCE POINTERS STRUCTURE
    SolutionVectors solutionvectors;
    solutionvectors.q = &q;
    solutionvectors.qn = &qn;
    solutionvectors.v = &v;

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
    printf("\nSaving starting configuration (iteration -1)...\n");
    WriteResults::CreateSaveDir(simu);
    Energy_fid = createOutputEnergyFile(simu); // done in inits
   // handleInitialEvents(eventlist,
   //              electrodes,
   //              alignment,
   //              simu,
   //              geometries,
   //              solutionvectors,
   //              Kpot,
   //              lc,
   //              settings);
    printf("OK\n");


    // TEMPORARY ARRAYS, USED IN POTENTIAL CONSISTENCY CALCULATIONS
    simu.setPotCons( Off );
    double* v_cons = NULL ;
    double* q_cons = NULL ;
    double* qn_cons= NULL ;
    if ( (simu.getPotCons() != Off) && (simu.getdt() > 0 ) )
    {
        v_cons  = (double*) malloc(v.getnDoF() * sizeof(double) );
        q_cons  = (double*) malloc(q.getnDoF() * q.getnDimensions() * sizeof(double) );
        qn_cons = (double*) malloc(q.getnDoF() * q.getnDimensions() * sizeof(double) );
    }

//===============================================
//	 Start simulation
//
//===============================================


        simu.setCurrentIteration( 1 );
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

        // mve this to event handling/result output
        CalculateFreeEnergy(Energy_fid, &simu, &lc, &geom1, &v, &q);

        // CALCULATES Q-TENSOR AND POTENTIAL
        maxdq = updateSolutions(v,q,qn,v_cons,q_cons,qn_cons,
                                    geom1, simu, lc, settings,
                                    alignment, electrodes,
                                    Kq, Kpot);
        /// UPDATE CURRENT TIME
        simu.setCurrentTime( simu.getCurrentTime() + simu.getdt() );
        /// CALCULATE NEW TIME STEP SIZE BASED ON maxdq
        adjustTimeStepSize( simu, maxdq );
        //EVENTS (ELECTRODES SWITCHING ETC.)
        handleEvents(eventlist,
                     electrodes,
                     alignment,
                     simu,
                     geometries,
                     solutionvectors,
                     Kpot,
                     lc,
                     settings);


        /// INCREMENT ITERATION COUNTER
        simu.IncrementCurrentIteration();
//* AUTOREFINEMENT ---- MOVE TO EVENTS
        /*
        if (( (meshrefinement.isRefinementIteration( simu.getCurrentIteration() ) ) ))  // if this is a refinement iteration
             //|| (simu.getCurrentIteration() == 1) ) &&                               // or if first iteration this should be done before start of simulation
            //(needsRefinement(geom1,q,meshrefinement) ) ) {
            {
            delete Kpot;
            delete Kq;
            autoref(geom_orig, geom_prev, geom1, q, qn, v,  meshrefinement, simu,alignment, electrodes, lc);
            Kpot = createSparseMatrix(geom1, v );
            Kq   = createSparseMatrix(geom1, q, MAT_DOMAIN1);
            calcpot3d(Kpot,&v,&q,&lc,geom1, &settings, &electrodes);

            if (simu.getdt() > 0 )// if not steady state
            {
                simu.setdt( simu.getMindt() );
            }

            WriteResults::WriteResult(&simu, &lc, &geom1, &v, &q, &meshrefinement);
        }
*/

	// if need end-refinement
        /*
    if ( (!simu.IsRunning() ) && (needsEndRefinement( geom1, q , meshrefinement ) ) )
	{
        printf("End refinement\n");

        delete Kpot;
        delete Kq;
                    
        autoref(geom_orig, geom_prev, geom1, q, qn, v,  meshrefinement, simu,alignment, electrodes, lc);
        Kpot = createSparseMatrix(geom1, v );
        Kq   = createSparseMatrix(geom1, q, MAT_DOMAIN1);
        calcpot3d(Kpot,&v,&q,&lc,geom1, &settings, &electrodes);

        FilesysFun::setCurrentDirectory( simu.getSaveDir() );
        WriteResults::WriteResult(&simu, &lc, &geom1, &v, &q, &meshrefinement);
        FilesysFun::setCurrentDirectory( simu.getCurrentDir() );

        if ( simu.getdt() > 0.0 )
        {
            simu.setdt( simu.getMindt() );
        }
    }// end if need end-refinement
    */
}while ( simu.IsRunning() ); // end MAIN LOOP - while simulation is runnning

    printf("\nSaving final result file...\n");
    simu.setCurrentIteration( SIMU_END_SIMULATION );
    FilesysFun::setCurrentDirectory( simu.getSaveDir() ); // GOTO OUPUT DIR
    WriteResults::WriteResult(&simu, &lc , &geom1, &v, &q);
    FilesysFun::setCurrentDirectory( simu.getCurrentDir() );// GO BACK TO EXECUTABLE DIR
    printf("OK\n");
    closeEnergyFile(Energy_fid , simu);



    delete Kpot;
    delete Kq;

    if (v_cons) free(v_cons);
    if (q_cons) free(q_cons);
    if (qn_cons) free(qn_cons);

    return (0);
}
// end main
