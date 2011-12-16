#include <refinement.h>
#include <meshrefinement.h>
#include <geometry.h>
#include <solutionvector.h>
#include <algorithm>
#include <alignment.h>
#include <lc.h>
#include <box.h>
#include <qlc3d.h>


double getMaxS(SolutionVector& q){

int npLC = q.getnDoF();
double * dir = tensortovector(q.Values, npLC );

double max = *max_element( dir+3*npLC, dir+4*npLC );

if (dir) delete [] dir;

return max;

}


double intepolate_scalar(double* loc , double* S){
/*! Interpolates scalar value S[4] to a single value using four local coordinates in loc[4]*/
	double val = 0;
	for (size_t i = 0 ; i < 4 ; i++){
		val+= loc[i]*S[i];
	}
    return val;
}

void interpolate(SolutionVector& qnew,
				 Geometry& geom_new,
				 SolutionVector& qold,
				 Geometry& geom_old){
/*! Interpolates Q-tensor qold from olde geometry geom_old to new geometry geom_new*/

	// MAKE COORDINATE TO CONTAINING ELEMENT IDEX pint
    vector < unsigned int> pint; // index from point to containing tet
	geom_old.genIndToTetsByCoords( pint, geom_new.getPtrTop(), geom_new.getnpLC() );

    unsigned int npLC_old = (unsigned int) geom_old.getnpLC();
    double* dir = tensortovector( qold.Values , npLC_old ); // director on old mesh


    // Loop over all new nodes and set Q-tensor
    for (size_t ind = 0 ; ind < (size_t) geom_new.getnpLC() ; ind++){ // for all new nodes

#ifdef DEBUG
        if (geom_old.t->getMaterialNumber(pint[ind]) != MAT_DOMAIN1 ){
            printf("error in intepolate (autorefinement.cpp)  \n");
            printf("new node %i, at coordinate %e,%e,%e in old structure\n", (int) ind, geom_new.getpX(ind) , geom_new.getpY(ind) , geom_new.getpZ(ind) );
            printf("the material of found tet is %i, which is not DOMAIN1 \n", geom_old.t->getMaterialNumber(pint[ind]) );
            exit(1);
        }
#endif


        double loc[4]; // LOCAL ELEMENT COORDS
		double* coord = geom_new.getPtrTop() + (3*ind); // pointer to this nodes coordinates

        // calculate local element coordinates loc of global coordinate coord.
		geom_old.t->CalcLocCoords( pint[ind], geom_old.getPtrTop(), coord, loc);

		int n[4];
		n[0] = geom_old.t->getNode( pint[ind], 0);
		n[1] = geom_old.t->getNode( pint[ind], 1);
		n[2] = geom_old.t->getNode( pint[ind], 2);
		n[3] = geom_old.t->getNode( pint[ind], 3);
        int max_node = *max_element(n, n+3);
#ifdef DEBUG
        if (max_node >= (int)   npLC_old){
            printf("error - node is larger than npLC(= %i)\n", npLC_old);
            printf(" interpolate in autorefinement.cpp\n");
            geom_old.t->PrintElement(pint[ind] );
            exit(1);
        }
#endif
        // IF MAXIMUM LOCAL COORDINATE VALUE is more or less 1 -> the node is an exiting one
        size_t ind_max = max_element(loc, loc+4) - loc ;
        if ( loc[ ind_max ] >= 0.99999){// EXISTING CORNER NODE, STRAIGHT COPY OF OLD VALUES
            qnew.setValue(ind, 0 , qold.getValue ( n[ind_max] , 0 ) );
		    qnew.setValue(ind, 1 , qold.getValue ( n[ind_max] , 1 ) );
		    qnew.setValue(ind, 2 , qold.getValue ( n[ind_max] , 2 ) );
		    qnew.setValue(ind, 3 , qold.getValue ( n[ind_max] , 3 ) );
		    qnew.setValue(ind, 4 , qold.getValue ( n[ind_max] , 4 ) );
		}
        else{

            // INTERPOLATE USING Q-TENSOR COMPONENTS AS SCALARS
            for (int i = 0 ; i  < 5 ; i++){ // loop over each tensor component
                double qo[4] = { qold.getValue( n[0], i) , // q tensor at four corner nodes of old tet
                                 qold.getValue( n[1], i) ,
                                 qold.getValue( n[2], i) ,
                                 qold.getValue( n[3], i) };
                double qn = intepolate_scalar( loc , qo );
                qnew.setValue( ind , i , qn ); // set i'th dimension of value at node ind to qn
            }// end for i
		}
	}// end for all new coords

    if (dir) delete [] dir;
}



double get_elem_maxdQ(const unsigned int& elem, Geometry& geom, SolutionVector& q){
/*! returns absolute maximum change in any of the 5 Q-tensor components in element elem*/
    double qe[4] = {0,0,0,0};
    double maxdq = 0;
    for (unsigned int dim = 0 ; dim < 5 ; dim ++){ // for q1 -> q5
        for ( int j = 0 ; j < geom.t->getnNodes() ; j++){ // for each node in this tet
            int nn = geom.t->getNode(elem, j); // node number
            qe[j] = q.getValue( nn , dim ); // get g
        } // end for each node

        double mxq = *max_element(qe , qe+4);
        double mnq = *min_element(qe , qe+4);

        if ( (mxq-mnq) > maxdq ) {maxdq = mxq-mnq;}

    }// end for each dimension

    return maxdq;
}

void get_index_to_tred(Geometry& geom_orig,
                       SolutionVector& q,
                       vector <unsigned int>& i_tet,
                       MeshRefinement& meshref,
                       const int& refiter,
                       bool isEndRefinement = false
                       ){
/*! determines which elements need to be refined. i_tet is set to RED_TET for those elements
    that need refinement. If isEndRefinement is 'true', uses EndRefinement, other wise uses AutoRefinment */

    i_tet.clear();
    i_tet.assign( geom_orig.t->getnElements() , 0 );

    if (! isEndRefinement)
    {

        // Loop over all elements
        for (unsigned int i = 0 ; i < (unsigned int) geom_orig.t->getnElements() ; i++){
            if (geom_orig.t->getMaterialNumber(i) <= MAT_DOMAIN7){ // if LC element
                vector <AutoRef> :: iterator ar = meshref.AutoRefinement.begin();
                vector <AutoRef> :: iterator end_itr = meshref.AutoRefinement.end();

                for ( ; ar != end_itr ; ar++ ){ // for each AUTOREF obj

                    double elem_sze = geom_orig.t->getDeterminant( i ) * 1e18; // <- a guesstimate of element side length obtained from element volume
                    elem_sze = pow( elem_sze , 1.0/3.0);
                    double maxdq = get_elem_maxdQ(i, geom_orig, q );

                    if (( ar->Type == Change ) &&
                     ( ar->getMinSize()  < elem_sze) &&
                     ( ar->getMaxValue(refiter) < maxdq  ))
                        {
                        i_tet[i] = RED_TET;
                        }
                }// end for every autorefinement object
            }// end if LC element
        }// end for all tets
    }
    else // it is End Refinement
    {
       // printf(" get_index_to_tred:EndRef\n");


        for (unsigned int i = 0 ; i < (unsigned int) geom_orig.t->getnElements() ; i++){

            if (geom_orig.t->getMaterialNumber(i) <= MAT_DOMAIN7){ // if LC element
                vector <EndRef> :: iterator ar = meshref.EndRefinement.begin();
                vector <EndRef> :: iterator end_itr = meshref.EndRefinement.end();

                for ( ; ar != end_itr ; ar++ ){ // for each ENDREF obj

                    double elem_sze = geom_orig.t->getDeterminant( i ) * 1e18; // <- a guesstimate of element side length obtained from element volume
                    elem_sze = pow( elem_sze , 1.0/3.0);
                    double maxdq = get_elem_maxdQ(i, geom_orig, q );

                    if (( ar->Type == Change ) &&
                     ( ar->getMinSize()  < elem_sze) &&
                     ( ar->getMaxValue(refiter) < maxdq  ))
                        {
                        i_tet[i] = RED_TET;

                        }
                }// end for every autorefinement object
            }// end if LC element
        }// end for all tets
    }



}

bool needsRefinement(Geometry &geom, SolutionVector &q, MeshRefinement &meshrefinement){
/*! checks whether mesh refinement is needed. Returns true if needed, false otherwise*/
    vector <AutoRef> :: iterator ar;
    // Loop over all elements
    for (unsigned int i = 0 ; i < (unsigned int) geom.t->getnElements() ; i++){
        if (geom.t->getMaterialNumber(i) <= MAT_DOMAIN7){ // if LC element
            for (ar = meshrefinement.AutoRefinement.begin(); ar != meshrefinement.AutoRefinement.end() ; ar++ ){ // for each AUTOREF obj

                for(int refiter = 0; refiter < (int) ar->getNumIterations(); refiter++){ // for each refiter in this autoref object

                     double elem_sze = geom.t->getDeterminant( i ) * 1e18; // <- a guesstimate of element side length obtained from element volume
                     elem_sze = pow( elem_sze , 1.0/3.0);
                     double maxdq = get_elem_maxdQ(i, geom, q );

                     // Check if all conditions for refinement are satisfied. return true and exit if yes
                     if (( ar->Type == Change ) &&
                         ( ar->getMinSize()  < elem_sze) &&
                         ( ar->getMaxValue(refiter) < maxdq  )){
                         return true;
                     }// end if
                }// end for refiter
            }// end for every autorefinement object
        }// end if LC element
    }// end for all tets

    return false; // if no tets that need refinement are found, return false and exit
}


bool needsEndRefinement(Geometry &geom, SolutionVector &q, MeshRefinement &meshrefinement){
/*! checks whether mesh refinement is needed. Returns true if needed, false otherwise*/
    vector <EndRef> :: iterator er;
    // Loop over all elements

    for (unsigned int i = 0 ; i < (unsigned int) geom.t->getnElements() ; i++){

        if (geom.t->getMaterialNumber(i) <= MAT_DOMAIN7){ // if LC element
            for (er = meshrefinement.EndRefinement.begin(); er != meshrefinement.EndRefinement.end() ; er++ ){ // for each AUTOREF obj

                //er->printAutoref();

                for(int refiter = 0; refiter < (int) er->getNumIterations(); refiter++){ // for each refiter in this autoref object

                     double elem_sze = geom.t->getDeterminant( i ) * 1e18; // <- a guesstimate of element side length obtained from element volume
                     elem_sze = pow( elem_sze , 1.0/3.0);
                     double maxdq = get_elem_maxdQ(i, geom, q );

                     // Check if all conditions for refinement are satisfied. return true and exit if yes
                     if (( er->Type == Change ) &&
                         ( er->getMinSize()  < elem_sze) &&
                         ( er->getMaxValue(refiter) < maxdq  )  //&&
                         //( meshrefinement.EndRefinement[0].getEndRefIteration() ) // check for max number of End-Refinements
                         ){
                         return true;
                     }// end if
                }// end for refiter
            }// end for every autorefinement object
        }// end if LC element
    }// end for all tets

    return false; // if no tets that need refinement are found, return false and exit
}






void autoref(   Geometry& geom_orig, Geometry& geom_prev, Geometry& geom_new,
                SolutionVector& q, SolutionVector& qn,
                SolutionVector& v,
                MeshRefinement& meshrefinement,
                Simu& simu, Alignment& alignment,Electrodes& electrodes, LC& lc){
    bool bRefined   = false; // indicates whether mesh is changed or not
    int refiter     = 0;     // refinement iteration counter
    int maxrefiter  = meshrefinement.getMaxNumAutoRefIterations(); // maximum possible number of refinement iterations

    bool IsEndRefinement = false; // switch between Auto-Refinement and End-Refinement
    if ( !simu.IsRunning() )
    {
         printf("\nEnd-Refinement\n");
         IsEndRefinement = true;
         maxrefiter = meshrefinement.getMaxNumEndRefIterations();

         refiter    = meshrefinement.EndRefinement[0].getEndRefIteration();
         int maxEndRef = meshrefinement.EndRefinement[0].getRefIter(); // maximum times end-refinement can be done
         if (refiter >= maxEndRef) // if maximum number of end refinements reached
         {
             maxrefiter = 0;
         }
         //maxrefiter = maxrefiter > meshrefinement.EndRefinement[0].getEndRefIteration() ? 0:maxrefiter;

         printf("\nEnd-Refinement, refiter = %i , maxrefiter = %i\n", refiter, maxrefiter );
    }



    Geometry geom_temp;     // temporary "working" geometry
    SolutionVector q_prev;  // INPUT Q-TENSOR

    if (maxrefiter == 0){   // IF NO REFIENMENT
        return;             // LEAVE REFINEMENT FUNCTION NOW
        simu.setMeshModified(false);
    }
    else{// ELSE YES REFINEMENT
        // TEMPORARY SOLUTION VECTOR USED TO HOLD PREVIOUS Q-TENSOR
        q_prev = q;

        // PREVIOUS Q IS INTERPOLATED TO ORIGINAL MESH - downsampling, this loses information
        // Q IS RESIZED TO MATCH THE ORIGIANL MESH
        q.Resize( geom_orig.getnpLC(), 5 ); // presumably this resets q
        interpolate(q, geom_orig, q_prev, geom_prev ); // overwrites with interplated values

        // REFINEMENT START FROM ORIGINAL GEOMETRY
        geom_new.setTo(  &geom_orig );
        geom_temp.setTo( &geom_orig );
    }

    //=====================================
    // DO REFINEMENT ITERATIONS, SPLIT TETS
    //=====================================
    printf("=========================================== \n");
    printf("Doing a maximum of %i refinement iterations \n", maxrefiter);
    printf("=========================================== \n");
    for ( refiter = 0 ; refiter < maxrefiter ; refiter ++){ // for max refiter
        printf("starting refiter = %i\n", refiter);
        // GET INDEX TO RED TETS IN geom_new
        vector <unsigned int> i_tet; // <- refinement type counter
        get_index_to_tred( geom_new, q, i_tet, meshrefinement, refiter, IsEndRefinement);



        // LEAVE REF LOOP IF NO REFINABLE TETS FOUND
        if (*max_element(i_tet.begin() , i_tet.end() ) < RED_TET ){
            printf("|------>  Refiter %i of %i <------| \n", refiter+1 , maxrefiter );
            printf("   No refinement this iteration\n");
            break;
        }

        printf("|------>  Refiter %i of %i <------| \n", refiter+1 , maxrefiter );
        Refine( geom_temp, geom_new , i_tet);
        printf("New node count: %i\n", geom_new.getnp());

        // UPDATE Q-TENSOR VALUES BY INTERPOLATING TO NEW MESH
        //SolutionVector qtemp( geom_new.getnpLC() , 5);	// TEMPORARY Q-TENSOR
        q.Resize( geom_new.getnpLC() , 5 ); // resets q-tensor
        interpolate(q, geom_new, q_prev, geom_prev); // interpolates from previous result

        // UPDATE geom_temp, TO USE IT IN NEXT REFINEMENT ITERATION
        geom_temp.setTo( &geom_new);

        bRefined = true; // YES, MESH HAS BEEN CHANGED
        printf("done refiter = %i\n", refiter);
    }// end for max refiters
    // REFIENEMENT IS DONE. DO STUFF (CALCULATE DETERMINANTS, TAKE CARE OF FIXED NODES ETC.) TO MAKE NEW GEOMETRY VALID.

    geom_prev.setTo( &geom_new ); // SET PREV TO NEW FOR FUTURE REFIENEMENTS

    // REALLOCATE POTENTIAL
    v.Allocate( geom_new.getnp() , 1);
    v.setFixedNodesPot(electrodes, geom_new.e, 0.0 ); // <-- sets current time to zero!!! this may be a problem
    v.setPeriodicEquNodes( &geom_new);
    //cout << "v reallocated " << endl;

    // SET BOUNDARY CONDITIONS
    //setSurfacesQ(&q, &alignment, &lc, &geom_new);
    setStrongSurfacesQ(&q, &alignment, &lc, &geom_new);
    q.setFixedNodesQ(&alignment, geom_new.e );

    q.setPeriodicEquNodes( & geom_new );


    qn = q; // USE CURRENT Q FOR PREVIOUS TIME STEP Q-TENSOR


    // NEW MESH FILE NEEDS TO BE WRITTEN WHEN RESULTS ARE OUTPUT
    // LET REST OF PROGRAM KNOW THAT GEOMETRY HAS BEEN MODIFIED
    if (bRefined)
    {
        cout << "=============done refining mesh=============" << endl;
        cout << "new nodecount = " << geom_new.getnp() << endl;
        //if (simu.getdt() > 0)
        //    simu.setdt( simu.getMindt() );
        simu.IncrementMeshNumber();     // output mesh name will be appended with this number
        simu.setMeshModified( true );   // this is a flag set to notify that a new output file needs to be written
		
		// if this was an end-refinement, need to make changes to simu, so that additional
		// simulation steps are taken
        if ( !simu.IsRunning() ) 
        {
            meshrefinement.EndRefinement[0].incrementEndRefIteration(); // increments end refinement counter
            simu.resetEndCriterion(); // resets end citerion = forces additional simulation steps on new mesh
        }
    }

}




