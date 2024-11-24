#include <refinement.h>
#include <meshrefinement.h>
#include <geometry.h>
#include <solutionvector.h>
#include <algorithm>
#include <alignment.h>
#include <simulation-state.h>
#include <util/logging.h>
#include <geom/vec3.h>
#include <geom/coordinates.h>
#include <inits.h>

double intepolate_scalar(double *loc , double *S) {
    /*! Interpolates scalar value S[4] to a single value using four local coordinates in loc[4]*/
    double val = 0;
    for (size_t i = 0 ; i < 4 ; i++) {
        val += loc[i] * S[i];
    }
    return val;
}

void interpolate(SolutionVector &qnew,
                 Geometry &geom_new,
                 SolutionVector &qold,
                 Geometry &geom_old) {
    // Interpolates Q-tensor qold from ye olde geometry geom_old to new geometry geom_new
    // MAKE COORDINATE TO CONTAINING ELEMENT IDEX pint
    vector <idx> pInTet; // index from point to containing tet
    geom_old.genIndToTetsByCoords(pInTet,
                                  geom_new.getCoordinates(),
                                  true,    // ERROR IF COORDINATE IS NOT FOUND
                                  true     // ONLY CONSIDER LC ELEMENTS
                                 );
    //unsigned int npLC_old = (unsigned int) geom_old.getnpLC();
    //double *dir = tensortovector(qold.Values , npLC_old);   // director on old mesh
    // Loop over all new nodes and set Q-tensor
    qnew.Resize(geom_new.getnpLC(), 5);
    auto &oldTets = geom_old.getTetrahedra();
    for (size_t ind = 0 ; ind < (size_t) geom_new.getnpLC() ; ind++) { // for all new nodes
        // only LC volumes are interpolated
        assert(geom_old.getTetrahedra().getMaterialNumber(pInTet[ind]) == MAT_DOMAIN1);

        double loc[4]; // LOCAL ELEMENT COORDS

        Vec3 targetPoint = geom_new.getCoordinates().getPoint(ind);
        // calculate local element coordinates loc of global coordinate coord.
        oldTets.calcLocCoords(pInTet[ind], geom_old.getCoordinates(), targetPoint, loc);
        size_t n[4];
        n[0] = oldTets.getNode(pInTet[ind], 0);
        n[1] = oldTets.getNode(pInTet[ind], 1);
        n[2] = oldTets.getNode(pInTet[ind], 2);
        n[3] = oldTets.getNode(pInTet[ind], 3);

        // IF MAXIMUM LOCAL COORDINATE VALUE is more or less 1 -> the node is an exsiting one
        size_t ind_max = max_element(loc, loc + 4) - loc ;
        if (loc[ ind_max ] >= 0.99999) {// EXISTING CORNER NODE, STRAIGHT COPY OF OLD VALUES
            qnew.setValue(ind, 0 , qold.getValue(n[ind_max] , 0));
            qnew.setValue(ind, 1 , qold.getValue(n[ind_max] , 1));
            qnew.setValue(ind, 2 , qold.getValue(n[ind_max] , 2));
            qnew.setValue(ind, 3 , qold.getValue(n[ind_max] , 3));
            qnew.setValue(ind, 4 , qold.getValue(n[ind_max] , 4));
        } else {
            // INTERPOLATE USING Q-TENSOR COMPONENTS AS SCALARS
            for (int i = 0 ; i  < 5 ; i++) { // loop over each tensor component
                double qo[4] = { qold.getValue(n[0], i) ,  // q tensor at four corner nodes of old tet
                                 qold.getValue(n[1], i) ,
                                 qold.getValue(n[2], i) ,
                                 qold.getValue(n[3], i)
                               };
                double qn = intepolate_scalar(loc , qo);
                qnew.setValue(ind , i , qn);   // set i'th dimension of value at node ind to qn
            }// end for i
        }
    }// end for all new coords
    //if (dir) delete [] dir;
}

bool needsInterpolatedQ(const list<RefInfo> &refInfos,
                        const unsigned int refiter) {
// CHECKS WHETHER IT IS NECESSARY TO INTERPOLATE Q-TENSOR ONTO A DIFFERENT MESH.
// SOME REFINEMENT TYPES NEED THIS WHILE OTHERS DON'T. INTERPOLATION CAN BE
// SLOW, SO AVOIDING IT CAN SPEED UP THINGS
    list<RefInfo>::const_iterator ritr = refInfos.begin();
    for (; ritr != refInfos.end() ; ++ritr) {
        // IF THIS REFITER IS NOT DEFINED FOR THIS REFINFO OBJECT, SKIP IT
        if (refiter > (*ritr).getRefIter())
            continue;
        switch ((*ritr).getType()) {
        case (RefInfo::Change):
            return true;
        default: {} // DO NOTHING
        }
    }
    return false;
}

void get_index_to_tred(Geometry &geom_curr, // CURRENT CALCULATION GEOMETRY
                       Geometry &geom_work, // REFINED WORKING GEOMETRY
                       SolutionVector &q,   // Q_TENSOR CORRESPONDING TO geom_curr
                       vector <idx> &i_tet,
                       const list<RefInfo> &refInfos,
                       const unsigned int refiter,
                       bool isEndRefinement = false
                      ) {
    // DETERMINES WHICH ELEMENTS NEED TO BE REFINED. VALUES IN i_tet
    // ARE SET TO RED_TET FOR THOSE TETRAHEDRA THAT NEED TO BE SPLIT.
    i_tet.clear();
    i_tet.assign(geom_work.getTetrahedra().getnElements() , 0);
    SolutionVector q_temp;
    // IF INTERPOLATION NEEDED. DO IT
    if (needsInterpolatedQ(refInfos, refiter)) {
        interpolate(q_temp, geom_work, q, geom_curr);
    }
    // LOOP OVER EACH REFINFO OBJECT AND PROCESS IT
    list<RefInfo>::const_iterator ritr = refInfos.begin();
    for (ritr = refInfos.begin() ; ritr != refInfos.end() ; ++ritr) {
        // CHECK WHETHER THIS RefInfo OBJECT IS DEFINED FOR THIS REFINEMENT ITERATION
        if (refiter > (*ritr).getRefIter())
            continue;
        // DIFERENT REFINFO TYPES REQUIRE DIFFERENT ELEMENT SEARCHES
        // ELEMENT SEARCH FUCNTIONS DEFINED IN findrefelems.cpp
        switch ((*ritr).getType()) {
            case (RefInfo::Change): {
                Log::info("Searching for tetrahedra by refinement type = CHANGE");
                findTets_Change((*ritr) , i_tet, refiter, geom_work , q_temp);
                break;
            }
            case (RefInfo::Sphere): {
                Log::info("Searching for tetrahedra by refinement type = SPHERE");
                findTets_Sphere((*ritr) , i_tet, refiter, geom_work);
                break;
            }
            case (RefInfo::Box): {
                Log::info("Searching for tetrahedra by refinement type = BOX");
                findTets_Box((*ritr), i_tet, refiter, geom_work);
                break;
            }
            default: {
                throw std::runtime_error(fmt::format("error in {}, {}, unhandled refinement type", __FILE__, __func__));
            }
        }
    }
}

idx getMaxRefiterCount(const list<RefInfo> &refInfos) {
// DETERMINES MAXIMUM NUMBER OF REFINEMENT ITERATIONS THAT MAY BE PERFORMED
    list<RefInfo>::const_iterator ri_itr = refInfos.begin();
    idx maxRefIter(0);
    for (; ri_itr != refInfos.end() ; ++ri_itr) {
        maxRefIter = maxRefIter > (*ri_itr).getRefIter() ? maxRefIter : (*ri_itr).getRefIter();
    }
    return maxRefIter;
}

bool autoref(Geometry &geom_orig,
             Geometry &geom,
             SolutionVector &q,
             SolutionVector &qn,
             SolutionVector &v,
             const list<RefInfo> &refInfos,
             Simu &simu,
             SimulationState &simulationState,
             Alignment &alignment,
             const Electrodes &electrodes,
             double S0) {
    bool bRefined{false};   // indicates whether mesh is changed or not
    unsigned int refiter{0};         // refinement iteration counter
    unsigned int maxrefiter = getMaxRefiterCount(refInfos);

    if (maxrefiter == 0) {   // IF NO REFINEMENT
        simulationState.meshModified(false);
        return false;             // LEAVE REFINEMENT FUNCTION NOW
    }
    // CREATE TEMPORARY WORKING COPY OF GEOMETRY
    // THIS WILL BE MODIFIED (REFINED/"DEREFINED")
    Geometry geom_temp;     // temporary "working" geometry
    geom_temp.setTo(&geom_orig);
    //=====================================
    // DO REFINEMENT ITERATIONS, SPLIT TETS
    //=====================================
    Log::info("Doing a maximum of {} refinement iterations.", maxrefiter);
    for (refiter = 0 ; refiter < maxrefiter ; refiter ++) { // for max refiter
        Log::info("Refinement iteration {} of {}.", refiter + 1, maxrefiter);
        // GET INDEX TO RED TETS IN geom
        vector <idx> i_tet(geom_temp.getTetrahedra().getnElements(), 0);  // REFINEMENT TYPE INDICATOR
        // SELECT RED TETS
        get_index_to_tred(geom ,
                          geom_temp,
                          q,
                          i_tet,
                          refInfos,
                          refiter,
                          false);//IsEndRefinement);
        // LEAVE REF LOOP IF NO REFINABLE TETS FOUND
        if (*max_element(i_tet.begin() , i_tet.end()) < RED_TET) {
            Log::info("No tetrahedra requiring refinement found during iteration {}.", refiter + 1);
            continue;
        }
        Refine(geom_temp  , i_tet);
        bRefined = true;                        // YES, MESH HAS BEEN CHANGED
        Log::info("Completed refinement iteration {} of {}. New node count is {}",
                  refiter + 1, maxrefiter, geom_temp.getnp());
    }
    //=============================================================
    //  DONE WITH REFINEMENT.
    //  DO CLEANUP AND INTERPOLATE RESULT ON NEW MESH.
    //=============================================================
    //geom_temp.t->CalculateDeterminants3D( geom_temp.getPtrTop() );
    //geom_temp.t->ScaleDeterminants( 1e-18);// scale to microns cubed
    geom_temp.calculateNodeNormals();
    //geom_temp.genIndWeakSurfaces(alignment);
    geom_temp.makeRegularGrid(simu.getRegularGridXCount(),
                              simu.getRegularGridYCount(),
                              simu.getRegularGridZCount());
    // RECREATE POTENTIAL SOLUTIONVECTOR FROM SCRATCH FOR THE NEW GEOMETRY.
    v.Resize(geom_temp.getnp() , 1);
  v.initialisePotentialBoundaries(geom_temp, electrodes.getCurrentPotentials(simulationState.currentTime().getTime()));

    // REALLOCATE Q-TENSOR
    qn = q; // temp swap
    q.Resize(geom_temp.getnpLC(), 5);       // ALLOCATE FOR NEW MESH SIZE
    interpolate(q, geom_temp, qn, geom);    // INTERPOLATE FROM PREVIOUS MESH
    // SET BOUNDARY CONDITIONS
    setStrongSurfacesQ(q, alignment, S0, geom_temp);
  q.initialiseLcBoundaries(geom_temp, alignment);

    q.EnforceEquNodes(geom_temp);
    qn = q;                                     // USE CURRENT Q FOR PREVIOUS TIME STEP Q-TENSOR
    geom.setTo(&geom_temp);
    // NEW MESH FILE NEEDS TO BE WRITTEN WHEN RESULTS ARE OUTPUT
    // LET REST OF PROGRAM KNOW THAT GEOMETRY HAS BEEN MODIFIED
    Log::info("Completed mesh refinement. New node count = {}, tetrahedron count = {}, triangle count = {}",
              geom.getnp(), geom.getTetrahedra().getnElements(), geom.getTriangles().getnElements());

    if (simu.simulationMode() == TimeStepping) {
        simulationState.dt(simu.getMindt());
        simulationState.restrictedTimeStep(true);
    }

    simulationState.incrementMeshNumber();
    simulationState.meshModified(true);

    return bRefined; // WHETHER MESH WAS REFINED
}




