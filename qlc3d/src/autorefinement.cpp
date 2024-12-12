#include <refinement.h>
#include <meshrefinement.h>
#include <refinement/refinement-spec.h>
#include <refinement/q-tensor-interpolator.h>
#include <geometry.h>
#include <solutionvector.h>
#include <algorithm>
#include <alignment.h>
#include <simulation-state.h>
#include <util/logging.h>
#include <inits.h>



bool needsInterpolatedQ(const std::vector<const RefinementSpec*> &specs,
                        const unsigned int refiter) {
// CHECKS WHETHER IT IS NECESSARY TO INTERPOLATE Q-TENSOR ONTO A DIFFERENT MESH.
// SOME REFINEMENT TYPES NEED THIS WHILE OTHERS DON'T. INTERPOLATION CAN BE
// SLOW, SO AVOIDING IT CAN SPEED UP THINGS
    for (const RefinementSpec* spec : specs) {
        if (refiter > spec->getRefIter())
            continue;
        if (spec->getType() == RefinementSpec::Type::Change)
            return true;
    }
    return false;
}

void get_index_to_tred(Geometry &geom_curr, // CURRENT CALCULATION GEOMETRY
                       Geometry &geom_work, // REFINED WORKING GEOMETRY
                       SolutionVector &q,   // Q_TENSOR CORRESPONDING TO geom_curr
                       vector <idx> &i_tet,
                       const std::vector<const RefinementSpec*> &specs,
                       const unsigned int refiter,
                       bool isEndRefinement = false
                      ) {
    // DETERMINES WHICH ELEMENTS NEED TO BE REFINED. VALUES IN i_tet
    // ARE SET TO RED_TET FOR THOSE TETRAHEDRA THAT NEED TO BE SPLIT.
    i_tet.clear();
    i_tet.assign(geom_work.getTetrahedra().getnElements() , 0);
    SolutionVector q_temp;
    // IF INTERPOLATION NEEDED. DO IT
    if (needsInterpolatedQ(specs, refiter)) {
        q_temp = interpolateQTensor(geom_work, geom_curr, q);
    }
    // LOOP OVER EACH SPEC AND PROCESS IT
    for (const RefinementSpec* spec : specs) {
        if (refiter > spec->getRefIter())
            continue;
        switch (spec->getType()) {
            case (RefinementSpec::Type::Change): {
                Log::info("Searching for tetrahedra by refinement type = CHANGE");
                findTets_Change(*spec, i_tet, refiter, geom_work, q_temp);
                break;
            }
            case (RefinementSpec::Type::Sphere): {
                Log::info("Searching for tetrahedra by refinement type = SPHERE");
                findTets_Sphere(*spec, i_tet, refiter, geom_work);
                break;
            }
            case (RefinementSpec::Type::Box): {
                Log::info("Searching for tetrahedra by refinement type = BOX");
                findTets_Box(*spec, i_tet, refiter, geom_work);
                break;
            }
            default: {
                throw std::runtime_error(fmt::format("error in {}, {}, unhandled refinement type", __FILE__, __func__));
            }
        }
    }
}

idx getMaxRefiterCount(const std::vector<const RefinementSpec*> &specs) {
// DETERMINES MAXIMUM NUMBER OF REFINEMENT ITERATIONS THAT MAY BE PERFORMED
    idx maxRefIter(0);
    for (const RefinementSpec* spec : specs) {
        maxRefIter = maxRefIter > spec->getRefIter() ? maxRefIter : spec->getRefIter();
    }
    return maxRefIter;
}

bool autoref(Geometry &geom_orig,
             Geometry &geom,
             SolutionVector &q,
             SolutionVector &v,
             const std::vector<const RefinementSpec*> &specs,
             Simu &simu,
             SimulationState &simulationState,
             Alignment &alignment,
             const Electrodes &electrodes,
             double S0) {

  if (geom_orig.getTetrahedra().getElementType() != ElementType::LINEAR_TETRAHEDRON) {
    throw NotYetImplementedException("Only linear tetrahedra are supported in mesh refinement. Found element type " +
      toString(geom_orig.getTetrahedra().getElementType()));
  }

    bool bRefined{false};   // indicates whether mesh is changed or not
    unsigned int refiter{0};         // refinement iteration counter
    unsigned int maxrefiter = getMaxRefiterCount(specs);

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
                          specs,
                          refiter,
                          false);
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
    v.initialisePotentialBoundaries(electrodes.getCurrentPotentials(simulationState.currentTime().getTime()), geom_temp);

    // REALLOCATE Q-TENSOR
    SolutionVector qTemp(q.getnDoF(), 5);
    qTemp = q; // temp swap
    q = interpolateQTensor(geom_temp, geom, qTemp);    // INTERPOLATE FROM PREVIOUS MESH
    // SET BOUNDARY CONDITIONS
    setStrongSurfacesQ(q, alignment, S0, geom_temp);
    q.initialiseLcBoundaries(geom_temp, alignment);
    geom_temp.clearPeriodicNodesMapping(); // release resources that are not needed anymore

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




