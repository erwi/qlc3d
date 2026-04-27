#include <refinement.h>
#include <meshrefinement.h>
#include <refinement/refinement-spec.h>
#include <refinement/q-tensor-interpolator.h>
#include <mesh/element-split-convert.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <solutionvector.h>
#include <algorithm>
#include <alignment.h>
#include <simulation-state.h>
#include <util/logging.h>
#include <inits.h>
#include <util/exception.h>
#include <regulargrid-factory.h>


namespace {

RawMeshData geometryToRawMeshData(const Geometry &geom) {
    const Mesh &tets = geom.getTetrahedra();
    const Mesh &tris = geom.getTriangles();

    const auto tetOrder = getElementOrder(tets.getElementType());
    const auto triOrder = getElementOrder(tris.getElementType());
    if (tetOrder == 0 || triOrder == 0 || tetOrder != triOrder) {
        RUNTIME_ERROR(fmt::format("Cannot convert geometry to raw mesh data: tetrahedra are {}, triangles are {}.",
                                  tets.getElementType(), tris.getElementType()));
    }

    std::vector<Vec3> points;
    points.reserve(geom.getnp());
    for (idx i = 0; i < geom.getnp(); ++i) {
        points.push_back(geom.getCoordinates().getPoint(i));
    }

    std::vector<idx> tetNodes;
    std::vector<idx> tetMaterials;
    tetNodes.reserve(tets.getnElements() * tets.getnNodes());
    tetMaterials.reserve(tets.getnElements());
    for (idx i = 0; i < tets.getnElements(); ++i) {
        for (idx n = 0; n < tets.getnNodes(); ++n) {
            tetNodes.push_back(tets.getNode(i, n));
        }
        tetMaterials.push_back(tets.getMaterialNumber(i));
    }

    std::vector<idx> triNodes;
    std::vector<idx> triMaterials;
    triNodes.reserve(tris.getnElements() * tris.getnNodes());
    triMaterials.reserve(tris.getnElements());
    for (idx i = 0; i < tris.getnElements(); ++i) {
        for (idx n = 0; n < tris.getnNodes(); ++n) {
            triNodes.push_back(tris.getNode(i, n));
        }
        triMaterials.push_back(tris.getMaterialNumber(i));
    }

    return {tetOrder, std::move(points), std::move(tetNodes), std::move(tetMaterials),
            std::move(triNodes), std::move(triMaterials)};
}

} // namespace


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
             double S0,
             std::unique_ptr<RegularGrid>& regGridOut) {

  const ElementType originalTetType = geom_orig.getTetrahedra().getElementType();
  const ElementType originalTriType = geom_orig.getTriangles().getElementType();
  const bool quadraticInput = originalTetType == ElementType::QUADRATIC_TETRAHEDRON;

  if (quadraticInput) {
    if (originalTriType != ElementType::QUADRATIC_TRIANGLE) {
      RUNTIME_ERROR("Quadratic mesh refinement expects quadratic triangles alongside quadratic tetrahedra. Found triangle element type " +
                    toString(originalTriType));
    }
  } else if (originalTetType != ElementType::LINEAR_TETRAHEDRON) {
    RUNTIME_ERROR("Only linear or quadratic tetrahedra are supported in mesh refinement. Found element type " +
                  toString(originalTetType));
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
    if (quadraticInput) {
        RawMeshData rawLinear = geometryToRawMeshData(geom_orig);
        splitQuadraticGeometryToLinear(geom_orig, rawLinear);
        auto coordinates = std::make_shared<Coordinates>(std::move(rawLinear.points));
        geom_temp.setCoordinates(coordinates);
        geom_temp.getTetrahedra().setElementData(ElementType::LINEAR_TETRAHEDRON,
                                                 std::move(rawLinear.tetNodes),
                                                 std::move(rawLinear.tetMaterials));
        geom_temp.getTriangles().setElementData(ElementType::LINEAR_TRIANGLE,
                                                 std::move(rawLinear.triNodes),
                                                 std::move(rawLinear.triMaterials));
        geom_temp.getTetrahedra().calculateDeterminants3D(geom_temp.getCoordinates());
        geom_temp.getTetrahedra().ScaleDeterminants(qlc3d::units::CUBIC_MICROMETER_TO_CUBIC_METER);
        geom_temp.getTriangles().setConnectedVolume(&geom_temp.getTetrahedra());
        geom_temp.getTriangles().calculateSurfaceNormals(geom_temp.getCoordinates(), &geom_temp.getTetrahedra());
        geom_temp.getTriangles().ScaleDeterminants(qlc3d::units::SQUARE_MICROMETER_TO_SQUARE_METER);
    } else {
        geom_temp.setTo(&geom_orig);
    }
    //=====================================
    // DO REFINEMENT ITERATIONS, SPLIT TETS
    //=====================================
    if (!quadraticInput) {
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
    } else {
        bRefined = true;
        Log::info("Quadratic input: using split–promote path without invoking the linear refinement loop.");
    }

    Geometry *geom_work = &geom_temp;
    if (quadraticInput) {
        RawMeshData rawRefined = geometryToRawMeshData(geom_temp);
        convertLinearMeshDataToQuadratic(rawRefined);
        auto coordinates = std::make_shared<Coordinates>(std::move(rawRefined.points));
        geom_temp.setCoordinates(coordinates);
        geom_temp.getTetrahedra().setElementData(ElementType::QUADRATIC_TETRAHEDRON,
                                                std::move(rawRefined.tetNodes),
                                                std::move(rawRefined.tetMaterials));
        geom_temp.getTriangles().setElementData(ElementType::QUADRATIC_TRIANGLE,
                                                std::move(rawRefined.triNodes),
                                                std::move(rawRefined.triMaterials));
    }
    //=============================================================
    //  DONE WITH REFINEMENT.
    //  DO CLEANUP AND INTERPOLATE RESULT ON NEW MESH.
    //=============================================================
    //geom_temp.t->CalculateDeterminants3D( geom_temp.getPtrTop() );
    //geom_temp.t->ScaleDeterminants(qlc3d::units::CUBIC_MICROMETER_TO_CUBIC_METER);
    geom_work->calculateNodeNormals();
    // Build regular grid for the refined geometry
    regGridOut = buildRegularGrid(simu.getRegularGridXCount(),
                                  simu.getRegularGridYCount(),
                                  simu.getRegularGridZCount(),
                                  *geom_work);

    // RECREATE POTENTIAL SOLUTIONVECTOR FROM SCRATCH FOR THE NEW GEOMETRY.
    v.Resize(geom_work->getnp() , 1);
    v.initialisePotentialBoundaries(electrodes.getCurrentPotentials(simulationState.currentTime().getTime()), *geom_work);

    // REALLOCATE Q-TENSOR
    SolutionVector qTemp(q.getnDoF(), 5);
    qTemp = q; // temp swap
    q = interpolateQTensor(*geom_work, geom, qTemp);    // INTERPOLATE FROM PREVIOUS MESH
    // SET BOUNDARY CONDITIONS
    setStrongSurfacesQ(q, alignment, S0, *geom_work);
    q.initialiseLcBoundaries(*geom_work, alignment);
    geom_work->clearPeriodicNodesMapping(); // release resources that are not needed anymore

    geom.setTo(geom_work);
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




