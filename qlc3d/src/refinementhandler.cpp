#include <eventlist.h>
#include <eventhandler.h>
#include <refinement.h> // declares autorefinement etc.
#include <refinement/refinement-spec.h>
#include <list>
#include <vector>
#include <qlc3d.h>
#include <util/logging.h>

#include <spamtrix_ircmatrix.hpp>

bool handleMeshRefinement(std::list<Event*>& refEvents,
                          Geometries& geometries,
                          SolutionVectors& solutionvectors,
                          Simu& simu,
                          SimulationState &simulationState,
                          Alignment& alignment,
                          Electrodes& electrodes,
                          double S0) {
    Log::info("Doing {} mesh refinements.", refEvents.size());

    // COLLECT ALL REFINEMENTSPEC POINTERS
    std::vector<const RefinementSpec*> specs;
    for (Event* ev : refEvents) {
        const RefinementSpec* spec = ev->getRefinementSpec();
        if (spec != nullptr) {
            specs.push_back(spec);
        }
    }
    Log::info("{} Refinement objects.", specs.size());

    // TRY TO DO REFINEMENT
    bool isRefined = autoref(*geometries.geom_orig,
            *geometries.geom,
            *solutionvectors.q,
            *solutionvectors.v,
            specs,
            simu,
            simulationState,
            alignment,
            electrodes,
            S0);

    // DELETE ALL REFINEMENT EVENTS. ALWAYS
    for (Event* ev : refEvents) {
        delete ev;
    }
    return isRefined;
}


void handlePreRefinement(std::list<Event*>& refEvents,
                         Geometries& geometries,
                         SolutionVectors& solutionvectors,
                         Simu& simu,
                         SimulationState &simulationState,
                         Alignment& alignment,
                         Electrodes& electrodes,
                         double S0)
{
// PRE REFINMENT MODIFICATIONS TO THE MESH CARRY THROUGH THE
// REST OF THE SIMULATION. THAT IS, THE INITIAL GEOMETRY IS
// MODIFIED TOO

    Log::info("Doing pre-refinement.");
    // MAKE LIST OF ALL REFINFO OBJECTS
    handleMeshRefinement(refEvents,
                         geometries,
                         solutionvectors,
                         simu,
                         simulationState,
                         alignment,
                         electrodes,
                         S0);
    // "ORIGINAL" MESH IS MODIFIED
    geometries.geom_orig->setTo( geometries.geom );

}
