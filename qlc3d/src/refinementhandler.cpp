#include <eventlist.h>
#include <eventhandler.h>
#include <refinement.h> // declares autorefinement etc.
#include <list>
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

    // MAKE LIST OF ALL REFINFO OBJECTS THAT ARE HANDLES NOW
    std::list<RefInfo> refInfos;
    std::list<Event*>::iterator evitr = refEvents.begin();
    for ( ; evitr!= refEvents.end() ; ++evitr){
        RefInfo* ref = static_cast<RefInfo*>( (*evitr)->getEventDataPtr() );
        refInfos.push_back( *ref );
    }
    Log::info("{} Refinement objects.", refInfos.size());

    // TRY TO DO REFINEMENT
    bool isRefined(false);
    isRefined = autoref(*geometries.geom_orig,
            *geometries.geom,
            *solutionvectors.q,
            *solutionvectors.v,
            refInfos,
            simu,
            simulationState,
            alignment,
            electrodes,
            S0);

    // DELETE ALL REFINEMENT EVENTS. ALWAYS
    for (evitr = refEvents.begin() ; evitr != refEvents.end() ; ++evitr){
        delete (*evitr);
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
