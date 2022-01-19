#include <eventlist.h>
#include <eventhandler.h>
#include <refinement.h> // declares autorefinement etc.
#include <list>
#include <qlc3d.h>
#include <util/logging.h>

#include <spamtrix_ircmatrix.hpp>

void handleMeshRefinement(std::list<Event*>& refEvents,
                          Geometries& geometries,
                          SolutionVectors& solutionvectors,
                          Simu& simu,
                          SimulationState &simulationState,
                          Alignment& alignment,
                          Electrodes& electrodes,
                          double S0,
                          SpaMtrix::IRCMatrix &Kpot,
                          SpaMtrix::IRCMatrix &Kq)
{
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
            *solutionvectors.qn,
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
    // IF MESH HAS BEEN REFINED NEED TO RECREATE MATRIXES
    // FOR Q-TENSOR AND POTENTIAL
    if (isRefined){
        Log::info("Creating new matrices for potential and Q-tensor solutions.");

        Kpot = createPotentialMatrix(*geometries.geom,
                                     *solutionvectors.v,
                                     0,
                                     electrodes);

        Kq = createQMatrix(*geometries.geom, *solutionvectors.q);
    }
}


void handlePreRefinement(std::list<Event*>& refEvents,
                         Geometries& geometries,
                         SolutionVectors& solutionvectors,
                         Simu& simu,
                         SimulationState &simulationState,
                         Alignment& alignment,
                         Electrodes& electrodes,
                         double S0,
                         SpaMtrix::IRCMatrix &Kpot,
                         SpaMtrix::IRCMatrix &Kq)
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
                         S0,
                         Kpot,
                         Kq);
    // "ORIGINAL" MESH IS MODIFIED
    geometries.geom_orig->setTo( geometries.geom );

}
