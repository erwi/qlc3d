#include <eventlist.h>
#include <eventhandler.h>
#include <refinfo.h>
#include <refinement.h> // declares autorefinement etc.
#include <list>

void handleMeshRefinement(std::list<Event*>& refEvents,
                          Geometries& geometries,
                          SolutionVectors& solutionvectors,
                          Simu& simu,
                          Alignment& alignment,
                          Electrodes& electrodes,
                          LC& lc
                          )
{
    printf("%u REFINEMENT EVENTS NOW\n", refEvents.size() );

    // MAKE LIST OF ALL REFINFO OBJECTS
    std::list<RefInfo> refInfos;
    std::list<Event*>::iterator evitr = refEvents.begin();
    for ( ; evitr!= refEvents.end() ; evitr++)
    {
        RefInfo* ref = static_cast<RefInfo*>( (*evitr)->getEventDataPtr() );
        refInfos.push_back( *ref );
    }
    printf("%u RefInfo objects\n", refInfos.size() );

    autoref(*geometries.geom_orig,
            *geometries.geom_prev,
            *geometries.geom,
            *solutionvectors.q,
            *solutionvectors.qn,
            *solutionvectors.v,
            refInfos,
            simu,
            alignment,
            electrodes,
            lc
            );



    // DELETE ALL REFINEMENT EVENTS
    for (evitr = refEvents.begin() ; evitr != refEvents.end() ; evitr++)
    {
        delete (*evitr);
    }
    return;

}
