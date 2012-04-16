#include <eventlist.h>
#include <eventhandler.h>
#include <refinfo.h>
#include <refinement.h> // declares autorefinement etc.
#include <list>
#include <qlc3d.h>
void handleMeshRefinement(std::list<Event*>& refEvents,
                          Geometries& geometries,
                          SolutionVectors& solutionvectors,
                          Simu& simu,
                          Alignment& alignment,
                          Electrodes& electrodes,
                          LC& lc,
                          SparseMatrix& Kpot,
                          SparseMatrix& Kq
                          )
{
    printf("%u REFINEMENT EVENTS NOW\n", refEvents.size() );

    // MAKE LIST OF ALL REFINFO OBJECTS THAT ARE HANDLES NOW
    std::list<RefInfo> refInfos;
    std::list<Event*>::iterator evitr = refEvents.begin();
    for ( ; evitr!= refEvents.end() ; ++evitr)
    {
        RefInfo* ref = static_cast<RefInfo*>( (*evitr)->getEventDataPtr() );
        refInfos.push_back( *ref );
    }
    printf("%u RefInfo objects\n", refInfos.size() );

    // TRY TO DO REFINEMENT
    bool isRefined(false);
    isRefined = autoref(*geometries.geom_orig,
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

    // DELETE ALL REFINEMENT EVENTS. ALWAYS
    for (evitr = refEvents.begin() ; evitr != refEvents.end() ; ++evitr)
    {
        delete (*evitr);
    }


    // IF MESH HAS BEEN REFINED NEED TO RECREATE MATRIXES
    // FOR Q-TENSOR AND POTENTIAL
    if (isRefined)
    {
        Kpot.PrintInfo();

        printf("recreating matrixes for: "); fflush(stdout);
        Kpot.~SparseMatrix();
        Kq.~SparseMatrix();
        printf("V..."); fflush(stdout);
        Kpot = *createSparseMatrix(*geometries.geom,
                                   *solutionvectors.v);
        printf("OK. Q...");fflush(stdout);
        Kq = *createSparseMatrix(*geometries.geom,
                                 *solutionvectors.q,
                                 MAT_DOMAIN1);
       Kpot.PrintInfo();

        printf("OK\n"); fflush(stdout);
    }
}


void handlePreRefinement(std::list<Event*>& refEvents,
                         Geometries& geometries,
                         SolutionVectors& solutionvectors,
                         Simu& simu,
                         Alignment& alignment,
                         Electrodes& electrodes,
                         LC& lc,
                         SparseMatrix& Kpot,
                         SparseMatrix& Kq)
{
// PRE REFINMENT MODIFICATIONS TO THE MESH CARRY THROUGH THE
// REST OF THE SIMULATION. THAT IS, THE INITIAL GEOMETRY IS
// MODIFIED TOO

    printf("PRE-REFINEMENT\n");
    // MAKE LIST OF ALL REFINFO OBJECTS
    handleMeshRefinement(refEvents,
                         geometries,
                         solutionvectors,
                         simu,
                         alignment,
                         electrodes,
                         lc,
                         Kpot,
                         Kq);
    // "ORIGINAL" MESH IS MODIFIED
    geometries.geom_orig->setTo( geometries.geom );

}
