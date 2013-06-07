#include <eventlist.h>
#include <eventhandler.h>
#include <refinfo.h>
#include <refinement.h> // declares autorefinement etc.
#include <list>
#include <qlc3d.h>

#include <spamtrix_ircmatrix.hpp>

void handleMeshRefinement(std::list<Event*>& refEvents,
                          Geometries& geometries,
                          SolutionVectors& solutionvectors,
                          Simu& simu,
                          Alignment& alignment,
                          Electrodes& electrodes,
                          LC& lc,
                          SpaMtrix::IRCMatrix &Kpot,
                          SpaMtrix::IRCMatrix &Kq
                          )
{

    cout << refEvents.size() << " REFINEMENTS NOW\n" << endl;
    // MAKE LIST OF ALL REFINFO OBJECTS THAT ARE HANDLES NOW
    std::list<RefInfo> refInfos;
    std::list<Event*>::iterator evitr = refEvents.begin();
    for ( ; evitr!= refEvents.end() ; ++evitr){
        RefInfo* ref = static_cast<RefInfo*>( (*evitr)->getEventDataPtr() );
        refInfos.push_back( *ref );
    }

    cout << refInfos.size() << "Refinement objects" << endl;
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
    for (evitr = refEvents.begin() ; evitr != refEvents.end() ; ++evitr){
        delete (*evitr);
    }


    // IF MESH HAS BEEN REFINED NEED TO RECREATE MATRIXES
    // FOR Q-TENSOR AND POTENTIAL
    if (isRefined)
    {
        std::cout << "creating new matrixes" << std::endl;


        std::cout << "old size: " << Kpot.getNumRows() <<", " << Kpot.getNumCols() << std::endl;
        Kpot.clear();
        Kpot = createPotentialMatrix(*geometries.geom,
                                     *solutionvectors.v,
                                     0,
                                     electrodes);

        std::cout << "new size: " << Kpot.getNumRows() <<", " << Kpot.getNumCols() << std::endl;
        std::cout << "exiting in " << __func__ << " need to make new matrix " << std::endl;
        exit(-1979);

        /*// DISABLED DURING SPAMTRIX MIGRATION
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
        */
    }
}


void handlePreRefinement(std::list<Event*>& refEvents,
                         Geometries& geometries,
                         SolutionVectors& solutionvectors,
                         Simu& simu,
                         Alignment& alignment,
                         Electrodes& electrodes,
                         LC& lc,
                         SpaMtrix::IRCMatrix &Kpot,
                         SpaMtrix::IRCMatrix &Kq)
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
