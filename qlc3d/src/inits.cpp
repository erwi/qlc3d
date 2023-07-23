#include <iostream>
#include <inits.h>
#include <simu.h>
#include <eventlist.h>
#include <refinement.h>
#include <io/meshreader.h>
#include <util/logging.h>
#include <util/exception.h>
/**
 * Goes through all triangle material numbers and tries to check that all is well.
 * @param mate triangle material numbers
 * @param ne number of triangles
 */
void validateTriangleMaterials(const idx* const mate, idx ne, const Electrodes &electrodes) {
    for (idx i = 0; i < ne; ++i) {
        idx m = mate[i];
        if ((m == MAT_PERIODIC) || (m == MAT_NEUMANN)) {
            continue;
        }
        size_t eNum = MATNUM_TO_ELECTRODE_NUMBER((size_t) m);
        size_t fNum = MATNUM_TO_FIXLC_NUMBER((size_t) m);
        if (m < MAT_ELECTRODE1) {
            RUNTIME_ERROR(fmt::format("Triangle {}, invalid material number {} < {}.", i, m, MAT_ELECTRODE1));
        }
        else if (eNum > 9) {
            RUNTIME_ERROR(fmt::format("Triangle {}, invalid Electrode number {} > 9", i, m));
        }
        else if (fNum > 9) {
            RUNTIME_ERROR(fmt::format("Triangle {}, invalid FixLC number {} > 9", i, m));
        }
        else if (eNum > electrodes.getnElectrodes()) {
            RUNTIME_ERROR(fmt::format("Triangle {} has electrode number {} but only {} electrode(s) have been defined.", i, eNum, electrodes.getnElectrodes()));
        }
    }
}

/**
 * Goes through each material number for tetrahedra. Throws runtime_exception if invalid ones are found.
 * @param matt tetrahdra material numbers
 * @param nt number of tetrahedra
 */
void validateTetrahedralMaterials(const idx* const matt, idx nt) {
    using namespace std;
    using namespace fmt;
    for (idx i = 0; i < nt; ++i) {
        const idx m = matt[i];
        if ((m == MAT_DOMAIN1) ||
                (m==MAT_DIELECTRIC1) || (m==MAT_DIELECTRIC2) || (m==MAT_DIELECTRIC3) ||
                (m==MAT_DIELECTRIC4) || (m==MAT_DIELECTRIC5) || (m==MAT_DIELECTRIC6) ||
                (m==MAT_DIELECTRIC7) ){
            continue;
        }
        else {
            RUNTIME_ERROR(fmt::format("Invalid material number {} found in mesh for tetrahedron {}. "
                                       "Valid material numbers are Domain1 and Dielectrics 1-7.", m, i));
        }
    }
}

void prepareGeometry(Geometry& geom,
                     const std::filesystem::path &meshFileName,
                     Simu& simu,
                     Alignment& alignment,
                     Electrodes& electrodes) {
    idx np,nt,ne;
    double *p;
    idx *t;
    idx *e;
    idx *emat;
    idx *tmat;

    // read mesh data from file. Allocates the data arrays.
    MeshReader::readMesh(meshFileName, &p, &np, &t, &nt, &e, &ne, &tmat, &emat);

    // Throw exception if invalid element material numbers are detected.
    validateTriangleMaterials(emat, ne, electrodes);
    validateTetrahedralMaterials(tmat, nt);

    // Count number of electrodes in mesh
    for (idx i = 0; i < np;   i++) {	// scale mesh
        p[3*i + 0] = p[3*i + 0]* simu.getStretchVectorX();
        p[3*i + 1] = p[3*i + 1]* simu.getStretchVectorY();
        p[3*i + 2] = p[3*i + 2]* simu.getStretchVectorZ();
    }

    geom.setCoordinates(p, (size_t) np);
    if (p) free(p);

    // set up volume mesh
    geom.t->setnElements(nt);		// set up tetrahedral mesh
    geom.t->setnNodes(4);		// number of nodes per element
    geom.t->setDimension(3);		// 3D element
    geom.t->AllocateMemory();		// allocate memory for arrays
    geom.t->setAllNodes(t);		// copy node numbers
    geom.t->setAllMaterials(tmat);	// copy material numbers
    geom.t->setMaxNodeNumber( np ); // total number of nodes in tet mesh
    free(t);
    free(tmat);

    // set up surface mesh
    geom.e->setnElements(ne);
    geom.e->setnNodes(3);
    geom.e->setDimension(2);
    geom.e->AllocateMemory();
    geom.e->setAllNodes(e);
    geom.e->setAllMaterials(emat);
    free(e);
    free(emat);

    geom.ReorderDielectricNodes(); // Dielectric nodes are moved last
    geom.e->setConnectedVolume(geom.t);		// neighbour index tri -> tet
    geom.t->CalculateDeterminants3D(geom.getPtrTop());		// calculate tetrahedral determinants
    geom.t->ScaleDeterminants(1e-18); // scale to microns

    geom.e->CalculateSurfaceNormals(geom.getPtrTop() , geom.t);		// calculate triangle determinants and surface normal vectors
    geom.e->ScaleDeterminants(1e-12); // scale to microns

    geom.setNodeNormals();
    geom.checkForPeriodicGeometry(); // also makes periodic node indexes

    geom.genIndWeakSurfaces(alignment);
    geom.makeRegularGrid(simu.getRegularGridXCount(),
                         simu.getRegularGridYCount(),
                         simu.getRegularGridZCount());
}

FILE* createOutputEnergyFile(Simu& simu) {
    FILE* fid = nullptr;
    if (simu.getOutputEnergy() == 1) {
      std::filesystem::path energyFilePath = simu.getSaveDirAbsolutePath() / "energy.m";
      fid = fopen( energyFilePath.string().c_str(), "w");
      if (fid == nullptr) {
        RUNTIME_ERROR(fmt::format("could not open output file for free energy by filename = {}.", energyFilePath));
      }
    }
    return fid;
}

void createMeshRefinementEvents(const MeshRefinement &meshRefinement,
                                EventList &eventListOut) {
    const unsigned int repRefIter = meshRefinement.getRepRefIter();
    const double repRefTime = meshRefinement.getRepRefTime();
    bool hasPeriodicRefinement = repRefIter > 0 || repRefTime > 0;
    Log::info("Creating periodic mesh refinement events at every {} iterations, {} seconds", repRefIter, repRefTime);
    eventListOut.setRepRefIter(repRefIter);
    eventListOut.setRepRefTime(repRefTime);

    if (hasPeriodicRefinement) {
        // create refinement event for periodically occurring mesh refinement.
        // These have no explicitly defined iterations or times
        for (auto &ref : meshRefinement.getRefinementConfig()) {
            if (!ref.occursPeriodically()) {
                continue;
            }
            RefInfo *info = RefInfo::ofPeriodicMeshRefinement(ref.type_, ref.values_, ref.x_, ref.y_, ref.z_);
            Event *event = Event::ofPeriodicMeshRefinement(info);
            eventListOut.addRepRefInfo(event);
        }

        // e.g. if user has defined RepRefIter > 0 or RepRefTime > 0, but no refinement objects with
        // empty explicit times/iterations lists
        if (eventListOut.getNumPeriodicRefinementObjects() == 0) {
            RUNTIME_ERROR("No refinement objects defined for periodic mesh refinement");
        }
    }

    // add refinement events occurring at explicitly defined iterations/times
    for (auto &ref : meshRefinement.getRefinementConfig()) {
        if (ref.occursPeriodically()) {
            continue;
        }

        // convert each explicitly defined refinement iteration to an *Event and add to event list
        for (unsigned int iter : ref.iterations_) {
            RefInfo * refInfo = RefInfo::make(ref.type_, iter, -1,
                                              ref.values_, ref.x_, ref.y_, ref.z_);
            auto *event = new Event(EVENT_REFINEMENT, iter, (void*) refInfo);
            eventListOut.insertIterEvent(event);
        }

        // convert each explicitly defined refinement time to an *Event adn add to event list
        for (double time : ref.times_) {
            RefInfo *refInfo = RefInfo::make(ref.type_, -1, time,
                                             ref.values_, ref.x_, ref.y_, ref.z_);
            auto *event = new Event(EVENT_REFINEMENT, time, (void*) refInfo);
            eventListOut.insertTimeEvent(event);
        }
    }
}