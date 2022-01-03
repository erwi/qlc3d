#include <iostream>
#include <inits.h>
#include <simu.h>
#include <globals.h>
#include <eventlist.h>
#include <refinement.h>
#include <io/meshreader.h>

void validateTriangleMaterials(const idx* const mate, idx ne){
    // goes through all triangle material numbers and tries to check that all is well.
    // if problems are found, error message is printed and program aborted.
    using std::cout;
    using std::endl;
    bool isBad = false;
    idx m;
    for (idx i = 0 ; i < ne ; ++i){
        m = mate[i];
        if ( (m == MAT_PERIODIC) || (m == MAT_NEUMANN)){
            continue;
        }
        size_t eNum = MATNUM_TO_ELECTRODE_NUMBER((size_t) m);
        size_t fNum = MATNUM_TO_FIXLC_NUMBER((size_t) m);
        if (m < MAT_ELECTRODE1){
            isBad = true;
            break;
        }
        else if (eNum > 9){
            cout << "\nbad electrode number detected " << eNum << endl;
            isBad = true;
            break;
        }
        else if (fNum > 9){
            cout << "\nerror, bad FixLC number detected " << fNum << endl;
            isBad = true;
            break;
        }
    }
    if (isBad){
        cout << "error, bad triangle material number " << m << " found in mesh - bye!" << endl;
        exit(1);
    }
}// end void validateTriangleMatrials

size_t countElectrodes(const idx* mate, idx ne) {
    /*!
    * Counts the number of electrodes defined in mesh. As in
    * Electrode 1, 2 ... not the number of electrode elements.
    *
    * Also checks that electrode surfaces are defined contiguously
    * in mesh file (i.e. 1, 2, 3,...  and not 1, 3, 4...
    *
    * mate is material numbers array and ne is array length.
    */
    // First set flags for those electrode numbers that are found
    std::vector<int> eFoundFlag(MAT_MAX_ELECTRODES_COUNT + 1);
    for (idx i = 0; i < ne; i++) {
        const size_t enumber = MATNUM_TO_ELECTRODE_NUMBER(mate[i]);
        if (enumber > 0) // if valid electrode number, set flag
            eFoundFlag[enumber] = 1;
    }
    //
    // Then find largest electrode whose flag is set
    size_t eCount = 0;
    for (idx i = eFoundFlag.size()-1; i > 1 ; i--) { // backwards loop
        if ((eFoundFlag[i] > 0) && (i > eCount))
            eCount = i;
        //
        // detect error, non-contiguous electrodes definitions
        if ((eFoundFlag[i] > 0) && (eFoundFlag[i-1] == 0)) {
            std::cerr << "\nerror in " << __func__ << std::endl;
            std::cerr << "Found error in mesh file: non-contiguously numbered electrode surfaces." << std::endl;
            std::cerr << "Electrode " << i << " (at least) is missing although "<< i+1 << " exists.\nbye!" << std::endl;
            std::exit(1);
        }
    }
   return eCount+1;
}

void validateTetrahedralMaterials(const idx* const matt, idx nt){
    // goes through each material number for tetrahedra. if bad ones are
    // found, error is printed and program terminated
    using std::cout;
    using std::endl;
    for (idx i = 0 ; i < nt ; ++i){
        const idx m = matt[i];
        if ((m == MAT_DOMAIN1) ||
                (m==MAT_DIELECTRIC1) || (m==MAT_DIELECTRIC2) || (m==MAT_DIELECTRIC3) ||
                (m==MAT_DIELECTRIC4) || (m==MAT_DIELECTRIC5) || (m==MAT_DIELECTRIC6) ||
                (m==MAT_DIELECTRIC7) ){
            continue;
        }
        else{
            cout << "error, bad tetrahedral (volume) material number " << m <<" found in mesh." << endl;
            cout << "valid volume material are Domain1 and Dielectrics 1-7. Goodbye!" << endl;
            exit(1);
        }
    }
}

void prepareGeometry(Geometry& geom,
                     const std::string &meshFileName,
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

    // PEOPLE DO STUPID THINGS IN GiD. NEED TO VALIDATE MATERIAL NUMBERS
    // IF VALIDATION FAILS, ERROR MESSAGE IS PRINTED AND PROGRAM ABORTED
    validateTriangleMaterials(emat, ne);
    validateTetrahedralMaterials(tmat, nt);

    // Count number of electrodes in mesh
    electrodes.setnElectrodes(countElectrodes(emat, ne));

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

FILE* createOutputEnergyFile(Simu& simu){

    FILE* fid = NULL;

    if (simu.getOutputEnergy() == 1){
        string energy_fn = simu.getSaveDir() + "/" + "energy.m";
        fid = fopen( energy_fn.c_str() , "w");
        if (fid == NULL)	{
            printf("error - could not open output file for free energy- bye ! \n");
            exit(0);
        }
    }
    return fid;
}

void createMeshRefinementEvents(const MeshRefinement &meshRefinement,
                                EventList &eventListOut) {
    const unsigned int repRefIter = meshRefinement.getRepRefIter();
    const double repRefTime = meshRefinement.getRepRefTime();
    bool hasPeriodicRefinement = repRefIter > 0 || repRefTime > 0;
    cout << "periodic mesh refinement at every " << repRefIter << " iterations, " << repRefTime << " seconds" << endl;
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
            throw runtime_error("no refinement objects defined for periodic mesh refinement");
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