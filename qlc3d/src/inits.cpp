#include <qlc3d.h>
#include <simu.h>
void prepareGeometry(Geometry& geom,
                     Simu& simu,
                     Alignment& alignment)
{

    int np,nt,ne;
    double *p;
    int *t;
    int *e;
    int *emat;
    int *tmat;

    // Check whether to read mesh from text file or binary 'geo'-file
    // determine by file extension
    std::string filename = simu.MeshName;
    size_t point = filename.find_last_of('.');  // index to separator point
    string extension = filename.substr(point + 1);

    if (extension.compare("geo") == 0 )
    {
        printf(" reading .geo file\n");
        readBinaryMesh(filename, p ,t, tmat, e, emat, &np, &nt, &ne);
    }
    else
    {
        printf(" reading .%s file\n", extension.c_str() );
        ReadGiDMesh3D(&simu,&p,&np,&t,&nt,&e,&ne,&tmat,&emat);
    }
    for (int i = 0; i < 4*nt; i++)  t[i]--;	// change GiD mesh to zero indexing
    for (int i = 0; i < 3*ne; i++)  e[i]--;

    for (int i = 0; i < np;   i++){	// scale mesh
        p[3*i + 0] = p[3*i + 0]* simu.getStretchVectorX();
        p[3*i + 1] = p[3*i + 1]* simu.getStretchVectorY();
        p[3*i + 2] = p[3*i + 2]* simu.getStretchVectorZ();
    }


    geom.setCoordinates(p, np);
    if (p) free(p);

    // set up volume mesh
    geom.t->setnElements(nt);		// set up tetrahedral mesh
    geom.t->setnNodes(4);		// number of nodes per element
    geom.t->setDimension(3);		// 3D element
    geom.t->AllocateMemory();		// allocate memory for arrays
    geom.t->setAllNodes(t);		// copy node numbers
    geom.t->setAllMaterials(tmat);	// copy material numbers
    geom.t->setMaxNodeNumber( (unsigned int) np ); // total number of nodes in tet mesh
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
    geom.checkForPeriodicGeometry();
    geom.genIndWeakSurfaces(alignment);

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
