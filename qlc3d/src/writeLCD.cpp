#include <geometry.h>
#include <stdio.h>
#include <qlc3d.h>
#include <string>
#include <iostream>
#include <simu.h>


void WriteMesh(Simu* simu, double* p, Mesh* t, Mesh* e, int np)
{
	int i;

	string fname = simu->getSaveDir() + simu->getMeshFileNameOnly();
	FILE *fid = fopen(fname.c_str(),"w");

	if (fid!=NULL)
	{
	fputs ("MESH    dimension 3 ElemType Tetrahedra  Nnode 4\nCoordinates\n",fid);
	for (i=0;i<np;i++)
		{
		fprintf(fid,"%i\t%f\t%f\t%f\n",i+1,p[3*i],p[3*i+1],p[3*i+2]);
		}
	fprintf(fid,"end coordinates\n\nElements\n");
	
	for (i = 0 ; i < t->getnElements() ; i++)
		fprintf(fid,"%i\t%i\t%i\t%i\t%i\t%i\n",i+1, t->getNode(i,0) +1, t->getNode(i,1)+1, t->getNode(i,2)+1, t->getNode(i,3)+1, t->getMaterialNumber(i));
	
	fprintf(fid,"end elements\nMESH    dimension 3 ElemType Triangle  Nnode 3\nCoordinates\nend coordinates\n\nElements\n");
	
        //int nNodes = e->getnNodes();
	for (i = 0 ; i < e->getnElements() ; i++)
		fprintf(fid,"%i\t%i\t%i\t%i\t%i\n",i+1, e->getNode(i,0)+1, e->getNode(i,1)+1, e->getNode(i,2)+1, e->getMaterialNumber(i) );
	fprintf(fid,"end elements\n");
	fclose(fid);
	}
	else
	{
	printf("WriteLCD - could not open file for output mesh: %s\n - bye!", fname.c_str());
	exit(1);
	}


}

void WriteLCD(double *p, Mesh *t, Mesh *e, SolutionVector *v, SolutionVector *q, Simu* simu)
{
// First write mesh
	int np = v->getnDoF();
    if ( simu->IsMeshModified() ){ // only output mesh file if it has changed
		WriteMesh(simu,p,t,e,np);
        simu->setMeshModified( false ); // no need to rewrite same mesh again
    }

// MESH FILE WRITTEN - THEN WRITE RESULT
	FILE* fid;
	int i;
	
    string resname = simu->getSaveDir();//"res/result";//
	char str[15];
	
	// check whether final result in sumlation, if yes, use special filename 
	if (simu->getCurrentIteration() != SIMU_END_SIMULATION )
		sprintf(str, "result%05i", simu->getCurrentIteration()-1 );
	else
		sprintf(str,"result_final");
	
	resname.append(str);
	resname.append(".dat");
	
	
	
	fid = fopen(resname.c_str(),"w");
    if (fid!=NULL){
		int npLC = q->getnDoF();
		double *n = tensortovector(q->Values,npLC); // get vector data 
		sprintf(str,"%f",simu->getCurrentTime());
		std::string text = "** Result Time :    ";
		text.append(str);
		text.append("\n** z Compression Ratio :  1.00000\n");

		// trickery needed for QT/no QT filenaming
		std::string meshname("");
		#ifdef NO_QT
		    std::string temp;
		    size_t pos = simu->getSaveDir().find_last_of("/");
		    temp = simu->getSaveDir().substr(pos+1, simu->getSaveDir().length() - pos);
		    meshname = temp;
		#endif
		meshname + simu->getMeshFileNameOnly() + "\n";
		//cout <<"meshname :" << mehsname << endl;
		text.append(meshname);
		fprintf(fid,"%s", text.c_str());//** Result Time :    0.00000000\n** z Compression Ratio :  1.00000\nmeshout.txt\n");
	
		for (i=0;i<np;i++)
			{
			if (i<npLC)
				fprintf(fid,"    %i   %f   %f   %f   %f   %f   %f\n",i+1,n[i],n[i+npLC],n[i+2*npLC],v->Values[i],n[i+3*npLC],n[i+4*npLC]);
			else
				fprintf(fid,"%i\t0\t0\t0\t%f\t1.0\t0.0\n",i+1,v->Values[i]);
			}
	
		fclose(fid);
		free(n);
	}
	else
	{
        printf("writeLCD - could not open result file %s!\n", resname.c_str() );
		exit(1);
	}
}
// end WriteLCD

void WriteLCD_B(double *p, Mesh *t, Mesh *e, SolutionVector *v, SolutionVector *q,
                Simu* simu, LC* lc)
{

    int npLC = q->getnDoF();
 /*

    double *n = tensortovector(q->Values,npLC); // get vector data
    for (int i = 0 ; i < npLC ; i++){

        if (n[i+3*npLC] > 0.6 ){
       // printf("%f, %f, %f, %f, %f \n ", n[i] , n[i+npLC], n[i+2*npLC], n[i+3*npLC] , n[i+4*npLC]);
        printf("node %i  = %f\n",i, n[i+3*npLC]);
      }

    }
    free( n );
//*/



    int np = v->getnDoF();

	// First write mesh
    if ( simu->IsMeshModified() ){ // only output mesh file if it has changed
		WriteMesh(simu,p,t,e,np);
        simu->setMeshModified( false ); // no need to rewrite same mesh again
    }
	
	// THEN WRITE RESULT FILE
	FILE* fid;
    int i;
    //npLC = q->getnDoF();
	

        string resname = simu->getSaveDir();
	char str[15];
	
	// check whether final result in simulation, if yes, use special filename
	if (simu->getCurrentIteration() != SIMU_END_SIMULATION )
            sprintf(str, "result%05i", simu->getCurrentIteration()-1 );
	else
            sprintf(str,"result_final");
	
	resname.append(str);
	resname.append(".dat");
		
	fid = fopen(resname.c_str(),"wb");

	char time[20];
		
	sprintf(time,"%1.9f\n", simu->getCurrentTime());
		
	string text ="** Result Time :   ";
	text.append(time);
		
	fprintf(fid,"%s\n",text.c_str());
	fprintf(fid,"** z Compression Ratio :  1.00000\n");

	// trickery needed for QT/no QT filenaming
	std::string meshname("");
	#ifdef NO_QT
	    std::string temp;
	    size_t pos = simu->getSaveDir().find_last_of("/");
	    temp = simu->getSaveDir().substr(pos+1, simu->getSaveDir().length() - pos);
	    meshname = temp;
	#endif
	meshname.append( simu->getMeshFileNameOnly() );// + "\n";
	meshname.append("\n");
	//cout << "meshname =" << meshname << endl;
	fprintf(fid,"%s",meshname.c_str());

	fprintf(fid,"RAW FLOAT TRI - S0, np, nsols\n");
	fprintf(fid,"%g %d %d\r\n",lc->getS0(),np,6);
	
		
	
	for (i=0;  i<np; i++){
	    if (i<npLC){ // WRITE LC REGIONS
		float value;

		value = (float)q->getValue(i,0);
		fwrite((void*) &value,sizeof(float),1, fid );
		value = (float)q->getValue(i,1);
		fwrite((void*) &value,sizeof(float),1, fid );
		value = (float)q->getValue(i,2);
		fwrite((void*) &value,sizeof(float),1, fid );
		value = (float)q->getValue(i,4);
		fwrite((void*) &value,sizeof(float),1, fid );
		value = (float)q->getValue(i,3);
		fwrite((void*) &value,sizeof(float),1, fid );
		value = v->getValue(i);
		fwrite((void*) &value , sizeof(float),1 ,fid );
    }
	    else{ // WRITE DIELECTRIC REGIONS
		float value = 0;
		for (int qc =0; qc<5; qc++)	// q-component
		fwrite((void*) &value,sizeof(float),1,fid);

		value = v->getValue(i);
		fwrite((void*) &value, sizeof(float),1,fid);
	    }
	}
	fclose(fid);
		
}



void ReadLCD_B(Simu* simu, SolutionVector *q)
{
// READS BINARY FORMATED RESULT FILE
    string filename = "res/" + simu->getLoadDir() + simu->getLoadQ();
    printf( "Loading Q-tensor from: %s\n",filename.c_str());
    FILE* fid = fopen( filename.c_str() , "rb" );
		
    char str[100];
	if (fid)
	{
		
		float S0;
		int	np, nsol;
		char* tempch;

		tempch = fgets (str , 100 , fid); // READS LINE x 5
		tempch = fgets (str , 100 , fid);
		tempch = fgets (str , 100 , fid);
		tempch = fgets (str , 100 , fid);
		tempch = fgets (str , 100 , fid);
		
		int i = fscanf (fid, "%f %i %i\n",&S0,&np,&nsol);
		i = 0; // no warnings...
		
		printf("S0 = %f, np = %i, nsol = %i \n", S0, np, nsol);
		if (np < q->getnDoF() )
        {
			printf("The file you're trying to load (\"%s\") doesn't seem to match the mesh - bye!\n",filename.c_str() );
			fclose(fid);
			exit(1);
		}
		
		float q1,q2,q3,q4,q5,temp;
	
                size_t size = 0;
		
		for (int i = 0 ; i < q->getnDoF() ; i ++){
			size  = fread ( &q1, sizeof(float), 1 , fid );
			size  = fread ( &q2, sizeof(float), 1 , fid );
			size  = fread ( &q3, sizeof(float), 1 , fid );
			size  = fread ( &q5, sizeof(float), 1 , fid );
			size  = fread ( &q4, sizeof(float), 1 , fid );
			q->setValue(i,0,q1);
			q->setValue(i,1,q2);
			q->setValue(i,2,q3);
			q->setValue(i,3,q4);
			q->setValue(i,4,q5);
									
			for (int j = 0 ; j < nsol-5 ; j ++)			// READ&DISCARD POTENTIAL AND FLOW
				size = fread ( (void *) &temp, sizeof(float), 1 , fid );
		}
		
		
		printf("OK\n");
				fclose(fid);
	}
	else
	{
		printf("could not open file: %s. It probably doesn't exist, but I could be wrong...\n", filename.c_str());
		exit(1);
	}

}

void WriteSettings(Simu* simu , LC* lc, Boxes* box, Alignment* alignment, Electrodes* electrodes){
	string fname = simu->getSaveDir() + "settings.qfg";

	FILE *fid = fopen(fname.c_str(),"wt");
	if (fid!=NULL)
	{
                //std::cout<< "writing ";
                simu->WriteSimu(fid);
                //std::cout << "simu";
		lc->WriteLC(fid);
                //std::cout<< "lc";
		box->WriteBoxes(fid);
		alignment->WriteAlignment(fid);
		electrodes->WriteElectrodes(fid);
		
		fclose(fid);
                //std::cout<<"done"<<std::endl;
	}
	else
	{
		printf("error in WriteSettings - could not open settings file %s - bye!\n",fname.c_str());
		exit(1);
	}

}

void WriteResult(
	Simu* simu, 
	LC* lc,
	Geometry* geom,
	SolutionVector* v,
    SolutionVector* q,
    MeshRefinement* meshref){

// DETERMINE IF SAVING OF RESULT IS NEEDED THIS ITERATION.
// IF MODULUS OF CURRENT ITERATION NUMBER AND SAVEITER IS 0 -> SAVE
// ALSO, SAVE INITIAL AND FINAL RESULTS

int mod = 1; // assume NO by default
// take care of undefined x%0
if (simu->getSaveIter() != 0)
	mod = simu->getCurrentIteration() % simu->getSaveIter() ;

// If a new mesh has been generated
if ( (meshref) && (meshref->isNeedsNewMesh() ) ){
    mod = 0;
}

if (// WRITE RESULT FILE IF...
	(mod == 0 )|| // ITERATION IS A SAVE ITERATION
	(simu->getCurrentIteration() == 0) ||	// ...OR ITERATION IS INITIAL CONFIGURATION
	(simu->getCurrentIteration() == SIMU_END_SIMULATION )) // ...OR FINAL RESULT

    {
        switch( simu->getOutputFormat() ){
        case SIMU_OUTPUT_FORMAT_BINARY:{
            WriteLCD_B(geom->getPtrTop(), geom->t, geom->e, v, q, simu, lc); // WRITES BINARY FILE
                break;
        }
        case SIMU_OUTPUT_FORMAT_TEXT:{
            WriteLCD(geom->getPtrTop(), geom->t, geom->e, v, q, simu); // WRITES TEXT FILE
            break;
        }
        default:{
            printf("error in WriteResult, Simu.OutputFormat is not recognised - bye!\n");
            fflush(stdout);
            exit(1);
        }
		} // end switch
	}// end if write result file
}
// end void WriteREsult


void CreateSaveDir(Simu* simu){
// CREATES NEW SAVE DIRECTORY IF COMILED WITH QT
// IF QT IS DISABLED, THE SAVEDIR STRING IS SET AS
// THE BEGINNING PART OF THE FILENAME
	using std::cout;
	using std::endl;
	std::string saveroot("./res/");
    std::string temp;
    size_t pos;


// NO QT OPTION
//#ifdef NO_QT
	saveroot.append(simu->getSaveDir());
	// MAKE SURE LAST CHARACTER OF SAVEDIR IS NOT A '/'

	   pos = saveroot.find_last_of("/");
	   if (pos == saveroot.length() )
	   saveroot = saveroot.substr(0 , pos-1);
	   simu->setSaveDir(saveroot);
	   std::cout << "save dir is set to " << simu->getSaveDir() << std::endl;
	   return;
//#endif





}
// end void CreateSaveDir
