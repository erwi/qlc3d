#include <geometry.h>
#include <stdio.h>
#include <qlc3d.h>
#include <string>
#include <iostream>
#include <simu.h>
#include <filesysfun.h>
#include <globals.h>
#include <resultoutput.h>


namespace LCviewIO{

void WriteMesh(Simu* simu,
               double* p,
               Mesh* t,
               Mesh* e,
               idx np)
{
    idx i;
    string fname = simu->getMeshFileNameOnly();
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


        for (i = 0 ; i < e->getnElements() ; i++)
            fprintf(fid,"%i\t%i\t%i\t%i\t%i\n",i+1, e->getNode(i,0)+1, e->getNode(i,1)+1, e->getNode(i,2)+1, e->getMaterialNumber(i) );

        fprintf(fid,"end elements\n");
        fclose(fid);
    }
    else
    {
        printf("error - could not open file for output mesh: %s\n - bye!", fname.c_str());
        exit(1);
    }


}

void WriteLCD(double *p, Mesh *t, Mesh *e, SolutionVector *v, SolutionVector *q, Simu* simu)
{
    /*!
      Writes LCview result in text format.
    */

    // First write mesh
    int np = v->getnDoF();
    if ( simu->IsMeshModified() ){ // only output mesh file if it has changed
        WriteMesh(simu,p,t,e,np);
        simu->setMeshModified( false ); // no need to rewrite same mesh again
    }

    // MESH FILE WRITTEN - THEN WRITE RESULT
    FILE* fid;
    int i;


    char str[15];
    // check whether final result in sumlation, if yes, use special filename
    if (simu->getCurrentIteration() != SIMU_END_SIMULATION )
        sprintf(str, "result_t%05i", simu->getCurrentIteration() );
    else
        sprintf(str,"result_final");

    string resname = str;
    resname.append(".dat");

    fid = fopen(resname.c_str(),"w");
    if (fid!=NULL){
        int npLC = q->getnDoF();
        double *n = tensortovector(q->Values,npLC); // get vector data
        sprintf(str,"%f",simu->getCurrentTime());
        std::string text = "** Result Time :    ";
        text.append(str);
        text.append("\n** z Compression Ratio :  1.00000\n");

        std::string meshname = simu->getMeshFileNameOnly() + "\n"; // meshfilename

        text.append(meshname);
        fprintf(fid,"%s", text.c_str());//** Result Time :    0.00000000\n** z Compression Ratio :  1.00000\nmeshout.txt\n");

        for (i=0;i<np;i++){
            if (i<npLC)
                fprintf(fid,LCVIEW_TEXT_FORMAT_STRING,i+1,n[i],n[i+npLC],n[i+2*npLC],v->Values[i],n[i+3*npLC],n[i+4*npLC]);
            else

                fprintf(fid,LCVIEW_TEXT_FORMAT_STRING,i+1,0.,0.,0.,v->Values[i],0.,0.);
        }

        fclose(fid);
        delete [] n;
    }
    else{
        printf("writeLCD - could not open result file %s!\n", resname.c_str() );
        exit(1);
    }
}
// end WriteLCD

void WriteLCD_B(double *p, Mesh *t, Mesh *e, SolutionVector *v, SolutionVector *q,
                Simu* simu, LC* lc)
{

    int npLC = q->getnDoF();
    int np = v->getnDoF();

    // First write mesh
    if ( simu->IsMeshModified() )// only output mesh file if it has changed
    {
        WriteMesh(simu,p,t,e,np);
        simu->setMeshModified( false ); // no need to rewrite same mesh again
    }

    // THEN WRITE RESULT FILE
    FILE* fid;
    int i;
    char str[15];

    // check whether final result in simulation, if yes, use special filename
    if (simu->getCurrentIteration() != SIMU_END_SIMULATION )
        sprintf(str, "result%05i", simu->getCurrentIteration() );
    else
        sprintf(str,"result_final");

    string resname = str;
    resname.append(".dat");

    fid = fopen(resname.c_str(),"wb");

    char time[20];

    sprintf(time,"%1.9f\n", simu->getCurrentTime());

    string text ="** Result Time :   ";
    text.append(time);

    fprintf(fid,"%s\n",text.c_str());
    fprintf(fid,"** z Compression Ratio :  1.00000\n");

    std::string meshname = simu->getMeshFileNameOnly();
    meshname.append("\n");
    //cout << "meshname =" << meshname << endl;
    fprintf(fid,"%s",meshname.c_str());
    fprintf(fid,"RAW FLOAT TRI - S0, np, nsols\n");
    fprintf(fid,"%g %d %d\r\n",lc->getS0(),np,6);

    for (i=0;  i<np; i++)
    {
        if (i<npLC)
        { // WRITE LC REGIONS
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
        else
        { // WRITE DIELECTRIC REGIONS
            float value = 0;
            for (int qc =0; qc<5; qc++)	// q-component
                fwrite((void*) &value,sizeof(float),1,fid);

            value = v->getValue(i);
            fwrite((void*) &value, sizeof(float),1,fid);
        }
    }
    fclose(fid);

}


void ReadResult(Simu& simu, SolutionVector& q)
{
    /*!
 * Tries to figure out whether a LCView result file on disk
 * is in text or binary format and then load the data using appropriate
 * loading function
 */
    string filename = simu.getLoadQ();
    cout << "\nReading result file: " << filename << endl;
    FILE* fid = fopen (filename.c_str() , "rt" );
    if (!fid)
    {
        cout << "error, could not open file: " << filename << " - bye!" << endl;
        exit(1);
    }

    /// RAD SOME LINES FROM THE FILE AND TRY TO FIND OUT WHICH TYPE IT IS
    bool isBinary = false;
    char line[100];
    // IF FILE CONTAINS BELOW MAGIC TEXT, IT IS IN BINARY MODE
    const char binaryMarker[] = "RAW FLOAT TRI";
    for (int i = 0 ; i < 5 ; i++)
    {
        fgets(line, 100, fid); // returns null pointer if fails to read

        string sline = line;
        if (sline.find(binaryMarker) < std::string::npos )
        {
            //cout << "on line " << i << " detected key-word " <<binaryMarker << endl;
            isBinary = true;
            break;
        }
    }

    fclose(fid);

    if (isBinary) {
        cout << "Result file format is binary" << endl;
        LCviewIO::ReadLCD_B(&simu,&q);
    }
    else  {
        cout << "Result file format is text" << endl;
        LCviewIO::ReadLCD_T(simu,q);
    }



}

void ReadLCD_T(Simu &simu, SolutionVector &q)
{
    /*!
     * loads Q-tensor from a result file, assuming the file is written as a text file
     */

    string filename = simu.getLoadQ();
    FILE* fid = fopen( filename.c_str(), "rt");
    if (!fid)
    {
        cout << "error, could not open " << filename << " - bye!" << endl;
        exit(1);
    }
    const int lineLen = 256;
    char line[lineLen];

    // READ 3  LINES OF HEADER DATA
    char* tempc = fgets(line, lineLen, fid);
    tempc = fgets(line, lineLen, fid);
    tempc = fgets(line, lineLen, fid);

    int id;
    float n[3], S[2], v;

    //fgets(line, lineLen, fid);
    std::vector<float> q1;
    std::vector<float> q2;
    std::vector<float> q3;
    std::vector<float> q4;
    std::vector<float> q5;

    double rt6 = sqrt(6.0);
    double rt2 = sqrt(2.0);

    // READ FROM FILE UNTIL EOF OR END OF LC REGION (|n| < 1)
    while ( fscanf(fid, LCVIEW_TEXT_FORMAT_STRING,
                   &id, &n[0], &n[1], &n[2], &v, &S[0], &S[1]) != EOF )
    {
        //cout << id << " n " << n[0]<< "," << n[1]<<"," << n[2] << endl;

        // MAKE SURE DIRECTOR LENGTH ~= 1.
        // ALL ZERO DIRECTOR MEANS DIELECTRIC REGION, WHICH WE DONT WANT TO READ
        float dirLen = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
        if (dirLen < 0.95 )
            break;


        // CONVERT DIRECTOR TO Q-TENSOR
        double a1, a2, a3, a4, a5; // "NORMAL" Q-TENSOR COMPONENTS
        a1 = S[0]*(3*n[0]*n[0]-1) / 2.0; // Qxx
        a2 = S[0]*(3*n[1]*n[1]-1) / 2.0; // Qyy
        a3 = S[0]*(3*n[0]*n[1]) / 2.0;   // Qxy
        a4 = S[0]*(3*n[1]*n[2]) / 2.0;   // Qyz
        a5 = S[0]*(3*n[0]*n[2]) / 2.0;   // Qxz

        // CONVERT TO TRACELESS BASIS
        q1.push_back( 0.5*(a1+a2) *rt6 );
        q2.push_back( (a1+(a1+a2)/2.0)*rt2 );
        q3.push_back( a3*rt2 );
        q4.push_back( a4*rt2 );
        q5.push_back( a5*rt2 );
    }
    fclose(fid);

    if ( (idx) q1.size() != q.getnDoF() )
    {
        cout << "error, the loaded result file does not math the mesh size - bye!" << endl;
        exit(1);
    }

    for (idx i = 0 ; i < q.getnDoF() ; ++i)
    {
        q.setValue(i,0,q1[i]);
        q.setValue(i,1,q2[i]);
        q.setValue(i,2,q3[i]);
        q.setValue(i,3,q4[i]);
        q.setValue(i,4,q5[i]);
    }

}

void ReadLCD_B(Simu* simu, SolutionVector *q)
{
    // READS BINARY FORMATED RESULT FILE
    string filename = simu->getLoadQ();
    // printf( "Loading Q-tensor from: %s\n",filename.c_str());
    FILE* fid = fopen( filename.c_str() , "rb" );

    char str[100];
    if (fid)
    {
        float S0;
        idx np, nsol;

        // READS 5 LINES DISCARIDING DATA
        char *tempc = fgets (str , 100 , fid);
        tempc = fgets (str , 100 , fid);
        tempc = fgets (str , 100 , fid);
        tempc = fgets (str , 100 , fid);
        tempc = fgets (str , 100 , fid);

        fscanf (fid, "%f %i %i\n",&S0,&np,&nsol);

        if (np < q->getnDoF() ){
            cout << "error, the loaded result file does not math the mesh size - bye!" << endl;
            //printf("The file you're trying to load (\"%s\") doesn't seem to match the mesh - bye!\n",filename.c_str() );
            fclose(fid);
            exit(1);
        }

        float q1,q2,q3,q4,q5,temp;

        for (idx i = 0 ; i < q->getnDoF() ; i ++){
            fread ( &q1, sizeof(float), 1 , fid );
            fread ( &q2, sizeof(float), 1 , fid );
            fread ( &q3, sizeof(float), 1 , fid );
            fread ( &q5, sizeof(float), 1 , fid );
            fread ( &q4, sizeof(float), 1 , fid );
            q->setValue(i,0,q1);
            q->setValue(i,1,q2);
            q->setValue(i,2,q3);
            q->setValue(i,3,q4);
            q->setValue(i,4,q5);
            for (idx j = 0 ; j < nsol-5 ; j ++)	// READ&DISCARD POTENTIAL AND FLOW
                fread ( (void *) &temp, sizeof(float), 1 , fid );
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



void WriteLCViewResult(
        Simu* simu,
        LC* lc,
        Geometry* geom,
        SolutionVector* v,
        SolutionVector* q)

{
 //
 //  Writes LCView result file(s) to disk. TODO: get rid of this function.
 //
    size_t sf = simu->getSaveFormat();
    
    if (sf & Simu::LCview){
        WriteLCD_B(geom->getPtrTop(), geom->t, geom->e, v, q, simu, lc); // WRITES BINARY FILE
    }
    
    if (sf & Simu::LCviewTXT){
        WriteLCD(geom->getPtrTop(), geom->t, geom->e, v, q, simu); // WRITES TEXT FILE
    }
 }
// end void WriteREsult


void CreateSaveDir(Simu& simu)
{
    // CREATES DIRECTORY FOR RESULTS, IF IT DOES NOT
    // ALREADY EXIST


    // first check if savedit already exists
    if ( FilesysFun::dirExists( simu.getSaveDir() ) )
        return;             // if yes, can leave

    // try to create savedir.
    if ( FilesysFun::createDirectory(simu.getSaveDir() ) )
        return;             // if succesful, can leave
    else

        // if not succesful:
        cout << "error - cold not create SaveDir:\n"<< simu.getSaveDir() <<"\nbye!"<<endl;
    exit(1);

}
// end void CreateSaveDir
} // end namespace // WriteResults
