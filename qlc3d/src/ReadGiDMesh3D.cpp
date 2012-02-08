#include <stdio.h>
#include <string>
#include <string.h>
#include <fstream>
#include <iostream>
#include <simu.h>
#include <geometry.h>
//using namespace std;


#define MAX_DIGITS 	1000000000
#define TOLER		1e-9



int CountNodes(ifstream* fin) // counts number of nodes
{
    string Coordinates = "Coordinates";
    string end_coordinates = "end coordinates";
    int np = 0;
    char *charray = (char*)malloc(200*sizeof(char));
    memset((void*)charray, 0 , 200*sizeof(char));
    // read until "Coordintes"
    while ( Coordinates.compare(0,Coordinates.size() , charray , 0 , Coordinates.size() ) != 0)
        fin->getline(charray,200,'\n');


    while ( end_coordinates.compare(0, end_coordinates.size() , charray , 0 , end_coordinates.size()  ) != 0 )
    {
        fin->getline(charray,200,'\n');
        np++;
    }
    np--;
    free(charray);
    return np;
} //end int CountNodes
int CountTetrahedra(ifstream* fin)// counts number of tetrahedral elements
{
    string end_elements = "end elements";
    string Elements 	= "Elements";


    int nt = 0;
    char *charray = (char*)malloc(200*sizeof(char));
    memset((void*)charray, 0 , 200*sizeof(char));

    // read file until Elements start
    while (Elements.compare(0, Elements.size() , charray , 0 , Elements.size()  ) != 0 )
        fin->getline(charray,200,'\n');

    while ( end_elements.compare(0 , end_elements.size() , charray , 0 , end_elements.size() ) )
    {
        fin->getline(charray,200,'\n');
        nt++;
    }
    nt--;
    free(charray);

    return nt;
}//end int CountTetrahedra
int CountTriangles(ifstream* fin)// counts number of triangles
{
    string Elements = "Elements";
    string end_elements = "end elements";

    char *charray = (char*)malloc(200*sizeof(char));
    memset((void*)charray, 0 , 200*sizeof(char));
    while ( Elements.compare(0 , Elements.size() , charray , 0 , Elements.size() ) != 0)
        fin->getline(charray,200,'\n');

    int ne = 0;
    while ( end_elements.compare(0 , end_elements.size() , charray , 0 , end_elements.size() ) != 0)
    {
        fin->getline(charray,200,'\n');
        ne++;
    }
    ne--;
    free(charray);
    return ne;
}//end int CountTriangles
int CountPrisms( ifstream* fin)// counts number of prisms (for periodic nodes)
{
    char* charray = (char*)malloc(200*sizeof(char));
    memset((void*)charray, 0 , 200*sizeof(char));
    string Elements = "Elements";
    string end_elements = "end elements";

    while (Elements.compare(0, Elements.size() , charray, 0, Elements.size() ) != 0 )
        fin->getline(charray,200,'\n');

    int npr = 0;
    while ( end_elements.compare(0 , end_elements.size() ,  charray , 0 ,end_elements.size() ) !=0 )
    {
        fin->getline(charray,200,'\n');
        npr++;
    }
    npr--;
    free(charray);
    return npr;
}//end int CountPrisms

void ReadNodes(ifstream* fin,int np, double* dp)
{
    char* charray = (char*)malloc(200*sizeof(char));
    memset((void*)charray, 0 , 200*sizeof(char));
    printf("\tReading %i nodes...",np); fflush(stdout);
    double  temp;

    string Coordinates = "Coordinates";

    while (Coordinates.compare(0,Coordinates.size() , charray , 0 , Coordinates.size() ) != 0 )
        fin->getline(charray,200,'\n');
    free(charray);
    for (int i = 0 ; i < np ; i++ )
    {
        *fin >> temp; //node number - not needed
        *fin >> dp[i*3+0];
        *fin >> dp[i*3+1];
        *fin >> dp[i*3+2];
        //if(i<10)
        //printf("[%f,%f,%f]\n",dp[i*3+0] , dp[i*3+1] , dp[i*3+2]);

    }

    cout << "OK\n";
}// end void ReadNodes
void ReadTetrahedra(ifstream* fin, idx nt, idx* dt, idx* dmatt)
{
    printf("\tReading %i tetrahedra...", nt); fflush(stdout);

    char* charray 	= (char*)malloc(200*sizeof(char));
    memset((void*)charray, 0 , 200*sizeof(char));
    string Elements = "Elements";

    while ( Elements.compare(0 , Elements.size() , charray, 0 , Elements.size() ) != 0)
        fin->getline(charray,200,'\n');	//find start of tet data
    free(charray);
    idx temp;

    for (idx i = 0 ; i < nt ; i++ )
    {
        *fin >> temp; //tetrahedra number - not needed
        *fin  >> dt[i*4+0];
        *fin  >> dt[i*4+1];
        *fin  >> dt[i*4+2];
        *fin  >> dt[i*4+3];
        *fin  >> dmatt[i];
    }

    printf("OK\n"); fflush(stdout);
}// end void ReadTetrahedra
void ReadPrisms(ifstream* fin, idx npr, idx* pr)
{
    printf("\tReading %i prisms (periodic boundaries)...",npr); fflush(stdout);
    char* charray = (char*)malloc(200*sizeof(char));
    memset((void*) charray , 0 , 22*sizeof(char) );
    idx temp;
    string Elements = "Elements";
    while ( Elements.compare(0 , Elements.size() , charray , 0 , Elements.size() ) != 0 )
        fin->getline(charray,200,'\n'); //find start of periodic nodes
    free(charray);

    for (idx i = 0 ; i < npr ; i++)
    {
        *fin >> temp;
        *fin >> pr[i*6+0];
        *fin >> pr[i*6+1];
        *fin >> pr[i*6+2];
        *fin >> pr[i*6+3];
        *fin >> pr[i*6+4];
        *fin >> pr[i*6+5];
    }
    printf("OK\n"); fflush(stdout);


}
void ReadTriangles(ifstream* fin, idx ne, idx* e, idx* emat)
{
    printf("\tReading %i triangles...",ne); fflush(stdout);
    char* charray 	= (char*)malloc(200*sizeof(char));
    memset((void*)charray, 0 , 200*sizeof(char));
    idx temp;

    string Elements = "Elements";
    while (Elements.compare(0,Elements.size() , charray,0,Elements.size() )!= 0)
        fin->getline(charray,200,'\n');		//find start of triangle data
    free(charray);


    for (idx i =0 ; i<ne ; i++)
    {
        *fin >> temp; //triangle number - not needed
        *fin >> e[i*3+0];
        *fin >> e[i*3+1];
        *fin >> e[i*3+2];

        if (fin->peek() == 10)
            emat[i] = 0; // if newline character ( = no material number assigned), convert to 0
        else
            *fin >> emat[i];
    }
    printf("OK\n"); fflush(stdout);

}

void ReadGiDMesh3D(Simu* simu,double **p, idx *np, idx **t, idx *nt,idx **e,
                   idx *ne, idx **matt, idx **mate)
{

    using namespace std;
    idx nperi = 0;
    ifstream fin;
    char *charray 	= (char*)malloc(200*sizeof(char));

    string filename = simu->MeshName;//"testcube.msh";
    cout <<"Attempting to open mesh file: " << filename  << endl;
    fin.open(filename.c_str() );//.c_str(),ifstream::in);

    if (!fin.good() )
    {
        cout << "failed opening: " << filename.c_str() << endl;
        exit(1);
    }

    else{
        printf("Reading GID mesh file: %s \n", filename.c_str()); fflush(stdout);
        //printf("a");
        np[0] = 0;
        nt[0] = 0;
        ne[0] = 0;
        bool tets_first = false;

        string Tets 	= "Nnode 4"; //MESH    dimension 3 ElemType Tetrahedra  Nnode 4";
        string Prisms  	= "Nnode 6";//MESH    dimension 3 ElemType Prisma  Nnode 6";
        string Tris 	= "Nnode 3";//MESH    dimension 3 ElemType Triangle  Nnode 3";



        while(!fin.eof()){ // count  - while loop
            fin.getline(charray,200);
            string line = charray;

            //if ( Tets.compare(0 , Tets.size() , charray, 0 , Tets.size() ) == 0 ){ // IF TETS ARE FIRST IN MESH FILE
            if (line.find(Tets) != std::string::npos ){
                if (np[0] == 0){	// if tets are defined first
                    np[0] = CountNodes(&fin);

                    tets_first = true;
                }

                nt[0] = CountTetrahedra(&fin);
            }// end if tets before prisms

            //else if ( Prisms.compare(0 , Prisms.size() , charray , 0 , Prisms.size() ) == 0 ){ // IF PERIODIC NODES
            else if (line.find(Prisms) != std::string::npos){
                // if prisms are defined before nodes
                if (np[0] == 0)
                    np[0] = CountNodes(&fin);
                nperi = CountPrisms(&fin);
            }

            //else if (  Tris.compare(0 , Tris.size() , charray , 0 , Tris.size() ) == 0 ){ // IF TRIANGLES
            else if (line.find(Tris) != std::string::npos){
                if (np[0] == 0)
                    np[0] = CountNodes(&fin);

                ne[0] = CountTriangles(&fin);
            }
        } // end count while loop

        printf("Tets = %u , Tris = %u , nodes = %u\n", nt[0], ne[0] , np[0]);


        //ALLOCATE MEMORY FOR MESH DATA ---- use of double pointers
        // ERROR CHECKING BLOCK
        {
            bool error = false;
            if (*np == 0 )
            {
                printf("\nerror in reading GiD mesh.\nYour mesh seems to be missing points. Bye!\n");
                error = true;
            }
            if (*ne == 0)
            {
                printf("\nerror in reading GiD mesh.\nYour mesh seems to be missing triangles. Bye!\n");
                error = true;
            }
            if (*nt==0)
            {
                printf("\nerror in reading GiD mesh.\nYour mesh seems to be missing tetrahedra. Bye!\n");
                error = true;
            }

            if (error)
            {
                fflush(stdout);
                exit(1);
            }
        }   // END ERROR CHECKING BLOCK

	double *dp 		= (double*) malloc(3*np[0]*sizeof(double));

        idx *dt			= (idx*)  malloc(4*nt[0]*sizeof(idx));
        idx *de 		= (idx*)  malloc(3*ne[0]*sizeof(idx));

        idx *dmatt 		= (idx*)malloc(nt[0]*sizeof(idx));
        idx *dmate 		= (idx*)malloc(ne[0]*sizeof(idx));

        idx* pr=NULL;
	if (nperi>0)
            pr 	= (idx*)malloc(nperi*6*sizeof(idx) );


        //REWIND TO START OF FILE
	fin.clear();              // forget we hit the end of file
	fin.seekg(0, ios::beg);   // move to the start of the file


	while(!fin.eof())	// read - while loop
	{
            fin.getline(charray,200,'\n');
            string line = charray;
            //if ( Tets.compare(0 , Tets.size() , charray , 0 , Tets.size()  ) == 0)
            if ( line.find(Tets) != std::string::npos ){
                if (tets_first)	// if tets are defined first
                    ReadNodes(&fin,np[0],dp);

                ReadTetrahedra(&fin,nt[0],dt,dmatt);
            }
            //else if ( Prisms.compare(0, Prisms.size() , charray , 0 , Prisms.size() )== 0 )
            else if ( line.find(Prisms) != std::string::npos ){
                if (!tets_first) // if prisms are defined first
                    ReadNodes(&fin,np[0],dp);

                ReadPrisms(&fin,nperi,pr);
            }
            //else if ( Tris.compare(0 , Tris.size() , charray , 0 , Tris.size() )== 0 ) // triangles
            else if ( line.find(Tris) != std::string::npos ){
                ReadTriangles(&fin,ne[0],de,dmate);
            }
	} // end read while loop



        // ONLY POSITIVE COORDINATES ALLOWED - MOVE IF NECESSARY
        double xmin,ymin,zmin,xmax,ymax,zmax;
	xmin = 1e9 ; ymin = 1e9 ; zmin = 1e9;
	xmax = -1e9; ymax = -1e9; zmax = -1e9;

        for (idx i=0 ; i<np[0] ; i++)
	{
            if(dp[i*3+0]<xmin) xmin = dp[i*3+0]; // find min
            if(dp[i*3+1]<ymin) ymin = dp[i*3+1];
            if(dp[i*3+2]<zmin) zmin = dp[i*3+2];

            if(dp[i*3+0]>xmax) xmax = dp[i*3+0]; // find max
            if(dp[i*3+1]>ymax) ymax = dp[i*3+1];
            if(dp[i*3+2]>zmax) zmax = dp[i*3+2];
	}

        for (idx i=0;i<np[0];i++)
	{
            if(xmin<0) dp[i*3+0]-=xmin; // move mesh
            if(ymin<0) dp[i*3+1]-=ymin;
            if(zmin<0) dp[i*3+2]-=zmin;
	}


        if (xmin<0) {xmax-=xmin;  cout << "\tshifting all x by :"<<-xmin<<endl;xmin=0;}
        if (ymin<0) {ymax-=ymin;  cout << "\tshifting all y by :"<<-ymin<<endl;ymin=0;}
        if (zmin<0) {zmax-=zmin;  cout << "\tshifting all z by :"<<-zmin<<endl;zmin=0;}

        //*
        // REMOVE NUMERICAL NOISE FROM BOUNDARY NODES
        for (idx i=0;i<np[0];i++)
	{
            if( dp[i*3+0]-xmin <= TOLER) dp[i*3+0]=xmin;
            if( xmax-dp[i*3+0] <= TOLER) dp[i*3+0]=xmax;

            if( dp[i*3+1]-ymin <= TOLER) dp[i*3+1]=ymin;
            if( ymax-dp[3*i+1] <= TOLER) dp[i*3+1]=ymax;

            if( dp[i*3+2]-zmin <= TOLER) dp[i*3+2]=zmin;
            if( zmax-dp[3*i+2] <= TOLER) dp[i*3+2]=zmax;
	}
        //*/


	if (pr!=NULL)
            free(pr);
	*p=dp;
	*t=dt;
	*e=de;
	*matt = dmatt;
	*mate = dmate;
	free(charray);
        /*
 */
	//exit(0);
    }//end if mesh file opened ok


}

void writeBinaryMesh(Simu& simu, Geometry& geom){
    /*! Writes binary representation of geometry geom to file*/

    // MAke output filename.
    std::string filename = simu.MeshName;
    size_t ind = filename.find_last_of('.');
    filename.erase(ind+1);
    filename.append("geo");
    printf("binary output filename = %s\n", filename.c_str()); fflush(stdout);

    fstream file ( filename.c_str() , ios::out | ios::binary);
    if (!file.good() ){
        printf("error - could not write %s\n", filename.c_str()); fflush(stdout);
        exit(1);
    }

    //file << (unsigned int) geom.getnp() << (unsigned int) geom.t->getnElements() << (unsigned int) geom.e->getnElements();
    unsigned int numbers[3] = { geom.getnp() , geom.t->getnElements() , geom.e->getnElements() };
    file.write( (char*) numbers , 3*sizeof(unsigned int) );
    // write coordinates
    file.write( (char*) geom.getPtrTop(), 3* geom.getnp() * sizeof( double ) );
    // write tets node indexes
    file.write( (char*) geom.t->getPtrToElement(0) , 4*numbers[1]*sizeof(int) );
    // write tris node indexes
    file.write( (char*) geom.e->getPtrToElement(0) , 3*numbers[2]*sizeof(int) );
    // write tets materials
    file.write( (char*) geom.t->getPtrToMaterialNumber(0), numbers[1]*sizeof(int) );
    // write tris materials
    file.write( (char*) geom.e->getPtrToMaterialNumber(0), numbers[2]*sizeof(int) );
    file.close();

    //exit(1);
}

void readBinaryMesh(std::string filename,
                    double *&p,
                    idx *&t, idx *&tmat,
                    idx *&e, idx *&emat,
                    idx *np, idx *nt, idx *ne){

    fstream file( filename.c_str() , ios::in | ios::binary);
    if ( !file.good() ){
        printf("error - could not read %s\n", filename.c_str() ); fflush(stdout);
    }
    //unsigned int np, nt, ne;

    file.read( (char*) np , sizeof(int) );
    file.read( (char*) nt , sizeof(int) );
    file.read( (char*) ne , sizeof(int) );

    printf(" np : %i\n nt : %i\n ne : %i", np[0], nt[0],ne[0]);

    p = (double*)  malloc( 3**np*sizeof(double) );
    t = (idx*)     malloc( 4**nt*sizeof(idx) );
    tmat = (idx*)  malloc( *nt*sizeof(idx) );
    e    = (idx*)  malloc( 3**ne*sizeof(idx) );
    emat = (idx*)  malloc( *ne*sizeof(idx) );

    if (!p)     {printf("error - could not allocate p\n"); exit(1);}
    if (!e)     {printf("error - could not allocate e\n"); exit(1);}
    if (!t)     {printf("error - could not allocate t\n"); exit(1);}
    if (!emat)  {printf("error - could not allocate emat\n"); exit(1);}
    if (!tmat)  {printf("error - could not allocate tmat\n"); exit(1);}

    file.read((char*) p , 3**np*sizeof(double) );
    file.read((char*) t , 4**nt*sizeof(int) );
    file.read((char*) e , 3**ne*sizeof(int) );
    file.read((char*) tmat, *nt*sizeof(int) );
    file.read((char*) emat, *ne*sizeof(int) );

    file.close();

    printf("\nfile read OK\n"); fflush(stdout);

}
