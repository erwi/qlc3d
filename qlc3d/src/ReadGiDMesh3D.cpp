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

void strToLower(std::string& str)
{
    // CONVERTS STRING TO ALL LOWER CASE LETTERS
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);

}
void charsToLower(char* str)
{
    int i = -1;
    while(true)
    {
        i++;
        char a = str[i];

        if (a == '\0')
            break;

        if (a < 65 )    // IF NUMBER OR OTHER SPECIAL CHARACTER
        {
            continue;
        }

        str[i] = tolower( a );

    }
}

void EOFError(const char* where)
{
    std::cout << "error reading mesh file.\nUnexpected end of file encountered when reading: "<< where<<std::endl;
    std::cout << "bye!" << std::endl;
    exit(1);
}


idx forwardToLine(ifstream* fin, const char* key)
{
    // FINDS LINE IN FILE fin CONTAINING TEXT key AND RETURNS NUMBER OF LINES READ
    // MAKE SURE key IS ALL LOWER CASE!!
    char charray[256];
    idx c = 0;
    while (true)
    {
        fin->getline( charray, 256,'\n');
        charsToLower(charray);

        // CONVERT READ LINE TO STRING AND SEARCH FOR SUBSTRINGS
        // THIS FIXES PROBLEM WITH DIFFERENCES BETWEEN WINDOWS/LINUX
        // END OF LINE CHARACTERS (Win = '\r\n', Linux = '\n')
        std::string sline = charray;
        if ( sline.find(key) != std::string::npos ) // IF LINE CONTAINS SUBSTRING.
            break;

        if ( fin->eof() )
            EOFError( key );

        c++;
    }

    return c;
}

idx CountNodes(ifstream* fin) // counts number of nodes
{
    idx l1 = forwardToLine(fin, "coordinates" );
    idx l2 = forwardToLine(fin, "end coordinates");

    return l2-l1;
} //end int CountNodes
idx CountTetrahedra(ifstream* fin)// counts number of tetrahedral elements
{
    forwardToLine(fin, "elements");
    idx l2 = forwardToLine(fin, "end elements");
    return l2;// - l1;
}//end int CountTetrahedra

idx CountTriangles(ifstream* fin)// counts number of triangles
{
    forwardToLine(fin, "elements");
    idx l2 = forwardToLine(fin, "end elements");

    return l2;
}//end int CountTriangles
idx CountPrisms( ifstream* fin)// counts number of prisms (for periodic nodes)
{

    forwardToLine(fin, "elements");
    idx l2 = forwardToLine(fin, "end elements");
    return l2;

}//end int CountPrisms

void ReadNodes(ifstream* fin,idx np, double* dp)
{


    forwardToLine(fin,"coordinates");

    for (idx i = 0 ; i < np ; i++ )
    {
        double temp;
        *fin >> temp; //node number - not needed
        *fin >> dp[i*3+0];
        *fin >> dp[i*3+1];
        *fin >> dp[i*3+2];
    }

    cout << "OK\n";
}// end void ReadNodes
void ReadTetrahedra(ifstream* fin, idx nt, idx* dt, idx* dmatt)
{
    printf("\tReading %i tetrahedra...", nt); fflush(stdout);
    forwardToLine(fin, "elements");
    for (idx i = 0 ; i < nt ; i++ )
    {
        idx temp;
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
    /*
    char* charray = (char*)malloc(200*sizeof(char));
    memset((void*) charray , 0 , 22*sizeof(char) );
    idx temp;
    string Elements = "Elements";
    while ( Elements.compare(0 , Elements.size() , charray , 0 , Elements.size() ) != 0 )
        fin->getline(charray,200,'\n'); //find start of periodic nodes
    free(charray);
    */
    forwardToLine(fin, "elements");
    for (idx i = 0 ; i < npr ; i++)
    {
        idx temp;
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

    forwardToLine(fin, "elements");
    for (idx i =0 ; i<ne ; i++)
    {
        idx temp;
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
    else // FILE OPENED OK
    {
        printf("Reading GID mesh file: %s \n", filename.c_str()); fflush(stdout);
        //printf("a");
        np[0] = 0;
        nt[0] = 0;
        ne[0] = 0;
        bool tets_first = false;

        string Tets 	= "Nnode 4"; //MESH    dimension 3 ElemType Tetrahedra  Nnode 4";
        string Prisms  	= "Nnode 6";//MESH    dimension 3 ElemType Prisma  Nnode 6";
        string Tris 	= "Nnode 3";//MESH    dimension 3 ElemType Triangle  Nnode 3";



        while(!fin.eof()) // count  - while loop
        {
            fin.getline(charray,200);
            string line = charray;


            if (line.find(Tets) != std::string::npos )
            {
                if (np[0] == 0)	// if tets are defined first
                {
                    np[0] = CountNodes(&fin);
                    tets_first = true;
                }

                nt[0] = CountTetrahedra(&fin);
            }// end if tets before prisms


            else if (line.find(Prisms) != std::string::npos)
            {
                // if prisms are defined before nodes
                if (np[0] == 0)
                    np[0] = CountNodes(&fin);
                nperi = CountPrisms(&fin);
            }


            else if (line.find(Tris) != std::string::npos)
            {
                if (np[0] == 0)
                    np[0] = CountNodes(&fin);

                ne[0] = CountTriangles(&fin);
            }
        } // end count while loop

        if( fin.eof() )
        {
            bool error = false;
            if ( *np == 0 )
            {
                printf("error, could not find coordinates - bye!\n");
                error = true;
            }
            if (*ne == 0)
            {
                printf("error, could not find triangle elements - bye\n");
                error = true;
            }
            if (*nt == 0)
            {
                printf("error, could not find tetrhedral elements - bye\n");
                error = true;
            }
            if (error)
                exit(1);

        }


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

        /*
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
