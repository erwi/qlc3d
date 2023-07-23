#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <filesystem>

#include <simu.h>
#include <geometry.h>
#include <util/exception.h>
#include <util/logging.h>

#define MAX_DIGITS 	1000000000
#define TOLER		1e-9

void strToLower(std::string& str)
{
    // CONVERTS STRING TO ALL LOWER CASE LETTERS
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);

}
void charsToLower(char* str) {
    int i = -1;
    while(true) {
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

idx forwardToLine(ifstream* fin, const char* key) {
    // FINDS LINE IN FILE fin CONTAINING TEXT key AND RETURNS NUMBER OF LINES READ
    // MAKE SURE key IS ALL LOWER CASE!!
    char charray[256];
    idx c = 0;
    while (true) {
        fin->getline( charray, 256,'\n');
        charsToLower(charray);

        // CONVERT READ LINE TO STRING AND SEARCH FOR SUBSTRINGS
        // THIS FIXES PROBLEM WITH DIFFERENCES BETWEEN WINDOWS/LINUX
        // END OF LINE CHARACTERS (Win = '\r\n', Linux = '\n')
        std::string sline = charray;
        if ( sline.find(key) != std::string::npos ) // IF LINE CONTAINS SUBSTRING.
            break;

        if ( fin->eof() ) {
            RUNTIME_ERROR(fmt::format("unexpected end of file when reading {}", key));
        }
        c++;
    }

    return c;
}

idx CountNodes(ifstream* fin) {
    idx l1 = forwardToLine(fin, "coordinates" );
    idx l2 = forwardToLine(fin, "end coordinates");
    return l2-l1;
}

idx CountTetrahedra(ifstream* fin) {
    forwardToLine(fin, "elements");
    idx l2 = forwardToLine(fin, "end elements");
    return l2;// - l1;
}

idx CountTriangles(ifstream* fin) {
    forwardToLine(fin, "elements");
    idx l2 = forwardToLine(fin, "end elements");
    return l2;
}

idx CountPrisms( ifstream* fin) {
    forwardToLine(fin, "elements");
    idx l2 = forwardToLine(fin, "end elements");
    return l2;
}

/**
 * Reads coordinate values into dp array.
 */
void ReadNodes(ifstream* fin,idx np, double* dp) {
    forwardToLine(fin,"coordinates");

    for (idx i = 0 ; i < np ; i++ ) {
        double temp;
        *fin >> temp; //node number - not needed
        *fin >> dp[i*3+0];
        *fin >> dp[i*3+1];
        *fin >> dp[i*3+2];
    }
}// end void ReadNodes

void ReadTetrahedra(ifstream* fin, idx nt, idx* dt, idx* dmatt) {
    Log::info("Reading {} tetrahedra", nt);
    forwardToLine(fin, "elements");
    char cbuff[256];
    for (idx i = 0; i < nt; i++) {
        fin->getline(cbuff,256);
        std::stringstream ss(cbuff);

        idx tmp=0;
        bool err = !( ss >> tmp >> dt[i*4] >> dt[i*4+1] >> dt[i*4+2] >> dt[i*4+3] >> dmatt[i]);

        if (err) {
            fin->close();
            std::string errorMsg = fmt::format("Bad file format reading tetrahedra. Expected 6 numbers per line in "
                                               "the format <tet number> <node 1> <node 2> <node 3> <node 4> <tet material>"
                                               ", but found {} instead.", cbuff);
            RUNTIME_ERROR(errorMsg);
         }

        // GiD mesh files use 1-based indexing, but we want 0-based, so decrement each new index by one.
        dt[i * 4]--;
        dt[i * 4 + 1]--;
        dt[i * 4 + 2]--;
        dt[i * 4 + 3]--;
    }
}

void ReadPrisms(ifstream* fin, idx npr, idx* pr) {
    Log::info("Reading {} prisms.", npr);
    forwardToLine(fin, "elements");
    for (idx i = 0; i < npr; i++) {
        idx temp;
        *fin >> temp;
        *fin >> pr[i*6+0];
        *fin >> pr[i*6+1];
        *fin >> pr[i*6+2];
        *fin >> pr[i*6+3];
        *fin >> pr[i*6+4];
        *fin >> pr[i*6+5];
    }
}

void ReadTriangles(ifstream* fin, idx ne, idx* e, idx* emat) {
    Log::info("Reading {} triangles", ne);
    forwardToLine(fin, "elements");
    for (idx i =0 ; i<ne ; i++) {
        idx temp;
        *fin >> temp; //triangle number - not needed
        *fin >> e[i*3+0];
        *fin >> e[i*3+1];
        *fin >> e[i*3+2];

        // convert from 1-based to 0-based indexing
        e[i * 3 + 0]--;
        e[i * 3 + 1]--;
        e[i * 3 + 2]--;
        if (fin->peek() == 10) { // 10 is ASCII FOR NEW LINE
            emat[i] = 0; // if newline character ( = no material number assigned), convert to 0
        } else {
            *fin >> emat[i];
        }
    }
}

void ReadGiDMesh3D(const std::filesystem::path &meshFileName,
                   double **p,
                   idx *np,
                   idx **t,
                   idx *nt,
                   idx **e,
                   idx *ne,
                   idx **matt,
                   idx **mate) {
    using namespace std;
    if (!filesystem::exists(meshFileName)) {
        RUNTIME_ERROR(fmt::format("Can not read mesh from path {}, it does not exist.", meshFileName));
    }
    idx nperi = 0;
    ifstream fin;
    char *charray 	= (char*)malloc(200*sizeof(char));

    Log::info("Attempting to open mesh file: {}", meshFileName.string());
    fin.open(meshFileName.string().c_str());
    if (!fin.good() ) {
        RUNTIME_ERROR("could not open mesh file " + meshFileName.string());
    }
    else { // File opened OK
        Log::info("Reading GID mesh file: {}", meshFileName.string());

        np[0] = 0;
        nt[0] = 0;
        ne[0] = 0;
        bool tets_first = false;

        string Tets 	= "Nnode 4";//MESH    dimension 3 ElemType Tetrahedra  Nnode 4";
        string Prisms  	= "Nnode 6";//MESH    dimension 3 ElemType Prisma  Nnode 6";
        string Tris 	= "Nnode 3";//MESH    dimension 3 ElemType Triangle  Nnode 3";

        while(!fin.eof()){ // count  - while loop
            fin.getline(charray,200);
            string line = charray;
            if (line.find(Tets) != std::string::npos ) {
                if (np[0] == 0) {	// if tets are defined first
                    np[0] = CountNodes(&fin);
                    tets_first = true;
                }
                nt[0] = CountTetrahedra(&fin);
            }// end if tets before prisms


            else if (line.find(Prisms) != std::string::npos) {
                // if prisms are defined before nodes
                if (np[0] == 0)
                    np[0] = CountNodes(&fin);
                nperi = CountPrisms(&fin);
            }
            else if (line.find(Tris) != std::string::npos) {
                if (np[0] == 0){
                    np[0] = CountNodes(&fin);
                }
                ne[0] = CountTriangles(&fin);
            }
        } // end count while loop

        if (fin.eof()) {
            bool error = false;
            if (*np == 0) {
                Log::error("Could not find coordinates.");
                error = true;
            }
            if (*ne == 0) {
                Log::error("Could not find triangle elements.");
                error = true;
            }
            if (*nt == 0) {
                Log::error("Could not find tetrahedral elements.");
                error = true;
            }
            if (error) {
                RUNTIME_ERROR("Could not find all required data in mesh file " + meshFileName.string());
            }
        }

        Log::info("Number of tetrahedra = {}, triangles = {}, prisms = {}, coordinates = {}",
                  nt[0], ne[0], nperi, np[0]);
        double *dp 		= (double*) malloc(3*np[0]*sizeof(double));
        idx *dt			= (idx*)  malloc(4*nt[0]*sizeof(idx));
        idx *de 		= (idx*)  malloc(3*ne[0]*sizeof(idx));
        idx *dmatt 		= (idx*)malloc(nt[0]*sizeof(idx));
        idx *dmate 		= (idx*)malloc(ne[0]*sizeof(idx));
        idx* pr=NULL;

        if (nperi > 0) {
            pr 	= (idx*) malloc(nperi*6*sizeof(idx) );
        }

        //REWIND TO START OF FILE
        fin.clear();              // forget we hit the end of file
        fin.seekg(0, ios::beg);   // move to the start of the file

        while(!fin.eof()){	// read - while loop
            fin.getline(charray,200,'\n');
            string line = charray;

            if ( line.find(Tets) != std::string::npos ) {
                if (tets_first)	// if tets are defined first
                    ReadNodes(&fin,np[0],dp);

                ReadTetrahedra(&fin,nt[0],dt,dmatt);
            }
            else if ( line.find(Prisms) != std::string::npos ) {
                if (!tets_first) // if prisms are defined first
                    ReadNodes(&fin,np[0],dp);

                ReadPrisms(&fin,nperi,pr);
            }
            else if ( line.find(Tris) != std::string::npos ) {
                ReadTriangles(&fin,ne[0],de,dmate);
            }
        } // end read while loop

        // ONLY POSITIVE COORDINATES ALLOWED - MOVE IF NECESSARY
        // TODO: is this really required??
        double xmin,ymin,zmin,xmax,ymax,zmax;
        xmin = 1e9 ; ymin = 1e9 ; zmin = 1e9;
        xmax = -1e9; ymax = -1e9; zmax = -1e9;

        for (idx i=0 ; i<np[0] ; i++) {
            if(dp[i*3+0]<xmin) xmin = dp[i*3+0]; // find min
            if(dp[i*3+1]<ymin) ymin = dp[i*3+1];
            if(dp[i*3+2]<zmin) zmin = dp[i*3+2];

            if(dp[i*3+0]>xmax) xmax = dp[i*3+0]; // find max
            if(dp[i*3+1]>ymax) ymax = dp[i*3+1];
            if(dp[i*3+2]>zmax) zmax = dp[i*3+2];
        }

        for (idx i=0; i<np[0]; i++) {
            if(xmin < 0) dp[i * 3 + 0] -= xmin; // move mesh
            if(ymin < 0) dp[i * 3 + 1] -= ymin;
            if(zmin < 0) dp[i * 3 + 2] -= zmin;
        }

        if (xmin < 0) {
            Log::info("Shifting all x by {}.", -xmin);
            xmax -= xmin;
            xmin = 0;
        }

        if (ymin < 0) {
            Log::info("Shifting all y by {}.", -ymin);
            ymax -= ymin;
            ymin=0;
        }

        if (zmin < 0) {
            Log::info("Shifting all z by {}.", -zmin);
            zmax-=zmin;
            zmin=0;
        }

        if (pr != nullptr) {
            free(pr);
        }
        *p=dp;
        *t=dt;
        *e=de;
        *matt = dmatt;
        *mate = dmate;
        free(charray);
     }
}
