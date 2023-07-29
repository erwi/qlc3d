#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <filesystem>

#include <simu.h>
#include <geometry.h>
#include <util/exception.h>
#include <util/logging.h>
#include <geom/vec3.h>

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
void readNodes(ifstream* fin, idx np, std::vector<Vec3> &points) {
    forwardToLine(fin,"coordinates");

    for (idx i = 0 ; i < np ; i++ ) {
        double temp, x, y, z;
        *fin >> temp; //node number - not needed
        *fin >> x;
        *fin >> y;
        *fin >> z;
        points.emplace_back(x, y, z);
    }
}// end void ReadNodes

void readTetrahedra(ifstream* fin, idx numTets, std::vector<idx> &tetNodes, std::vector<idx> &tetMaterials) {
    Log::info("Reading {} tetrahedra", numTets);
    forwardToLine(fin, "elements");
    char cbuff[256];

    tetNodes.resize(numTets * 4, NOT_AN_INDEX);
    tetMaterials.resize(numTets, NOT_AN_INDEX);
    for (idx i = 0; i < numTets; i++) {
        fin->getline(cbuff,256);
        std::stringstream ss(cbuff);

        idx tmp=0;
        bool err = !( ss >> tmp >> tetNodes[i * 4] >> tetNodes[i * 4 + 1] >> tetNodes[i * 4 + 2] >> tetNodes[i * 4 + 3] >> tetMaterials[i]);

        if (err) {
            fin->close();
            std::string errorMsg = fmt::format("Bad file format reading tetrahedra. Expected 6 numbers per line in "
                                               "the format <tet number> <node 1> <node 2> <node 3> <node 4> <tet material>"
                                               ", but found {} instead.", cbuff);
            RUNTIME_ERROR(errorMsg);
         }

        // GiD mesh files use 1-based indexing, but we want 0-based, so decrement each new index by one.
        tetNodes[i * 4]--;
        tetNodes[i * 4 + 1]--;
        tetNodes[i * 4 + 2]--;
        tetNodes[i * 4 + 3]--;
    }
}

void readTriangles(ifstream* fin, idx numTris, std::vector<idx> &triNodes, std::vector<idx> &triMaterials) {
    Log::info("Reading {} triangles", numTris);
    forwardToLine(fin, "elements");
    triNodes.clear();
    triMaterials.clear();
    triNodes.resize(numTris * 3, NOT_AN_INDEX);
    triMaterials.resize(numTris, NOT_AN_INDEX);

    for (idx i =0 ; i < numTris ; i++) {
        idx temp;
        *fin >> temp; //triangle number - not needed
        *fin >> triNodes[i * 3 + 0];
        *fin >> triNodes[i * 3 + 1];
        *fin >> triNodes[i * 3 + 2];

        // convert from 1-based to 0-based indexing
        triNodes[i * 3 + 0]--;
        triNodes[i * 3 + 1]--;
        triNodes[i * 3 + 2]--;
        if (fin->peek() == 10) { // 10 is ASCII FOR NEW LINE
            triMaterials[i] = 0; // if newline character ( = no material number assigned), convert to 0
        } else {
            *fin >> triMaterials[i];
        }
    }
}

void ReadGiDMesh3D(const std::filesystem::path &meshFileName,
                   std::vector<Vec3> &pointsOut,
                   std::vector<idx> &tetNodes,
                   std::vector<idx> &tetMaterials,
                   std::vector<idx> &triNodes,
                    std::vector<idx> &triMaterials) {
    using namespace std;
    if (!filesystem::exists(meshFileName)) {
        RUNTIME_ERROR(fmt::format("Can not read mesh from path {}, it does not exist.", meshFileName));
    }
    idx nperi = 0;
    idx numPoints = 0;
    ifstream fin;
    char charray[200];

    Log::info("Attempting to open mesh file: {}", meshFileName.string());
    fin.open(meshFileName.string().c_str());
    if (!fin.good() ) {
        RUNTIME_ERROR("could not open mesh file " + meshFileName.string());
    }
    else { // File opened OK
        Log::info("Reading GID mesh file: {}", meshFileName.string());

        idx numTets = 0;
        idx numTris = 0;
        bool tets_first = false;

        string Tets 	= "Nnode 4";//MESH    dimension 3 ElemType Tetrahedra  Nnode 4";
        string Prisms  	= "Nnode 6";//MESH    dimension 3 ElemType Prisma  Nnode 6";
        string Tris 	= "Nnode 3";//MESH    dimension 3 ElemType Triangle  Nnode 3";

        while(!fin.eof()) { // count  - while loop
            fin.getline(charray,200);
            string line = charray;
            if (line.find(Tets) != std::string::npos ) {
                if (numPoints == 0) {	// if tets are defined first
                    numPoints = CountNodes(&fin);
                    tets_first = true;
                }
                numTets = CountTetrahedra(&fin);
            }// end if tets before prisms


            else if (line.find(Prisms) != std::string::npos) {
                // if prisms are defined before nodes
                if (numPoints== 0)
                    numPoints = CountNodes(&fin);
                nperi = CountPrisms(&fin);
            }
            else if (line.find(Tris) != std::string::npos) {
                if (numPoints == 0) {
                    numPoints = CountNodes(&fin);
                }
                numTris = CountTriangles(&fin);
            }
        } // end count while loop

        if (fin.eof()) {
            bool error = false;
            if (numPoints == 0) {
                Log::error("Could not find coordinates.");
                error = true;
            }
            if (numTris == 0) {
                Log::error("Could not find triangle elements.");
                error = true;
            }
            if (numTets == 0) {
                Log::error("Could not find tetrahedral elements.");
                error = true;
            }
            if (error) {
                RUNTIME_ERROR("Could not find all required data in mesh file " + meshFileName.string());
            }
        }

        Log::info("Number of tetrahedra = {}, triangles = {}, prisms = {}, coordinates = {}", numTets, numTris, nperi, numPoints);
        pointsOut.clear();
        pointsOut.reserve(numPoints);

        //REWIND TO START OF FILE
        fin.clear();              // forget we hit the end of file
        fin.seekg(0, ios::beg);   // move to the start of the file

        while(!fin.eof()) {	// read - while loop
            fin.getline(charray,200,'\n');
            string line = charray;

            if ( line.find(Tets) != std::string::npos ) {
                if (tets_first)	// if tets are defined first
                  readNodes(&fin, numPoints, pointsOut);

              readTetrahedra(&fin, numTets, tetNodes, tetMaterials);
            }
            else if ( line.find(Prisms) != std::string::npos ) {
                if (!tets_first) // if prisms are defined first
                  readNodes(&fin, numPoints, pointsOut);
            }
            else if ( line.find(Tris) != std::string::npos ) {
              readTriangles(&fin, numTris, triNodes, triMaterials);
            }
        } // end read while loop

        // ONLY POSITIVE COORDINATES ALLOWED - MOVE IF NECESSARY
        // TODO: is this really required??
        double xmin = pointsOut[0].x();
        double ymin = pointsOut[0].y();
        double zmin = pointsOut[0].z();

        for (idx i=0; i < numPoints; i++) {
          xmin = std::min(xmin, pointsOut[i].x());
          ymin = std::min(ymin, pointsOut[i].y());
          zmin = std::min(zmin, pointsOut[i].z());
        }

        Vec3 delta = {xmin < 0 ? -xmin : 0,
                      ymin < 0 ? -ymin : 0,
                      zmin < 0 ? -zmin : 0};

        if (delta.norm() > 0) {
          Log::info("Shifting all coordinates by {}.", delta);
          for (Vec3 &p: pointsOut) {
            p += delta;
          }
        }
     }
}
