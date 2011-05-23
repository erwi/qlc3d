
#ifndef QLC3D_H
#define QLC3D_H

#include <electrodes.h>
#include <lc.h>
#include <box.h>
#include <alignment.h>
#include <simu.h>
#include <mesh.h>
//#include "nodes.h"
#include <solutionvector.h>
#include <material_numbers.h>
#include <sparsematrix.h>
#include <eventlist.h>
#include <settings.h>
#include <geometry.h>
#include <energy.h>
#include <line.h>
#include <meshrefinement.h>
//#include <meshrefinement.h"
//#include "line.h"

#include <vector>
#include <list>
#include <iostream>
using std::vector;
using std::string;

//#define N_THREADS 4
#define PI 3.14159
#define eps0 8.8541878176e-12

// struct spm{				// define sparse matrix link
//		int row;
//		spm *next;
//		spm *prev;
//	};

//typedef struct llnode spm;

void ReadGiDMesh3D(Simu *simu,double **p, int *np, int **t, int *nt,int **e, int *ne, int **matt, int **mate);


//void writeBinaryMesh(Simu& simu, Geometry& geom); // writes mesh in binary format. in definition in "ReadGiDMesh3D.cpp"
void readBinaryMesh(std::string filename ,  // same as above
                    double *&p,
                    int *&t, int *&tmat,
                    int *&e, int *&emat,
                    int *np, int *nt, int *ne);
void calcpot3d(SparseMatrix* Kpot,SolutionVector *v, SolutionVector *q,LC* lc, Mesh *mesh, Mesh* surf_mesh, double *p, Settings* settings, Electrodes* electrodes);

void solve_pcg(SparseMatrix *K, double *b, double *x ,Settings* settings);
void solve_gmres(SparseMatrix *K, double *b, double *x ,Settings* settings);

// Assembles previous time step part of RHS when doing non-linear Crank-Nicholson
void assemble_prev_rhs(double* Ln,
		       SolutionVector& qn,
		       SolutionVector& v,
		       Mesh& t,
		       Mesh& e,
		       double* p ,
		       LC&mat_par,
		       Simu& simu);

void assembleQ(SparseMatrix* K,
	       double* L,	// current RHS
           SolutionVector *q,
           SolutionVector* v,
	       Mesh* t,
	       Mesh* e,
	       double* p,
	       LC* mat_par,
	       Simu* simu,
	       Settings* settings,
	       Alignment* alignment,
	       double* NodeNormals);

double calcQ3d(SolutionVector *q,
	       SolutionVector* qn,
		   SolutionVector *v,
	       Mesh *t,
	       Mesh *e,
	       double *p,
	       LC* mat_par,
	       Simu* simu,
	       SparseMatrix* Kq,
           Settings* settings,
	       Alignment* alignment,
	       double* NodeNormals);


//void WriteLCD(double *p, Mesh *t, Mesh *e, SolutionVector *v, SolutionVector *q,Simu* simu);	// writes result as text file
//void WriteLCD_B(double *p, Mesh *t, Mesh *e, SolutionVector *v, SolutionVector *q,Simu* simu, LC* lc); // writes result file in binary format
void WriteResult(
	Simu* simu, 		// Simulation settings
	LC* lc,				// LC material paramters
	Geometry* geom,		// mesh geometry data 
	SolutionVector* v,  // potential solution
    SolutionVector* q,  // Q-tensor solution
    MeshRefinement* meshref = NULL); // meshrefinement info. including whether a new mesh has been generated

void ReadLCD_B(Simu* simu, SolutionVector* q);

void ReadSettings(
        string settings_filename,
        Simu* simu,
        LC* lc,
        Boxes* boxes,
        Alignment* alignment,
	Electrodes* electrodes,
	MeshRefinement* meshrefinement
	);

void ReadSolverSettings(const char* filename, Settings* settings);
void WriteSettings(Simu* simu, LC* lc, Boxes* box, Alignment* alignment , Electrodes* electrodes);
void CreateSaveDir(Simu* simu); //creates new save dir, if needed

// -----------------------------
//
// INITIALISATION FUCTIONS
//
// -----------------------------

void prepareGeometry(Geometry& geom, Simu& simu);// defined in inits.cpp
FILE* createOutputEnergyFile(Simu& simu); // defined in inits.cpp


void SetVolumeQ(SolutionVector *q, LC* lc, Boxes* boxes, double* p);

// Sets all Q-tensor values to those specified in Alignment
void setSurfacesQ(
		SolutionVector* q,
		Alignment* alignment,
		LC* lc, 
		Geometry* geom);
		
// Sets Q-tensor values only for fixed nodes (strong, homeotropic anchoring), as specified in Alignment
// defined in file : setsurfacesq.cpp
void setStrongSurfacesQ(
		SolutionVector* q,
		Alignment* alignment,
		LC* lc, 
		Geometry* geom);
		
		
int ReorderNodes(double *p, int np,int *t, int nt, int *e, int ne,int *tmat,int *emat);
double* tensortovector(double *a, int npLC);
void tensorToEigs(double* a,		// input Q-tensor, traceless basis
				  double* eigVal,	// result eigenvalues
				  double* eigVec    // result eigenvectors
				  );
// THESE SHOULD BE DEFINED AS FRIEND FUNCTIONS TO CLASS SparseMatrix ?
//SparseMatrix* createSparseMatrix(Geometry* geom, SolutionVector* u); // sparse matrix for all domains
SparseMatrix* createSparseMatrix(Geometry& geom, SolutionVector& sol, const int& MatNum = 0); // sparse matrix for selected domain = LC only
//SparseMatrix* createSparseMatrix(Mesh* m);
SparseMatrix* createSparseMatrix( vector<Line>& lines);
#endif

