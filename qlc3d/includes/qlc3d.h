
#ifndef QLC3D_H
#define QLC3D_H

//#include <electrodes.h>
#include <lc.h>
#include <box.h>
#include <alignment.h>
#include <simu.h>
#include <mesh.h>

#include <solutionvector.h>
#include <material_numbers.h>
//#include <sparsematrix.h>
#include <eventlist.h>
#include <settings.h>
#include <geometry.h>
#include <energy.h>
#include <line.h>
#include <meshrefinement.h>
#include <eventlist.h>
#include <vector>
#include <list>
#include <iostream>

// SPAMTRIX INCLUDES
//#include <spamtrix_ircmatrix.hpp>

using std::vector;
using std::string;
namespace SpaMtrix{
    class IRCMatrix;
    class Vector;
}
#ifndef PI
    #define PI 3.14159265358979323846264338327950288419716939937510
#endif

#define eps0 8.8541878176e-12
#define COMPLEX std::complex<double>


void ReadGiDMesh3D(Simu *simu,
                   double **p,
                   idx *np,
                   idx **t,
                   idx *nt,
                   idx **e,
                   idx *ne,
                   idx **matt,
                   idx **mate);


//void writeBinaryMesh(Simu& simu, Geometry& geom); // writes mesh in binary format. in definition in "ReadGiDMesh3D.cpp"
void readBinaryMesh(std::string filename ,  // same as above
                    double *&p,
                    idx *&t, idx *&tmat,
                    idx *&e, idx *&emat,
                    idx *np, idx *nt, idx *ne);

// Assembles previous time step part of RHS when doing non-linear Crank-Nicholson
void assemble_prev_rhs(SpaMtrix::Vector &Ln,
		       SolutionVector& qn,
		       SolutionVector& v,
                       LC&mat_par,
                       Simu& simu,
                       Geometry& geom);



void assembleQ(SpaMtrix::IRCMatrix &K,
           SpaMtrix::Vector &L,	// current RHS
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

// UPDATE Q-TENSOR USING IMPLICIT METHODS
double calcQ3d(SolutionVector *q,
               SolutionVector* qn,
               SolutionVector *v,
               Geometry& geom,
               LC* mat_par,
               Simu* simu,
               SpaMtrix::IRCMatrix &Kq,
               Settings* settings,
               Alignment* alignment);//
               //double* NodeNormals);

class SparseMatrix;
double calcQExplicit(SolutionVector &q,
                     SolutionVector &v,
                     SparseMatrix &K,
                     Geometry &geom,
                     LC &lc,
                     Alignment &alignment,
                     Simu &simu,
                     Settings &settings
                     );
class Electrodes;
class EventList;
void ReadSettings(
        string settings_filename,
        Simu* simu,
        LC& lc,
        Boxes* boxes,
        Alignment* alignment,
        Electrodes* electrodes,
        MeshRefinement* meshrefinement,
        EventList& eventlist
	);

void ReadSolverSettings(const char* filename, Settings* settings);
void CreateSaveDir(Simu* simu); //creates new save dir, if needed

// -----------------------------
//
// INITIALISATION FUCTIONS
//
// -----------------------------

void prepareGeometry(Geometry& geom,    // defined in inits.cpp
                     Simu& simu,
                     Alignment& ali);
FILE* createOutputEnergyFile(Simu& simu); // defined in inits.cpp
void selectQMatrixSolver(Simu &simu, const LC &lc); // defined in ints.cpp



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
SparseMatrix* createSparseMatrix(const Geometry& geom, const SolutionVector& sol, const int& MatNum = 0); // sparse matrix for selected domain = LC only
SparseMatrix* createSparseMatrix( vector<Line>& lines);
// CREATES SPARSE MATRIX OF CORRECT TYPE FOR SOLUTION OF Q-TENSOR
// DEPENDING ON SOLVER METHOD USED (IMPLICIT / EXPLICIT)
//SparseMatrix* createSparseMatrixQ(const Geometry &geom,
//                                 const SolutionVector &q,
//                                 const Settings &set);





// CREATES SPAMTRIX SPARSE MATRIC FOR POTENTIAL
SpaMtrix::IRCMatrix createPotentialMatrix(Geometry &geom,
                              SolutionVector &sol,
                              const int &MatNum,
                              const Electrodes &electrodes);
// CREATES SPAMTRIX SPARSE MATRIX FOR Q-TENSOR
SpaMtrix::IRCMatrix createQMatrix(Geometry &geom,
                        SolutionVector &q,
                        const int& MatNum = MAT_DOMAIN1);


#endif

