
#ifndef QLC3D_H
#define QLC3D_H


#include <lc.h>
#include <box.h>
#include <alignment.h>
#include <simu.h>
#include <mesh.h>

#include <solutionvector.h>
#include <material_numbers.h>
#include <eventlist.h>
#include <solver-settings.h>
#include <geometry.h>
#include <energy.h>
#include <line.h>
#include <meshrefinement.h>
#include <eventlist.h>
#include <vector>
#include <list>
#include <iostream>

namespace SpaMtrix{
    class IRCMatrix;
    class Vector;
}
#ifndef PI
    #define PI 3.14159265358979323846264338327950288419716939937510
#endif

class Configuration;

const static double TIME_RESOLUTION = 1e-15;

// Assembles previous time step part of RHS when doing non-linear Crank-Nicholson
void assemble_prev_rhs(SpaMtrix::Vector &Ln,
		       SolutionVector& qn,
		       SolutionVector& v,
                       LC&mat_par,
                       double dt,
                       Geometry& geom);



void assembleQ(SpaMtrix::IRCMatrix &K,
           SpaMtrix::Vector &L,	// current RHS
           SolutionVector *q,
           SolutionVector* v,
           Mesh* t,
           Mesh* e,
           const Coordinates &coordinates,
           LC* mat_par,
           double dt,
           Alignment* alignment,
           const std::vector<Vec3> &nodeNormals);

// UPDATE Q-TENSOR USING IMPLICIT METHODS
double calcQ3d(SolutionVector *q,
               SolutionVector* qn,
               SolutionVector *v,
               Geometry& geom,
               LC* mat_par,
               Simu* simu,
               SimulationState &simulationState,
               SpaMtrix::IRCMatrix &Kq,
               SolverSettings* settings,
               Alignment* alignment);

class Electrodes;
class EventList;
void ReadSettings(const std::filesystem::path &settingsFilePath, Boxes &boxes, EventList &eventlist);

void CreateSaveDir(Simu* simu); //creates new save dir, if needed

/*! tensortovector returns the director and two order parameters in an array.
 *  The returned vector must be freed outside this function.
 *  TODO: make the return a shared_ptr
 *  @param a Q-tensor in traceless base
 *  @param npLC number of LC nodes (= 5x length of a)
 *  @return the vector representation of the LC. The first npLC elements are nx, then ny, nz, S1, S2.
 */
double* tensortovector(double *a, int npLC);

// THESE SHOULD BE DEFINED AS FRIEND FUNCTIONS TO CLASS SparseMatrix ?
//SparseMatrix* createSparseMatrix(const Geometry& geom, const SolutionVector& sol, const int& MatNum = 0); // sparse matrix for selected domain = LC only
//SparseMatrix* createSparseMatrix( vector<Line>& lines);
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
// TODO: delete this declaration whne sparsematrix.h is used everywhere instead
SpaMtrix::IRCMatrix createQMatrix(Geometry &geom,
                        SolutionVector &q,
                        const int& MatNum = MAT_DOMAIN1);

struct Qlc3dInfo {
    const std::string buildDate = __DATE__;
    const std::string buildTime = __TIME__;

    // this may not be reliable!
#ifdef NDEBUG
    const bool isDebug = false;
#else
    const bool isDebug = true;
#endif

#ifdef QLC3D_SHA
  const char* gitCommitSha = QLC3D_SHA;
#else
  const char* gitCommitSha = "unknown";
#endif
};

#endif

