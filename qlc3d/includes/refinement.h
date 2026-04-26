#ifndef REFINEMENT_H
#define REFINEMENT_H


#include <line.h>
#include <meshrefinement.h>
#include <solutionvector.h>
#include <alignment.h>
#include <simu.h>
#include <vector>
#include <globals.h>
#include <geometry.h>
#include <refinement/refinement-spec.h>
#include <regulargrid.h>
#include <memory>
#define MIN(X,Y) ((X) < (Y) ? : (X) : (Y))


#define RED_TET	    6
#define GREEN3_TET  3
#define GREEN2_TET  2
#define GREEN1_TET  1


struct Num_Ref_Tet{
    // COUNTER STRUCT
    unsigned int red;
    unsigned int green3;
    unsigned int green2;
    unsigned int green1;
    Num_Ref_Tet(): red(0), green3(0), green2(0), green1(0){}
};



void Refine(Geometry& geom,             // refined geometry
            vector <unsigned int>& i_tet);  // index to refinable tets

/**
 * Create refined tetrahedra and triangles for a geometry using the supplied
 * refinement classifications and bisected edges.
 *
 * @param geom  Geometry to refine.
 * @param i_tet  Tetrahedron refinement types.
 * @param i_tri  Triangle refinement types.
 * @param lines  Bisected edges used to create midpoint nodes.
 * @param t_to_l  Tet-to-bisected-edge index mapping.
 * @param e_to_l  Triangle-to-bisected-edge index mapping.
 * @param new_p  Output coordinates for newly created midpoint nodes.
 * @param new_t  Output tetrahedron connectivity data.
 * @param new_mat_t  Output tetrahedron material numbers.
 * @param new_e  Output triangle connectivity data.
 * @param new_mat_e  Output triangle material numbers.
 */
void create_new_elements(   Geometry& geom,    // refinable geom
                            vector <idx>& i_tet,               // tet refinement types, 0=none, 1=green1 ...
                            vector <idx>& i_tri,               // tri refinement types, 0=none
                            const vector <Line>& lines,        // bisectable lines
                            const vector < set<idx> >& t_to_l, // index from tet to its bisectable lines
                            const vector < set<idx> >& e_to_l, // index from tri to its bisectable lines
                            vector <double>& new_p,            // return value, new coordinates created here
                            vector <idx>& new_t,               // return value, new tet elements created here
                            vector <idx>& new_mat_t,           // return value, new material numbers created here
                            vector <idx>& new_e,               // return value, new tria elements
                            vector <idx>& new_mat_e            // return value, new tri materials
                            );

bool autoref(Geometry& geom_orig,
             Geometry& geom_new,
             SolutionVector& q,
             SolutionVector& v,
             const std::vector<const RefinementSpec*>& specs,
             Simu& simu,
             SimulationState &simulationState,
             Alignment& alignment,
             const Electrodes& electrodes,
             double S0,
             std::unique_ptr<RegularGrid>& regGridOut);

// Checks for maximum dq within an element, as specified in meshrefinement.
// returns true is dq within any element is too large.
bool needsRefinement(Geometry& geom, SolutionVector& q, MeshRefinement& meshrefinement);

// returns true if End-Refinement is needed
bool needsEndRefinement(Geometry& geom,
                        SolutionVector& q,
                        MeshRefinement& meshrefinement);


// ELEMENT SEARCH FUNCTIONS - DEFINED IN findrefelems.cpp

// SELECTS ELEMENTS WHERE Q-TENSOR CHANGE IS ABOVE ALLOWED THRESHOLD,
// AS SPECIFIED IN REFINEMENTSPEC OBJECT
void findTets_Change(const RefinementSpec& refinfo,
                     vector <idx>& i_tet,
                     const int refiter,
                     const Geometry& geom,
                     const SolutionVector& q
                     );
// SELECTS ELEMENTS THAT CONTAIN NODES WITHIN A SPHERICAL REGION
// CENTRED AT X Y Z COORDINATES SPECIFIED IN spec
void findTets_Sphere(const RefinementSpec& refInfo,
                     vector<idx> &i_tet,
                     const int refIter,
                     const Geometry& geom);

void findTets_Box(const RefinementSpec& refInfo,
                  vector<idx> &i_tet,
                  const int refIter,
                  const Geometry& geom);


#endif
