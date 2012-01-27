#ifndef REFINEMENT_H
#define REFINEMENT_H


#include <line.h>
#include <meshrefinement.h>
#include <solutionvector.h>
#include <alignment.h>
#include <lc.h>
#include <simu.h>
#include <box.h>
#include <vector>
#define MIN(X,Y) ((X) < (Y) ? : (X) : (Y))


// this is probably bad... but makes things som much easier in tet-splitting
	#define nA	no[0] // shortcuts to old nodes
	#define nB	no[1]
	#define nC	no[2]
	#define nD	no[3]

	#define nAB nn[0] // shortcut to new nodes
	#define nAC nn[1]
	#define nAD nn[2]
	#define nBC nn[3]
	#define nBD nn[4]
	#define nCD nn[5]


#define RED_TET	    6
#define GREEN3_TET  3
#define GREEN2_TET  2
#define GREEN1_TET  1


struct Num_Ref_Tet{
    unsigned int red;
    unsigned int green3;
    unsigned int green2;
    unsigned int green1;
};

struct peri_lines{
	// faces lines
	vector <Line> lfront;
	vector <Line> lback;
	vector <Line> lleft;
	vector <Line> lright;
	vector <Line> ltop;
	vector <Line> lbottom;
	// corners lines
	vector <Line> lc0;
	vector <Line> lc1;
	vector <Line> lc2;
	vector <Line> lc3;
	
	vector <Line> lca;
	vector <Line> lcb;
	vector <Line> lcc;
	vector <Line> lcd;
	
	vector <Line> lcA;
	vector <Line> lcB;
	vector <Line> lcC;
	vector <Line> lcD;
};


void Refine(Geometry& geom_orig, // original GiD(?) mesh
	    Geometry& geom_prev, // mesh from previous iteration, used for interpolating results from
	    Geometry& geom_new,  // new geometry created in here
	    MeshRefinement* = NULL);

void Refine(Geometry& srce,					// source geometry
			Geometry& dest,					// destination geometry
			vector <unsigned int>& i_tet);	// index to refinable tets

void create_new_elements(   Geometry& geom,		    // refinable geom
				vector <unsigned int>& i_tet,		// tet refinement types, 0=none, 1=green1 ...
				vector <unsigned int>& i_tri,		// tri refinement types, 0=none
				vector <Line>& lines,				// bisectable lines
				vector < set<unsigned int> > t_to_l,// index from tet to its bisectable lines
				vector < set<unsigned int> > e_to_l,// index from tri to its bisectable lines
				vector <double>& new_p,				// return value, new coordinates created here
				vector <unsigned int>& new_t,       // return value, new tet elements created here
				vector <int>& new_mat_t,			// return value, new material numbers created here
				vector <unsigned int>& new_e,		// return value, new tria elements
				vector <int>& new_mat_e             // return value, new tri materials
			 );

void autoref(Geometry& geom_orig, Geometry& geom_prev, Geometry& geom_new,
             SolutionVector& q, SolutionVector& qn,
             SolutionVector& v,
             //MeshRefinement& meshrefinement,
             const list<RefInfo>& refInfos,
             Simu& simu, Alignment& alignment,  Electrodes& electrodes, LC& lc);

// Checks for maximum dq within an element, as specified in meshrefinement.
// returns true is dq within any element is too large.
bool needsRefinement(Geometry& geom, SolutionVector& q, MeshRefinement& meshrefinement);

// returns true if End-Refinement is needed
bool needsEndRefinement(Geometry& geom,
                       SolutionVector& q,
                       MeshRefinement& meshrefinement);


#endif
