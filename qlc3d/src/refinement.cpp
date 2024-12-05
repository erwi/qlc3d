
#include <refinement.h>
#include <material_numbers.h>
#include <meshrefinement.h>
#include <line.h>
#include <set>
#include <globals.h>
#include <geom/coordinates.h>
#include <util/exception.h>
#include <util/logging.h>
using std::set;

void find_all_core_lines( vector<Line>& lines, const vector <unsigned int>& i_tet, Mesh* m){
    // creates unique set of lines from all bisectable lines of red tets
    // only works for tets!
    vector <unsigned int>::const_iterator itr; // iterator to red tet indexes

    for (unsigned int i = 0 ; i < i_tet.size() ; i++){ // loop over all tets
	if  ( i_tet[i] == 6 ){ // if this is a red tet
        idx tet[4] = {  m->getNode(i, 0), // shortcut to tet nodes
			    m->getNode(i, 1),
			    m->getNode(i, 2),
			    m->getNode(i, 3)};

	    // EACH TET CONTAINS 6 LINES - INSERT THESE
	    lines.push_back( Line( tet[0] , tet[1] ) ); // 0->1
	    lines.push_back( Line( tet[0] , tet[2] ) ); // 0->2
	    lines.push_back( Line( tet[0] , tet[3] ) ); // 0->3
	    lines.push_back( Line( tet[1] , tet[2] ) ); // 1->2
	    lines.push_back( Line( tet[1] , tet[3] ) ); // 1->3
	    lines.push_back( Line( tet[2] , tet[3] ) ); // 2->3
	}//end if re tet
    }

    // MAKE LINES LIST UNIQUE

    sort( lines.begin(), lines.end() ); // SORT
    vector <Line>::iterator u;
    u = unique( lines.begin() , lines.end() );                // REORDER
    lines.resize( u - lines.begin() );  // RESIZE
}
//*/
void count_lines(   Mesh* m,                                // mesh elements
                    vector< set < idx > >& p_to_m ,         // point to elements indices
                    const vector<Line>& lines ,             // list of lines
                    vector < idx >& i_elem,                 // counter array
                    vector < set < idx> >& m_to_l           // elements to lines indices
                    ){
    // COUNT NUMBER OF OCCURRENCIES OF LINES IN ELEMENTS.
    // LOOP OVER EACH LINE. USE p_to_m INDICES TO FIND ELEMENTS
    // THAT CONTAIN *BOTH* NODES OF THE LINE.

    i_elem.clear();
    i_elem.reserve( m->getnElements() );
    i_elem.assign( m->getnElements() , 0 );

    // ALLOCATE MEMORY FOR m_to_l
    m_to_l.clear();
    m_to_l.reserve( m->getnElements() );
    set <unsigned int> empty;
    m_to_l.assign( m->getnElements() , empty );

    // pragma omp loopable?
    int num_points = (int) p_to_m.size();
    for ( idx l = 0 ; l< (idx) lines.size() ; l++)// for each line
    {
        int n1 = lines[l].L[0]; // FIRST AND SECOND NODE NOBERS OF LINE l
        int n2 = lines[l].L[1];

        // IF ELEMENT TYPE IS TRIANGLE, NOT ALL NODES HAVE CONNECTIONS TO A TRIANGLE
        // CHECK THAT NODENUMBERS ARE VALID, THESE
        // SHOULD ALWAYS BE FOR TETS, BUT NOT NOT ALWAYS FOR TRIS
        if ( (n1 < num_points ) && ( n2 < num_points ) )
        {
            // INDEXES TO ALL TETS CONNECTED TO FIRST  AND SECOND NODES OF THIS LINE
            set <unsigned int> elems1;
            set <unsigned int> elems2;

            elems1.insert( p_to_m[ n1 ].begin() , p_to_m[ n1 ].end() );
            elems2.insert( p_to_m[ n2 ].begin() , p_to_m[ n2 ].end() );

            vector <unsigned int> elems_res;

            set_intersection( elems1.begin() , elems1.end() , elems2.begin() , elems2.end() ,back_inserter(elems_res) );
            // elems_res now contains list of all elements that share this line
            // increase counter for those elements that contain lines in 'line;
            for (int i = 0; i < (int) elems_res.size() ; i++) // loop over all elements sharing the line
            {
                // increment element line count
                i_elem[ elems_res[i] ]++;
                // adds index from element to this line
                m_to_l[ elems_res[i] ].insert( l ); // can be made faster by inserting last using iterators
            }// end for
        } // end check for valid nodenumber
    }// end for all lines
}// end count lines
void find_tet_refinement_types(vector < idx >& i_tet, // tet line counts
                               Num_Ref_Tet& nrt )
{ // tet to lines index
    nrt.green1 = 0; // RESET
    nrt.green2 = 0;
    nrt.green3 = 0;
    nrt.red    = 0;

    for (idx i = 0 ; i < (idx) i_tet.size() ; i++) // loop over tet bisection line count
    {
        if ( i_tet[i]==0 ){} // do nothing
        else if ( i_tet[i] ==GREEN1_TET )
        {
            nrt.green1++;
        }
        else if ( i_tet[i] == GREEN2_TET )
        {
            nrt.green2++;
        }
        else if ( i_tet[i] == GREEN3_TET )
        {
            nrt.green3++;
        }
        else if ( i_tet[i] > GREEN3_TET)
        {
            nrt.red ++;
            i_tet[i] = RED_TET ; // enforces anything above and including 4 is actually 6
        }
    }// end for
}


void fix_green3_red_confusions( vector <idx>& i_tet,        // tets line counter
                                vector<Line>& lines,        // bisected lines
                                vector< set<idx> > t_to_l)  // index from tets to  lines
{
    // SOME ELEMENTS MARKED AS GREEN 3 ARE IN FACT RED TETS
    // GREEN 3 CONTAIN 3 BISECTABLE LINES FORMING A CLOSED LOOP -> 3 UNIQUE NODES
    // IN RED TETS, SOME OTHER CONFIGURATION OF 3 LINES EXIST -> 4 UNIQUE NODES

    // LOOP OVER EACH TET-LINE COUNTER PRAGMA OMP'ABLE
    //#pragma omp parallel for
    for (idx i = 0 ; i < i_tet.size() ; i++) // for all tets
    {
        // IF PREVIOUSLY MARKED AS A GREEN 3 TET
        if ( i_tet[i] == 3)
        {
            // LOOP OVER ALL ITS BISECABLE LINES AND CHECK
            // NUMBER OF UNIQUE NODES
            set < idx > ::iterator itr;
            set < idx > nodes; // nodes set, size of this determines type
            for (itr = t_to_l[i].begin() ; itr != t_to_l[i].end() ; ++itr) // for index to lines
            {
                nodes.insert( lines[*itr].L[0]); // insert both nodes of line
                nodes.insert( lines[*itr].L[1]);
            }// end for index to lines

            if (nodes.size() == 4 )
            {
                // THIS ELEMENT IS ACTUALLY RED. MAKE IT SO!
                i_tet[i] = 6; // 6 is magic for red
            }

            else if ( nodes.size() != 3) {
                RUNTIME_ERROR(fmt::format("Expected 3 nodes, got {}.", nodes.size()));
            }

        }
    }// end for all tets
}// end fix_green3_red_confusions



void gen_peri_lines_vec(peri_lines& plines,
                        Geometry& geom)
{/*! creates vector of lines that are on periodic surfaces of geometry geom.
 The lines are unique, i.e. no repetitions. */

    plines.lfront.clear(); plines.lback.clear();
    plines.lleft.clear();  plines.lright.clear();
    plines.lc0.clear(); plines.lc1.clear(); plines.lc2.clear(); plines.lc3.clear();
    plines.lca.clear(); plines.lcb.clear(); plines.lcc.clear(); plines.lcd.clear();
    plines.lcA.clear(); plines.lcB.clear(); plines.lcC.clear(); plines.lcD.clear();

    // this is a silly shortcut to save typing on detrmining which
    // surfaces are periodic
    bool peritype[3] = {geom.getfront_back_is_periodic() ,
                        geom.getleft_right_is_periodic() ,
                        geom.gettop_bottom_is_periodic() };


    // loop over each surface triangle
    for (idx i = 0 ; i < geom.getTriangles().getnElements() ; i++ ) {
        // if triangle is periodic
        if ( geom.getTriangles().getMaterialNumber(i) == MAT_PERIODIC ) {
          idx n[3];
          geom.getTriangles().loadNodes(i, n);
          //idx* n = geom.e->getPtrToElement( i ); // pointer to element node numbers
            Line lines[3] = {Line(n[0], n[1]), Line(n[0], n[2]) , Line(n[1], n[2]) };

            for ( idx l = 0 ; l < 3 ; l++)
            { // loop over each line and determine its location
                // Consider 3 cases of periodicity:
                // 1. FRONT/BACK ONLY
                // 2. FRONT/BACK and LEFT/RIGHT
                // 3. FRONY/BACK, LEFT/RIGHT and TOP/BOTTOM

                //CASE 1, front back only
                if ( peritype[0] && (!peritype[1]) && (!peritype[2]) )
                {
                    if ( lines[l].isOnFrontSurface( &geom ) )
                        plines.lfront.push_back( lines[l] );
                    else
                        if ( lines[l].isOnBackSurface( &geom ) )
                            plines.lback.push_back( lines[l] );
                }
                else
                    // CASE 2, front/back, left/right
                    if ( (peritype[0] ) && (peritype[1]) && (!peritype[2] ) )
                    {
                        // each line may be in multiple lists when along a corner
                        if ( lines[l].isOnFrontSurface( &geom ) )
                            plines.lfront.push_back( lines[l] );
                        else
                            if ( lines[l].isOnBackSurface( &geom) )
                                plines.lback.push_back( lines[l] );
			
                        if ( lines[l].isOnLeftSurface( &geom ) )
                            plines.lleft.push_back( lines[l] );
                        else
                            if (lines[l].isOnRightSurface(&geom) )
                                plines.lright.push_back( lines[l] );

                        // corner lines
                        if ( lines[l].isCorn0( &geom ) )
                            plines.lc0.push_back( lines[l] );
                        else
                            if ( lines[l].isCorn1(&geom) )
                                plines.lc1.push_back( lines[l] );
                            else
                                if ( lines[l].isCorn2(&geom ) )
                                {
                                    plines.lc2.push_back( lines[l] );
                                    //lines[l].PrintLine();
                                }
                                else
                                    if ( lines[l].isCorn3(&geom) )
                                        plines.lc3.push_back( lines[l] );


			
                    } // end if peritype
            }// end for l
        }// end if periodic triangle
    }// end of each triangle i

    // remove repeated lines. function defined as inline in line.h
    uniquefy_line_vector( plines.lleft );
    uniquefy_line_vector( plines.lright );
    uniquefy_line_vector( plines.lfront );
    uniquefy_line_vector( plines.lback );

    uniquefy_line_vector( plines.lc0 );
    uniquefy_line_vector( plines.lc1 );
    uniquefy_line_vector( plines.lc2 );
    uniquefy_line_vector( plines.lc3 );

}

bool find_transl_line( Line& l1,
                       vector<Line>& lines,
                       Geometry& geom,
                       double* dir,
                       vector<Line>:: iterator& other
                       )
{
// compares line l1 with lines in vector lines and checks whether it is
// a tranlation described by dir. sets other to found line and returns true
// if successfull

    // loop over all possible lines and compare
    for (other = lines.begin() ; other!= lines.end() ; other++)
    {
        if ( l1.isTranslationOf(*other, &geom, dir ) )
            return true;
    }
    return false;
}

void expand_periodic_boundaries(vector <Line>& lines, // lines to split
                                peri_lines& plines,  // lines on periodic boundary
                                Geometry& geom )
{
    /*! This function checks each line in vector lines for its periodic
  *  equivalencies and adds them to the list of bisectable lines.
  *  This is necessary to maintain periodicity of a mesh when it is
  *  refined near a periodic boundary.
  *
  *  plines is a convenience structure of vectors of lines along periodic surfaces
  *  and corners.
  */

    // minimum requirement for periodicity is that at least front/back
    // surfaces are periodic. if not, can return
    if ( !geom.getfront_back_is_periodic() )
    {
        return;
    }
    vector <Line> newlines;
    vector <Line> :: iterator litr;
    vector <Line> :: iterator lotr;
    //vector <Line> :: iterator loend; // end iterator
    for (litr = lines.begin() ; litr != lines.end() ; ++litr ) // loop over each line
    {

        if ( geom.getfront_back_is_periodic() )
        {


            bool found = true;
            double dir[3] = {1,0,1};// compare X and Z coorinates (shift in Y is allowed)

            if ( litr->isOnFrontSurface(&geom) )
            {
                found = false;
                if ( find_transl_line( *litr, plines.lback, geom, dir, lotr) )
                {
                    newlines.push_back(*lotr);
                    found = true;
                }
            }
            else
                if ( litr->isOnBackSurface(&geom ) )
                {
                    found = false;
                    if ( find_transl_line( *litr, plines.lfront, geom, dir, lotr) )
                    {
                        newlines.push_back(*lotr);
                        found = true;
                    }
                }

            if (!found) {
                RUNTIME_ERROR("Could not find periodic front/back face line.")
            }

        }// end if front/back

        if (geom.getleft_right_is_periodic() )
        {
            double dir[3] = {0,1,1}; // compare Y,Z. Allow X-shift
            bool found = true;
            if ( litr->isOnLeftSurface(&geom) )
            {
                found = false;
                if ( find_transl_line(*litr, plines.lright, geom, dir, lotr ) )
                {
                    newlines.push_back(*lotr);
                    found = true;
                }
            }
            else
                if (litr->isOnRightSurface(&geom) )
                {
                    found = false;
                    if ( find_transl_line(*litr, plines.lleft, geom, dir, lotr) )
                    {
                        newlines.push_back(*lotr);
                        found = true;
                    }
                }
            if (!found) {
                RUNTIME_ERROR("Could not find periodic left/right face line.")
            }


            // in case of l/r periodicity, also vertical corners must be
            // taken into account. For each corner line, must find 3 other
            // corner lines
            if ( litr->isTopBottomCornerLine(&geom) )
            {
                // Search all 4 possible corners and make sure exactly
                // 3 corresponding lines are found
                int cc = 0; //corner count
                double dir[3] = {0,0,1}; // make sure Z matches, allow X and Y shifts
                if ( find_transl_line(*litr, plines.lc0, geom, dir, lotr ) )
                {
                    newlines.push_back(*lotr);
                    cc++;
                }
                if ( find_transl_line(*litr, plines.lc1, geom, dir, lotr ) )
                {
                    newlines.push_back(*lotr);
                    cc++;
                }
                if ( find_transl_line(*litr, plines.lc2, geom, dir, lotr ) )
                {
                    newlines.push_back(*lotr);
                    cc++;
                }
                if ( find_transl_line(*litr, plines.lc3, geom, dir, lotr ) )
                {
                    newlines.push_back(*lotr);
                    cc++;
                }

                if ( cc != 3) {
                    RUNTIME_ERROR("Problem finding periodic top-bottom corner lines.")
                }
            }// end if vertical corner
        }// end if left/right

        if (geom.gettop_bottom_is_periodic() )
        {

            bool found = true;
            double dir[3] = {1,1,0};// compare X and Y, shift in Z

            if ( litr->isOnTopSurface(&geom) )
            {
                found = false;
                if ( find_transl_line(*litr, plines.lbottom, geom, dir, lotr) )
                {
                    newlines.push_back(*lotr);
                    found = true;
                }
            }
            else
                if ( litr->isOnBottomSurface( &geom ) )
                {
                    found = false;
                    if ( find_transl_line(*litr, plines.ltop, geom , dir, lotr) )
                    {
                        newlines.push_back(*lotr);
                        found = true;
                    }
                }
            if (!found) {
                RUNTIME_ERROR("Could not find periodic top/bottom line.")
            }

            // in case of top/bottom periodicity (=fully periodic),
            // horizontal corners must be taken into account. Two cases exist:
            // horizontal along X and horizontal along Y

            // corners A,B,C,D
            if (litr->isFrontBackCornerLine( &geom )) {
                // same strategy as for vertical corners c0,c1,c2,c3 above
                int cc = 0; //corner count
                double dir[3] = {0,1,0}; // make sure Y matches, allow X and Z shifts
                if (find_transl_line(*litr, plines.lcA, geom, dir, lotr )) {
                    newlines.push_back(*lotr);
                    cc++;
                }
                if (find_transl_line(*litr, plines.lcB, geom, dir, lotr )) {
                    newlines.push_back(*lotr);
                    cc++;
                }
                if (find_transl_line(*litr, plines.lcC, geom, dir, lotr )) {
                    newlines.push_back(*lotr);
                    cc++;
                }
                if (find_transl_line(*litr, plines.lcD, geom, dir, lotr )) {
                    newlines.push_back(*lotr);
                    cc++;
                }
                if (cc != 3) {
                    RUNTIME_ERROR("Problem finding periodic front-back corner lines.")
                }

            }// end horizontal corner along y
            else if (litr->isLeftRightCornerLine( &geom )) {
                int cc = 0; //corner count
                double dir[3] = {1,0,0}; // make sure X matches, allow Y and Z shifts
                if ( find_transl_line(*litr, plines.lca, geom, dir, lotr )) {
                    newlines.push_back(*lotr);
                    cc++;
                }
                if ( find_transl_line(*litr, plines.lcb, geom, dir, lotr )) {
                    newlines.push_back(*lotr);
                    cc++;
                }
                if ( find_transl_line(*litr, plines.lcc, geom, dir, lotr )) {
                    newlines.push_back(*lotr);
                    cc++;
                }
                if ( find_transl_line(*litr, plines.lcd, geom, dir, lotr )) {
                    newlines.push_back(*lotr);
                    cc++;
                }
                if ( cc != 3) {
                    RUNTIME_ERROR("Problem finding periodic left-right corner lines.")
                }
            }// end horizontal corner along X
            RUNTIME_ERROR("Top-bottom periodicity not implemented in mesh refinement yet.");
        }// end if top/bottom surface is periodic
    }//end for i loop over each line

    // Add newly found periodic lines to bisectable lines and remove repetitions
    lines.insert(lines.end() , newlines.begin() , newlines.end() );
    uniquefy_line_vector( lines );
}

void expand_refinement_region(vector <unsigned int>& i_tet,	// index to tet bisectable lines counts
			      Num_Ref_Tet& nrt,
			      vector <Line>& lines,
			      Geometry& geom_prev,		// geometry
			      vector < set <unsigned int> >& t_to_l        // tets to lines index
			      ){	// index to core nodes
    // REFINEMENT REGION MAY NEED TO BE EXPANDED IN ORDER TO MAINTAIN
    // PERIODIC BOUNDARY CONDITIONS OR IN SOME CASES TO
    // AVOID SPECIALLY AWKWARD CONCAVE REFINEMENT REGIONS
    if ( nrt.red == 0 )
	return; // exit if nothing to do

    unsigned int n_tred_old = 0;        // NUMBER OF RED TETS FROM PREVIOUS EXPANSION LOOP
    unsigned int d_n_tred   = nrt.red;	// CHANGE IN NUMBER OF RED TETS PER EXPANSION LOOP

    vector < set <unsigned int> > p_to_t;   // POINTS TO TETS INDEX
    geom_prev.getTetrahedra().gen_p_to_elem( p_to_t );
    
    peri_lines plines;
    // Construct list of periodic line elements iff structure has periodic boundaries
    if (  geom_prev.getleft_right_is_periodic()  ||
          geom_prev.getfront_back_is_periodic()  ||
          geom_prev.gettop_bottom_is_periodic()  )
    {
        gen_peri_lines_vec( plines,	geom_prev);
    }
    
    while( d_n_tred > 0) // loop while refinement region grows
    {

        find_all_core_lines(lines, i_tet, &geom_prev.getTetrahedra()); // creates all bisectable lines
        expand_periodic_boundaries( lines , plines, geom_prev );

        count_lines(&geom_prev.getTetrahedra(), p_to_t, lines, i_tet , t_to_l);

        fix_green3_red_confusions(i_tet, lines, t_to_l);
        n_tred_old = nrt.red;

        find_tet_refinement_types( i_tet, nrt);


        // CHECK WHETHER NEW ELEMENTS HAVE BEEN ADDED
        d_n_tred = nrt.red - n_tred_old;

    }// end while refinement region grows
}// end expand_refinement_region





void find_triangle_reftypes(Geometry& geom,
                            vector <unsigned int>& i_tri,
                            vector <Line>& lines,
                            vector < set <unsigned int> >& e_to_l
                            ){

    // MAKE LINKS FROM NODES TO TRIS
    vector < set <unsigned int> > p_to_e;

    geom.getTriangles().gen_p_to_elem( p_to_e );
    // MAKE INDEXES FROM TRIS TO LINES

    count_lines(&geom.getTriangles(), p_to_e, lines, i_tri, e_to_l);

}

void modify_geometry(Geometry& geom,
                     vector <idx>& i_tet, vector <idx>& i_tri,
                     vector <double>& new_p,
                     vector <idx>& new_t, vector <idx>& new_e,
                     vector <idx> new_mat_t ,
                     vector <idx> new_mat_e
                     ){

    if ( i_tet.empty() ) return;

    geom.addCoordinates( new_p ); // APPENDS NEW COORDINATE DATA
    // MAKE LIST OF TETS AND TRIS TO BE REMOVED. THESE ARE REPLACED BY THE NEWLY
    // CREATED ONES
    set <idx> ind_remove_tets;
    for (idx i = 0 ; i < (idx) i_tet.size() ; i++)
    {
        if ( i_tet[i] >0 )
        {
            ind_remove_tets.insert(i);
        }
    }
    set <idx> ind_remove_tris;
    for( idx i = 0 ; i < (idx) i_tri.size(); i++ )
    {
        if ( i_tri[i] )
        {
            ind_remove_tris.insert(i);
        }
    }

    auto &tets = geom.getTetrahedra();
    auto &tris = geom.getTriangles();
    tets.removeElements( ind_remove_tets );
    tets.appendElements(new_t, new_mat_t);
    tris.removeElements( ind_remove_tris );
    tris.appendElements(new_e, new_mat_e);

    geom.ReorderDielectricNodes(); // Dielectric nodes are moved last
    //geom.t->setMaxNodeNumber( (unsigned int) geom.getnp() );
    tris.setConnectedVolume(&tets);
    tets.calculateDeterminants3D(geom.getCoordinates());
    tets.ScaleDeterminants( 1e-18);// scale to microns cubed
    tris.calculateSurfaceNormals(geom.getCoordinates(), &tets);
    tris.ScaleDeterminants( 1e-12); // scale to microns squared
    geom.initialisePeriodicity();
}

void Refine(Geometry& geom,                 // SOURCE (OLD) GEOMETRY
            vector <idx> & i_tet)        // REFINEMENT TYPES VECTOR
{

    Num_Ref_Tet nrt;    // TYPE COUNTERS
    find_tet_refinement_types( i_tet, nrt); // COUNTS HOW MANY OF EACH EXIST

    if (nrt.red ==  0) { // EXIT IF NO RED TETS
        RUNTIME_ERROR("No red tetrahedra selected.")
    }
    Log::info("{} Red tetrahedra selected.", nrt.red);

    //=================================
    //	2.  EXPAND REFINEMENT REGION TO SATISFY
    //      PERIODICITY OF STRUCTURE
    //=================================
    vector <Line> lines;
    lines.clear();
    vector < set <idx> > t_to_l; // index from tets to lines
    t_to_l.clear();
    expand_refinement_region(i_tet, nrt, lines, geom, t_to_l );
    //======================================
    //  3.  FIND TRIANGLE REFINEMENT TYPES
    //      THESE DO NOT AFFECT TET TYPES AND
    //      CAN BE DONE LAST
    //======================================
    vector <idx> i_tris;
    vector < set <idx> > e_to_l;
    find_triangle_reftypes( geom, i_tris, lines , e_to_l);

    //====================================================
    //  5.  CREATE NEW ELEMENTS. BOTH TETS AND TRIANGLES
    //====================================================
    vector <double> new_p;
    vector <idx>  new_t;
    vector <idx>  new_e;
    vector <idx> new_mat_t;
    vector <idx> new_mat_e;

    // in refinement2.cpp
    create_new_elements( geom,
                         i_tet, i_tris,
                         lines,
                         t_to_l, e_to_l,
                         new_p,
                         new_t, new_mat_t,
                         new_e, new_mat_e);

    //=================================================
    //  6.  MODIFY GEOMETRY BY REPLACING RED ELEMENTS
    //      WITH NEW ONES CREATED IN STEP 5
    //=================================================
    modify_geometry(geom,
                    i_tet, i_tris,
                    new_p, new_t, new_e,
                    new_mat_t, new_mat_e );
}// end void Refine
