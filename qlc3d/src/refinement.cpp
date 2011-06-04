
#include <refinement.h>
#include <material_numbers.h>
#include <meshrefinement.h>
#include <line.h>
#include <set>
using std::set;

struct Sphere{
    double centre[3];
    double radius;
};


void tets_in_sphere( vector<unsigned int>& itet, Geometry& geom , RefReg& sphere, const int& refiter){

    if ( refiter >= sphere.getNumIterations() ) // EXIT IF THIS REFITERATION IS NOT DEFINED FOR THIS SPHERE
	return;


    printf("REFREG%i=Sphere\n", refiter);
	for ( size_t i = 0 ; i < itet.size() ; i++){ // loop over tets
            if (geom.t->getMaterialNumber( i ) <= MAT_DOMAIN7 ){
                for (int j = 0 ; j < geom.t->getnNodes() ; j++){
                    double d[3] ={ geom.getAbsXDist(geom.t->getNode(i,j) , sphere.X[0]),
				   geom.getAbsYDist(geom.t->getNode(i,j) , sphere.Y[0]),
				   geom.getAbsZDist(geom.t->getNode(i,j) , sphere.Z[0])};
                    double dist = sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );

                    if ( dist<= sphere.Distance[refiter] )
                        itet[i] = RED_TET;// 6 = red tet
                    }// end for nodes
	    }// end if LC element
    }// end for tets

}

void tets_near_surface( vector<unsigned int>& itet, Geometry& geom, RefReg& region, const int& refiter){
//    /*! finds all tets that are close to a surface*/
    // A brute-force approach

    // Form a list of all triangles with this material number
    std::vector <unsigned int> surf_elems;
    geom.e->listElementsOfMaterial( surf_elems , 2048 );

    printf("num sur_elems: %i\n", (int) surf_elems.size() );




    // convert surface triangle list to unique surface nodes list
    std::set< unsigned int > surf_nodes;
    for ( unsigned int i = 0; i < surf_elems.size() ; i++ ){
        for (int j = 0 ; j < geom.e->getnNodes() ; j++){
            surf_nodes.insert( (unsigned int) geom.e->getNode( surf_elems[i] , j ) );
        }// end for nodes per element
    }
printf("num surfnodes: %u\n", surf_nodes.size());

    //std::set <unsigned int> ::iterator itr;
    //for (itr = surf_nodes.begin() ; itr != surf_nodes.end() ; itr++){
    //    printf("node %u = [%f,%f,%f]\n", *itr, geom.getpX(*itr), geom.getpY(*itr), geom.getpZ(*itr) );
    //}


    // perform comparison between EVERY tetrahedron and EVERY surface node and calculate
    // distances. This is probably very slow
//*
    for (unsigned int elem = 0 ; elem < (unsigned int) geom.t->getnElements() ; elem++){

        double bary[3] = {0,0,0};
        geom.t->CalcElemBary( elem , geom.getPtrTop() , bary); // gets barycenre of element elem

        std::set< unsigned int > :: iterator node;
        for ( node = surf_nodes.begin() ; node != surf_nodes.end() ; node++ ){
            double pn[3] = {0,0,0};// = geom.getPtrTop() + 3 * (*node); // pointer to coordinates of node
            pn[0] = geom.getpX( *node );
            pn[1] = geom.getpY( *node );
            pn[2] = geom.getpZ( *node );
            // squared distance between barycentre and surface node
            double dsqr =   (bary[0] - pn[0])*(bary[0] - pn[0]) +
                            (bary[1] - pn[1])*(bary[1] - pn[1]) +
                            (bary[2] - pn[2])*(bary[2] - pn[2]) ;

            double mindist = region.getDistance(refiter);

            if (dsqr < ( mindist*mindist) ){ // close enough, mark this element as RED

                itet[elem] = RED_TET;
                break;
            }

        }
    }// end for elements
//*/

}


// GENERATES LIST (vector) OF INDEXES TO TETS THAT ARE CHOSEN FOR REFINEMENT

void get_index_to_tred(vector <unsigned int>& i_tet,
		       Geometry& g_prev,
		       MeshRefinement* meshrefinement,
		       const int& refiter ){

    i_tet.clear();
    i_tet.reserve( g_prev.t->getnElements() );
    i_tet.assign( g_prev.t->getnElements() , 0); // number of elements * 0

    vector<RefReg> :: iterator ritr;
    cout << "num refregs : " <<meshrefinement->RefinementRegion.size() << endl;

    for (ritr = meshrefinement->RefinementRegion.begin(); ritr != meshrefinement->RefinementRegion.end(); ritr++){
	switch (ritr->Type){
	case (RefReg_Sphere):
	    tets_in_sphere( i_tet, g_prev, *ritr, refiter);
	    break;
        case (RefReg_Surface):
            printf("Surface REFREG\n");
            tets_near_surface(i_tet, g_prev, *ritr, refiter);
            break;
        default:
	    cout << "Unknown REFREG - bye!" <<endl;
	    exit(1);
	}// end switch-case

    }// end for refregs
}// end get_index_to_tred



// REPOULATES LIST i_tred WITH IDEXES TO TETS IN m THAT CONTAIN AT LEAST num NODES FROM LIST i_p
void contains_nodes( Mesh* m,				// mesh whose elements are checked
			 set <unsigned int>& i_p,		// index to all internal nodes
			 const int& num,				// minimum number of nodes needed to be included in
			 set< unsigned int>& i_tred){	// index to elements that contain at least 'num' nodes from 'i_p'
    //
    // THIS COULD BE MADE MUCH FASTER IF AN INDEX FROM NODES->ELEMENTS WAS MADE FIRST
    //
    i_tred.clear();
    set <unsigned int> :: iterator p_itr;
    for (unsigned int t = 0 ; t < (unsigned int) m->getnElements() ; t++){ // for all elements
	int count = 0;
	for (unsigned int i = 0 ; i < (unsigned int) m->getnNodes() ; i ++){ // for all nodes in this element
	    for ( p_itr = i_p.begin() ; p_itr != i_p.end() ; p_itr ++ ){
		if ( *p_itr == (unsigned int) m->getNode((int) t , (int) i) ) // if this element contains node *p_itr
		    count ++;

	    }// end for
	}// end for nodes
	if (count >= num )
	    i_tred.insert( t );
    }// end for elements
}// end contains_nodes

//
void find_unique_nodes( Mesh* m , set <unsigned int>& ind_m , set <unsigned int>& i_p){
    i_p.clear();
    set <unsigned int> :: iterator m_itr; // iterator to index of chosen tets
    int nnodes = m->getnNodes();
    for ( m_itr = ind_m.begin() ; m_itr != ind_m.end() ; m_itr++ ) { // loop over each element in mesh
	for ( int n = 0 ; n < nnodes ; n++){ // loop over each node in this element. assumes same number of nodes per element
	    i_p.insert( m->getNode(*m_itr , n) );
	}// end for loop over nodes
    }// end for loop over elements
}

void find_all_core_lines( vector<Line>& lines, const vector <unsigned int>& i_tet, Mesh* m){
// creates unique set of lines from all bisectable lines of red tets
// only works for tets!
    vector <unsigned int>::const_iterator itr; // iterator to red tet indexes

    for (unsigned int i = 0 ; i < i_tet.size() ; i++){ // loop over all tets
	if  ( i_tet[i] == 6 ){ // if this is a red tet
	    int tet[4] = {  m->getNode(i, 0), // shortcut to tet nodes
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

void count_lines(   Mesh* m,						// mesh elements
			vector< set <unsigned int> >& p_to_m ,	// point to elements indices
			const vector<Line>& lines ,				// list of lines
			vector <unsigned int>& i_elem,			// counter array
			vector < set < unsigned int> >& m_to_l	// elements to lines indices
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
    for ( unsigned int l = 0 ; l< (unsigned int) lines.size() ; l++){// for each line
       //if ( (l == 15 ) && (lines.size() == 3653) )
       //    int b = 0;

        int n1 = lines[l].L[0];
        int n2 = lines[l].L[1];

        // IF ELEMENT TYPE IS TRIANGLE, NOT ALL NODES HAVE CONNECTIONS TO A TRIANGLE
        // CHECK THAT NODENUMBERS ARE VALID, THESE
        // SHOULD ALWAYS BE FOR TETS, BUT NOT NOT ALWAYS FOR TRIS

        if ( (n1 < num_points ) && ( n2 < num_points ) ){
            // INDEXES TO ALL TETS CONNECTED TO FIRST  AND SECOND NODES OF THIS LINE
            set <unsigned int> elems1;
            set <unsigned int> elems2;

            elems1.insert( p_to_m[ n1 ].begin() , p_to_m[ n1 ].end() );
            elems2.insert( p_to_m[ n2 ].begin() , p_to_m[ n2 ].end() );

            vector <unsigned int> elems_res;

            set_intersection( elems1.begin() , elems1.end() , elems2.begin() , elems2.end() ,back_inserter(elems_res) );
            // elems_res now contains list of all elements that share this line
            // increase counter for those elements that contain lines in 'line;
            for (int i = 0; i < (int) elems_res.size() ; i++){ // loop over all elements sharing the line
                // increment element line count
                i_elem[ elems_res[i] ]++;
                // adds index from element to this line
                m_to_l[ elems_res[i] ].insert( l ); // can be made faster by inserting last using iterators
            }// end for
        } // end check for valid nodenumber



	}// end for all lines
}// end count lines
void find_tet_refinement_types(vector <unsigned int>& i_tet, // tet line counts
							   Num_Ref_Tet& nrt ){ // tet to lines index
    nrt.green1 = 0; // RESET
    nrt.green2 = 0;
    nrt.green3 = 0;
    nrt.red    = 0;

	for (int i = 0 ; i < (int) i_tet.size() ; i++){ // loop over tet bisection line count
		//cout << " i = " << i << endl;
		if ( i_tet[i]==0 ){} // do nothing
		else
		if ( i_tet[i] ==GREEN1_TET ){
			//i_tet[i] = 0;
			nrt.green1++;
		}
		else
		if ( i_tet[i] == 2 ){
			//i_tet[i] = 0;
			nrt.green2++;
		}
		else
		if ( i_tet[i] == 3 ){
			//i_tet[i] = 0;
			nrt.green3++;
		}
		else
		if ( i_tet[i] > GREEN3_TET){
			nrt.red ++;
			i_tet[i] = 6 ; // enforces anything above and including 4 is actually 6
	}

	}// end for

}


void fix_green3_red_confusions( vector <unsigned int>& i_tet,		// tets line counter
							   vector<Line>& lines,					// bisected lines
							   vector< set<unsigned int> > t_to_l){ // index from tets to  lines
// SOME ELEMENTS MARKED AS GREEN 3 ARE IN FACT RED TETS
// GREEN 3 CONTAIN 3 BISECTABLE LINES FORMING A CLOSED LOOP -> 3 UNIQUE NODES
// IN RED TETS, SOME OTHER CONFIGURATION OF 3 LINES EXIST -> 4 UNIQUE NODES


	// LOOP OVER EACH TET-LINE COUNTER PRAGMA OMP'ABLE
	//#pragma omp parallel for
	for (unsigned int i = 0 ; i < i_tet.size() ; i++){ // for all tets
		// IF PREVIOUSLY MARKED AS A GREEN 3 TET
		if ( i_tet[i] == 3){
			// LOOP OVER ALL ITS BISECABLE LINES AND CHECK
			// NUMBER OF UNIQUE NODES
			set <unsigned int> ::iterator itr;
			set < int > nodes; // nodes set, size of this determines type
			for (itr = t_to_l[i].begin() ; itr != t_to_l[i].end() ; itr++){ // for index to lines
				nodes.insert( lines[*itr].L[0]); // insert both nodes of line
				nodes.insert( lines[*itr].L[1]);
			}// end for index to lines
			//cout << "total number of nodes " << nodes.size() << endl;
			if (nodes.size() == 4 ){
				// THIS ELEMENT IS ACTUALLY RED. MAKE IT SO!
				i_tet[i] = 6; // 6 is magic for red
			}
			#ifdef DEBUG
			else if ( nodes.size() != 3){ // debig
				cout << " error in fix_green_red_confusions - this should not happen "<< endl;
				cout << " something has gone badly wrong - bye!" << endl;
				exit(1);
			}
			# endif
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
	for (int i = 0 ; i < geom.e->getnElements() ; i++ ) 
	{
		// if triangle is periodic
		if ( geom.e->getMaterialNumber(i) == MAT_PERIODIC )
		{
			int* n = geom.e->getPtrToElement( i ); // pointer to element node numbers
						
						
			Line lines[3] = {Line(n[0], n[1]), Line(n[0], n[2]) , Line(n[1], n[2]) };
			
			for ( int l = 0 ; l < 3 ; l++)
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

void expand_periodic_boundaries( 	vector <Line> lines, // lines to split
									peri_lines& plines,  // lines on periodic boundary
									Geometry& geom )
{
	vector <Line> newlines; 
	vector <Line> :: iterator litr;
	vector <Line> :: iterator lotr;
	vector <Line> :: iterator loend; // end iterator
	for (litr = lines.begin() ; litr != lines.end() ; litr++ ) // loop over each line
	{
		
		if ( geom.getfront_back_is_periodic() )
		{
			
			lotr = loend; // stops comparison
			bool found = true;	
			if ( litr->isOnFrontSurface(&geom) )
			{
				lotr = plines.lback.begin();
				loend = plines.lback.end();
				printf("front...");
				found = false;
			}
			else
			if ( litr->isOnBackSurface(&geom ) )
			{
				lotr = plines.lfront.begin();
				loend = plines.lfront.end();
				printf("back...");
				found = false;
			}
							
			for (; lotr != loend ; lotr++)
			{
				if ( litr->isTranslationOf(*lotr, geom) )
				{
					newlines.push_back(*lotr);
					printf("found\n");
					found = true;
					break;
				}
			}
			if (!found)
			{
				printf("error - could not find periodic front/back line - bye!\n");
				litr->PrintLine();
				exit(1);
			}
		}// end if front/back
		
		if (geom.getleft_right_is_periodic() )
		{
			lotr = loend;
			if ( litr->isOnLeftSurface(&geom) )
			{
				lotr = plines.lright.begin();
				loend = plines.lright.end();
				printf("left...");
			}
			else
			if (litr->isOnRightSurface(&geom) )
			{
				lotr = plines.lleft.begin();
				loend = plines.lleft.end();
				printf("right...");
			}
			
			for(; lotr!=loend; lotr++)
			{
				if (litr->isTranslationOf(*lotr, geom ) )
				{
					newlines.push_back(*lotr);
					printf("found\n");
					break;
				}
			}
		
			// in case of l/r periodicity, also vertical corners must be 
			// taken into account. For each corner line, must find 3 other
			// corner lines
			if ( litr->isTopBottomCornerLine(&geom) )
			{
				int cc = 0; //corner count
				vector <Line>::iterator cs1;
				vector <Line>::iterator ce1;
				vector <Line>::iterator cs2;
				vector <Line>::iterator ce2;
				vector <Line>::iterator cs3;
				vector <Line>::iterator ce3;
			
				cs1 = ce1; 
				cs2 = ce2; 
				cs3 = ce3;
			
				if (litr->isCorn0(&geom) )
				{
					printf("corn0...");
					cs1 = plines.lc1.begin(); ce1 = plines.lc1.end();
					cs2 = plines.lc2.begin(); ce2 = plines.lc2.end();
					cs3 = plines.lc3.begin(); ce3 = plines.lc3.end();
				}
				else
				if (litr->isCorn1(&geom) )
				{
					printf("corn1...");
					cs1 = plines.lc0.begin(); ce1 = plines.lc0.end();
					cs2 = plines.lc2.begin(); ce2 = plines.lc2.end();
					cs3 = plines.lc3.begin(); ce3 = plines.lc3.end();
				}
				else
				if (litr->isCorn2(&geom) )
				{
					printf("corn2...");
					cs1 = plines.lc0.begin(); ce1 = plines.lc0.end();
					cs2 = plines.lc1.begin(); ce2 = plines.lc1.end();
					cs3 = plines.lc3.begin(); ce3 = plines.lc3.end();
				}
				else
				if (litr->isCorn3(&geom) )
				{
					printf("corn3...");
					cs1 = plines.lc0.begin(); ce1 = plines.lc0.end();
					cs2 = plines.lc1.begin(); ce2 = plines.lc1.end();
					cs3 = plines.lc2.begin(); ce3 = plines.lc2.end();
				}
					
				for (;cs1 != ce1 ; cs1++)
				{
					if (litr->isTranslationOf(*cs1 , geom ) )
					{
						cc++;
						newlines.push_back(*cs1);
						printf("c1 ");
						break;
					}
				}
		
				for (;cs2 != ce2 ; cs2++)
				{
					if(litr->isTranslationOf(*cs2 , geom ) )
					{
						cc++;
						newlines.push_back(*cs2);
						printf("c2 ");
						break;
					}
				}
			
				for (;cs3 != ce3; cs3++)
				{
					if (litr->isTranslationOf(*cs3 , geom) )
					{
						cc++;
						newlines.push_back(*cs3);
						printf("c3");
						break;
					}
				}
			
		
			printf("cc = %i\n",cc);
			}// end if tvertical corner		
		}// end if left/right
		
		
	}//end for i loop over each line



}

void expand_refinement_region(vector <unsigned int>& i_tet,	// index to tet bisectable lines counts
			      Num_Ref_Tet& nrt,
			      vector <Line>& lines,
			      Geometry& geom_prev,		// geometry
			      vector < set <unsigned int> >& t_to_l        // tets to lines index
			      ){	// index to core nodes
    if ( nrt.red == 0 )
	return; // exit if nothing to do

	unsigned int n_tred_old = 0;//nrt.red ;					// old number of red tets
	unsigned int d_n_tred   = nrt.red;				// change in number of tred tets

    vector < set <unsigned int> > p_to_t;			// p to tets index
    geom_prev.t->gen_p_to_elem( p_to_t );
    
    peri_lines plines;
    // Construct list of periodic line elements iff structure has periodic boundaries
    if (  geom_prev.getleft_right_is_periodic()  ||
		  geom_prev.getfront_back_is_periodic()  ||
		  geom_prev.gettop_bottom_is_periodic()  )
	{
		gen_peri_lines_vec( plines,	geom_prev);
		printf(" f,b,l,r = %i,%i,%i,%i\n", plines.lfront.size(), plines.lback.size() , plines.lleft.size() , plines.lright.size() );
		printf(" c 0 1 2 3 = %i, %i, %i,%i\n", plines.lc0.size() , plines.lc1.size(), plines.lc2.size() , plines.lc3.size() );
		
	}
    
    while( d_n_tred > 0){// loop while refinement region grows

		find_all_core_lines(lines, i_tet, geom_prev.t); // creates all bisectable lines
		
		// periodic lines added here
		expand_periodic_boundaries( lines , plines, geom_prev );
		exit(1);
		count_lines(geom_prev.t, p_to_t, lines, i_tet , t_to_l);

		fix_green3_red_confusions(i_tet, lines, t_to_l);
		n_tred_old = nrt.red;

		find_tet_refinement_types( i_tet, nrt);


		// CHECK WHETHER NEW ELEMENTS HAVE BEEN ADDED
		d_n_tred = nrt.red - n_tred_old;
		//cout << "n_old " << n_tred_old << " n_now " << nrt.red << "change = " << d_n_tred<< endl;

		//break;
    }// end while refinement region grows

}// end expand_refinement_region





void find_triangle_reftypes(Geometry& geom,
							vector <unsigned int>& i_tri,
							vector <Line>& lines,
							vector < set <unsigned int> >& e_to_l
							){

	// MAKE LINKS FROM NODES TO TRIS
	vector < set <unsigned int> > p_to_e;

	geom.e->gen_p_to_elem( p_to_e );
	// MAKE INDEXES FROM TRIS TO LINES

	count_lines( geom.e , p_to_e, lines, i_tri, e_to_l);

}

void make_new_geometry(Geometry& geom_new, Geometry& geom_orig,
		       vector <unsigned int>& i_tet, vector <unsigned int>& i_tri,
		       vector <double>& new_p,
		       vector <unsigned int>& new_t, vector <unsigned int>& new_e,
		       vector <int> new_mat_t , vector <int> new_mat_e
		       ){
    geom_new.ClearGeometry();
    geom_new.setTo( &geom_orig );

    if (i_tet.size() == 0) return;

    geom_new.addCoordinates( new_p );
    // MAKE LIST OF TETS AND TRIS TO BE REMOVED. THESE ARE REPLACED BY THE NEWLY
    // CREATED ONES
    set <unsigned int> ind_remove_tets;
    for (unsigned int i = 0 ; i < (unsigned int) i_tet.size() ; i++){
	    if ( i_tet[i] ) {
		    ind_remove_tets.insert(i);
	    }
    }
    set <unsigned int> ind_remove_tris;
    for( unsigned int i = 0 ; i < (unsigned int) i_tri.size(); i++ ){
	    if ( i_tri[i] ){
		    ind_remove_tris.insert(i);
	    }
    }
    geom_new.t->removeElements( ind_remove_tets );
    geom_new.t->addElements(new_t , new_mat_t );
    geom_new.t->setMaxNodeNumber( (unsigned int) geom_new.getnp() );

    geom_new.t->CalculateDeterminants3D( geom_new.getPtrTop() );
    geom_new.t->CalculateDeterminants3D( geom_new.getPtrTop() );
    geom_new.t->ScaleDeterminants( 1e-18);// scale to microns cubed

    // IF TOTAL VOLUME HAS CHANGED, SOMETHING HAS GONE WRONG
    // 1e-20 IS ARBITRARY ACCURACY, SOME NUMERICAL NOISE IS EXPECTED, BUT NOT MUCH
    if ( fabs(geom_orig.t->getTotalSize() - geom_new.t->getTotalSize()) >= 1e-20 ){
	cout << "error - total volume changed in refinement" << endl;
	cout << "volume before ="<< geom_orig.t->getTotalSize() << endl;
	cout << "volume after = "<< geom_new.t->getTotalSize() << endl;
	exit(1);
    }
    geom_new.e->removeElements( ind_remove_tris );
    geom_new.e->addElements( new_e, new_mat_e );
	geom_new.ReorderDielectricNodes(); // Dielectric nodes are moved last
	geom_new.e->setConnectedVolume( geom_new.t );
    geom_new.e->CalculateSurfaceNormals( geom_new.getPtrTop(), geom_new.t);
    geom_new.e->ScaleDeterminants( 1e-12); // scale to microns squared


   geom_new.setNodeNormals();

    geom_new.checkForPeriodicGeometry();


}

void Refine(Geometry& srce, Geometry& dest, vector <unsigned int>& i_tet){


	Num_Ref_Tet nrt;
	find_tet_refinement_types( i_tet, nrt);



	if (nrt.red ==  0){ // EXIT NOW IF NO REFIENEMENT
		dest.ClearGeometry();
		dest.setTo( &srce );
//		cout << "refinement iteration: " << refiter << " no refinable tets. Done" << endl;
//		return;
	}

	//=================================
	//	2. EXPAND REFINEMENT REGION
	//=================================
	vector <Line> lines;
	lines.clear();
	vector < set <unsigned int> > t_to_l; // index from tets to lines
	t_to_l.clear();

	expand_refinement_region(i_tet, nrt, lines, srce, t_to_l );

	//======================================
	//  3.	FIND TRIANGLE REFINEMENT TYPES
	//	THESE DO NOT AFFECT TET TYPES AND
	//	CAN BE DONE LAST
	//======================================
	vector <unsigned int> i_tris;
	vector < set <unsigned int> > e_to_l;

	find_triangle_reftypes( srce, i_tris, lines , e_to_l);

	//======================================
	//  5. CREATE NEW ELEMENTS
	//======================================
	vector <double> new_p;
	vector <unsigned int>  new_t;
	vector <unsigned int>  new_e;
	vector <int> new_mat_t;
	vector <int> new_mat_e;

	create_new_elements( srce,
				i_tet, i_tris,
				lines,
				t_to_l, e_to_l,
				new_p,
				new_t, new_mat_t,
				new_e, new_mat_e);

	//======================================
	//  6. CREATE NEW GEOMETRY
	//======================================
	make_new_geometry(dest , srce, i_tet, i_tris, new_p, new_t, new_e, new_mat_t, new_mat_e );


}


void Refine(Geometry& geom_orig, Geometry& geom_prev, Geometry& geom_new ,
            MeshRefinement* meshrefinement){
/*! Pre - Refinement */

	geom_orig.getnp(); // NO WARNINGS

    int MaxRefIter = 0;
    if (meshrefinement){
        MaxRefIter = meshrefinement->getMaxNumRefIterations();
        if (MaxRefIter == 0){
            geom_new.setTo( &geom_prev );
            return;
        }
    }

    cout << "---------doing " << MaxRefIter << " refinement iterations-------" << endl;
    for (int refiter = 0 ; refiter < MaxRefIter ; refiter++){ // DO REFINEMENT ITERATIONS
//	cout << "\tIteration " << refiter << endl;

	//=================================
	//	1. SELECT RED TETRAHEDRA
	//=================================
	vector <unsigned int> i_tet; // each value correponds to number of bisected edges per element. size of this equals number of elements
	i_tet.clear();
	Num_Ref_Tet nrt;
	if (meshrefinement){
	    get_index_to_tred(i_tet, geom_prev, meshrefinement, refiter);
	}
	else{
		cout <<"No refinement regions => nothing to refine" << endl;
		return;
	}


	find_tet_refinement_types( i_tet, nrt);

	if (nrt.red ==  0){ // EXIT NOW IF NO REFIENEMENT
	    geom_new.ClearGeometry();
	    geom_new.setTo( &geom_prev );
	    cout << "refinement iteration: " << refiter << " no refinable tets. Done" << endl;
	    return;
	}

	//=================================
	//	2. EXPAND REFINEMENT REGION
	//=================================
	vector <Line> lines;
	lines.clear();
	vector < set <unsigned int> > t_to_l; // index from tets to lines
	t_to_l.clear();
	expand_refinement_region(i_tet, nrt, lines, geom_prev, t_to_l );

	//======================================
	//  3.	FIND TRIANGLE REFINEMENT TYPES
	//	THESE DO NOT AFFECT TET TYPES AND
	//	CAN BE DONE LAST
	//======================================
	vector <unsigned int> i_tris;
	vector < set <unsigned int> > e_to_l;
	find_triangle_reftypes( geom_prev, i_tris, lines , e_to_l);

	//======================================
	//  5. CREATE NEW ELEMENTS
	//======================================
	vector <double> new_p;
	vector <unsigned int>  new_t;
	vector <unsigned int>  new_e;
	vector <int> new_mat_t;
        vector <int> new_mat_e;



	create_new_elements( geom_prev,
				i_tet, i_tris,
				lines,
				t_to_l, e_to_l,
				new_p,
				new_t, new_mat_t,
				new_e, new_mat_e);

	//======================================
	//  6. CREATE NEW GEOMETRY
	//======================================
	make_new_geometry(geom_new , geom_prev, i_tet, i_tris, new_p, new_t, new_e, new_mat_t, new_mat_e );

	geom_prev.setTo( &geom_new ); // sets previous to new for next refinement iteration
	cout << "\tNumber of nodes: "<<geom_prev.getnp() << endl;
	cout << "\tNumber of Tets : "<<geom_prev.t->getnElements() << endl;
	cout << "\tNumber of Tris : "<<geom_prev.e->getnElements() << endl;
    }// end for refinement iterations


    printf("::::::::::::::::::::::::::::::::::::::::::::::\n");
    printf(":              Done refining                 :\n");
    printf(":     tets %i, tris %i, nodes %i             :\n", geom_new.t->getnElements(), geom_new.e->getnElements(), geom_new.getnp() );
    printf("::::::::::::::::::::::::::::::::::::::::::::::\n");

 //   geom_new.t->PrintElements();
    /*
    geom_new.checkForOverlapingNodes();
vector <int> refc;
int minc[2] = {100000,0};
int maxc[2] = {0,0};
    geom_new.countNodeReferences(refc, *geom_new.t);
    for (int i = 0 ; i < geom_new.getnp(); i++ ){
        if (refc[i] < minc[0] ) { minc[0] = refc[i]; minc[1] = i;}
        if (refc[i] > maxc[0] ) { maxc[0] = refc[i]; maxc[1] = i;}
        printf("node [%i] = %i\n ", i, refc[i]);
    }

    printf(" minc = %i, %i\n", minc[0], minc[1] );
    printf(" maxc = %i, %i\n", maxc[0], maxc[1] );
    geom_new.PrintNode( minc[1]);
    geom_new.PrintNode( maxc[1]);
*/
}
