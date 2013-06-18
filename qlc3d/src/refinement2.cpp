#include <refinement.h>
#include <algorithm>
#include <geometry.h>
#include <line.h>
#include <vector>
#include <qlc3d.h>
#include <globals.h>

#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_matrixmaker.hpp>
/*
void create_node_number_matrix(SparseMatrix*& nnumbers,
			       Geometry& geom,
			       vector <Line>& lines){
// COMMENTED DURING SPAMTRIX MIGRATION
    nnumbers = createSparseMatrix( lines ); // this creates a matrix for the mesh. should make one for the lines, but this will do for now

    // LOOP OVER ALL LINES AND ADD INDEXES TO NEW NODES
    unsigned int nold = geom.getnp(); // number of old nodes

    for (unsigned int i = 0 ; i < lines.size() ; i++){
	int n1 = lines[i].L[0];
	int n2 = lines[i].L[1];
	nnumbers->sparse_set( n1 , n2 , nold + i ); // adding same node on both diagonals is wasteful
	nnumbers->sparse_set( n2 , n1 , nold + i ); // but may make it easier to access matrix as M(i,j) or M(j,i)
    }
}
*/
SpaMtrix::IRCMatrix createNodeNumbersMatrix(Geometry &geom, vector<Line> &lines){
/*!Creates sparse matrix with node numbers for new nodes*/

    // USING A SPARSE MATRIX OF DOUBLE TYPE INSTED OF UNSIGNED INT MAY BE BAD IDEA...
    unsigned int nold = geom.getnp();
    SpaMtrix::MatrixMaker mm(nold,nold);
    unsigned int i = 0;
    for (const auto &l : lines ){
        mm.addNonZero(l.L[0], l.L[1], nold+i);
        mm.addNonZero(l.L[1], l.L[0], nold+i);
        i++;
    }
    return mm.getIRCMatrix();
}

void create_new_coordinates( Geometry& geom,
                             vector <Line>& lines,
                             vector <double>& new_p
                             ){
    new_p.clear();
    new_p.reserve( lines.size()*3 ); // reserve space for 3 coordinates per new node
    // LOOP OVER LINES, CALCULATE MID-POINT LOCATION AND ADD TO NEW COORDINATES
    for (unsigned int i = 0 ; i < lines.size() ; i++){

        double x1 = geom.getpX( lines[i].L[0] );
        double y1 = geom.getpY( lines[i].L[0] );
        double z1 = geom.getpZ( lines[i].L[0] );

        double x2 = geom.getpX( lines[i].L[1] );
        double y2 = geom.getpY( lines[i].L[1] );
        double z2 = geom.getpZ( lines[i].L[1] );

        double xn = ( x1 + x2 ) / 2.0;
        double yn = ( y1 + y2 ) / 2.0;
        double zn = ( z1 + z2 ) / 2.0;
        //printf( "new %f %f %f\n", xn, yn, zn);
        new_p.push_back( xn );
        new_p.push_back( yn );
        new_p.push_back( zn );
    }
}// end create_new_coordinates

void make_new_green1_tet( vector <idx>& new_t,
                          vector <idx>& new_mat_t,
                          Geometry& geom,
                          const idx& elem ,
                          vector< Line>& lines,
                          vector< set <unsigned int> >& t_to_l,
                          const SpaMtrix::IRCMatrix &nnodes){
                          //SparseMatrix* nnodes){ <-

    idx ln = *( t_to_l[elem].begin() ); // index to the only bisect line

    // GENERTE LISTS OF OLD NODES AND ONE NEW NODE
    vector <idx> no; // old nodes
    vector <idx> nn; // new nodes
    no.push_back( (idx) lines[ln].L[0] );
    no.push_back( (idx) lines[ln].L[1] );
    geom.t->CompleteNodesSet( elem , no );

    nn.push_back( (unsigned int) nnodes.sparse_get(nA, nB));//nnodes->sparse_get(nA , nB) ); <-

    // MAKE 2 NEW TETS
    // TET1 A,C,D,AB
    // TET2 B,C,D,AB
    unsigned int tet[8] = {	nA,	nC,	nD,	nAB,
				nAB,	nC,	nD,	nB};

    new_t.insert( new_t.end() , tet , tet + 8 );
    idx mat[2] = {geom.t->getMaterialNumber(elem), geom.t->getMaterialNumber(elem)};
    new_mat_t.insert ( new_mat_t.end() , mat , mat + 2 );

}// end make_new_green1_tet

void make_new_green2_tet( vector <unsigned int>& new_t,
                          vector <idx>& new_mat_t,
                          Geometry& geom,
                          const idx& elem ,
                          vector< Line>& lines,
                          vector< set <idx> >& t_to_l,
                          const SpaMtrix::IRCMatrix &nnodes)
                          //SparseMatrix* nnodes)
{
    //GENERATE LIST OF OLD AND NEW NODES
    vector <idx> no; // old;
    vector <idx> nn; // new;

    set <idx> :: iterator itr;
    set <idx> unodes ; //unique nodes
    itr = t_to_l[elem].begin();
    unodes.insert( lines[*itr].L[0] );
    unodes.insert( lines[*itr].L[1] );
    ++itr;
    unodes.insert( lines[*itr].L[0] );
    unodes.insert( lines[*itr].L[1] );
    no.insert( no.end() , unodes.begin() , unodes.end() ); // unique nodes


    // CHECK FOR NUMBER OF UNIQUE NODES RECOVERED FROM LINES
    // THIS DETERMINES WHETHER THIS IS green2a OR green2b
    idx nunique = (idx) unodes.size();
    idx* tet	= new idx[nunique * 4];
    idx* mat	= new idx[nunique];

    if ( nunique == 4 ){
        //cout << "TET 2 B" << endl;

        // RESET OLD NODES FROM LINES. THIS PRESERVES A-B, C-D NODE NUMBERING
        // OTHERWISE NONEXISTENT BC, AD etc. LINES MIGHT BE ATTEMPTED
        itr = t_to_l[elem].begin();
        nA = lines[*itr].L[0];
        nB = lines[*itr].L[1];
        ++itr;
        nC = lines[*itr].L[0];
        nD = lines[*itr].L[1];

        nn.push_back((unsigned int) nnodes.sparse_get(nA, nB));// nnodes->sparse_get( nA, nB) ); // AB
        unsigned int temp[4] = {0,0,0,0}; // AC, AD, BC, BD
        nn.insert( nn.end() ,temp, temp+4 ); // PADDING NODES, NOT USED
        nn.push_back((unsigned int) nnodes.sparse_get(nC, nD));//nnodes->sparse_get( nC, nD) ); // CD

        tet[0] = nA;
        tet[1] = nAB;
        tet[2] = nCD;
        tet[3] = nC;

        tet[4] = nA;
        tet[5] = nAB;
        tet[6] = nD;
        tet[7] = nCD;

        tet[8]  = nAB;
        tet[9]  = nB;
        tet[10] = nCD;
        tet[11] = nC;

        tet[12] = nAB;
        tet[13] = nB;
        tet[14] = nD;
        tet[15] = nCD;

        int m = geom.t->getMaterialNumber(elem);
        mat[0] = m;
        mat[1] = m;
        mat[2] = m;
        mat[3] = m;

    }
    else
	if ( nunique == 3 ){
            //cout << "TET2A"<< endl;
            //geom.t->CompleteNodesSet( elem , no ); // ADD MISSING FOURTH NODE nD from TET

            // CONDITIONS
            // nA is shared node
            // nB < nC
            // THESE HAVE TO BE ENFORCED IN TRIANGLES TOO

            // 1. FIND SHARED NODE NUMBER
            vector <unsigned int> nodes;
            itr = t_to_l[elem].begin();

            nodes.push_back( lines[*itr].L[0] );
            nodes.push_back( lines[*itr].L[1] );
            ++itr;
            nodes.push_back( lines[*itr].L[0] );
            nodes.push_back( lines[*itr].L[1] );

            // NODES NOW CONTAINS 4 VALUES, 3 UNIQUE
            // FIND REPEATING VALUE
            sort( nodes.begin(), nodes.end() );

            //printf("nodes = %u,%u,%u,%u\n", nodes[0], nodes[1], nodes[2], nodes[3] );
            vector< unsigned int> ::iterator rep = adjacent_find( nodes.begin(), nodes.end() );
            nA = *rep; // SHARED NOD = nA

            remove(nodes.begin(), nodes.end(), nA); // REMOVE nA OCURRENCES FROM nodes

            // 2. SET nB and nC
            nB = nodes[0] < nodes[1]? nodes[0]:nodes[1]; // return smaller
            nC = nodes[0] > nodes[1]? nodes[0]:nodes[1]; // return larger

            geom.t->CompleteNodesSet( elem, no);
            //printf("n = %u, %u, %u, %u\n", nA, nB, nC, nD);

            nn.push_back((unsigned int) nnodes.sparse_get(nA, nB) ); // <-
            nn.push_back((unsigned int) nnodes.sparse_get(nA, nC) ); // <-

            tet[ 0 ] = nA;
            tet[ 1 ] = nAB;
            tet[ 2 ] = nD;
            tet[ 3 ] = nAC;

            tet[ 4 ] = nAC;
            tet[ 5 ] = nAB;
            tet[ 6 ] = nD;
            tet[ 7 ] = nC;

            tet[ 8 ] = nAB;
            tet[ 9 ] = nB;
            tet[ 10] = nD;
            tet[ 11] = nC;

            int m = geom.t->getMaterialNumber(elem);
            mat[0] = m;
            mat[1] = m;
            mat[2] = m;
	}
	else{
            cout << "tet green2 has " << no.size() << "nodes - bye!" << endl;
            exit(1);
	}

    // add to 'global' new element and material lists
    new_t.insert(new_t.end() , tet, tet + (nunique*4) );
    new_mat_t.insert( new_mat_t.end() , mat, mat + nunique );

    delete [] tet;
    delete [] mat;
    //cout << "tet2, unique nodes :" << no.size() << endl;


}

void make_new_red_tet( vector <idx>& new_t,
                       vector <idx>& new_mat_t,
                       Geometry& geom,
                       const unsigned int& elem ,
                       vector< Line>& lines,
                       vector< set <idx> >& t_to_l,
                       const SpaMtrix::IRCMatrix &nnodes){
                       //SparseMatrix* nnodes){

    lines.begin(); // silence compiler warnings
    t_to_l.begin(); // NO WARNINGS


    // GENERATE LIST OF OLD AND NEW NODES
    vector <unsigned int> no;// old
    vector <unsigned int> nn;// new
    // make old nodes list
    for (int i = 0 ; i < 4 ; i++) no.push_back( geom.t->getNode( elem , i ) );

    // make new nodes list
    nn.push_back((unsigned int) nnodes.sparse_get(nA,nB)); // AB <-
    nn.push_back((unsigned int) nnodes.sparse_get(nA,nC)); // AC <-
    nn.push_back((unsigned int) nnodes.sparse_get(nA,nD)); // AD <-
    nn.push_back((unsigned int) nnodes.sparse_get(nB,nC)); // BC <-
    nn.push_back((unsigned int) nnodes.sparse_get(nB,nD)); // BD <-
    nn.push_back((unsigned int) nnodes.sparse_get(nC,nD)); // CD <-


    // CREATE 8 NEW ELEMENTS
    // the ordering could be improved to ensure positive determinants
    unsigned int tet [ 4 * 8] = { nA, nAB, nAC, nAD ,
                                  nB, nBC, nAB, nBD ,
                                  nC, nAC, nBC, nCD ,
                                  nD, nAD, nCD, nBD ,
                                  nAB, nAC, nAD, nBD,
                                  nAB, nAC, nBD, nBC,
                                  nAC, nAD, nBD, nCD,
                                  nAC, nBC, nCD, nBD};

    new_t.insert( new_t.end() , tet , tet + (4*8) );
    int m = geom.t->getMaterialNumber( elem );
    int mat[8] = {m,m,m,m,m,m,m,m};
    new_mat_t.insert ( new_mat_t.end() , mat, mat+ 8);
}// end void make_new_red_tet

void make_new_green3_tet( vector <idx>& new_t,
                          vector <idx>& new_mat_t,
                          Geometry& geom,
                          const idx& elem ,
                          vector< Line>& lines,
                          vector< set <idx> >& t_to_l,
                          const SpaMtrix::IRCMatrix &nnodes){
                          //SparseMatrix* nnodes){
    // GENERATE LIST OF OLD AND NEW NODES
    vector <unsigned int> no ; // old
    vector <unsigned int> nn ; // new
    // GET 3 UNIQUE NODES FROM 3 LINES, THESE ARE NODES A, B, C
    set <unsigned int > un;

    set <unsigned int> ::iterator itr;
    for ( itr = t_to_l[elem].begin() ; itr!= t_to_l[elem].end() ; ++itr){ // loop over line indexes
        un.insert( lines[*itr].L[0] );
        un.insert( lines[*itr].L[1] );
    }
    if (un.size() != 3) {
        cout << "error, three nodes expected - bye !" << endl;
        exit(1);
    }
    no.insert( no.end() , un.begin() , un.end() );

    geom.t->CompleteNodesSet( elem , no ); // get remaining node

    // POPULATE NEW NODES LIST
    nn.push_back( (unsigned int) nnodes.sparse_get(nA, nB) ); // AB <-
    nn.push_back( (unsigned int) nnodes.sparse_get(nA, nC) ); // AC <-
    nn.push_back( 0 ); // dummy AD
    nn.push_back( (unsigned int) nnodes.sparse_get(nB, nC) ); // BC <-
    //printf("newn = %u, %u, %u, %u", nAB, nAC, nAD, nBC);
    // MAKE 4 NEW TETS
    unsigned int tet[4*4] = { nA, nD, nAB, nAC,
                              nB, nD, nAB, nBC,
                              nC, nD, nAC, nBC,
                              nD, nAB,nAC, nBC};
    new_t.insert( new_t.end() , tet, tet+ 4*4);
    // MATERIAL NUMBERS FOR 4 NEW ELEMENTS

    int m = geom.t->getMaterialNumber(elem);
    int mat[4] = {m,m,m,m};


    new_mat_t.insert( new_mat_t.end() , mat , mat+4);


} // end void make_new_green3_tet

void make_new_tri1(vector <idx>& new_e,
                   vector <idx>& new_mat_e,
                   Geometry& geom,
                   const idx& elem,
                   vector <Line>& lines,
                   vector <set <idx> >& e_to_l,
                   const SpaMtrix::IRCMatrix &nnodes
                   //SparseMatrix* nnodes
                   ){
    // GENERTE LIST OF OLD AND NEW NODES
    vector <unsigned int> no; //old
    vector <unsigned int> nn; // new

    // MAKE OLD NODES LIST
    set <unsigned int> li_ind; // a single line index
    li_ind = e_to_l[elem];
    no.push_back( lines[*(li_ind.begin() ) ].L[0] );
    no.push_back( lines[*(li_ind.begin() ) ].L[1] );
    geom.e->CompleteNodesSet( elem , no );

    // MAKE NEW NODES LIST
    nn.push_back( (unsigned int) nnodes.sparse_get( nA, nB ) ); // <-


    // CREATE 2 NEW TRIANGLES
    unsigned int tri[3 * 2] = {nA, nAB, nC,
                               nB, nC , nAB};



    new_e.insert( new_e.end() , tri , tri + (3*2) );

    // ADD 2 NEW MATERIAL NUMBERS
    int m = geom.e->getMaterialNumber( elem );
    int mat[2] = {m,m};
    new_mat_e.insert( new_mat_e.end() , mat , mat + 2 );
} // end make new tri1
void make_new_tri2(vector <idx>& new_e,
                   vector <idx>& new_mat_e,
                   Geometry& geom,
                   const unsigned int& elem,
                   vector <Line>& lines,
                   vector <set <idx> >& e_to_l,
                   const SpaMtrix::IRCMatrix &nnodes
                   //SparseMatrix* nnodes
                   ){
    // GENERTE LIST OF OLD AND NEW NODES
    vector <unsigned int> no; //old
    vector <unsigned int> nn; // new

    set <unsigned int> nu; // unique, sorted old nodes
    set <unsigned int> ::iterator iter;
    iter = e_to_l[ elem ].begin();
    //printf("num lines %i\n", e_to_l[elem].size() );


    nu.insert( lines[*iter].L[0] );
    nu.insert( lines[*iter].L[1] );
    ++iter;
    nu.insert( lines[*iter].L[0] );
    nu.insert( lines[*iter].L[1] );



#ifdef DEBUG
    if (nu.size() != 3 ){
        printf(" error - tri2 has %i unique nodes\n", (int) nu.size() );
        exit(1);
    }
#endif
    // 1. FIND SHARED NODE NUMBER
    vector <unsigned int> nodes;
    iter = e_to_l[ elem ].begin();
    nodes.push_back( lines[*iter].L[0] ) ;
    nodes.push_back( lines[*iter].L[1] ) ;
    ++iter;
    nodes.push_back( lines[*iter].L[0] ) ;
    nodes.push_back( lines[*iter].L[1] ) ;

    // NODES NOW CONTAINS 4 VALUES, 3 UNIQUE
    // FIND REPEATING VALUE
    no.assign(4 , 0);
    sort( nodes.begin(), nodes.end() );
    vector< unsigned int> ::iterator rep = adjacent_find( nodes.begin(), nodes.end() );
    nA = *rep; // SHARED NOD = nA
    remove(nodes.begin(), nodes.end(), nA); // REMOVE nA OCURRENCES FROM nodes

    // 2. SET nB and nC
    nB = nodes[0] < nodes[1]? nodes[0]:nodes[1]; // return smaller
    nC = nodes[0] > nodes[1]? nodes[0]:nodes[1]; // return larger


    nn.push_back( nnodes.sparse_get( nA, nB ) ); // <-
    nn.push_back( nnodes.sparse_get( nA, nC ) ); // <-

    printf("new nodes = %u, %u\n", nn[0] , nn[1] );
    // MAKE 3 NEW TRIANGLES

    unsigned int tri[3*3] = {	nA, nAB, nAC ,
                                nAC, nAB, nC ,
                                nAB, nB, nC};

    int m = geom.e->getMaterialNumber( elem );
    int mat[3] = {m,m,m};

    new_e.insert( new_e.end() , tri , tri + (3*3) );
    new_mat_e.insert(new_mat_e.end() , mat , mat+3);

} // end make new tri2


void make_new_tri3(vector <idx>& new_e,
                   vector <idx>& new_mat_e,
                   Geometry& geom,
                   const idx& elem,
                   vector <Line>& lines,
                   vector <set <idx> >& e_to_l,
                   const SpaMtrix::IRCMatrix &nnodes
                   //SparseMatrix* nnodes
                   ){
    lines.begin(); e_to_l.begin(); // NO COMPILER WARNINGS
    // GENERATE LIST OF OLD AND NEW NODES
    vector <unsigned int> no; // old
    vector <unsigned int> nn; // new

    // MAKE OLD NODES LIST
    no.push_back( geom.e->getNode( elem, 0) );
    no.push_back( geom.e->getNode( elem, 1) );
    no.push_back( geom.e->getNode( elem, 2) );
    // MAKE NEW NODES LIST
    nn.push_back( nnodes.sparse_get(nA, nB) ); // AB <-
    nn.push_back( nnodes.sparse_get(nA, nC) ); // AC <-
    nn.push_back( 0 );							// AD ,dummy, not used , D does not exist in tris
    nn.push_back( nnodes.sparse_get(nB, nC) ); // BC <-

    // CREATE 4 NEW TRIANGLES

    unsigned int tri[3*4] = {nA, nAB, nAC,
                             nAB, nB, nBC,
                             nAC, nBC, nC,
                             nAC, nAB, nBC};
    new_e.insert( new_e.end() , tri , tri + (3*4) );

    // CREATE 4 NEW MATERIAL NUMBERS
    int m = geom.e->getMaterialNumber( elem );
    int mat[4] = {m,m,m,m};
    new_mat_e.insert( new_mat_e.end() , mat , mat + 4 );

}


void count_refinement_types(int& nred, int& ngreen1, int& ngreen2, int& ngreen3, vector<unsigned int>& i_tet){
    /*! Counts the number of tets of each refinement type */
    nred = 0;
    ngreen1 = 0;
    ngreen2 = 0;
    ngreen3 = 0;
    for (size_t i = 0 ; i < i_tet.size() ; i++){
        if (i_tet[i] == RED_TET ) nred++;
        if (i_tet[i] == GREEN1_TET ) ngreen1++;
        if (i_tet[i] == GREEN2_TET ) ngreen2++;
        if (i_tet[i] == GREEN3_TET ) ngreen3++;
    }
}

void create_new_elements(Geometry& geom,
                         vector <idx>& i_tet,
                         vector <idx>& i_tri,
                         vector <Line>& lines,
                         vector <set<idx> > t_to_l,
                         vector <set<idx> > e_to_l,
                         vector <double>& new_p,
                         vector <idx>& new_t,
                         vector <idx>& new_mat_t,
                         vector <idx>& new_e,
                         vector <idx>& new_mat_e
                         )
{
    // CREATES NEW NODES, ELEMENTS (TETS + TRIS ) AND MATERIAL
    // VALUES ARRAYS.

    // MAKE NEW NODES FIRST. A SPARSE MATRIX IS USED TO MAP TWO
    // OLD NODES TO A NEW ONE
    //SparseMatrix* nnumbers; <-

    SpaMtrix::IRCMatrix nnumbers = createNodeNumbersMatrix(geom,lines);
    //create_node_number_matrix(nnumbers, geom, lines); <-

    create_new_coordinates( geom, lines, new_p);

    // CREATE NEW TETRAHEDRA ELEMENTS
    for (idx i = 0 ; i < (idx) i_tet.size() ; i++){
        switch (i_tet[i]){
        case (0):   // NON-REFINABLE TET. DO NOTHING
            break;
        case (GREEN1_TET):
            make_new_green1_tet(new_t, new_mat_t, geom, i , lines, t_to_l, nnumbers);
            break;
        case (GREEN2_TET):
            make_new_green2_tet(new_t, new_mat_t, geom, i, lines,t_to_l, nnumbers);
            break;
        case(GREEN3_TET):
            make_new_green3_tet(new_t, new_mat_t , geom, i , lines , t_to_l , nnumbers);
            break;
        case(RED_TET):
            make_new_red_tet(new_t, new_mat_t, geom, i , lines, t_to_l, nnumbers);
            break;
        default:
            printf("error in %s, i_tet[%u] = %u - bye\n",__func__, i , i_tet[i]);
            exit(1);
        }
    }// end for i

    // CREATE NEW TRIANGLE ELEMENTS
    for ( idx i = 0 ; i < (idx) i_tri.size() ; i++){
	if (i_tri[i] == 0){}// do nothing
        else if ( i_tri[i] == 1){
            make_new_tri1( new_e, new_mat_e , geom, i, lines, e_to_l, nnumbers );
        }
        else if ( i_tri[i]== 2){
            make_new_tri2( new_e, new_mat_e , geom, i, lines, e_to_l, nnumbers );
        }
        else if ( i_tri[i] == 3){
            make_new_tri3 ( new_e, new_mat_e, geom , i , lines, e_to_l , nnumbers );
        }
        else{
            printf("error - i_tri[%u] = %u - bye!\n", i , i_tri[i] );
            exit(1);
        }
    }
    //delete nnumbers;
}


