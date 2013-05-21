
#include <geometry.h>

const idx Geometry::NOT_AN_INDEX = std::numeric_limits<idx>::max();



Geometry::Geometry():
    numWeakSurf(0),
    indWeakSurf(NULL),
    regularGrid(NULL)
{
    np 			= 0;
    npLC			= 0;
    p 			= NULL;
    NodeNormals		= NULL;
    t 			= new Mesh();
    e			= new Mesh();
    Xmin			= 0;
    Xmax    		= 0;
    Ymin			= 0;
    Ymax			= 0;
    Zmin			= 0;
    Zmax			= 0;
    t->setDimension(3u);
    e->setDimension(2);

    left_right_is_periodic = false;
    front_back_is_periodic = false;
    top_bottom_is_periodic = false;

}
Geometry::~Geometry()
{
    if (p != NULL){
        free(p);
    }
    if (NodeNormals != NULL){
        free (NodeNormals);
    }
    delete t;
    delete e;

    if ( indWeakSurf != NULL )
        free(indWeakSurf);

    delete regularGrid;

}
void Geometry::setTo(Geometry* geom)
{

    this->ClearGeometry();
    //printf("is equal %i\n" , geom->getnp());
    np 		= geom->getnp();						// number of nodes
    npLC 	= geom->getnpLC();					// number of LC nodes
    Xmin 	= geom->getXmin();
    Xmax	= geom->getXmax();
    Ymin	= geom->getYmin();
    Ymax	= geom->getYmax();
    Zmin	= geom->getZmin();
    Zmax	= geom->getZmax();


    t->CopyMesh(geom->t);
    e->CopyMesh(geom->e);


    if (p!=NULL) free(p);
    p = (double*) malloc(3*getnp()*sizeof(double));

    if (NodeNormals!=NULL) free(NodeNormals);
    NodeNormals = (double*) malloc(3*getnp()*sizeof(double));

    // copy coordinates
    double* temp = geom->getPtrTop();
    for (size_t i = 0 ; i < 3*getnp(); i++)
        p[i] = temp[i];

    // copy node normals
    temp = geom->getPtrToNodeNormals();
    for (size_t i = 0 ; i < 3*getnp(); i ++)
        NodeNormals[i] = temp[i];

    left_right_is_periodic = geom->getleft_right_is_periodic();
    front_back_is_periodic = geom->getfront_back_is_periodic();
    top_bottom_is_periodic = geom->gettop_bottom_is_periodic();

    numWeakSurf = geom->numWeakSurf;
    if (numWeakSurf)
    {
        if (indWeakSurf==NULL) free(indWeakSurf);
        indWeakSurf = (size_t*) malloc( numWeakSurf*sizeof(size_t) );
        memcpy(indWeakSurf , geom->indWeakSurf, numWeakSurf*sizeof(size_t) );

    }

    this->periNodes_.clear();
    periNodes_.insert(periNodes_.end(), geom->periNodes_.begin(), geom->periNodes_.end() );



    if (this->regularGrid) delete regularGrid;
    if ( geom->regularGrid )
        regularGrid = new RegularGrid( *geom->regularGrid );

}
void Geometry::ClearGeometry()
{
    np = 0;
    npLC = 0;
    if (p!=NULL) free(p);
    p = NULL;

    if (NodeNormals != NULL ) free( NodeNormals );
    NodeNormals = NULL;

    Xmin = 0;
    Xmax = 0;
    Ymin = 0;
    Ymax = 0;
    Zmin = 0;
    Zmax = 0;

    left_right_is_periodic = false;
    front_back_is_periodic = false;
    top_bottom_is_periodic = false;

    peri_equ_nodes.clear();	// this is never used anyways ??

    t->ClearMesh();
    e->ClearMesh();


}

void Geometry::setCoordinates(double* coords, const size_t& np)
{
    if (p!=NULL) free(p);
    if (NodeNormals!=NULL) free(NodeNormals);

    p = (double*)malloc(3*np*sizeof(double));
    NodeNormals = (double*)malloc(3*np*sizeof(double));

    if ( (p == NULL) || (NodeNormals == NULL))
    {
        printf("error - Geometry::setCoordintes - could not allocate memtory - bye!\n");
        exit(1);
    }
    memset(NodeNormals, 0, 3 * np * sizeof(double)); //reset node normals to all zero

    Xmin = 1e9; Xmax = -1e9;
    Ymin = 1e9; Ymax = -1e9;
    Zmin = 1e9; Zmax = -1e9;

    for (size_t i = 0 ; i < np ; i ++) // copy coordinates
    {
        p[i*3+0] = coords[i*3+0];
        p[i*3+1] = coords[i*3+1];
        p[i*3+2] = coords[i*3+2];

        if (p[i*3+0] < Xmin) Xmin = p[i*3+0]; 	// find xmin
        if (p[i*3+0] > Xmax) Xmax = p[i*3+0];
        if (p[i*3+1] < Ymin) Ymin = p[i*3+1];	// find ymin
        if (p[i*3+1] > Ymax) Ymax = p[i*3+1];
        if (p[i*3+2] < Zmin) Zmin = p[i*3+2];	// find zmin
        if (p[i*3+2] > Zmax) Zmax = p[i*3+2];

    }

    setnp(np);
    setnpLC(np);
}
void Geometry::addCoordinates(double* coords, const size_t& npn){
    // APPENDS npn NEW COORDIANTES IN coords TO this->p
    if ( ( npn > 0 ) && (coords) ){ // don't do anything if no new coordinates
        size_t np_new = np + npn;
	double* pnew = (double*) malloc( 3*np_new*sizeof(double));
	
	// first copy old coordinates
        for(size_t i = 0 ; i < 3*np ; i++)
            pnew[i] = p[i];
	// then add new ones
        for(size_t i = 0 ; i < 3* npn ; i ++)
            pnew[3*np+i] = coords[i];
	
	// set p to point to new extended cordinates
	if (p!=NULL) free(p); 
	p = pnew;
	
	np = np_new;
	setnpLC(np_new); /// needs reordering after this
    }// end if no new coords
}
void Geometry::addCoordinates(vector<double> &coords){
    if (coords.size()>0){
        int np_new = np + coords.size() / 3 ;
        double* pnew = (double*) malloc(3*np_new*sizeof(double) ); // allocate space for new coodinates

        // first copy odl coordinates
        memcpy(pnew , p , 3*np*sizeof(double) );

        // then add new ones. would memcpy work here?
        for (size_t i = 0; i < coords.size() ; i++)
            pnew[3*np+i] = coords[i];

        // update pointers
        if (p) free(p);
        p = pnew;

        np = np_new;

        setnpLC( np_new ); // needs reordering after this.

    }
}


void Geometry::updateMaxNodeNumbers()
{
// UPDATE MaxNodeNumber FOR SURFACE AND VOLUME MESHES
// THIS NEEDS TO BE DONE AFTER NODE REOREDERING

// VOLUMES
    idx* ielem = this->t->getPtrToElement(0);
    idx maxnn = *max_element(ielem, ielem+ t->getnElements()*t->getnNodes() );
    t->setMaxNodeNumber( maxnn );

// SURFACES
    ielem = e->getPtrToElement(0);
    maxnn = *max_element(ielem, ielem+e->getnElements()*e->getnNodes() );
    e->setMaxNodeNumber( maxnn );
}

void Geometry::setnp(int n)		{np = n;}
void Geometry::setnpLC(int n)	{npLC = n;}
void Geometry::setNodeNormals()
{
    // CALCULATES NODE NORMALS BASED ON CONNECTED TRIANGLE NORMALS
    if (e->getnElements() == 0)
    {
        printf("error - Geometry::setNodeNormals() - surface mesh not defined, bye!\n");
        exit(1);
    }

    if (NodeNormals!=NULL) free(NodeNormals);						// ALLOCATE
    NodeNormals = (double*)malloc(3 * np * sizeof(double));

    memset(NodeNormals, 0, 3 * np * sizeof(double));

    double tempn[3] ={0,0,0};

    for (idx i = 0 ; i < e->getnElements() ; i ++) // add neighbouring surface normals
    {
        int m = e->getMaterialNumber(i);
        //if ((m != MAT_PERIODIC) && (m!= MAT_NEUMANN))//if element is not periodic or neuman surface <---SHOULD THIS BE FIXLC ONLY?? ARE NODE NORMALS USED BESIDES IN ANCHORING?
        if ( MATNUM_TO_FIXLC_NUMBER(m) ) // IF FIXLC SURFACE
        {
            e->CopySurfaceNormal(i,&tempn[0]);	// copy surface triangle normal to temp normal

            for (idx j = 0; j < e->getnNodes() ; j++)
            {
                NodeNormals[e->getNode(i,j)*3 + 0] += tempn[0]; // add x, y, z components for each node in triangle i
                NodeNormals[e->getNode(i,j)*3 + 1] += tempn[1];
                NodeNormals[e->getNode(i,j)*3 + 2] += tempn[2];
            }// end for j

            tempn[0] = 0; tempn[1] = 0; tempn[2] = 0;
        }// end if not periodic
    }// end for i

    for (size_t i = 0 ; i < np ; i ++) // normalise length
    {
        double len = sqrt( NodeNormals[i*3+0]*NodeNormals[i*3+0] +
                           NodeNormals[i*3+1]*NodeNormals[i*3+1] +
                           NodeNormals[i*3+2]*NodeNormals[i*3+2] );
        if (len>0)
        {
            NodeNormals[i*3+0] = NodeNormals[i*3+0] / len;
            NodeNormals[i*3+1] = NodeNormals[i*3+1] / len;
            NodeNormals[i*3+2] = NodeNormals[i*3+2] / len;
        }
	
    } // end normalise loop, i


}
void Geometry::PrintNodeNormals()
{
    printf("going to print %i node normals\n",np);
    fflush(stdout);
    for (size_t i = 0 ; i < np ; i ++)
    {
        printf("NodeNormals[%u] = [%f,%f,%f]\n",(idx) i, NodeNormals[i*3], NodeNormals[i*3+1], NodeNormals[i*3+2]);
        fflush(stdout);
    }

}

void Geometry::ReorderDielectricNodes()
{
    if ((t->getnElements() <= 0 ) || ( getnp() <=0) || (t->getPtrToMaterialNumber(0) == NULL))
    {
        printf("error - Geometry::ReorderDielectricNodes - geometry not fully defined yet - bye!\n");
        exit(1);
    }


    // check if dielectric elements exist
    bool DE_exist = false;
    //int npLC = np;

    idx* tmat = t->getPtrToMaterialNumber(0);
    for (idx i = 0 ; i < t->getnElements() ; i++)
        if (tmat[i]>= MAT_DIELECTRIC1)	// if material nuber > LC number
        {
            DE_exist = true;
            break;
        }

    setnpLC( getnp() );
    if (!DE_exist) // if no dielectric materials -> no need to reorder = exit
    {
        this->updateMaxNodeNumbers();
        return;
    }


    //
    //GENERATE LIST OF MATERIAL NUMBERS FOR NODES. IN CASE OF DUAL VALUE LC PRECEDES
    //

    // mark all LC nodes as 1 and others as 0
    idx* lcde = (idx*) malloc(getnp() * sizeof(idx) );
    memset(lcde , 0 , getnp() * sizeof(int) ); // start with everything 0

    for (idx i = 0 ; i < t->getnElements() ; i ++ )
    {
        if (t->getMaterialNumber(i) == MAT_DOMAIN1)
        {
            for (idx j = 0 ; j < t->getnNodes() ; j ++ ) // loop over all nodes
                lcde[t->getNode(i , j)] = 1; // LC --> 1
        }
    }

    npLC = 0;
    vector <int> v_mat_index;


    for (size_t i=0; i< getnp() ;i++)// first add all nodes marked as LC
        if (lcde[i] == 1)
        {
            v_mat_index.push_back(i);
            npLC++;
        }

    for (size_t i=0; i< getnp() ;i++) // then add all non-LC nodes
        if (lcde[i] == 0 )
            v_mat_index.push_back(i);

    //printf("np = %i, npLC = %i\n",np, npLC);
    free(lcde);

    //make inverse map
    vector <int> v_invmap;
    v_invmap.resize( v_mat_index.size() , -1);
    for (size_t i = 0 ; i < getnp() ; i ++)
        v_invmap[v_mat_index[i]]= i;//v_mat_index[i];//= i;


    //for( int i  = 0 ; i < getnp() ; i ++ )
    //printf("%i <->  %i\n", v_invmap[i] , v_mat_index[i] );


    //REORDER NODES
    double* newp = (double*)malloc( 3*getnp() *sizeof(double)); // memory for reordered node coordinates
    for (size_t i = 0 ; i < getnp() ; i++)
    {
        newp[i*3+0]=getpX( v_mat_index[i] ); //x-coord
        newp[i*3+1]=getpY( v_mat_index[i] ); //y-coord
        newp[i*3+2]=getpZ( v_mat_index[i] ); //z-coord
    }

    free(p); // make p = new reordered p
    p=newp;

    //REORDER TETRAHEDRA
    idx* newt = (idx*)malloc(t->getnElements() * t->getnNodes() * sizeof(idx));
    for (idx i = 0 ; i < t->getnElements() ; i++ )
    {
        for (idx j = 0 ; j < t->getnNodes() ; j ++ ) // loop over all nodes of element i
            newt[i*t->getnNodes() + j ] = v_invmap[ t->getNode(i , j) ];
    }
    t->setAllNodes(newt ); // copy node numbers
    free(newt);
    //exit(1);

    //REORDER TRIANGLES
    idx* newe = (idx*)malloc(e->getnElements() * e->getnNodes() * sizeof(idx));

    for(idx i = 0 ; i < e->getnElements() ; i++)
    {
        for (idx j = 0 ; j < e->getnNodes() ; j++ )
            newe[i * e->getnNodes() + j ] = v_invmap[ e->getNode(i,j) ];
    }
    e->setAllNodes(newe ); // copy node numbers
    free(newe);


    this->updateMaxNodeNumbers();

}




void Geometry::setFacePeriNodes(list<size_t> &face0,
                                list<size_t> &face1,
                                const int &norm)
{
    // SETS INDEX VALUES IN periNodes_ FOR FACES
    // norm is face normal, need coordinates to vectors parallel to it
    // norm = 0 -> COMPARE Y,Z
    // norm = 1 -> COMPARE X,Z
    // notm = 2 -> COMPARE X,Y
    int ind1 = ( norm + 1 ) % 3; // ind is a pre-calculated offeset to coordinate comparison in p.
    int ind2 = ( norm + 2 ) % 3; // Selects comparison of x,y or z coordinate. e.g. norm = 0 -> ind1 = 1, ind2 = 2
    double eps = 1e-5; // accuracy of coordinate comparison
    // SEARCH FOR NODE EQUIVALENCIES BY COMPARING COORDINATES
    // THAT ARE PERPENDICULAR TO FACE NORMAL
    list <size_t>:: iterator F0;    // FACE 0 NODES
    list <size_t>:: iterator F1;    // FACE 1 NODES
    int fc, bc; // debug counters
    for (F0 = face0.begin(), fc = 0; F0!=face0.end() ; ++F0, ++fc ) // LOOP OVER FACE 0
    {
        bool found = false;
#ifdef DEBUG
        int ind_n  = 0;     // index to neares (debug)
        double dist = 1000000;
#endif
        double f1 = p[3*(*F0) + ind1 ]; // coordinates of node F2 in face0
        double f2 = p[3*(*F0) + ind2 ];

        for (F1 = face1.begin() , bc = 0; F1!= face1.end() ; ++F1, ++bc)
        {
            double fa = p[3*(*F1) + ind1]; // coordinates of node F1 in face 1
            double fb = p[3*(*F1) + ind2];

            //compare coordinates
            double dist1 = fabs(f1 - fa); // distances in plane
            double dist2 = fabs(f2 - fb);
            double tdist = dist1*dist1 + dist2*dist2;
#ifdef DEBUG
            if (tdist < dist) // debug info only, keep track of nearest found node
            {
                dist = tdist; // nearest distance
                ind_n = *F1;  // index to nearest distance
            }
#endif
            if (tdist < eps*eps) // compare squared distances
            {
                this->periNodes_[*F1] =*F0;
                found = true;
                break;
            }
        }// end for B
        if (!found)
        {
            printf("error, periodic boundaries do not seem to match - bye!\n");
#ifdef DEBUG

            double* c0 = getPtrTop() + *F0;  // coordinates of face 0 coords
            double* c1 = getPtrTop() + ind_n;  // face 1 coords

            printf("error - in file %s, function %s.\n", __FILE__, __func__);
            cout << "normal = " << norm << endl;
            cout << "coordinate "<< *F0 <<" FACE0: "<< *c0 <<","<< *(c0+1) <<","<< *(c0+2) << endl;
            cout << "nearest match in FACE1: index "<<ind_n<<" , dist = " <<dist << endl;
            cout << "coordinate: "<< *c1 << "," << *(c1+1) << "," << *(c1+2) << endl;
#endif



            exit(1);
        }
    }//end for F
}


void Geometry::setEdgePeriNodes(list<size_t> &edge0,
                                list<size_t> &edge1,
                                const int &dim)
{

    const double BIGNUM = 1e99;
    list<size_t> :: iterator e0;
    list<size_t> :: iterator e1;
    double eps = 1e-5;

    // LOOP OVER EACH NODE IN EDGE 0 AND FIND ITS EQUIVALENT IN OTHER EDGES
    for (e0 = edge0.begin() ; e0 != edge0.end() ; ++e0 )
    { // outer corner loop - corn0
        bool found = false;
        double dist = 0;

        double p0 = p[(*e0) * 3 + dim ]; // dim coordinate of node c0
        int iNearest = -1;
        double minDist = BIGNUM;

        // COMPARISON WITH C1
        for (e1 = edge1.begin() ; e1 != edge1.end() ; ++e1)
        { // inner corner loop - corn1
            double p1 = p[(*e1) * 3 + dim];
            dist = fabs( p0 - p1 );

            if ( dist <= minDist )
            {
                minDist = dist;
                iNearest = *e1;
            }


            if ( dist <= eps ) // if same coordinate -> match!
            {

                periNodes_[*e1] = *e0;
                found = true;
                break;
            }
        }
        if (!found)
        {

            printf("error - edge node 1 not found - bye\n");
            printf("edge node 0 = %u,\tat [%e,%e,%e]\n", (unsigned int)*e0, p[*e0*3 + 0], p[*e0*3 + 1], p[*e0*3 + 2] );
            printf("nearest node  = %u,\tat [%e,%e,%e]\n", (unsigned int)iNearest, p[iNearest*3 + 0], p[iNearest*3 + 1], p[iNearest*3 + 2]);
            printf("distance = %e\n", minDist);
            exit(1);
        }
    }// end loop over all nodes in edge0
}


void Geometry::makePeriEquNodes()
{
    // GENERATES INDEX OF PERIODIC NODE EQUIVALENCIES

    // LIST OF INDEXES TO ALL PERIODIC NODES
    vector <unsigned int> nodes;
    e->listNodesOfMaterial( nodes , MAT_PERIODIC );

    size_t numPeri = nodes.size();
    if ( !numPeri )
    {
        periNodes_.clear();
        return;
    }
    periNodes_.clear();
    periNodes_.reserve( getnp() );
    for ( size_t i = 0 ; i < getnp() ; i++ )
        periNodes_.push_back(i);


    // NEED TO CONSIDER 3 DIFFERENT CASES, DEPENDING ON TYPE OF PERIODICITY OF MODELLING WINDOW
    double eps = 1e-5; // accuracy for coordinate comparisons
    double xmin = getXmin(); // convenience shortcuts to geometry min and max dimensions
    double xmax = getXmax();
    double ymin = getYmin();
    double ymax = getYmax();
    double zmin = getZmin();
    double zmax = getZmax();
    /// PROBABLY EVIL, BUT SO CONVENIENT...
#define LEFT    ( getAbsXDist(n, xmin) <= eps )
#define RIGHT   ( getAbsXDist(n, xmax) <= eps )
#define FRONT   ( getAbsYDist(n, ymin) <= eps )
#define BACK    ( getAbsYDist(n, ymax) <= eps )
#define BOTTOM  ( getAbsZDist(n, zmin) <= eps )
#define TOP     ( getAbsZDist(n, zmax) <= eps )

    //CASE 1 FRONT BACK IS PERIODIC ONLY
    if (getfront_back_is_periodic() &&
            !getleft_right_is_periodic() &&
            !gettop_bottom_is_periodic() )
    {

        //SEPARATE NODES INTO TWO LISTS FOR FRONT AND BACK SURFACES
        list <size_t> front;
        list <size_t> back;

        for (size_t i = 0 ; i < nodes.size() ; i ++ )
        {
            unsigned int n = nodes[i];

            if ( FRONT ) // check if node i is on front surface
            { front.push_back(n) ; }
            else
                if ( BACK) // check if node i is on back surface
                { back.push_back(n); }

                else // ERROR
                {
                    printf("error in %s, CASE 1 - bye \n", __func__ );
                    exit(1);
                }
        }//end for i


        /// MAKE SURE EQUAL NUMBER OF NODES HAVE BEEN FOUND ON BOTH SURFACES
        if (front.size() != back.size() ){
            printf("error - SolutionVector::setPeriosdicEquNodes\n");
            printf("front and back surfaces do not have same number of nodes\n");
            printf("front = %i , back = %i - bye!\n " , (int) front.size() , (int) back.size() );
            exit(1);
        }

        // SEARCH FOR NODE EQUIVALENCIES BY COMPARING X AND Z COORDINATES
        // BACK NODES MAP TO FRONT NODES
        this->setFacePeriNodes( front, back, 1 );
        // setFaceElim( front, back, Elim, 1, geom->getPtrTop() );
    }

    // CASE 2 FRONT-BACK AND LEFT-RIGHT ARE PERIODIC
    else if (    getfront_back_is_periodic() &&
                 getleft_right_is_periodic() &&
                 !gettop_bottom_is_periodic() )
    {

        // separate nodes into 8 lists, 4 x corners left/right and front/back planes
        list <size_t> edge0; //x = 0, y = 0
        list <size_t> edge1; //x = 0, y = max
        list <size_t> edge2; //x = max, y = max
        list <size_t> edge3; //x = max, y = 0
        list <size_t> front; //x = 0
        list <size_t> back;  //x = max
        list <size_t> right; //y = max
        list <size_t> left;  //y = 0

        for (size_t i = 0 ; i < nodes.size() ; i++) // loop over all nodes and insert to correct list
        {
            int n = nodes[i];
            if ( LEFT && FRONT ) // edge0
            {edge0.push_back(n);}
            else if ( LEFT && BACK ) // edge1
            {edge1.push_back(n);}
            else if ( RIGHT  &&  BACK ) // edge2
            {edge2.push_back(n);}
            else if ( RIGHT && FRONT ) // edge3
            {edge3.push_back(n);}
            else if ( FRONT ) // front surface
            {front.push_back(n);}
            else if ( BACK ) // back surface
            {back.push_back(n);}
            else if ( LEFT ) // left surface
            {left.push_back(n);}
            else if ( RIGHT ) // right surface
            {right.push_back(n);}
        }

        if ( ( edge0.size() != edge1.size() ) ||
             ( edge0.size() != edge2.size() )  ||
             ( edge0.size() != edge3.size() )  ||
             ( left.size() != right.size() ) ||
             ( front.size() != back.size() ) )
        {
            printf("error, different number of periodic boundary nodes\n");
            printf("corners 0,1,2,3 = %i,%i,%i,%i\n", (int) edge0.size() , (int) edge1.size() , (int) edge2.size() , (int) edge3.size() );
            printf("front/back , left/right = %i/%i , %i/%i \n", (int) front.size() , (int) back.size() , (int) left.size() , (int) right.size() );
            exit(1);
        }
        setEdgePeriNodes(edge0, edge1, 2 );
        setEdgePeriNodes(edge0, edge2, 2 );
        setEdgePeriNodes(edge0, edge3, 2 );
        setFacePeriNodes(left  , right, 0); // COMPARE Y, Z
        setFacePeriNodes(front , back , 1 ); // COMPARE X, Z


    }// END CASE 2

    // CASE 3 FRONT-BACK, LEFT-RIGHT AND TOP-BOTTOM ARE PERIODIC
    else if ( getfront_back_is_periodic() &&
              getleft_right_is_periodic() &&
              gettop_bottom_is_periodic() )

    {
        // separate nodes into lists, 12 x edges left/right, front/back and top/bottom planes
        // Vertical corners along Z

        // Additionally, 7 corner nodes must point to bottom left (origin xmin,ymin,zmin) corner

        list <size_t> edge0; //x = 0, y = 0
        list <size_t> edge1; //x = 0, y = max
        list <size_t> edge2; //x = max, y = max
        list <size_t> edge3; //x = max, y = 0

        // Horizontal edges along X
        list <size_t> edgea; // y = 0, z = 0
        list <size_t> edgeb; // y = max, z = 0
        list <size_t> edgec; // y = max, z = max
        list <size_t> edged; // y = 0, z = max

        // Horizontal edges along Y
        list <size_t> edgeA; // x = 0, z = 0
        list <size_t> edgeB; // x = max, z = 0
        list <size_t> edgeC; // x = max, z = max
        list <size_t> edgeD; // x = 0, z = max

        list <size_t> front; //x = 0
        list <size_t> back;  //x = max
        list <size_t> right; //y = max
        list <size_t> left;  //y = 0
        list <size_t> top;   //z = max
        list <size_t> bottom;//z = 0;

        // LOOP OVER ALL NODES AND INSERT TO CORRECT LIST
        int corner_nodes[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
        for (size_t i = 0 ; i < nodes.size() ; i++)
        {
            int n = nodes[i];

            // CORNER NODES TAKE PRECEDENCE OVER OTHER NODES
            // FRONT LEFT BOTTOM
            if ( FRONT && LEFT && BOTTOM )
                corner_nodes[0] = n;
            // FRONT RIGHT BOTTOM
            else if ( FRONT && RIGHT && BOTTOM )
                corner_nodes[1] = n;
            else if ( FRONT && LEFT && TOP )  /// FRONT LEFT TOP
                corner_nodes[2] = n;
            else if (FRONT && RIGHT && TOP)   /// FRONT RIGHT TOP
                corner_nodes[3] = n;
            else  if (BACK && LEFT && BOTTOM)
                corner_nodes[4] = n;
            else  if (BACK && RIGHT && BOTTOM)
                corner_nodes[5] = n;
            else  if (BACK && LEFT && TOP)
                corner_nodes[6] = n;
            else  if (BACK && RIGHT && TOP)
                corner_nodes[7] = n;

            // EDGE NODES
            // 4 x Vertical Corners
            else if ( LEFT && FRONT ) // edge0
            {edge0.push_back(n);}
            else if ( LEFT  && BACK ) // edge1
            {edge1.push_back(n);}
            else if (RIGHT  && BACK ) // edge2
            {edge2.push_back(n);}
            else if ( RIGHT && FRONT ) // edge3
            {edge3.push_back(n);}
            // 4 x Horizontal along X
            else if ( FRONT && BOTTOM)	// ymin and zmin
            {edgea.push_back(n);}
            else if (BACK && BOTTOM) // ymax and zmin
            {edgeb.push_back(n);}
            else if (BACK && TOP) // ymax and zmax
            {edgec.push_back(n);}
            else if (FRONT && TOP)  // ymin and zmax
            {edged.push_back(n);}
            // 4 x Horizontal along Y
            else if( LEFT && BOTTOM)
            {edgeA.push_back(n);} // xmin and zmin
            else if( RIGHT && BOTTOM)
            {edgeB.push_back(n);} // xmax and zmin
            else if( RIGHT && TOP )
            {edgeC.push_back(n);} // xmax and zmax
            else if( LEFT && TOP)
            {edgeD.push_back(n);} // xmin and zmax
            // FRONT/BACK, LEFT/RIGHT, TOP/BOTTOM FACES
            else if ( FRONT ) // front surface
            {front.push_back(n);}
            else if ( BACK ) // back surface
            {back.push_back(n);}
            else if ( LEFT ) // left surface
            {left.push_back(n);}
            else if ( RIGHT ) // right surface
            {right.push_back(n);}
            else if ( BOTTOM ) // bottom surface
            {bottom.push_back(n); }
            else if ( TOP ) // top surface
            {top.push_back(n); }

            else
            {
                printf("error in %s, periodic node not on external surface - bye!\n",__func__);
                exit(1);
            }

        }// end for i, loop over all nodes

        // CHECK THAT OPPOSITE FACES HAVE EQUAL NUMBER OF NODES
        {// start dummy scope

            int mincorner = *min_element( corner_nodes, corner_nodes+8);
            if (mincorner < 0)
            {
                printf(" error - corner nodes not found - bye!\n");
                printf("indexes are = [%i,%i,%i,%i,%i,%i,%i,%i]\n", corner_nodes[0],corner_nodes[1],corner_nodes[2],corner_nodes[3],
                       corner_nodes[4],corner_nodes[5],corner_nodes[6],corner_nodes[7]);
                exit(1);

            }


            if (top.size() != bottom.size() )
            {
                printf("error - top and bottom surfaces do not match - bye!\n");
                printf("sizes are top,bottom = %i,%i\n", (int) top.size(), (int) bottom.size() );
                exit(1);
            }
            if (left.size() != right.size() )
            {
                printf("error - left and right surfaces do not match - bye!\n");
                exit(1);
            }
            if (front.size() != back.size() )
            {
                printf("error - front and back surfaces do not match - bye!\n");
                exit(1);
            }
            // CHECK ALL CORNERS HAVE CORRECT NUMBER OF NODES
            size_t s0, s1, s2, s3;
            s0 = edge0.size(); s1 = edge1.size(); s2 = edge2.size() ; s3 = edge3.size();
            if ( (s1!=s0) || (s2 != s0) || (s3!=s0) )
            {
                printf("error - vertical corner node counts do not match\n");
                exit(1);
            }
            s0 = edgea.size(); s1 = edgeb.size(); s2 = edgec.size(); s3 = edged.size();
            if ( (s1!=s0) || (s2 != s0) || (s3!=s0) )
            {
                printf("error - horizontal corner (along x) node counts do not match\n");
                exit(1);
            }
            s0 = edgeA.size(); s1 = edgeB.size(); s2 = edgeC.size(); s3 = edgeD.size();
            if ( (s1!=s0) || (s2 != s0) || (s3!=s0) )
            {
                printf("error - horizontal corner (along y) node counts do not match\n");
                exit(1);
            }
        }// end dummy scope

        // set corner nodes
        periNodes_[corner_nodes[1] ] = corner_nodes[0];
        periNodes_[corner_nodes[2] ] = corner_nodes[0];
        periNodes_[corner_nodes[3] ] = corner_nodes[0];
        periNodes_[corner_nodes[4] ] = corner_nodes[0];
        periNodes_[corner_nodes[5] ] = corner_nodes[0];
        periNodes_[corner_nodes[6] ] = corner_nodes[0];
        periNodes_[corner_nodes[7] ] = corner_nodes[0];

        // match edge nodes
        setEdgePeriNodes(edge0, edge1 , 2 );  // VERTICAL EDGES
        setEdgePeriNodes(edge0, edge2 , 2 );
        setEdgePeriNodes(edge0, edge3 , 2 );

        setEdgePeriNodes(edgea, edgeb, 0); // HORIZONTAL ALONG X
        setEdgePeriNodes(edgea, edgec, 0);
        setEdgePeriNodes(edgea, edged, 0);

        setEdgePeriNodes(edgeA, edgeB, 1);  // HORIZONTAL ALONG Y
        setEdgePeriNodes(edgeA, edgeC, 1);
        setEdgePeriNodes(edgeA, edgeD, 1);

        setFacePeriNodes(left, right, 0);   // COMPARE Y,Z
        setFacePeriNodes(front, back, 1);   // COMPARE X,Z
        setFacePeriNodes(bottom, top, 2);   // COMPARE X,Y

    }// end if 3 different periodicity cases


}// end void MakePEriEquNodes()

void Geometry::checkForPeriodicGeometry()
{
    // CHECKS FOR TYPES OF PERIODICITY PRESENT IN CURRENT STRUCTURE.
    // POSSIBLE PERIODIC SURFACES ARE:
    //      LEFT/RIGHT
    //      FRONT/BACK
    //      TOP/BOTTOM

    for (idx i = 0 ; i < e->getnElements() ; i++)
    {
        if (e->getMaterialNumber(i) == MAT_PERIODIC) // if surface is periodic
        {
            // check to see which side surface is on by looking at the surface normal
            double* snorm = e->getPtrToSurfaceNormal(i);

            // IF SURFACE NORMAL X-COMPONENT = 1
            if (fabs( fabs(snorm[0]) - 1.0 ) < EPS)
            {
                left_right_is_periodic = true;
            }
            else// IF SURFACE NORMAL Y-COMPONENT = 1
                if (fabs( fabs(snorm[1]) - 1.0 ) < EPS)
                {
                    front_back_is_periodic = true;
                }
                else // IF SURFACE NORMAL Z-COMPONENT = 1
                    if (fabs( fabs(snorm[2]) - 1.0 ) < EPS)
                    {
                        top_bottom_is_periodic = true;
                    }
                    else
                    {
                        printf("error - checkForPeriodicGeometry() - periodic surface element %i has invalid normal:\n [%f, %f, %f] - bye!\n",i, snorm[0] , snorm[1], snorm[2]);
                        exit(1);
                    }

            // IF ALL SURFACES HAVE ALREADY BEEN IDENTIFIED AS PERIODIC
            // NO NEED TO CHECK FURTHER TRIANGLES
            if ((getleft_right_is_periodic()) &&
                    (getfront_back_is_periodic()) &&
                    (gettop_bottom_is_periodic() ))
                break;
        }// end if periodic surface
    }// end for i
    // IF ANY PERIODIC TRIANGLES WERE DETECTED
    if (    getleft_right_is_periodic()
            ||  getfront_back_is_periodic()
            ||  gettop_bottom_is_periodic() )
    {
        makePeriEquNodes();
    }

}


void Geometry::makeRegularGrid(const size_t &nx,
                               const size_t &ny,
                               const size_t &nz)
{
    // CREATES REUGLAR GRID OBJECT USED FOR INTERPOLATING VALUES FROM
    // TETRAHEDRAL MESH ONTO REGULARLY SPACED GRID

    printf("generating regular grid lookup..."); fflush(stdout);
    if ( regularGrid ) delete regularGrid;
    regularGrid = new RegularGrid();

    if ( !nx || !ny || !nz )    // IF NO REGULAR GRID SIZE IS DEFINED. DONT MAKE ONE
        return;

    regularGrid->createFromTetMesh(nx, ny, nz, this );

    printf("OK\n"); fflush(stdout);
}

bool Geometry::brute_force_search( unsigned int &ind,            // return index
                                   double* coord,                // search coordinate
                                   const bool& terminateOnError, // terminate if not found>
                                   const bool& requireLCEelement // only LC element index may be returned
                                   )
{
    // BRUTE FORCE DEBUG SEARCH FOR TETRAHEFRON THAN CONTAINS POINT WITH COORDINATES IN coord
    // coord IS ASSUMED TO BE OF LENGTH 3, FOR x, y, z
    // loop over each element
    for (idx i = 0 ; i < t->getnElements() ; i++)
    {
        if  ( t->ContainsCoordinate( i , getPtrTop(), coord ) ) // If coord is in tet i
        {
            if (requireLCEelement) // WANT LC
            {
                if (t->getMaterialNumber( i ) <= MAT_DOMAIN7 ) // IF LC
                {
                    ind = i;
                    return true;
                }
            }
            else // DON'T CARE WHETHER LC OR DE
            {
                ind = i ;
                return true; // exit function when found
            }
        }
    }// end for loop over all elems

    // IF COORDINATE WAS NOT FOUND APPLICATION MAY NEED TO BE TERMINATED
    if (terminateOnError)
    {
        printf("error - brute_force_search could not find coord %f,%f,%f - bye ! \n", coord[0], coord[1], coord[2] );
        exit(1);
    }
    // SIGNAL A NON-FOUND COORDINATE BY RETURNING FALSE
    return false;
}

bool Geometry::getContainingTet(vector< set < unsigned int> >& p_to_t,
                                double crd[3],  // coords to search
                                unsigned int& t0) // start tet
{
    // TRIENS TO FIND INDEX TO TETRAHEDRON THAT CONTAINS THE POINT
    // WITH CORRDINATES DEFINED IN crd. USING A NEIGHBOUR SEARCH.
    // STARTS FROM INITIAL TET t0, AND ADVANCES TO ITS NEIGHBOUR WHOSE
    // BARYCENTRE IS CLOSEST TO crd.
    // SOMETIMES THIS FAILS AND A BRUTE FORCE SEARCH MAY BE NEEDED



    bool found = false;
    unsigned int t_prev = t->getnElements();
    for (size_t i = 0 ; i < (size_t) t->getnElements() ; i++) // for all elements
    {

        if  (  t->ContainsCoordinate(t0, p , crd)  ) // if element t0 contains coordinate
        {
            found = true;
            break;
        }
        set<unsigned int> :: iterator nitr;
        set<unsigned int> setneighs;
        vector <unsigned int> vecneighs;
        // make UNIQUE list of all tets sharing nodes with t0
        setneighs.insert( p_to_t[ t->getNode(t0, 0 ) ].begin(), p_to_t[ t->getNode(t0, 0 ) ].end() );
        setneighs.insert( p_to_t[ t->getNode(t0, 1 ) ].begin(), p_to_t[ t->getNode(t0, 1 ) ].end() );
        setneighs.insert( p_to_t[ t->getNode(t0, 2 ) ].begin(), p_to_t[ t->getNode(t0, 2 ) ].end() );
        setneighs.insert( p_to_t[ t->getNode(t0, 3 ) ].begin(), p_to_t[ t->getNode(t0, 3 ) ].end() );

        // MAKE RANDOM ACCESS COPY OF UNIQUE NEIGHBOUR LIST
        vecneighs.insert( vecneighs.end(), setneighs.begin(), setneighs.end() );

        // CALCULATE DISTANCES FOR EACH NEIGBOUR
        vector <double > dst;
        for (unsigned int d = 0 ; d<vecneighs.size() ; d++)
        {
            dst.push_back( t->CalcBaryDistSqr(p , vecneighs[d], crd));
        }

        // FIND INDEX TO MINIMUM DISTANCE AND SET t0 TO THAT
        unsigned int ind_min = min_element(dst.begin(), dst.end() ) - dst.begin();
        t_prev = t0;
        t0 = vecneighs[ind_min];

        // SEARCH MAY GET STUCK WHEN POINT IS CLOSE TO A CORNER NODE OF A LONG THIN
        // TET CONTAINING THE COORD. IN THAT CASE, THE BARYCENTRE OF A NEIGHBORING
        // TET MAY BE CLOSER -> SEARCH WILL NOT FIND THE ACTUAL TET
        if (t_prev == t0) // IF STUCK, TRY ALL NEIGHBOURS IN ORDER OF DISTANCE
        {
            dst[ind_min] = DBL_MAX; // SET CURRENT t0 TO MAXIMUM POSSIBLE DISTANCE AS DEFINED IN <float.h>

            for (size_t ind = 0 ; ind < vecneighs.size() ; ind++) // LOOP OVER ALL NEIGHBOURS
            {
                ind_min = min_element(dst.begin(), dst.end() ) - dst.begin(); // refind new closest value
                if ( t->ContainsCoordinate( vecneighs[ind_min], p, crd ) )    // If contains coordinate
                {
                    found = true;
                    t0 = vecneighs[ind_min];
                    break;
                }
                else
                {
                    dst[ind_min] = DBL_MAX; // STILL NOT FOUND?, TRY AGAIN WITH TET THAT IS NEXT CLOSEST
                }
            }// end for all neighbours
            // IF REALLY STUCK EXIT AND TEST USING BRUTE FORCE
            break;
        } // end if stuck

    }// end for all tets
    return found;
}


size_t Geometry::recursive_neighbour_search(double crd[3],
                                            const vector<set<unsigned int> > &p_to_t,
                                            const size_t &currentTet,
                                            std::set<size_t>& tetHistory,
                                            const bool& requireLCElement    // only LC element index can be returned
                                            )
{
    // TRIES TO FIND TETRAHEDRON CONTAINING POINT crd bY RECURSIVELY
    // SELECTING NEIGHBOUR TET WHOSE BARYHENTRE IS NEARES TO crd

    if ( t->ContainsCoordinate(currentTet, p, crd) )
    {
        if (!requireLCElement)  // if not worried about whether LC or DE element
        {
            return currentTet;
        }
        else // LC element is required
        {
            if (  this->t->getMaterialNumber(currentTet) <= MAT_DOMAIN7 ) // if LC element
            {
                return currentTet;
            }
        }

    }
    tetHistory.insert( currentTet ); // history should be used to avoid visiting same element multiple times (this isn't implemented yet)

// RECURSION LIMIT. UGLY SOLUTION TO STOP OUT OF MEMORY (?) CRASH FOR LARGE GRIDS AND/OR MESHES
    int RECURSION_LIMIT = 1000;
    if (tetHistory.size() > RECURSION_LIMIT ) return NOT_AN_INDEX;



    // CREATE LIST OF NEIGHBOURING ELEMENTS
    std::vector <unsigned int> neighs;
    for (idx i = 0 ; i < t->getnNodes() ; i++)
    {
        int n = t->getNode(currentTet, i);
        neighs.insert(neighs.end(), p_to_t[n].begin(), p_to_t[n].end() );
    }

    sort(neighs.begin(), neighs.end() );
    vector<unsigned int>::iterator itr = unique( neighs.begin(),
                                                 neighs.end() );
    neighs.erase( itr, neighs.end() );
    // REMOVE REFERENCE TO SELF TOO....

    vector <double> dists;
    for(size_t i = 0 ; i < neighs.size() ; i++)
    {
        dists.push_back( t->CalcBaryDistSqr(p, neighs[i], crd) );
    }

    // FIND INDEX TO NEIGHBOUR THAT IS NEAREST
    size_t indn = min_element( dists.begin(), dists.end() ) - dists.begin();


    while ( dists[indn] < DBL_MAX )
    {
        size_t indt = neighs[indn]; // actual element number
        //double d = dists[indn];

        size_t ifound = NOT_AN_INDEX;

        if ( tetHistory.find(indt) == tetHistory.end() )
        {
            ifound =  recursive_neighbour_search( crd,
                                                  p_to_t,
                                                  indt,
                                                  tetHistory,
                                                  requireLCElement // WHETHER ONLY LC ELEMENTS ARE ACCEPTABLE
                                                  );

            if (ifound != NOT_AN_INDEX)
            {
                return ifound;
            }

        }


        // SET DISTANCE TO MAX TO INDICATE THAT THIS TET HAS ALREADY
        // BEEN TRIED
        dists[indn] = DBL_MAX;
        indn = min_element( dists.begin() , dists.end() ) - dists.begin();

    }

    // IF ALL NEIGHBOURS FAIL
    return NOT_AN_INDEX;


}



void Geometry::genIndToTetsByCoords(vector<unsigned int> &ind,   // return index
                                    double *coord,               // search cordinate values
                                    const unsigned int &nc,      // number of coordinate values
                                    const bool& terminateOnError,// whther to terminate app. if coordinate not found. default = true;
                                    const bool& requireLCElement)// only LC element can be re returned
{
    /*!
    Generates index to tetrahedron that contain coordinate coord.

The 'terminateOnError' flag is used to spcify whether to terminate app. if a coord
    is not found, or to mark it as NOT_AN_INDEX. This may occur e.g. when
    interpolating between two different meshes.

'requireLCElement' determines whether only LC elements can be considered. if this is false,
also dielectrinc elements indexes may be returned. This is often problematic when searching
for an LC node on the boundary between LC and DE regions, i.e. it exists in both regions, but
is only properly defined in the LC element.
*/

    ind.clear();
    unsigned int nt = (unsigned int) this->t->getnElements();
    ind.assign( nc, nt ); // assing with a value that is one too much initially

    vector < set <unsigned int> > p_to_t;
    t->gen_p_to_elem(p_to_t);

    // FIND STARTING TET
    //---------------------------
    double mid[3] = {   ( getXmax() - getXmin() )/2.0 , // CENTRE COORDINATE OF STRUCTURE
                        ( getYmax() - getYmin() )/2.0 ,
                        ( getZmax() - getZmin() )/2.0 };

    unsigned int mt = 0;
    //if ( !getContainingTet( p_to_t, mid, mt) )
    //{
    //    mt = t->getnElements() / 2; // if mid tet not found, set to numtets/2
    //}
    std::set<size_t> tetHistory;
    mt = recursive_neighbour_search( mid,
                                     p_to_t,
                                     0,
                                     tetHistory
                                     );

    if ( mt == NOT_AN_INDEX ) // starting index at centre of structure not found (probably a hole)
    {
        mt = 0;
    }
    //---------------------------

    unsigned int n;
#ifndef DEBUG
#pragma omp parallel for
#endif
    for ( n = 0 ; n < nc ; n++) // for each coord
    {

        std::set<size_t> searchHistory; // keeps track of tested tets to avoid repeating work
        double crd[3];
        //crd = coord + n*3; // n'th coordinate
        crd[0] = *(coord + n*3 + 0);
        crd[1] = *(coord + n*3 + 1);
        crd[2] = *(coord + n*3 + 2);
        // nearest neigbour search

        size_t t0 = recursive_neighbour_search( crd,
                                                p_to_t,
                                                mt,
                                                searchHistory,
                                                requireLCElement    // WHETHER TO ONLY ACCEPT LC ELEMENTS ARE RETURN VALUE
                                                );

        if (t0 != NOT_AN_INDEX )
        {
            ind[n] = t0;
        }
        // if neares neighbour search fails, use brute force
        else
        {
            if (n==152)
            {
                int debug = 9;
            }

            unsigned int bfind = 0;
            // TRY BRUTE FORCE. THIS MAY TERMINATE APP., DEPENDING ON BOOL FLAG
            if ( brute_force_search( bfind,
                                     crd,
                                     terminateOnError,  // TEMINATE PROGRAM IF NODE NOT FOUND
                                     requireLCElement   // WHETHER ONLY LC ELEMENTS ARE ACCEPTED
                                     ) )
            {

                ind[n] = bfind;
            }
            else    // BRUTE FORCE FAIL IS ALLOWED (NODE MAY BE OUTSIDE MESH)
            {
                ind[n] = Geometry::NOT_AN_INDEX; // MARK INDEX AS INVALID
                cout << "could not find regular grid point " << n << "in tetrahedral mesh. Continuing" << endl;
                cout << crd[0] << ", " << crd[1] << ", " << crd[2] << endl;
            }
        }

    }// end for each coord

}

double* Geometry::getPtrTop() 	{ return p;}
double Geometry::getpX(int i) const
{ 
#ifdef DEBUG
    isValidNodeIndex(i);
#endif
    return p[i*3 + 0];
}
double Geometry::getpY(int i) const
{ 
#ifdef DEBUG
    isValidNodeIndex(i);
#endif
    return p[i*3 + 1];
}
double Geometry::getpZ(int i)	const
{ 
#ifdef DEBUG
    isValidNodeIndex(i);
#endif
    return p[i*3 + 2];
}

double Geometry::getAbsXDist(int i , double x)
{
    return fabs( getpX(i) - x );
}
double Geometry::getAbsYDist(int i , double y)
{
    return fabs( getpY(i) - y );
}
double Geometry::getAbsZDist(int i, double z)
{
    return fabs( getpZ(i) - z );
}

double Geometry::getAbsDistSqr(const unsigned int i, const double *const coord) const
{
    double* pp = p + (3*i); //shortcut to ith coordinate

    return  ( ( pp[0]-coord[0] )*( pp[0]-coord[0] )+
	      ( pp[1]-coord[1] )*( pp[1]-coord[1] )+
	      ( pp[2]-coord[2] )*( pp[2]-coord[2] ) ) ;
}

size_t Geometry::getTotalSize(){
    size_t size_p = np*3*sizeof(double);
    size_t size_NodeNormals = np*3*sizeof(double);

    // Tetrahedra mesh
    size_t sze_t_Elem = this->t->getnElements()*this->t->getnNodes()*sizeof(int);
    size_t sze_t_Mat  = this->t->getnElements()*sizeof(int);
    size_t sze_t_Det  = this->t->getnElements()*sizeof(double);


    // Triangle mesh
    size_t sze_e_Elem = this->e->getnElements()*this->e->getnNodes()*sizeof(int);
    size_t sze_e_Mat  = this->e->getnElements()*sizeof(int);
    size_t sze_e_Det  = this->e->getnElements()*sizeof(double);
    size_t sze_e_Con  = this->e->getnElements()*sizeof(int);
    size_t sze_e_Norm = this->e->getnElements()*3*sizeof(double);


    int MB = 1048576;
    printf("size/MB : \np = %i\nNodeNormals=%i\n", (int) size_p / MB, (int) size_NodeNormals/MB);
    printf("t Elements %i\n Mat %i\n Det %i\n",     (int) sze_t_Elem / MB, (int) sze_t_Mat / MB, (int)sze_t_Det / MB);
    printf("e Elements %i\n Mat %i\n Det %i\n Con %i\n Norm %i\n", (int)sze_e_Elem / MB ,
           (int) sze_e_Mat/MB, (int) sze_e_Det/MB, (int)sze_e_Con/MB, (int)sze_e_Norm/MB );

    size_t sze_total =  size_p + size_NodeNormals +
            sze_t_Elem + sze_t_Mat + sze_t_Det +
            sze_e_Elem + sze_e_Mat + sze_e_Det + sze_e_Con + sze_e_Norm;
    printf("Total Size = %i\n", (int) sze_total / MB);
    return sze_total;


}

bool Geometry::getleft_right_is_periodic() {return left_right_is_periodic;}
bool Geometry::getfront_back_is_periodic() {return front_back_is_periodic;}
bool Geometry::gettop_bottom_is_periodic() {return top_bottom_is_periodic;}

double* Geometry::getPtrToNodeNormals() {return NodeNormals;}
double Geometry::getNodeNormalsX(int i)		{return NodeNormals[i*3+0];}
double Geometry::getNodeNormalsY(int i)		{return NodeNormals[i*3+1];}
double Geometry::getNodeNormalsZ(int i)		{return NodeNormals[i*3+2];}
double Geometry::getXmin()	{return Xmin;}
double Geometry::getXmax() 	{return Xmax;}
double Geometry::getYmin() 	{return Ymin;}
double Geometry::getYmax()	{return Ymax;}
double Geometry::getZmin()	{return Zmin;}
double Geometry::getZmax() 	{return Zmax;}

void Geometry::PrintNodes()
{
    printf("printting %i nodes:\n",getnp());
    for (unsigned int i = 0 ; i < getnp() ; i ++)
    {
        printf("\t%u = [%f,%f,%f]\n", i , getpX(i) , getpY(i) , getpZ(i) );
    }


}
void Geometry::PrintNode(unsigned int i )
{
    printf("node %u = [%f,%f,%f]\n",i , getpX(i) , getpY(i) , getpZ(i));

}
void Geometry::PrintPeriodicNodes()
{

    for (unsigned int i = 0 ; i < periNodes_.size() ; i++)
        cout <<"periNodes_["<<i<<"] = " << periNodes_[i] << endl;

}

void Geometry::getTetBaryCentre(double* x, const unsigned int &it){
    x[0] = 0; x[1] = 0; x[2] = 0;
    unsigned int nn = this->t->getnNodes();

    for (unsigned int i = 0 ; i < nn ; i++){   // sums all coords or element it
        unsigned int n = this->t->getNode(it, i);
        x[0] += this->getpX( n );
        x[1] += this->getpY( n );
        x[2] += this->getpZ( n );



    }
    x[0]/= (double)nn; // = average x,y,z-coords
    x[1]/= (double)nn;
    x[2]/= (double)nn;
}
void Geometry::isValidNodeIndex(const unsigned int &i) const
{
// DEBUG ASSERTION TO MAKE SURE A VALID NODE IS ACCESSED
    if (i>(unsigned int) np){
        std::cout << "error, requesting node "<<i<<" , when np is: " << np << " - bye!" << std::endl;
        exit(1);
    }
}
void Geometry::genIndWeakSurfaces(Alignment& alignment)
{
    // GENERATES INDEX TO WEAK SURFACE ELEMENTS

    // FIRST MAKE TEMPORARY VECTOR OF WEAK SURFACES
    std::vector<size_t> ws;
    for (idx i = 0 ; i < (idx) e->getnElements() ; i++)
    {
        int FixLCNum = e->getFixLCNumber((int) i); // MATERIAL NUMBER / FIXLC1

        if ( (FixLCNum > 0) &&
             (!alignment.IsStrong( FixLCNum- 1) ) )  // IS WEAK
            ws.push_back( i );
    }

    numWeakSurf = ws.size();

    // IF WEAK ANCHORING SURFACES EXIST, COPY TO PERMANENT ARRAY
    if (numWeakSurf)
    {
        if (indWeakSurf) free( indWeakSurf );
        indWeakSurf = (size_t*) malloc(numWeakSurf * sizeof(size_t) );

        for (size_t i = 0 ; i < numWeakSurf ; i++)
            indWeakSurf[i] = ws[i];

    }

}


bool Geometry::checkForOverlapingNodes(){
    //bool isOK = true;
    double mindist = 100000000;
    unsigned int mini = 0 , minj = 0;
    for (unsigned int i = 0 ; i < (unsigned int) this->getnp() ; i++){

        for (unsigned int j = i ; j < (unsigned int) this->getnp() ; j++){
            double *pi = getPtrTop() + j*3;
            if (i != j){
                double dist = getAbsDistSqr(i, pi);
                if ( dist < mindist ) {
                    mindist = dist;
                    mini = i;
                    minj = j;
                }
                //printf("dist %i - %i = %e\n", i,j,dist);
            }
        }
    }

    printf("mindist = %e, i,j = %u,%u \n", mindist, mini,minj);
    this->PrintNode(mini);
    this->PrintNode(minj);
    return false;
}
void Geometry::countNodeReferences(vector <int>& refc, Mesh& mesh){
    refc.clear();
    refc.reserve( (size_t) this->getnp() );
    refc.assign( (size_t) this->getnp() , 0 ); // set all counts to zero

    for (idx i = 0 ; i < mesh.getnElements() ; i++)
    { // for each element
        for (idx n = 0 ; n < mesh.getnNodes() ; n++)
        { // for each node in element i
            refc[mesh.getNode(i,n)] ++ ;
        }
    }


}
