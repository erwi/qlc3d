
#include <geometry.h>
Geometry::Geometry()
{
	np 					= 0;
	npLC				= 0;
	p 					= NULL;
	NodeNormals			= NULL;
	t 					= new Mesh();
	e					= new Mesh();
	//oct				= NULL;
	//this_oct			= NULL;
	Xmin 				= 0;
	Xmax				= 0;
	Ymin				= 0;
	Ymax				= 0;
	Zmin				= 0;
	Zmax				= 0;
	t->setDimension(3);
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
		//t->~Mesh();
		//e->~Mesh();
        //printf("NOTE: ~Oct_Box() -> error, fix it!\n");
	//oct->~Oct_Box();
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
	for (int i = 0 ; i < 3*getnp(); i++)
		p[i] = temp[i];
		
	// copy node normals
	temp = geom->getPtrToNodeNormals();
	for (int i = 0 ; i < 3*getnp(); i ++)
		NodeNormals[i] = temp[i];
	
	left_right_is_periodic = geom->getleft_right_is_periodic();
	front_back_is_periodic = geom->getfront_back_is_periodic();
	top_bottom_is_periodic = geom->gettop_bottom_is_periodic();
		
		//Oct_Box* oct;				// octree mesh index for fast search
		//Oct_Box* this_oct;
		
	//exit(0);
	

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

void Geometry::setCoordinates(double* coords, int np)
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
	
	for (int i = 0 ; i < np ; i ++) // copy coordinates
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
void Geometry::addCoordinates(double* coords, int npn){
// APPENDS npn NEW COORDIANTES IN coords TO this->p
    if ( ( npn > 0 ) && (coords) ){ // don't do anything if no new coordinates
	int np_new = np + npn;
	double* pnew = (double*) malloc( 3*np_new*sizeof(double));
	
	// first copy old coordinates
	for(int i = 0 ; i < 3*np ; i++)
		pnew[i] = p[i];
	// then add new ones
	for(int i = 0 ; i < 3* npn ; i ++)
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




void Geometry::setnp(int n)		{np = n;}
void Geometry::setnpLC(int n)	{npLC = n;}
void Geometry::setNodeNormals(){
	
	if (e->getnElements() == 0){
			printf("error - Geometry::setNodeNormals() - surface mesh not defined, bye!\n");
			exit(1);
	}
	
	if (NodeNormals!=NULL) free(NodeNormals);						// ALLOCATE
		NodeNormals = (double*)malloc(3 * np * sizeof(double));

	memset(NodeNormals, 0, 3 * np * sizeof(double));
		
	double tempn[3] ={0,0,0};	
	//printf("nelemenst = %i, np = %i\n",e->getnElements(),np);
	
	for (int i = 0 ; i < e->getnElements() ; i ++){ // add neighbouring surface normals
		int m = e->getMaterialNumber(i);
		if ((m != MAT_PERIODIC) && (m!= MAT_NEUMANN)){ //if element is not periodic or neuman surface
			e->CopySurfaceNormal(i,&tempn[0]);	// copy surface triangle normal to temp normal
			
			for (int j = 0; j < e->getnNodes() ; j++){
				NodeNormals[e->getNode(i,j)*3 + 0] += tempn[0]; // add x, y, z components for each node in triangle i
				NodeNormals[e->getNode(i,j)*3 + 1] += tempn[1];
				NodeNormals[e->getNode(i,j)*3 + 2] += tempn[2];
			}// end for j
			tempn[0] = 0; tempn[1] = 0; tempn[2] = 0;
		}// end if not periodic 
	}// end for i

	for (int i = 0 ; i < np ; i ++){ // normalise length
		double len = sqrt( NodeNormals[i*3+0]*NodeNormals[i*3+0] + NodeNormals[i*3+1]*NodeNormals[i*3+1] + NodeNormals[i*3+2]*NodeNormals[i*3+2] );
		if (len>0){
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
for (int i = 0 ; i < np ; i ++)	{
		printf("NodeNormals[%i] = [%f,%f,%f]\n", i, NodeNormals[i*3], NodeNormals[i*3+1], NodeNormals[i*3+2]);
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
	
	int* tmat = t->getPtrToMaterialNumber(0);
	for (int i = 0 ; i < t->getnElements() ; i++)
		if (tmat[i]>= MAT_DIELECTRIC1)	// if material nuber > LC number
		{
			DE_exist = true;
			break;
		}
		
	setnpLC( getnp() );	
	if (!DE_exist) // if no dielectric materials -> no need to reorder = exit
		return;
		

	//
	//GENERATE LIST OF MATERIAL NUMBERS FOR NODES. IN CASE OF DUAL VALUE LC PRECEDES	
	//

	// mark all LC nodes as 1 and others as 0
	int* lcde = (int*) malloc(getnp() * sizeof(int) );
	memset(lcde , 0 , getnp() * sizeof(int) ); // start with everything 0
	
	for (int i = 0 ; i < t->getnElements() ; i ++ )
	{
		if (t->getMaterialNumber(i) == MAT_DOMAIN1)
		{
			for (int j = 0 ; j < t->getnNodes() ; j ++ ) // loop over all nodes
				lcde[t->getNode(i , j)] = 1; // LC --> 1
		}
	}
	
	//for (int i = 0 ; i < getnp() ; i ++ )
	//printf("lcde[%i] = %i\n", i, lcde[i]);
		
	
	//list <int> l_mat;	 // nodes
	npLC = 0;
	vector <int> v_mat_index;
	
	
	for (int i=0; i< getnp() ;i++)// first add all nodes marked as LC
		if (lcde[i] == 1)
		{
			v_mat_index.push_back(i);
			npLC++;
		}
	
	for (int i=0; i< getnp() ;i++) // then add all non-LC nodes
		if (lcde[i] == 0 )
			v_mat_index.push_back(i);
	
	//printf("np = %i, npLC = %i\n",np, npLC);
	free(lcde);
	
	//make inverse map
	vector <int> v_invmap;
	v_invmap.resize( v_mat_index.size() , -1);
	for (int i = 0 ; i < getnp() ; i ++)
		v_invmap[v_mat_index[i]]= i;//v_mat_index[i];//= i;
	
	
	//for( int i  = 0 ; i < getnp() ; i ++ )
	//printf("%i <->  %i\n", v_invmap[i] , v_mat_index[i] );
	

//REORDER NODES 
	double* newp = (double*)malloc( 3*getnp() *sizeof(double)); // memory for reordered node coordinates
	for (int i = 0 ; i < getnp() ; i++)
	{
		
		newp[i*3+0]=getpX( v_mat_index[i] ); //x-coord
		
		newp[i*3+1]=getpY( v_mat_index[i] ); //y-coord
		
		newp[i*3+2]=getpZ( v_mat_index[i] ); //z-coord
		
	}
	
	free(p); // make p = new reordered p
	p=newp;
	
//REORDER TETRAHEDRA
	int* newt = (int*)malloc(t->getnElements() * t->getnNodes() * sizeof(int));
	for (int i = 0 ; i < t->getnElements() ; i++ )
	{
		for (int j = 0 ; j < t->getnNodes() ; j ++ ) // loop over all nodes of element i
			newt[i*t->getnNodes() + j ] = v_invmap[ t->getNode(i , j) ];
	}
	t->setAllNodes(newt ); // copy node numbers
	free(newt);
	//exit(1);
	
//REORDER TRIANGLES
	int* newe = (int*)malloc(e->getnElements() * e->getnNodes() * sizeof(int));

	for(int i = 0 ; i < e->getnElements() ; i++)
	{
		for (int j = 0 ; j < e->getnNodes() ; j++ )
			newe[i * e->getnNodes() + j ] = v_invmap[ e->getNode(i,j) ];
	}
	e->setAllNodes(newe ); // copy node numbers
	free(newe);
	
	
}


void Geometry::MakePeriEquNodes(){
        //printf("in MakePeriEquNodes()...\n");
        if (p == NULL){
		printf("error - can't find periodic equ nodes - bye\n");
		exit(1);
	}
	peri_equ_nodes.resize(getnp());
	
}// end void MakePEriEquNodes()

void Geometry::checkForPeriodicGeometry()
{
// CHECKS FOR TYPES OF PERIODICITY PRESENT IN CURRENT STRUCTURE.
// POSSIBLE PERIODIC SURFACES ARE:
//      LEFT/RIGHT
//      FRONT/BACK
//      TOP/BOTTOM

    for (int i = 0 ; i < e->getnElements() ; i++)
    {
        if (e->getMaterialNumber(i) == MAT_PERIODIC) // if surface is periodic
        {
            // check to see which side surface is on by looking at the surface normal
            double* snorm = e->getPtrToSurfaceNormal(i);
			
            // IF SURFACE NORMAL X-COMPONENT = 1
            if (fabs( fabs(snorm[0]) - 1.0 ) < EPS)
                        left_right_is_periodic = true;
            else
            // IF SURFACE NORMAL Y-COMPONENT = 1
            if (fabs( fabs(snorm[1]) - 1.0 ) < EPS)
                front_back_is_periodic = true;
            else
            // IF SURFACE NORMAL Z-COMPONENT = 1
            if (fabs( fabs(snorm[2]) - 1.0 ) < EPS)
                top_bottom_is_periodic = true;
            else
            {
                printf("error - checkForPeriodicGeometry() - periodic surface element %i has invalid normal:\n [%f, %f, %f] - bye!\n",i, snorm[0] , snorm[1], snorm[2]);

                this->e->PrintNormals();

                exit(1);
            }
	
            // IF ALL SURFACES HAVE ALREADY BEEN IDENTIFIED AS PERIODIC
            // NO NEED TO CHECK FURTHER TRIANGLES
            if ((getleft_right_is_periodic()) &&
                    (getfront_back_is_periodic()) &&
                    (gettop_bottom_is_periodic() ))
                break;
	
        }
    }
    // IF ANY PERIODIC TRIANGLES WERE DETECTED
    if (    getleft_right_is_periodic()
        ||  getfront_back_is_periodic()
        ||  gettop_bottom_is_periodic() )
    {
        MakePeriEquNodes();
    }
	
}




void Geometry::brute_force_search( unsigned int &ind, double* coord){
// BRUTE FORCE DEBUG SEARCH FOR TETRAHEFRON THAN CONTAINS POINT WITH COORDINATES IN coord
// coord IS ASSUMED TO BE OF LENGTH 3, FOR x, y, z
	// loop over each element
	for (unsigned int i = 0 ; i < (unsigned int) t->getnElements() ; i++){
        if ( ( t->ContainsCoordinate( i , getPtrTop(), coord ) ) // If coord is in tet i AND
                && (t->getMaterialNumber(i) <= MAT_DOMAIN7) ){   // tet i material is LC
			ind = i ;
			return; // exit function when found
		}
	}// end for loop over all elems

    printf("error - brute_force_search could not find coord %f,%f,%f - bye ! \n", coord[0], coord[1], coord[2] );
	exit(1);

}

bool Geometry::getContainingTet(vector< set < unsigned int> >& p_to_t,
								double *crd,  // coords to search
								unsigned int& t0){ // start tet

	bool found = false;
	unsigned int t_prev = t->getnElements();
    for (size_t i = 0 ; i < (size_t) t->getnElements() ; i++){ // for all elements

        if  ( ( t->ContainsCoordinate(t0, p , crd) ) &&   // if element t0 contains coordinate AND
            ( t->getMaterialNumber(t0) <= MAT_DOMAIN7) ){ // is an LC element
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
			for (unsigned int d = 0 ; d<vecneighs.size() ; d++){
				dst.push_back( t->CalcBaryDistSqr(p , vecneighs[d], crd));
			}
			// FIND INDEX TO MINIMUM DISTANCE AND SET t0 TO THAT
			unsigned int ind_min = min_element(dst.begin(), dst.end() ) - dst.begin();
			t_prev = t0;
			t0 = vecneighs[ind_min];

			// SEARCH MAY GET STUCK WHEN POINT IS CLOSE TO A CORNER NODE OF A LONG THIN
			// TET CONTAINING THE COORD. IN THAT CASE, THE BARYCENTRE OF A NEIGHBORING
			// TET MAY BE CLOSER -> SEARCH WILL NOT FIND THE ACTUAL TET
			if (t_prev == t0){ // IF STUCK, TRY ALL NEIGHBOURS IN ORDER OF DISTANCE

                dst[ind_min] = DBL_MAX; // SET CURRENT t0 TO MAXIMUM POSSIBLE DISTANCE AS DEFINED IN <float.h>

				for (size_t ind = 0 ; ind < vecneighs.size() ; ind++){ // LOOP OVER ALL NEIGHBOURS
					ind_min = min_element(dst.begin(), dst.end() ) - dst.begin(); // refind new closest value

                    if ( ( t->ContainsCoordinate( vecneighs[ind_min], p, crd ) )&& // If contains coordinate AND
                    (t->getMaterialNumber(vecneighs[ind_min]) <= MAT_DOMAIN7 ) )   // If is LC element
                    {
						found = true;
						t0 = vecneighs[ind_min];
						break;
					}
					else{
                        dst[ind_min] = DBL_MAX; // STILL NOT FOUND?, TRY AGAIN WITH TET THAT IS NEXT CLOSEST
					}
				}// end for all neighbours
				// IF REALLY STUCK EXIT AND TEST USING BRUTE FORCE
				break;
			} // end if stuck

	}// end for all tets
	return found;
}


void Geometry::genIndToTetsByCoords(vector<unsigned int> &ind, double *coord, const unsigned int &nc){
    /*! Generates index to tetrahedron that contain coordinate coord*/

    ind.clear();
    ind.assign( nc, this->t->getnElements()); // assing with a value that is one too much initially

	vector < set <unsigned int> > p_to_t;
	t->gen_p_to_elem(p_to_t);

	// FIND STARTING TET
	//---------------------------
	double mid[3] = { ( getXmax() - getXmin() )/2.0 , // CENTRE COORDINATE OF STRUCTURE
					  ( getYmax() - getYmin() )/2.0 ,
					  ( getZmax() - getZmin() )/2.0 };

	unsigned int mt = 0;
	if ( !getContainingTet( p_to_t, mid, mt) ){
		mt = t->getnElements() / 2; // if mid tet not found, set to numtets/2
	}
	//---------------------------

	unsigned int n;

	#pragma omp parallel for
	for ( n = 0 ; n < nc ; n++){ // for each coord

		unsigned int t0 = mt ; // SET START ELEM

		double* crd;
		crd = coord + n*3; // n'th coordinate

        // nearest neigbour search
		if ( getContainingTet( p_to_t, crd, t0) ){
			//cout << "found :" << t0<< endl;
			ind[n] = t0;
		}
        // if neares neighbour search fails, use brute force
        else{
			cout << "warning - Geometry::genIndToTetsByCoords \n\tcould not find coord using fast method." << endl;
			cout << "\tTrying slower brute force instead...";
			fflush(stdout);
			unsigned int bfind = 0;
            brute_force_search( bfind, crd); // IF BRUTE FORCE DOES NOT FIND -> TERMINATE PROGRAM
			cout << "OK, found it now!" << endl;
			ind[n] = bfind;
		}

	}// end for each coord




}

int 	Geometry::getnp()		{ return np;}
int 	Geometry::getnpLC()		{ return npLC;}
double* Geometry::getPtrTop() 	{ return p;}
double Geometry::getpX(int i)	
{ 
	#ifdef DEBUG
	isValidNodeIndex(i);
	#endif
	return p[i*3 + 0];
}
double Geometry::getpY(int i)   
{ 
	#ifdef DEBUG
	isValidNodeIndex(i);
	#endif
	return p[i*3 + 1];
}
double Geometry::getpZ(int i)	
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

double Geometry::getAbsDistSqr(const unsigned int &i, double *coord){
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

void Geometry::PrintNodes(){
	printf("printting %i nodes:\n",getnp());
	for (int i = 0 ; i < getnp() ; i ++)
	{
		printf("\t%i = [%f,%f,%f]\n", i , getpX(i) , getpY(i) , getpZ(i) );
	}


}
void Geometry::PrintNode(int i )
{
	printf("node %i = [%f,%f,%f]\n",i , getpX(i) , getpY(i) , getpZ(i));

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
void Geometry::isValidNodeIndex(const unsigned int &i) const{
	if (i>(unsigned int) np){
		std::cout << "error, requesting node "<<i<<" , when np is: " << np << " - bye!" << std::endl;
		exit(1);
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

    for (unsigned int i = 0 ; i < (unsigned int) mesh.getnElements() ; i++){ // for each element
        for (int n = 0 ; n < mesh.getnNodes() ; n++){ // for each node in element i
           refc[mesh.getNode(i,n)] ++ ;
        }
    }


}
