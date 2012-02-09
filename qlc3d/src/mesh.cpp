#include <qlc3d.h>
#include <mesh.h>
#include <material_numbers.h>
#include <set>
Mesh::Mesh()
{
    nElements 		= 0;
    nNodes 			= 0;
    MaxNodeNumber           = 0;
    TotalSize		= 0;
    Dimension		= 0;
    Elem			= NULL;
    Mat			= NULL;
    Determinant		= NULL;
    SurfaceNormal	= NULL;
    ConnectedVolume = NULL;
}
Mesh::~Mesh()
{
    if (Elem!=NULL)	{
        free(Elem);
        Elem=NULL;
    }
    if (Mat !=NULL)	{
        free(Mat);
        Mat=NULL;
    }
    if (Determinant != NULL)	{
        free(Determinant);
        Determinant = NULL;
    }
    if (SurfaceNormal != NULL)	{
        free(SurfaceNormal);
        SurfaceNormal = NULL;
    }
    if (ConnectedVolume != NULL)	{
        free(ConnectedVolume);
        ConnectedVolume = NULL;
    }
}
Mesh::Mesh(idx n,idx nNod)
{

    nElements = n;
    nNodes = nNod;
    // allocate memory for mesh
    Elem = (idx*)malloc(n * nNodes * sizeof(idx));
    Mat  = (idx*)malloc(n * sizeof(idx));
    Determinant = (double*)malloc(n * sizeof(double));
    SurfaceNormal = NULL;
    ConnectedVolume = NULL;
}
void Mesh::AllocateMemory()
{
    if ( (nElements>0) && (nNodes>0) )
    {
        if (Elem != NULL) free(Elem);
        if (Mat != NULL) free(Mat);
        if (Determinant != NULL) free(Determinant);
        Elem = (idx*) malloc(nElements * nNodes * sizeof(idx) );
        Mat  = (idx*) malloc(nElements * sizeof(idx) );
        Determinant = (double*) malloc(nElements * sizeof(double) );
	
    }
    else
    {
        printf("error - Mesh::AllocateMemory - bye!\n");
        printf("nElements = %u , nNodes = %u\n", nElements , nNodes);
        exit(1);
    }
}

idx Mesh::getMaterialNumber(const idx e) const
{
#ifdef DEBUG
    if (Mat == NULL)	{
        printf("error - Mesh::getMaterialNumber - trying to acces NULL pointer, bye!\n");
        exit(1);
    }
    if ( (e<0) || (e>=nElements))
    {
        printf("error - Mesh::getMaterialNumber - index out of bounds, bye!\n");
        printf("trying to access material[%i], when nElements is %i\n", e, getnElements() );
        exit(1);
    }
#endif
    return Mat[e];
}
idx Mesh::getFixLCNumber(const idx e) const
{
    return (idx) getMaterialNumber(e) / MAT_FIXLC1;
}
idx Mesh::getDielectricNumber(const idx e) const
{
    return (int) getMaterialNumber(e) / MAT_DIELECTRIC1;
}
double Mesh::getDeterminant(const idx i) const
{
    // RETURNS THE PRE-CALCULATED DETERMINANT FOR ELEMENT i
#ifdef DEBUG
    if (Determinant == NULL)
    {
        printf("error - Mesh::getDeterminant - Determinants not initialised, bye!");
        exit(1);
    }
    if ((i<0) ||  (i > nElements))
    {
        printf("error - Mesh::getDeterminant - index out of bounds, bye!");
        exit(1);
    }
#endif
    return Determinant[i];
}


idx Mesh::getConnectedVolume(const idx e) const
{
    // RETURNS INDEX TO CONNECTED LC TETRAHEDRON, WHEN e IS INDEX TO A TRIANGLE
#ifdef DEBUG
    if ( ConnectedVolume == NULL)
    {
        printf("error - Mesh::getConnectedVolume - NULL pointer - bye!\n");
        exit(1);
    }
    if ( ( e > getnElements() ) || (e < 0) )
    {
        printf("error - Mesh::getConnectedVolume(int e) - e = %i - bye!\n",e);
        exit(1);
    }
    if (getDimension() != 2)
    {
        printf("error - Mesh::getConnectedVolume(int e) only works for surface meshes... sorry - bye\n");
        exit(1);
    }
#endif
    return ConnectedVolume[e];
}
idx* Mesh::getPtrToElement(const idx e) const
{
    return &Elem[e*nNodes];
}
idx* Mesh::getPtrToMaterialNumber(const idx e)const
{
    return &Mat[e*nNodes];
}
idx* Mesh::getPtrToConnectedVolume(const idx e)const
{
    return &ConnectedVolume[e];
}
double* Mesh::getPtrToDeterminant(const idx e)const
{
    return &Determinant[e];
}
double* Mesh::getPtrToSurfaceNormal(const idx e)const
{
    return &SurfaceNormal[e*3];
}


void Mesh::setDeterminant(const idx i, double det)
{
#ifdef DEBUG
    if (Determinant == NULL)
    {
        printf("error - Mesh::setDeterminant - Determinants not initialised, bye!");
        exit(1);
    }
    // ADD CHECKING OF i TOO LARGE OR SMALL

#endif
    Determinant[i] = det;
}
void Mesh::setSurfaceNormal(idx i, double norm[3])
{
    if ((i > getnElements()) || (SurfaceNormal==NULL) )
    {
        printf("error - Mesh::setSurfaceNormal(int i, double norm[3]) - bye!\n");
        exit(1);
    }

    SurfaceNormal[i*3 + 0] = norm[0];
    SurfaceNormal[i*3 + 1] = norm[1];
    SurfaceNormal[i*3 + 2] = norm[2];

}
void Mesh::setAllNodes(idx *nodes) // copies all values from array nodes to this->Elem
{
    if ((nElements>0) && (Elem!=NULL) ) // make sure mesh is initialised
    {
        for (idx i = 0 ; i < nElements*nNodes ; i++)
	    Elem[i] = nodes[i];
    }
    else{
	printf("error in Mesh::SetAllNodes - mesh not initilised\n");
	exit(1);
    }
}// end void SetAllNodes

void Mesh::setAllMaterials(idx *mat)
{
#ifdef DEBUG
    if (!Mat)
    {
        std::cout<< "error in" << __func__ <<", mesh not initialised - bye!"<< std::endl;
        exit(1);
    }
#endif
    memcpy( Mat, mat, nElements*sizeof(idx) );
}// end void SetAllMaterials
void Mesh::setDimension(const idx dim)
{
    Dimension = dim;
}

void Mesh::setnElements(idx nelem)
{
    nElements = nelem;
}
void Mesh::setnNodes(idx nnodes)
{
    nNodes = nnodes;
}



void dumbsearch( Mesh* vol, Mesh* surf, int ind){

    idx* nodes = surf->getPtrToElement( ind );
    bool found = false;
    for (idx i = 0 ; i < vol->getnElements() ; i++){
        if ( vol->ContainsAllNodes(i, 3, nodes) )
        {
            found = true;
            break;
        }
    }
    if (found )
        cout << "FOUND" << endl;
    else
        cout << "REALLY NOT FOUND !" << endl;


    exit(1);
}


void Mesh::setConnectedVolume(Mesh* vol)
{
    // finds indexes from triangles to tetrahedra in 'vol' that share 3 points.
    // these are needed for surface integrals, e.g. Neumann boundaries.
    // counts number of times volumes are connected to each node
    // then compares with nodes in surface mesh

    // Check that surface and volume mesh dimensions are OK
    if ( (getDimension() != 2 ) ||
         (vol->getDimension() != 3) )
    {
        printf("error - Mesh::setConnectedVolumes(Mesh* vol) - bye!\n");
        printf("this.Dimension = %u, t.Dimension = %u\n", this->getDimension(), vol->getDimension());
        exit(1);
    }
    // release old index, if exists
    if (ConnectedVolume != NULL){
        free(ConnectedVolume);
        ConnectedVolume = NULL;
    }

    // set up index and initialise all to -1
    ConnectedVolume = (idx*) malloc( getnElements() * sizeof(idx) );
    for (idx i = 0 ; i < getnElements() ; i ++)
        ConnectedVolume[i] = -1;

    vector <set <idx> > v_in_p ; // vector of sets containing volume element numbers connected to each node
    vector < set < idx> > p_to_t;
    vol->gen_p_to_elem( p_to_t);

    for (idx i = 0 ; i < this->getnElements() ; i++) // loop over surface elements
    {
        if (( getMaterialNumber( i ) < MAT_ELECTRODE1) ||   // exclude eletrode only surfaces as
                (  getMaterialNumber( i ) > MAT_ELECTRODE9) )    // these are fixed, and do not need to know
            // about neighbouring LC tets
        {
            // NOW HAVE 3 LISTS OF TETS CONNECTED TO THIS SURFACE
            // ONE PER CORNER NODE. INTERSECTION OF THESE LISTS
            // SHOULD LEAVE A SINGLE INDEX TO A TET THAT HAS THIS
            // TRIANGLE AS ONE OF ITS FACES
            set <idx> tets1;
            set <idx> tets2;
            set <idx> tets3;
            idx n1 = this->getNode( i , 0 );
            idx n2 = this->getNode( i , 1 );
            idx n3 = this->getNode( i , 2 );
            tets1.insert( p_to_t[ n1 ].begin() , p_to_t[ n1 ].end() ); // tets connected to node 1
            tets2.insert( p_to_t[ n2 ].begin() , p_to_t[ n2 ].end() ); // tets connected to node 2
            tets3.insert( p_to_t[ n3 ].begin() , p_to_t[ n3 ].end() ); // tets connected to node 3

            // FIND INTESRECTIONS OF THE LISTS
            set <idx> intr1;
            set <idx> intr2;
            set <idx> final;

            set_intersection( tets1.begin() , tets1.end() , tets2.begin() , tets2.end() , inserter( intr1 , intr1.begin() ) );
            set_intersection( tets2.begin() , tets2.end() , tets3.begin() , tets3.end() , inserter( intr2 , intr2.begin() ) );
            set_intersection( intr1.begin() , intr1.end() , intr2.begin() , intr2.end() , inserter( final , final.begin() ) );

            if (final.size() == 0)
            {
                printf("error - Mesh::setConnectedVolume - tri %u does not seem to be connected to any tets!\n", i);
                this->PrintElement( i );
                printf("%i,%i,%i\n", (int) tets1.size(), (int) tets2.size() , (int) tets3.size());
                printf("%i,%i\n", (int) intr1.size(), (int) intr2.size() );
                dumbsearch( vol, this, i);
                exit(1);
            }

            // SURFACE CAN BE BETWEEN LC AND DIELECTRIC LAYERS
            // THEN final.size() IS 2. MUST CHECK WHICH ONE IS LC
            // AND WHICH IS DIELECTRIC
            set <idx> ::iterator itr;
            itr = final.begin();
            if ( final.size() == 1) // IF SINGLE TET FOUND
            {
                if (vol->getMaterialNumber( *itr ) <= MAT_DOMAIN7 ) // only connect to an LC element
                {
                    ConnectedVolume[i] = *itr;
                }
                else
                    if ( ( vol->getMaterialNumber( *itr) >= MAT_DIELECTRIC1 ) &&
                         (vol->getMaterialNumber(*itr) <= MAT_DIELECTRIC7) )
                    {
                        ConnectedVolume[i] = -2; // set to -2 if connected to dielectric
                    }
            }
            else
                if ( final.size() == 2)
                {
                    // CHECK BOTH FIRST AND SECOND INDICES
                    if (vol->getMaterialNumber( *itr ) <= MAT_DOMAIN7 )// only connect to an LC element
                        ConnectedVolume[i] = *itr;
                    else{
                        itr++;
                        if (vol->getMaterialNumber( *itr ) <= MAT_DOMAIN7 )// only connect to an LC element
                            ConnectedVolume[i] = *itr;
                    }
                }
                else
                {
                    printf("oops. somethings gone wrong in Mesh::setConnectedVolume(Mesh* vol, int np) \n");
                    printf("tri %u seems to be connected to %i tets \n", i , (int) final.size() );
                    printf("triangle:\n");
                    this->PrintElement( i );
                    exit(1);
                }

        }// end exclude electrode only surfaces
    }


}
double Mesh::Calculate4x4Determinant(double* M) const
{
    // CALCULATES DETERMINANT OF THE 4x4 MATRIX M, WHERE M
    // MOST LIKELY HOLDS COORDINATE VALUES. THIS COULD BE MOVED ELSEWHERE
    // (GENERAL UTULITY FUNCTIONS?) AS IT DOES NOT TOUCH ANY OF THE Mesh DATA.
    double x[4], y[4], z[4];

    for (int i = 0 ; i < 4 ; i++){
        x[i] = M[i*3 + 0];
        y[i] = M[i*3 + 1];
        z[i] = M[i*3 + 2];
    }

    double A[3],B[3],C[3];
    A[0] = x[1] - x[0]; A[1] = y[1] - y[0]; A[2] = z[1] - z[0];
    B[0] = x[2] - x[0]; B[1] = y[2] - y[0]; B[2] = z[2] - z[0];
    C[0] = x[3] - x[0]; C[1] = y[3] - y[0]; C[2] = z[3] - z[0];

    // B x C
    double cross[3];
    cross[0] =   B[1]*C[2] - C[1]*B[2];
    cross[1] = - B[0]*C[2] + C[0]*B[2];
    cross[2] =   B[0]*C[1] - C[0]*B[1];

    // determinant = A . (B x C)
    double Determinant = A[0]*cross[0] + A[1]*cross[1] + A[2]*cross[2];
    return Determinant;
}
bool Mesh::ContainsAllNodes(const idx elem, const idx n, idx* nodes) const
{
// CHECKS WHETHER ELEMENT elem CONTAINS ALL n NODES IN ARRAY nodes
    // assumes that all values in nodes are unique,
    // i.e. elements are valid, with no repeated node numbers

    idx counter = 0;
    for (idx i = 0 ; i < n ; i++) // loop over all nodes in array nodes
    {
        for (idx j = 0 ; j < getnNodes() ; j ++) // loop over all nodes in element
        {
            if ( getNode(elem,j) == nodes[i] ) // if found
            {
                counter++;
                break;
            }

        }// end foor j
    }// end for i

    if (counter == n )
        return true;
    else
        return false;

}
void Mesh::ContainsNodes(list <idx>* elems, list <idx>* points)
{
// CREATES INDEX OF ELEMENTS THAT CONTAIN AT LEAST ONE OF THE NODES IN LIST points

    list <idx>::iterator int_itr;
    for (idx i = 0 ; i < getnElements() ; i ++ ) // loop over all elements
    {
        for (int_itr = points->begin() ; int_itr != points->end() ; int_itr++) // loop over all points in list
        {
            for (idx j = 0 ; j < getnNodes() ; j++) // loop over all nodes in element i
            {
                if ( getNode(i,j) == *int_itr)
                {
                    elems->push_back(i);
                }
            }
        }
    }
    // sort and remove duplicated entries from neighbouring elements
    elems->sort();
    elems->unique();
}
bool Mesh::ContainsCoordinate(const idx elem, const double *p, const double *coord) const
{
    // CHECKS TO SEE WHETHER ELEMENT elem CONTAINS COORDINATE coord. THAT IS, COORD
    // CAN BE ANYWHERE INSIDE elem, NOT JUST CORNER NODES
    // coord is expected to be of length 3 [x,y,z]. is most likely a pointer to somewhere in array *p
#ifdef DEBUG
    if (this->Dimension == 2){
	printf("Mesh::ContainsCoordinate only works on tets - bye!\n");
	exit(1);
    }
#endif

    // CALCULATE 4 DETERMINANTS FOR TETS FORMED USING NODES FROM THIS ELEMENT
    // AND THE TESTED COORDINATE AND ONE USING THIS TET ONLY. IF SUM OF DETS OF TETS FORMED USING coord
    // EQUALS THAT OF FOR THIS ELEMENT ONLY, THEN THE COORDINATE IS WITHIN THIS TET
    const double eps = 1e-15; // accuracy used for coordinate comparions
    double x[4] , y[4] , z[4];
    idx *n = &Elem[elem*nNodes]; // shortcut to all node numbers in this tet
    double Det = 0.0;    // cumulative determinant for tets formed using coord
    double Det_this = 0; // determinant for the whole element

    // FIRST CHECK WHETHER coord IS ONE OF THE CORNER NODES (IT OFTEN IS)
    // THIS IS A QUICK TEST
    x[0] = p[3*n[0]];   y[0] = p[3*n[0]+1]; z[0] = p[3*n[0]+2];
    x[1] = p[3*n[1]];   y[1] = p[3*n[1]+1]; z[1] = p[3*n[1]+2];
    x[2] = p[3*n[2]];   y[2] = p[3*n[2]+1]; z[2] = p[3*n[2]+2];
    x[3] = p[3*n[3]];   y[3] = p[3*n[3]+1]; z[3] = p[3*n[3]+2];

    if  ( (	  ( fabs( coord[0]- x[0]) < eps ) && (fabs(coord[1] - y[0]) < eps) && (fabs(coord[2] - z[0])<eps) )|| // node 1
          ( ( fabs( coord[0]- x[1]) < eps ) && (fabs(coord[1] - y[1]) < eps) && (fabs(coord[2] - z[1])<eps) )|| // node 2
          ( ( fabs( coord[0]- x[2]) < eps ) && (fabs(coord[1] - y[2]) < eps) && (fabs(coord[2] - z[2])<eps) )|| // node 3
          ( ( fabs( coord[0]- x[3]) < eps ) && (fabs(coord[1] - y[3]) < eps) && (fabs(coord[2] - z[3])<eps) )){ // node 4
	//cout << "coerner node found" << endl;
	return true;
    }


    // THEN CHECK IF coord IS AN INTERNAL NODE. THIS IS SLOWER
    for (unsigned int i = 0 ; i < (unsigned int) getnNodes()+1 ; i++ ){

	x[0] = p[3*n[0]];   y[0] = p[3*n[0]+1]; z[0] = p[3*n[0]+2];
	x[1] = p[3*n[1]];   y[1] = p[3*n[1]+1]; z[1] = p[3*n[1]+2];
	x[2] = p[3*n[2]];   y[2] = p[3*n[2]+1]; z[2] = p[3*n[2]+2];
	x[3] = p[3*n[3]];   y[3] = p[3*n[3]+1]; z[3] = p[3*n[3]+2];

        // OVERWRITE NODE i WITH INPUT COORDINATE VALUE
	if (i < (unsigned int) getnNodes() ){ // use coord as one of the values
	    x[i] = coord[0] ; y[i] = coord[1] ; z[i] = coord[2];
	}

	double A[3],B[3],C[3];
	A[0] = x[1] - x[0]; A[1] = y[1] - y[0]; A[2] = z[1] - z[0];
	B[0] = x[2] - x[0]; B[1] = y[2] - y[0]; B[2] = z[2] - z[0];
	C[0] = x[3] - x[0]; C[1] = y[3] - y[0]; C[2] = z[3] - z[0];

        // B x C
	double cross[3];
	cross[0] =   B[1]*C[2] - C[1]*B[2];
	cross[1] = - B[0]*C[2] + C[0]*B[2];
	cross[2] =   B[0]*C[1] - C[0]*B[1];

        // determinant = A . (B x C)
	if (i < (unsigned int) getnNodes() )
	    Det += fabs( A[0]*cross[0] + A[1]*cross[1] + A[2]*cross[2] ) ; // det for tets formed using coord
	else
	    Det_this = fabs( A[0]*cross[0] + A[1]*cross[1] + A[2]*cross[2] ) ; // determinant for the whole tet

    }

    if ( fabs( Det_this - Det) <= eps )
	return true;
    else
	return false;



}

void Mesh::CompleteNodesSet(const idx elem, std::vector <idx>& nodes) const
{
    // IF NODES IS EMPTY, COPY ALL NODE NUMBERS
    if (nodes.size() == 0 ){
        for (idx i = 0 ; i < this->getnNodes() ; i++)
        {
	    nodes.push_back( this->getNode( elem, i ));
	}
    }
    // OTHERWISE ONLY ADD THOSE NODES THAT ARE MISSING FROM nodes VECTOR
    // (WHY NOT JUST ALWAYS MAKE A COPY OF ELEMENT NODES? WHERE IS THIS USED?)
    else
    {
        for ( idx i = 0 ; i < this->getnNodes() ; i++)
        {
	    bool found = false;
            idx node = this->getNode(elem, i);
            for (int j = 0 ; j < (int) nodes.size() ; j++ )
            {
                if ( nodes[j] == (unsigned int) node )
                {
		    found = true;
		    break;
		}
	    }
	    if (found == false){
		nodes.push_back( node );
	    }
	}
    }
#ifdef DEBUG
    // MAKE SURE CORRECT NUMBER OF NODES IS RETURNED
    if (nodes.size() != (unsigned int) this->getnNodes() ){
        printf( "error in Mesh::CompleteNodesSet - return size does not mathc number of nodes in element\n");
        printf( "element %u has %u nodes, return vector has %u - bye!\n", elem, this->getnNodes(), (idx) nodes.size() );
        exit(1);
    }
#endif
}


void Mesh::CalculateDeterminants3D(double *p)
{

    if (Determinant !=NULL) free(Determinant);
    Determinant = (double*) malloc( getnElements()*sizeof(double));

    if ((Determinant != NULL) && (Dimension == 3) && (getnElements() >0) )
    {

        TotalSize = 0.0; // reset total mesh volume
        printf("calculating %i 3D determinants",getnElements());

        int fraction = (int) getnElements() /10;
        //#pragma omp parallel for // doesn't do progress dots... .. .
        for (idx i = 0 ; i < getnElements() ; i++)
        {
            idx* n;
            if (fraction!=0){ if (i%fraction==0) printf(".");}

            n = &Elem[i*nNodes];

            double x[4], y[4], z[4];
            x[0] = p[3*n[0]];   y[0] = p[3*n[0]+1]; z[0] = p[3*n[0]+2];
            x[1] = p[3*n[1]];   y[1] = p[3*n[1]+1]; z[1] = p[3*n[1]+2];
            x[2] = p[3*n[2]];   y[2] = p[3*n[2]+1]; z[2] = p[3*n[2]+2];
            x[3] = p[3*n[3]];   y[3] = p[3*n[3]+1]; z[3] = p[3*n[3]+2];

            double A[3],B[3],C[3];
            A[0] = x[1] - x[0]; A[1] = y[1] - y[0]; A[2] = z[1] - z[0];
            B[0] = x[2] - x[0]; B[1] = y[2] - y[0]; B[2] = z[2] - z[0];
            C[0] = x[3] - x[0]; C[1] = y[3] - y[0]; C[2] = z[3] - z[0];

            // B x C
            double cross[3];
            cross[0] =   B[1]*C[2] - C[1]*B[2];
            cross[1] = - B[0]*C[2] + C[0]*B[2];
            cross[2] =   B[0]*C[1] - C[0]*B[1];

            // determinant = A . (B x C)
            Determinant[i] = A[0]*cross[0] + A[1]*cross[1] + A[2]*cross[2];
            if (Determinant[i] != Determinant[i] )             printf("det[%i] = %e\n", i, Determinant[i]);
            if (Determinant[i] < 0 ){ // reorder nodes to ensure positive determinant
                Determinant[i] = -1.0 * Determinant[i];
            }

            TotalSize+= Determinant[i] ; // add determinat of this element to total

            if (Determinant[i] == 0 )
            {
                printf("warning - Mesh::CalculateDeterminants3D, Zero determinant in  element %i \n",i);
                printf("Nodes:\n");
                printf("%u = [%f,%f,%f]\n",n[0], x[0], y[0], z[0]);
                printf("%u = [%f,%f,%f]\n",n[1], x[1], y[1], z[1]);
                printf("%u = [%f,%f,%f]\n",n[2], x[2], y[2], z[2]);
                printf("%u = [%f,%f,%f]\n",n[3], x[3], y[3], z[3]);
                printf("have a nice day ! \n");
                exit(1);
            }

        }
        TotalSize = TotalSize / 6.0; // scale tet determinant to volume
        printf("OK\n");
    }
    else
    {
        printf("error - Mesh::CalculateDeterminants3D - Determinats have already been calculated\n");
        exit(1);
    }
}
void Mesh::CalculateSurfaceNormals(double *p, Mesh* tets)
{
    // CALCULATES THE SURFACE NORMAL FOR A TRIANGULAR MESH
    if ((Dimension != 2) ||
            (getnElements() <1) ||
            (getnNodes() <1))
    {
        printf("error - Mesh::CalculateSurfaceNormals(double *p) - cannot calculate surface normals, bye!");
        exit(1);
    }

    printf("calculating %i surface normals and 2D determinants ",getnElements());
    if (SurfaceNormal != NULL) free(SurfaceNormal);
    if (Determinant != NULL) free(Determinant);

    SurfaceNormal = (double*) malloc(3 * getnElements() * sizeof(double));
    Determinant = (double*) malloc( getnElements() * sizeof(double) );

    if ((SurfaceNormal==NULL) || (Determinant==NULL))
    {
        printf("error - MeshCalculateSurfaceNormals(dobule*) - could not allocate memory- bye!\b");
        exit(1);
    }


    // Loop over each element and calculate determinant and surface normal
    double Ax,Ay,Az, Bx,By,Bz;
    double cross[3];
    int fraction = (int) getnElements() / 10;
    TotalSize = 0;

    for ( idx i = 0 ; i < getnElements() ; i ++)
    {
        if (  i%fraction  == 0 ) // progress dots printing
            printf(".");

        Ax = p[getNode(i,1)*3+0] - p[getNode(i,0)*3 + 0];
        Ay = p[getNode(i,1)*3+1] - p[getNode(i,0)*3 + 1];
        Az = p[getNode(i,1)*3+2] - p[getNode(i,0)*3 + 2];
	
	
        Bx = p[getNode(i,2)*3 + 0] - p[getNode(i,0)*3 + 0];
        By = p[getNode(i,2)*3 + 1] - p[getNode(i,0)*3 + 1];
        Bz = p[getNode(i,2)*3 + 2] - p[getNode(i,0)*3 + 2];
	
        cross[0] = Ay*Bz - Az*By;
        cross[1] = Az*Bx - Ax*Bz;
        cross[2] = Ax*By - Ay*Bx;
        double det = sqrt( cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

        if (det < 0 )
        { // if det is negative make it not...
#ifdef DEBUG
            printf("triangle det[%i] < 0 , det = %f\n", i , det);
#endif
            det = -1.0*det;
        }

        setDeterminant(i,  det );
        TotalSize+= det / 2.0;
        //normalise cross
        cross[0] = cross[0] / det;
        cross[1] = cross[1] / det;
        cross[2] = cross[2] / det;

        // CHEK THAT SURFACE NORMAL POINTS TOWARDS LC REGION. I.E. INWARDS
        // IF NOT REVERSE ITS ORIENTATION

        if ((ConnectedVolume) && (tets) ) // only if connected volumes have been set and tets are provided
        {
            int t = getConnectedVolume( i ); // index to neighbouring tet
            if ( t > -1 )
            { // TRIANGLE IS NOT CONNECTED TO AN LC. PROBABLY A DIELECTRIC NRIGHBOUT ONLY TRIANGLE
                // triangle barycentre
                double tri_bary[3] = { ( p[getNode(i,0)*3+0] + p[getNode(i,1)*3+0] +p[getNode(i,2)*3+0] ) / 3.0 ,  // x
                                       ( p[getNode(i,0)*3+1] + p[getNode(i,1)*3+1] +p[getNode(i,2)*3+1] ) / 3.0 ,  // y
                                       ( p[getNode(i,0)*3+2] + p[getNode(i,1)*3+2] +p[getNode(i,2)*3+2] ) / 3.0 }; // z
                // tet barycentre
                int n[4] = { tets->getNode(t,0) , tets->getNode(t,1) , tets->getNode(t,2), tets->getNode(t,3) };
                double tet_bary[3] = {  (p[n[0]*3+0] + p[n[1]*3+0] + p[n[2]*3+0] + p[n[3]*3+0] )/ 4.0 ,	// x
                                        (p[n[0]*3+1] + p[n[1]*3+1] + p[n[2]*3+1] + p[n[3]*3+1] )/ 4.0 ,	// y
                                        (p[n[0]*3+2] + p[n[1]*3+2] + p[n[2]*3+2] + p[n[3]*3+2] )/ 4.0 };	// z

                // vector from tri barycentre to tet barycentre
                double v1[3] = {tet_bary[0] - tri_bary[0],
                                tet_bary[1] - tri_bary[1],
                                tet_bary[2] - tri_bary[2]};

                // dot product between surface normal and v1 determines orientation of surface normal
                double dot = cross[0]*v1[0] + cross[1]*v1[1] + cross[2]*v1[2];

                if (dot < 0.0 )
                {
                    // IF SURFACE NORMAL IS IN WRONG DIRECTION,
                    // REVERSE NODE ORDER AND INVERT VECTOR DIRECTION
                    idx* e = &Elem[i * getnNodes() ]; // pointer to first node in this element
                    idx temp = e[2];
                    e[2] = e[1];
                    e[1] = temp;

                    cross[0]*=-1.0;
                    cross[1]*=-1.0;
                    cross[2]*=-1.0;
                }
            }// end if index to neighbouring tet exists

        } // if connected volume found
        else
        {
            printf("error in Mesh::CalculateSurfaceNormals - Connected volumes or tets do not exist - bye\n");
            exit(1);
        }
        // if (getMaterialNumber(i)>4000)
        //     printf("surf[%i] = %e,%e,%e, mat = %i\n", i, cross[0], cross[1], cross[2], getMaterialNumber(i) );


        setSurfaceNormal(i,cross);
    }
    printf("OK\n");



}

void Mesh::CalcLocCoords(const idx elem, double *p, double *coord, double *loc)
{
    // calculates 4 local coordinates of coordinate coord in element elem
    // returned in loc
    idx* n = this->getPtrToElement( elem );

    double M[12]= {0,0,0,0,0,0,0,0,0,0,0,0};
    //cout << "loc: ";
    for (int c = 0 ; c < 4 ; c++){
        for (int i = 0 ; i < 4 ; i++)
        { // node number
            M[i*3 + 0] = p[ 3*n[i] + 0]; // x
            M[i*3 + 1] = p[ 3*n[i] + 1]; // y
            M[i*3 + 2] = p[ 3*n[i] + 2]; // z
        }

        M[c*3 + 0 ] = coord[0];
        M[c*3 + 1 ] = coord[1];
        M[c*3 + 2 ] = coord[2];

        double det = fabs( this->Calculate4x4Determinant(M) ) ;
        loc[c] = det / (this->getDeterminant(elem)*1e18 );

        //cout << loc[c] << " " ;
    }
    //cout << endl;

}

void Mesh::CalcElemBary(const idx elem, double *p, double *bary)const
{
#ifdef DEBUG
    if (Dimension == 2)
    {
        cout << "error - CalcElemBary currently only works on tets - bye!" << endl;
        exit(1);
    }
#endif
    idx* n = Elem + elem*4; // shortcut

    double* p1 , *p2, *p3, *p4;
    p1 = p + n[0]*3;
    p2 = p + n[1]*3;
    p3 = p + n[2]*3;
    p4 = p + n[3]*3;

    bary[0] = ( p1[0] + p2[0] + p3[0] + p4[0] ) / 4.0;
    bary[1] = ( p1[1] + p2[1] + p3[1] + p4[1] ) / 4.0;
    bary[2] = ( p1[2] + p2[2] + p3[2] + p4[2] ) / 4.0;

}

double Mesh::CalcBaryDistSqr(double *p, const idx elem, double *coord) const
{
// returns the squared distance between coordinate coord[3] and the barycentre of
// element elem. p is pointer to all node coordinates.

    double bary[3] = {0,0,0};
    CalcElemBary( elem , p, bary);


    return	( (bary[0]-coord[0])*(bary[0]-coord[0]) +
                  (bary[1]-coord[1])*(bary[1]-coord[1]) +
                  (bary[2]-coord[2])*(bary[2]-coord[2]) );
}

void Mesh::listElementsOfMaterial(std::vector<idx> &elems, const idx mat) const
{
    elems.clear();
    elems.reserve( this->nElements );

    for (idx i = 0 ; i < nElements ; i++)
    {
        if ( ( this->Mat[i] &  mat) == mat) // bitwise comparison FIXLC1+ELECTRODE1 & FIXLC1 -> FIXLC1
        {
            elems.push_back( i );
        }
    }

}

void Mesh::listNodesOfMaterial(std::vector<idx> &nodes, const idx mat) const
{
    // MAKES A LIST OF NODES BLEONGING TO ELEMENTS OF MATERIAL mat
    nodes.clear();

    for (idx i = 0 ; i < getnElements() ; i++)
    {
        if ( (this->getMaterialNumber(i) & mat ) == mat )
        {
            for (idx n = 0 ; n < this->getnNodes() ; n++)
                nodes.push_back( (unsigned int) this->getNode(i, n) );
        }
    }

    // REMOVE REPEATED ENTRIES
    sort(nodes.begin() , nodes.end() );
    std::vector< unsigned int> :: iterator u;
    u = unique( nodes.begin() , nodes.end() );
    nodes.erase( u , nodes.end() );

}

/*
void Mesh::FindIndexToMaterialNodes(idx mat_num, vector<idx> *index) const
{
    index->clear();
    vector <idx> ind_t;

    if (mat_num == 0)
    {
        printf("error - Mesh::FindIndexToMaterialNodes - a zero material number will result in inf, bye!");
        exit(1);
    }


    for (idx i = 0; i < nElements; i++)
    {// loop over all elements
        if (mat_num >= MAT_FIXLC1){ // check which FIXLC number
            int fix = mat_num / MAT_FIXLC1;
            if ((Mat[i]&31*MAT_FIXLC1)/MAT_FIXLC1 == fix)
                ind_t.push_back(i);
        }
        else
            if (mat_num >= MAT_ELECTRODE1){ // check which ELECTRODE number
                int el = mat_num / MAT_ELECTRODE1;
                if  ( (Mat[i]&31*MAT_ELECTRODE1)/MAT_ELECTRODE1 == el )
                    ind_t.push_back(i);
            }
            else{
                printf("Mesh::FindIndexToMaterialNodes cannot at the moment search for this material :( \n - bye\n");
                exit(1);
            }
    }// end loop over all elements

    if (ind_t.size()>0){// elements found?
        // then extract all nodes corresponding to elements in ind_t and make list unique
        list<int> ind_p;
        vector <int>::iterator i;

        for (i = ind_t.begin(); i != ind_t.end(); i++){// loop over elements
            for (int j= 0; j < nNodes ; j++){ // loop over nodes per element
                ind_p.push_back(Elem[(*i)*nNodes + j] );
            }
        }//end loop over elements
	
        //sort and remove repeated nodes
        ind_p.sort();
        ind_p.unique();
	
        list <int>::iterator liter;
        for (liter = ind_p.begin(); liter != ind_p.end(); liter++)
            index->push_back(*liter);

    }// if nodes found
}
*/

void Mesh::CopySurfaceNormal(idx i, double* norm)
{
#ifdef DEBUG
    if ( (i >= nElements) || (i < 0) || (SurfaceNormal==NULL) || (Dimension != 2))
    {
        printf("error - Mesh::CopySurfaceNormal(int, double*) - bad, bye!\n");
        printf("i = %i, nElements = %i, Dimension = %i",i,nElements,Dimension);

        if (SurfaceNormal == NULL)printf(" SurfaceNormal is unfortunately NULL!\n");
        exit(1);
    }
#endif

    norm[0] = SurfaceNormal[i*3+0];
    norm[1] = SurfaceNormal[i*3+1];
    norm[2] = SurfaceNormal[i*3+2];
}

void Mesh::removeElements(std::set <idx>& index)
{
    // REMOVES ELEMENTS IN INDEX.
    if ( !index.size() )
    {
        return ;
    }

    idx num_new_elements = getnElements() - (idx) index.size();

    // make sure index does not contain elements that do not exist
    idx maxindex = * max_element( index.begin() , index.end() );

    if ( maxindex >= getnElements() )
    {
        printf("error - Mesh::removeElements \ntrying to remove elements that do not exist\n");
        printf("size of mesh is %u, trying to remove element %u - bye! ",getnElements() , *(index.rbegin()));
        exit(0);
    }

    idx* tnew 		= (idx*) malloc( num_new_elements * getnNodes() * sizeof(idx) ); // new elements
    idx* mnew 		= (idx*) malloc( num_new_elements * sizeof(idx) );		// new materials
    double* dnew	= (double*) malloc( num_new_elements * sizeof(double) ); 	// new determinants

    double* SurfaceNormal_new 	= NULL;  // only needed for triangle meshes
    idx*    ConnectedVolume_new	= NULL;

    if (getDimension() == 2) // if a triangle mesh
    {
        SurfaceNormal_new   = (double*) malloc( num_new_elements * 3 * sizeof(double) );
        ConnectedVolume_new = (idx*) malloc( num_new_elements * sizeof(idx) );
    }

    if ( (tnew==NULL) || (mnew == NULL) || (dnew == NULL) )
    {
        printf("error - Mesh::removeElements - received null pointer - bye!\n");
        exit(1);
    }

    std::set <idx> :: iterator itr;
    itr = index.begin();

    size_t c = 0; // counter to new elements

    for (idx i = 0 ; i < this->getnElements() ; i ++) // for all elements
    {

        if (*itr != i) // KEEP ELEMENT I
        {
            mnew[c] = getMaterialNumber(i);
            dnew[c] = getDeterminant(i);

            if( getDimension() == 2) // if triangles
            {
                ConnectedVolume_new[c] = getConnectedVolume(i);
                SurfaceNormal_new[c*3+0] = *(getPtrToSurfaceNormal(i) 	); // normal X
                SurfaceNormal_new[c*3+1] = *(getPtrToSurfaceNormal(i)+1	); // normal Y
                SurfaceNormal_new[c*3+2] = *(getPtrToSurfaceNormal(i)+2 ); // normal Z
            }

            // COPY NODES
            for (idx j = 0 ; j < getnNodes() ; j ++ )
            {
                tnew[c*getnNodes() + j ] = getNode(i,j);
            }
            c++;

        }// end if
        else if (*itr != *index.rbegin() ) // if not end of list, increment (ugly comparison of values instead of pointers)
        {
            itr++;
        }
    }// end for all elements


    if (Elem!=NULL)	{
        free(Elem);
        Elem = tnew;
    }
    if (Mat != NULL){
        free(Mat);
        Mat = mnew;
    }
    if (Determinant != NULL ){
        free(Determinant);
        Determinant = dnew;
    }

    if (getDimension() == 2){ // if tris , set new surface normals and connected to volume links
        if (SurfaceNormal != NULL) {
            free(SurfaceNormal);
            SurfaceNormal = SurfaceNormal_new;
        }
        if (ConnectedVolume != NULL){
            free(ConnectedVolume);
            ConnectedVolume = ConnectedVolume_new;
        }
    }// end if tris

    setnElements(num_new_elements);


}

void Mesh::ClearMesh(){
    Dimension 	= 0;
    nElements	= 0;
    nNodes		= 0;

    if (Elem != NULL) free(Elem);
    Elem = NULL;

    if (Mat != NULL) free(Mat);
    Mat = NULL;

    if (Determinant != NULL ) free(Determinant);
    Determinant = NULL;

    if (ConnectedVolume != NULL) free(ConnectedVolume);
    ConnectedVolume = NULL;

    if (SurfaceNormal != NULL ) free(SurfaceNormal);
    SurfaceNormal = NULL;
}

void Mesh::addElements(vector<idx> &m_new, vector<idx> &mat_new)
{

    // CHECK THAT EQUAL SIZE OF MATERIAL NUMBERS AND ELEMENT NODES
    idx n_p_e = (unsigned int) this->getnNodes();
    idx num_new = m_new.size() / n_p_e;

    if ( num_new != mat_new.size() )
    {
        cout << "error - Mesh::addElements(vector<idx>&,vector<idx>&) " << endl;
        cout << "number of elements does not match number of material numbers - bye!" << endl;
        exit(1);
    }

    // ALLOCATE MEMORY FOR NEW ELEMENTS AND MAT NUMS
    idx* nelem = (idx*) malloc( (getnElements() + m_new.size() ) * sizeof (idx) * n_p_e );
    idx* nmat  = (idx*) malloc( (getnElements() + mat_new.size() ) * sizeof(idx) );

    // COPY OLD VALUES
    memcpy( nelem, Elem, n_p_e*getnElements() * sizeof(idx) );
    memcpy( nmat , Mat , getnElements() * sizeof(idx) );

    // ADD NEW VALUES
    for (idx i = 0 ; i < (idx) mat_new.size()  ; i++)
    {
        nmat [ i+ getnElements()] = mat_new[i];
        for (idx j = 0 ; j < n_p_e ; j++)
        {
            nelem[ getnElements()*n_p_e + i*n_p_e + j] = m_new[ i*n_p_e + j ];
        }
    }

    // REPLACE OLD
    if (Elem) free (Elem);
    Elem = nelem;

    if (Mat) free (Mat);
    Mat = nmat;

    // UDPDATE TOTAL NUMBER OF ELEMENTS
    this->setnElements( getnElements() + num_new );

}



void Mesh::CopyMesh(Mesh* rhs)
{
    setDimension(rhs->getDimension() );
    setnElements(rhs->getnElements() );		// set number of elements
    setnNodes(	rhs->getnNodes() );             // number of nodes per element
    TotalSize = rhs->TotalSize;

    setMaxNodeNumber( rhs->getMaxNodeNumber() );

    AllocateMemory();				// allocate memory for arrays
    setAllNodes(rhs->getPtrToElement(0) );	// copy node numbers
    setAllMaterials(rhs->getPtrToMaterialNumber(0) );// copy material numbers

    // COPY DETERMINANTS
    if (rhs->getPtrToDeterminant(0) != 0)
    {
        if (Determinant!=NULL) free(Determinant);
        Determinant = (double*) malloc( getnElements() * sizeof(double) );
        memcpy( Determinant, rhs->Determinant, getnElements() * sizeof(double) );
    }

    this->Dimension = rhs->Dimension;

    if (Dimension==2) // if triangle mesh
    {
        // COPY TRIANGLE/TET CONNECTIONS
        if(rhs->getPtrToConnectedVolume(0) != NULL)
        {
            if (ConnectedVolume != NULL) free(ConnectedVolume);
            ConnectedVolume = (idx*) malloc( getnElements() * sizeof(idx) );

            memcpy(ConnectedVolume, rhs->ConnectedVolume , getnElements() * sizeof(idx) );

        }
	
        // COPY SURFACE NORMALS
        if (rhs->getPtrToSurfaceNormal(0) != NULL)
        {
            if (SurfaceNormal != NULL) free(SurfaceNormal);
            SurfaceNormal = (double*) malloc(3 * getnElements() * sizeof(double));

            memcpy( SurfaceNormal, rhs->SurfaceNormal, 3*getnElements() * sizeof(double) );
        }
    } // end if triangle mesh
}

void Mesh::ScaleDeterminants( const double& s)
{
    if (Determinant!= NULL){
        for (idx i = 0 ; i < getnElements() ; i ++)
            Determinant[i] = Determinant[i]*s;

        TotalSize = TotalSize * s;
    }
    else
    {
        printf("error - Mesh::ScaleDeterminants - Determinants not initialised - bye !\n");
        exit(1);
    }

}

void Mesh::PrintElements()const
{
    printf("printing %i elements - with %i nodes:\n", nElements, nNodes);
    if (nElements>0)
    {
        idx i,j;
	
        for (i=0;i<nElements;i++)
        {
            printf("%u -- ",i);

            for(j = 0 ; j < nNodes ; j++)
                printf(" %u\t",Elem[i*nNodes +j]);

            printf("- mat = %u",Mat[i]);

            if (SurfaceNormal != NULL)
                printf("- normal = [%1.2f,%1.2f,%1.2f]", SurfaceNormal[3*i] , SurfaceNormal[3*i+1] , SurfaceNormal[3*i+2]);

            if (this->Determinant != NULL)
                printf(" det =  %e", this->getDeterminant(i) );

            printf("\n");
        }
    }
}

void Mesh::PrintElement(idx e) const
{
    if (e < nElements)
    {
        if (getDimension()==2){
            printf("tri[%u] = [", e);
        }
        else if (getDimension()==3){
            printf("tet[%u] = [", e);
        }
        for (idx j = 0 ; j < nNodes ; j ++)
            printf(" %u",Elem[e*nNodes + j]);
        printf("]");
        if (Mat){
            printf(" Mat = %u", getMaterialNumber( e ) );
        }
        if (ConnectedVolume){
            printf(" Connected to tet[%u]" , this->getConnectedVolume( e));
        }
        if (this->Determinant){
            printf(" Determinant %e", this->getDeterminant( e ));
        }
        printf("\n");
    }
    else
    {
        printf("cannot print element %u , nElements is %u\n",e,nElements);
    }
}

void Mesh::PrintNormals() const
{
    // PRINTS ALL SURFACE NORMALS, IF PRESENT
    if (!this->SurfaceNormal)
        return;
    for (idx i = 0 ; i < getnElements() ; i++)
    {
        double* nn = &SurfaceNormal[i*3];
        printf("normal[%u] = [%e,%e,%e]\n", i , nn[0], nn[1], nn[2] );
    }
}

bool Mesh::isOnXPlane(idx e, double X, double* p) const
{
    for (idx i = 0 ; i < getnNodes() ; i ++)
    {
        if (p[3*getNode(e,i)+0] != X)
            return false;
    }
    return true;
}

bool Mesh::isOnYPlane(idx e, double Y, double* p) const
{
    for (idx i = 0 ; i < getnNodes() ; i ++)
    {
        if (p[3*getNode(e,i)+1] != Y)
            return false;
    }
    return true;
}

bool Mesh::isOnZPlane(idx e, double Z, double* p) const
{ //returns true if all nodes are on ZPlane z, else false
    for (idx i = 0 ; i < getnNodes() ; i ++)
    {
        if (p[3*getNode(e,i)+2] != Z)
            return false;
    }
    return true;
}

bool Mesh::isNeighbours(const idx el1, const idx el2) const
{
    set < idx > nodes;

    // make set of combined nodes of both elements.
    for (idx i = 0 ; i < getnNodes() ; i++)
    {
	nodes.insert( getNode( el1, i) );
	nodes.insert( getNode( el2, i) );
    }
    // COMBINED NODE NUMBERS LIST SIZE WILL BE nNodes + 1, IF
    // BOTH ELEMENTS SHARE A FACE (tets) OR A LINE (tri)
    if ( (idx) nodes.size() == getnNodes() + 1)
	return true;
    else
	return false;
}

void Mesh::gen_p_to_elem(vector < set < idx > > &p_to_elem) const
{
    p_to_elem.clear();
    p_to_elem.reserve( this->getMaxNodeNumber() ); // pre-allocate space

    // initialise vector to allow [index] access later
    for ( idx i = 0 ; i < this->getMaxNodeNumber()+1 ; i++)
    {
        set <idx> empty;
        p_to_elem.push_back( empty );
    }

    for (idx i = 0 ; i <  this->getnElements() ; i++) // for every element
    {
        for (idx j = 0 ; j < this->getnNodes() ; j++) // for every node
        {
            idx n = this->getNode( i , j );
            //printf("inserting element %u to node %u\n",i, n);
            p_to_elem[n].insert( i );
        }
    }

}


void Mesh::gen_neighbour_list( vector<vector<idx> >& neigh) const
{
    vector < set <idx> > p_to_t;
    this->gen_p_to_elem( p_to_t );

    // LOOP OVER EACH ELEMENT AND FIND ITS NEIGHBOURS
    neigh.clear();
    neigh.reserve( this->getnElements() );

    for (idx i = 0 ; i < this->getnElements() ; i++)
    { // pre-allocate
        vector <idx> empty;
	neigh.push_back( empty );
    }
#ifndef DEBUG
#pragma omp parallel for
#endif
    for (idx elem = 0 ; elem < (unsigned int) getnElements() ; elem ++)// for all elements
    {
	// GET ALL NODE NEIGHBOURS, INCLUDIGN THOSE THAT ONLY SHARE A SINGLE NODE
        set <idx> tets;
        for (idx j = 0 ; j < getnNodes() ; j++) // for all nodes j
        {
            idx n = this->getNode( elem , j );
	    tets.insert(p_to_t[n].begin() , p_to_t[n].end() );
	}// end for all nodes j


	// CHECK WHICK OF THE POSSIBLE ELEMENTS ARE FACE/EDGE NEIGHBOURS
        set <idx> ::iterator itr;
        vector <idx> neighs;
        for (itr = tets.begin() ; itr != tets.end() ; itr++) // for all node neighbours
        {
            if (isNeighbours( elem , *itr) )
            {
		neighs.push_back( *itr );
                if ((idx) neighs.size() == getnElements() )
		    break;
	    }
	}// end for all node neighbours

        // REORDER NEIGHBOURS SO THAT FIRST NEIGHBOUR OF TRI A IS TRI B
        //      2   AND THIRD NEIGHBOUR OF TRI B IS TRI A ( WHEN A= [1,2,3] , B=[2,3,4])
        //     /|\  I.E., NEIGHBOURS ARE ORDERED W.R.T. TO OPPOSITE 'LONELY' NODE.
        //    / | \																			*/
        //  1/  |  \4																		*/
        //   \A | B/																		*/
        //    \ | /																			*/
        //     \|/																			*/
        //*      3																			*/

        vector <idx > no; // ordered neighbours
	no.assign(getnNodes() , getnElements() ); // initial value for all nodes is 1 too much allowing for checking of non-existent neighbours at boundaries

        for (idx i = 0 ; i < getnNodes() ; i++)// loop over each node of this element
        {
            idx n = getNode(elem, i);
            for (idx i2 = 0 ; i2 < (idx) neighs.size(); i2++)
            {
		// if neighs[i2] doesn't contain node i of this element, it is
		// the i'th neighbour
                if (!ContainsAllNodes( neighs[i2] , 1 , &n ) )
                {
		    no[i] = neighs[i2];
		    break;
		}
	    }// end for i2
	}// end for i, all nodes of this tet

	neigh[elem] =  no ;
    }// end for all elements

}

