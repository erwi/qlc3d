#include <qlc3d.h>
#include <mesh.h>
#include <material_numbers.h>
#include <fmt/format.h>
#include <util/logging.h>
#include <util/exception.h>
#include <set>

using fmt::format;
Mesh::Mesh() {
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

Mesh::~Mesh() {
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
Mesh::Mesh(idx numElements, idx numNodes) {
    nElements = numElements;
    nNodes = numNodes;
    // allocate memory for mesh
    Elem = (idx*)malloc(nElements * nNodes * sizeof(idx));
    Mat  = (idx*)malloc(nElements * sizeof(idx));
    Determinant = (double*)malloc(nElements * sizeof(double));
    SurfaceNormal = NULL;
    ConnectedVolume = NULL;
}
void Mesh::AllocateMemory() {
    if ((nElements > 0) && (nNodes > 0)) {
        if (Elem != NULL) free(Elem);
        if (Mat != NULL) free(Mat);
        if (Determinant != NULL) free(Determinant);
        Elem = (idx*) malloc(nElements * nNodes * sizeof(idx) );
        Mat  = (idx*) malloc(nElements * sizeof(idx) );
        Determinant = (double*) malloc(nElements * sizeof(double) );
    }
    else {
        RUNTIME_ERROR(format("Expected non-zero counts, got nElements={}, nNodes={}.", nElements, nNodes));
    }
}

idx Mesh::getMaterialNumber(const idx e) const {
#ifdef DEBUG
    assert(Mat);
    assert(e < nElements);
#endif
    return Mat[e];
}

idx Mesh::getFixLCNumber(const idx e) const {
    return (idx) getMaterialNumber(e) / MAT_FIXLC1;
}

idx Mesh::getDielectricNumber(const idx e) const {
    return (int) getMaterialNumber(e) / MAT_DIELECTRIC1;
}

double Mesh::getDeterminant(const idx i) const {
    // RETURNS THE PRE-CALCULATED DETERMINANT FOR ELEMENT i
#ifdef DEBUG
    assert(Determinant);
    assert(i < nElements);
#endif
    return Determinant[i];
}


idx Mesh::getConnectedVolume(const idx e) const {
    // RETURNS INDEX TO CONNECTED LC TETRAHEDRON, WHEN e IS INDEX TO A TRIANGLE
#ifdef DEBUG
    assert(ConnectedVolume);    // CONNECTIONS HAVE BEEN INITIALISED
    assert(e<getnElements() );  // IN BOUNDS
#endif
    return ConnectedVolume[e];
}
idx* Mesh::getPtrToElement(const idx e) const {
    return &Elem[e*nNodes];
}
idx* Mesh::getPtrToMaterialNumber(const idx e) const {
    return &Mat[e*nNodes];
}

idx* Mesh::getPtrToConnectedVolume(const idx e) const {
    return &ConnectedVolume[e];
}

double* Mesh::getPtrToDeterminant(const idx e) const {
    return &Determinant[e];
}

double* Mesh::getPtrToSurfaceNormal(const idx e) const {
    return &SurfaceNormal[e*3];
}

void Mesh::setDeterminant(const idx i, double det) {
#ifdef DEBUG
    assert(Determinant);
    // ADD CHECKING OF i TOO LARGE OR SMALL
#endif
    Determinant[i] = det;
}
void Mesh::setSurfaceNormal(idx i, double norm[3]) {
#ifdef DEBUG
    assert(i < getnElements());
    assert(SurfaceNormal != nullptr);
#endif
    SurfaceNormal[i*3 + 0] = norm[0];
    SurfaceNormal[i*3 + 1] = norm[1];
    SurfaceNormal[i*3 + 2] = norm[2];
}

/**
 * copies all values from array nodes to this->Elem
 * @param nodes
 */
void Mesh::setAllNodes(idx *nodes) {
#ifdef DEBUG
    assert(nElements > 0);
    assert(Elem != nullptr);
#endif
    memcpy(Elem, nodes, nElements * nNodes * sizeof(idx));
}

void Mesh::setAllMaterials(idx *mat) {
#ifdef DEBUG
    assert(Mat);
#endif
    memcpy( Mat, mat, nElements*sizeof(idx) );
}

void Mesh::setDimension(const idx dim) {
    Dimension = dim;
}

void Mesh::setnElements(idx nelem) {
    nElements = nelem;
}

void Mesh::setnNodes(idx nnodes) {
    nNodes = nnodes;
}

void Mesh::setConnectedVolume(Mesh* vol) {
    // finds indexes from triangles to tetrahedra in 'vol' that share 3 points.
    // these are needed for surface integrals, e.g. Neumann boundaries.
    // counts number of times volumes are connected to each node
    // then compares with nodes in surface mesh

    // Check that surface and volume mesh dimensions are OK
    if ((getDimension() != 2) ||
         (vol->getDimension() != 3)) {
        RUNTIME_ERROR(format("Unable to associate triangles with tetrahedra. Triangles dimension = {}, "
                             "tetrahedra dimension = {}", this->getDimension(), vol->getDimension()));
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

    for (idx i = 0 ; i < this->getnElements() ; i++){ // loop over surface elements
        if (( getMaterialNumber( i ) < MAT_ELECTRODE1) ||     // exclude eletrode only surfaces as
            (  getMaterialNumber( i ) > MAT_ELECTRODE9) ){    // these are fixed, and do not need to know
                                                                  // about neighbouring LC tets
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

            if (final.empty()) {
                auto msg = format("Triangle {} with nodes [{}, {}, {}] is not connected to any tetrahedron in {}, {}",
                                  i, getNode(i, 0), getNode(i, 1), getNode(i, 2), __FILE__, __func__);
                throw std::runtime_error(msg);
            }

            // SURFACE CAN BE BETWEEN LC AND DIELECTRIC LAYERS
            // THEN final.size() IS 2. MUST CHECK WHICH ONE IS LC
            // AND WHICH IS DIELECTRIC
            set <idx> ::iterator itr;
            itr = final.begin();
            if ( final.size() == 1){ // IF SINGLE TET FOUND
                if (vol->getMaterialNumber( *itr ) <= MAT_DOMAIN7 ){ // only connect to an LC element
                    ConnectedVolume[i] = *itr;
                }
                else if ( ( vol->getMaterialNumber( *itr) >= MAT_DIELECTRIC1 ) &&
                          (vol->getMaterialNumber(*itr) <= MAT_DIELECTRIC7) ){
                    // MAKE SURE THAT WE DONT HAVE A FIXLC SURFACE THAT IS ONLY CONNECTED TO A DIELECTRIC TET
                    idx surfMat = this->getMaterialNumber(i);
                    surfMat = MATNUM_TO_FIXLC_NUMBER((size_t) surfMat);
                    if ((surfMat > 0) && (surfMat < 10)) {
                        throw std::runtime_error(format("Found a FxLC{} surface only connected to a dielectric "
                                                        "volume, but no LC volume in {}, {}.", surfMat, __FILE__, __func__));
                    }

                    ConnectedVolume[i] = -2; // set to -2 if connected to dielectric
                }
            }
            else
                if ( final.size() == 2){
                    // CHECK BOTH FIRST AND SECOND INDICES
                    if (vol->getMaterialNumber( *itr ) <= MAT_DOMAIN7 )// only connect to an LC element
                        ConnectedVolume[i] = *itr;
                    else{
                        ++itr;
                        if (vol->getMaterialNumber( *itr ) <= MAT_DOMAIN7 )// only connect to an LC element
                            ConnectedVolume[i] = *itr;
                    }
                }
                else {
                    throw std::runtime_error(format("Triangle {} is connected to {} tetrahedra, but expected 2, "
                                                    "in {}, {}", i, final.size(), __FILE__, __func__));
                }
        } // end exclude electrode only surfaces
    }
}
double Mesh::Calculate4x4Determinant(double* M) const {
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

void Mesh::ContainsNodes(list <idx>* elems, list <idx>* points) {
    // CREATES INDEX OF ELEMENTS THAT CONTAIN AT LEAST ONE OF THE NODES IN LIST points

    list <idx>::iterator int_itr;
    for (idx i = 0 ; i < getnElements() ; i ++ ) // loop over all elements
    {
        for (int_itr = points->begin() ; int_itr != points->end() ; ++int_itr) // loop over all points in list
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

bool Mesh::ContainsCoordinate(const idx elem, const double *p, const double *coord ) const {
    // CHECKS TO SEE WHETHER ELEMENT elem CONTAINS COORDINATE coord. THAT IS, COORD
    // CAN BE ANYWHERE INSIDE elem, NOT JUST CORNER NODES
    // coord is expected to be of length 3 [x,y,z]. is most likely a pointer to somewhere in array *p
#ifdef DEBUG
    if (this->Dimension == 2) {
        RUNTIME_ERROR("This method only works on tetrahedra.");
    }
#endif

    // CALCULATE 4 DETERMINANTS FOR TETS FORMED USING NODES FROM THIS ELEMENT
    // AND THE TESTED COORDINATE AND ONE USING THIS TET ONLY. IF SUM OF DETS OF TETS FORMED USING coord
    // EQUALS THAT OF FOR THIS ELEMENT ONLY, THEN THE COORDINATE IS WITHIN THIS TET
    const double eps = qlc3d_GLOBALS::GLOBAL_COORD_EPS; // accuracy used for coordinate comparions
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
    if (nodes.size() != (unsigned int) this->getnNodes() ) {
        RUNTIME_ERROR(format("Return size does not match the number of nodes in element. "
                      "Element {} has {} nodes, return vector has {}", elem, getnNodes(), nodes.size()));
    }
#endif
}


void Mesh::CalculateDeterminants3D(double *p) {

    if (Determinant !=NULL) free(Determinant);
    Determinant = (double*) malloc( getnElements()*sizeof(double));

    if ((Determinant != NULL) && (Dimension == 3) && (getnElements() >0) ) {

        TotalSize = 0.0; // reset total mesh volume
        Log::info("Calculating {} 3D determinants.", getnElements());

        for (idx i = 0 ; i < getnElements() ; i++) {
            idx* n;

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

            if (Determinant[i] < 0 ) { // should reorder nodes to ensure positive determinant
                Log::warn("Negative determinant for element {}, negating sign and continuing.", i);
                Determinant[i] = -1.0 * Determinant[i];
            }

            TotalSize += Determinant[i] ; // add determinat of this element to total

            if (Determinant[i] == 0) {
                RUNTIME_ERROR(format("Zero determinant in element {} with "
                                     "nodes: {} = ({}, {}, {}), {} = ({}, {}, {}), {} = ({}, {}, {}), {} = ({}, {}, {}).", i,
                                     n[0], x[0], y[0], z[0],
                                     n[1], x[1], y[1], z[1],
                                     n[2], x[2], y[2], z[2],
                                     n[3], x[3], y[3], z[3]));
            }

        }
        TotalSize = TotalSize / 6.0; // scale tet determinant to volume
    }
    else {
        RUNTIME_ERROR("Determinants have already been calculated.");
    }
}

void Mesh::CalculateSurfaceNormals(double *p, Mesh* tets) {
    // CALCULATES THE SURFACE NORMAL FOR A TRIANGULAR MESH
    if ((Dimension != 2) ||
            (getnElements() <1) ||
            (getnNodes() <1)) {
        RUNTIME_ERROR("Can not calculate surface normals.");
    }

    Log::info("Calculating {} surface normals and 2D determinants.", getnElements());

    if (SurfaceNormal != NULL) free(SurfaceNormal);
    if (Determinant != NULL) free(Determinant);

    SurfaceNormal = (double*) malloc(3 * getnElements() * sizeof(double));
    Determinant = (double*) malloc( getnElements() * sizeof(double) );

    if ((SurfaceNormal==NULL) || (Determinant==NULL)) {
        RUNTIME_ERROR("Failed to allocate memory.");
    }

    // Loop over each element and calculate determinant and surface normal
    double Ax,Ay,Az, Bx,By,Bz;
    double cross[3];
    TotalSize = 0;

    for ( idx i = 0 ; i < getnElements() ; i ++) {
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

        if (det < 0 ) { // if det is negative make it not...
            Log::warn("Negative determinant for triangle {}, negating it and continuing.", i);
            det = -1.0 * det;
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
                idx n[4] = { tets->getNode(t,0) , tets->getNode(t,1) , tets->getNode(t,2), tets->getNode(t,3) };
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
        else {
            RUNTIME_ERROR("Connected volumes or tetrahedra do not exist.");
        }
        setSurfaceNormal(i,cross);
    }
}

void Mesh::CalcLocCoords(const idx elem, double *p, double *coord, double *loc) {
    // calculates 4 local coordinates of coordinate coord in element elem
    // returned in loc
    idx* n = this->getPtrToElement( elem );

    double M[12]= {0,0,0,0,0,0,0,0,0,0,0,0};

    for (int c = 0 ; c < 4 ; c++) {
        for (int i = 0 ; i < 4 ; i++) {
            M[i*3 + 0] = p[ 3*n[i] + 0]; // x
            M[i*3 + 1] = p[ 3*n[i] + 1]; // y
            M[i*3 + 2] = p[ 3*n[i] + 2]; // z
        }

        M[c*3 + 0 ] = coord[0];
        M[c*3 + 1 ] = coord[1];
        M[c*3 + 2 ] = coord[2];

        double det = fabs( this->Calculate4x4Determinant(M) ) ;
        loc[c] = det / (this->getDeterminant(elem)*1e18 );
    }
}

void Mesh::CalcElemBary(const idx elem, const double *p, double *bary) const {
#ifdef DEBUG
    assert(Dimension == 3); // currently only works for tetrahedra
#endif
    idx* n = Elem + elem*4; // shortcut

    const double *p1 , *p2, *p3, *p4;
    p1 = p + n[0]*3;
    p2 = p + n[1]*3;
    p3 = p + n[2]*3;
    p4 = p + n[3]*3;

    bary[0] = ( p1[0] + p2[0] + p3[0] + p4[0] ) / 4.0;
    bary[1] = ( p1[1] + p2[1] + p3[1] + p4[1] ) / 4.0;
    bary[2] = ( p1[2] + p2[2] + p3[2] + p4[2] ) / 4.0;

}

double Mesh::CalcBaryDistSqr(const double *p, const idx elem, const double *coord) const {
    // returns the squared distance between coordinate coord[3] and the barycentre of
    // element elem. p is pointer to all node coordinates.

    double bary[3] = {0,0,0};
    CalcElemBary( elem , p, bary);

    return	( (bary[0]-coord[0])*(bary[0]-coord[0]) +
              (bary[1]-coord[1])*(bary[1]-coord[1]) +
              (bary[2]-coord[2])*(bary[2]-coord[2]) );
}

void Mesh::listNodesOfMaterial(std::vector<idx> &nodes, const idx mat) const {
    // MAKES A LIST OF NODES BLEONGING TO ELEMENTS OF MATERIAL mat
    nodes.clear();
    for (idx i = 0 ; i < getnElements() ; i++){
        const idx curMat = this->getMaterialNumber(i);

        // THERE MAY BE A BUG HERE IN SOME CASES
        const idx bitRes = curMat & mat;
        if ( bitRes == mat ){
            for (idx n = 0 ; n < this->getnNodes() ; n++){
                const idx nodeNum= this->getNode(i,n); // nth node in ith element
                nodes.push_back( nodeNum );
            }
        }
    }

    // REMOVE REPEATED ENTRIES
    sort(nodes.begin() , nodes.end() );
    std::vector< unsigned int> :: iterator u;
    u = unique( nodes.begin() , nodes.end() );
    nodes.erase( u , nodes.end() );
}

void Mesh::listFixLCSurfaces(std::vector<idx> &nodes, const idx FixLCNumber) const {

#ifdef DEBUG
    assert((FixLCNumber <=9) &&(FixLCNumber >0)); // valid FixLC surfaces are in range 1-9
#endif

    nodes.clear();
    for(idx i = 0; i < getnElements(); i++){
        const idx curMat = this->getMaterialNumber(i);
        const idx curFixLCNum = curMat / MAT_FIXLC1;
        if (curFixLCNum == FixLCNumber){ // if matching FixLC number
            for (idx n = 0 ; n < this->getnNodes() ; n++){
                const idx nodeNum = this->getNode(i,n);
                nodes.push_back(nodeNum);
            }

        }
    }
}

void Mesh::CopySurfaceNormal(idx i, double* norm) {
#ifdef DEBUG
    assert(i < nElements);
    assert(SurfaceNormal);
#endif

    norm[0] = SurfaceNormal[i*3+0];
    norm[1] = SurfaceNormal[i*3+1];
    norm[2] = SurfaceNormal[i*3+2];
}

void Mesh::removeElements(std::set <idx>& index) {
    // REMOVES ELEMENTS IN INDEX.
    if (index.empty()) {
        return ;
    }

    idx num_new_elements = getnElements() - (idx) index.size();

    // make sure index does not contain elements that do not exist
    idx maxindex = * max_element( index.begin() , index.end() );

    if (maxindex >= getnElements()) {
        throw std::runtime_error(format("Element {} does not exist. Number of elements = {}.", maxindex, getnElements()));
    }

    idx* tnew 		= (idx*) malloc( num_new_elements * getnNodes() * sizeof(idx) ); // new elements
    idx* mnew 		= (idx*) malloc( num_new_elements * sizeof(idx) );		// new materials
    double* dnew	= (double*) malloc( num_new_elements * sizeof(double) ); 	// new determinants

    if ((tnew == nullptr) || (mnew == nullptr) || (dnew == nullptr)) {
        throw std::runtime_error(format("Allocation failed in {}, {}.", __FILE__, __func__));
    }

    double* SurfaceNormal_new 	= NULL;  // only needed for triangle meshes
    idx*    ConnectedVolume_new	= NULL;

    if (getDimension() == 2) { // if a triangle mesh
        SurfaceNormal_new   = (double*) malloc( num_new_elements * 3 * sizeof(double) );
        ConnectedVolume_new = (idx*) malloc( num_new_elements * sizeof(idx) );
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
            ++itr;
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

void Mesh::addElements(vector<idx> &m_new, vector<idx> &mat_new) {

    // CHECK THAT EQUAL SIZE OF MATERIAL NUMBERS AND ELEMENT NODES
    idx n_p_e = (unsigned int) this->getnNodes();
    idx num_new = m_new.size() / n_p_e;

    if ( num_new != mat_new.size() ) {
        throw std::runtime_error(format("Number of elements {} does not match number of materials {} in {}, {}",
                                        num_new, mat_new.size(), __FILE__, __func__ ));
    }

    // ALLOCATE MEMORY FOR NEW ELEMENTS AND MAT NUMS
    idx* nelem = (idx*) malloc( (getnElements() + m_new.size() ) * sizeof (idx) * n_p_e );
    idx* nmat  = (idx*) malloc( (getnElements() + mat_new.size() ) * sizeof(idx) );

    // COPY OLD VALUES
    memcpy( nelem, Elem, n_p_e*getnElements() * sizeof(idx) );
    memcpy( nmat , Mat , getnElements() * sizeof(idx) );

    // ADD NEW VALUES
    for (idx i = 0 ; i < (idx) mat_new.size()  ; i++) {
        nmat [ i+ getnElements()] = mat_new[i];
        for (idx j = 0 ; j < n_p_e ; j++) {
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

void Mesh::CopyMesh(Mesh* rhs) {
    setDimension(rhs->getDimension() );
    setnElements(rhs->getnElements() );		// set number of elements
    setnNodes(	rhs->getnNodes() );             // number of nodes per element
    TotalSize = rhs->TotalSize;

    setMaxNodeNumber( rhs->getMaxNodeNumber() );

    AllocateMemory();				// allocate memory for arrays
    setAllNodes(rhs->getPtrToElement(0) );	// copy node numbers
    setAllMaterials(rhs->getPtrToMaterialNumber(0) );// copy material numbers

    // COPY DETERMINANTS
    if (rhs->getPtrToDeterminant(0) != 0) {
        if (Determinant!=NULL) free(Determinant);
        Determinant = (double*) malloc( getnElements() * sizeof(double) );
        memcpy( Determinant, rhs->Determinant, getnElements() * sizeof(double) );
    }

    this->Dimension = rhs->Dimension;

    if (Dimension == 2) // if triangle mesh
    {
        // COPY TRIANGLE/TET CONNECTIONS
        if(rhs->getPtrToConnectedVolume(0) != NULL)
        {
            if (ConnectedVolume != NULL) free(ConnectedVolume);
            ConnectedVolume = (idx*) malloc( getnElements() * sizeof(idx) );

            memcpy(ConnectedVolume, rhs->ConnectedVolume , getnElements() * sizeof(idx) );

        }

        // COPY SURFACE NORMALS
        if (rhs->getPtrToSurfaceNormal(0) != NULL) {
            if (SurfaceNormal != NULL) free(SurfaceNormal);
            SurfaceNormal = (double*) malloc(3 * getnElements() * sizeof(double));

            memcpy( SurfaceNormal, rhs->SurfaceNormal, 3*getnElements() * sizeof(double) );
        }
    } // end if triangle mesh
}

void Mesh::ScaleDeterminants( const double& s) {
#ifdef DEBUG
    assert(Determinant != nullptr);
#endif

    for (idx i = 0; i < getnElements(); i ++) {
        Determinant[i] = Determinant[i] * s;
    }

    TotalSize = TotalSize * s;
}

void Mesh::gen_p_to_elem(vector<set<idx>> &p_to_elem) const {
    p_to_elem.clear();
    p_to_elem.reserve( this->getMaxNodeNumber() ); // pre-allocate space

    // initialise vector to allow [index] access later
    for ( idx i = 0 ; i < this->getMaxNodeNumber()+1 ; i++) {
        set <idx> empty;
        p_to_elem.push_back( empty );
    }

    for (idx i = 0 ; i <  this->getnElements() ; i++) // for every element
    {
        for (idx j = 0 ; j < this->getnNodes() ; j++) // for every node
        {
            idx n = this->getNode( i , j );
            p_to_elem[n].insert( i );
        }
    }
}
