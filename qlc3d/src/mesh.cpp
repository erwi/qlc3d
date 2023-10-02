#include <qlc3d.h>
#include <mesh.h>
#include <material_numbers.h>
#include <fmt/format.h>
#include <util/logging.h>
#include <util/exception.h>
#include <set>
#include <geom/coordinates.h>
#include <geom/vec3.h>

#include <vector>

using fmt::format;

Mesh::Mesh(unsigned int dimension, unsigned int nodesPerElement) :
        Dimension{dimension}, nNodes{nodesPerElement}, nElements{0}, TotalSize{0} {}

Mesh::~Mesh() { }

idx Mesh::getMaterialNumber(const idx e) const {
  return materials[e];
}

idx Mesh::getFixLCNumber(const idx e) const {
    return (idx) getMaterialNumber(e) / MAT_FIXLC1;
}

idx Mesh::getDielectricNumber(const idx e) const {
    return (int) getMaterialNumber(e) / MAT_DIELECTRIC1;
}

double Mesh::getDeterminant(const idx i) const {
#ifdef DEBUG
  assert(i < nElements);
#endif
  return determinants[i];
}

idx Mesh::getConnectedVolume(const idx e) const {
  // RETURNS INDEX TO CONNECTED LC TETRAHEDRON, WHEN e IS INDEX TO A TRIANGLE
#ifdef DEBUG
  assert(e < connectedVolumes.size());
#endif
  return connectedVolumes[e];
}

void Mesh::setSurfaceNormal(idx i, const Vec3 &normal) {
#ifdef DEBUG
  assert(i < getnElements());
  assert(i < surfaceNormals.size());
#endif
  surfaceNormals[i] = normal;
}

Vec3 Mesh::getSurfaceNormal(unsigned int i) const {
#ifdef DEBUG
  assert(i < getnElements());
  assert(getDimension() == 2);
#endif
  return surfaceNormals[i];
}


void Mesh::setElementData(std::vector<unsigned int> &&nodes, std::vector<unsigned int> &&materials) {
  unsigned int numElements = nodes.size() / getnNodes();

  if (nodes.size() % getnNodes() != 0) {
    RUNTIME_ERROR(format("Number of nodes ({}) is not a multiple of number of nodes per element ({}).", nodes.size(), getnNodes()));
  }

  if (numElements != materials.size()) {
    RUNTIME_ERROR(format("Number of elements ({}) does not match number of materials ({}).", numElements, materials.size()));
  }

  this->nodes = nodes;
  this->materials = materials;
  this->nElements = numElements;
}

/**
 * copies all values from array nodes to this->Elem
 * @param nodes
 * @deprecated use setElementData instead
 */
void Mesh::setAllNodes(idx *nodes) {
#ifdef DEBUG
    assert(nElements > 0);
    assert(this->nodes.size() == nElements * nNodes); // only used in init. TODO fix this
#endif
    for (unsigned int i = 0; i < this->nodes.size(); i++) { // assumes already correct size TODO: don't assume
        this->nodes[i] = nodes[i];
    }
}

void Mesh::setnElements(idx nelem) {
    nElements = nelem;
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
    connectedVolumes.resize(getnElements(), NOT_AN_INDEX);

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
            if ( final.size() == 1) { // IF SINGLE TET FOUND
                if (vol->getMaterialNumber( *itr ) <= MAT_DOMAIN7 ) { // only connect to an LC element
                    connectedVolumes[i] = *itr;
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

                    connectedVolumes[i] = -2; // set to -2 if connected to dielectric
                }
            }
            else
                if ( final.size() == 2) {
                    // CHECK BOTH FIRST AND SECOND INDICES
                    if (vol->getMaterialNumber( *itr ) <= MAT_DOMAIN7 ) {// only connect to an LC element
                      connectedVolumes[i] = *itr;
                    }
                    else {
                        ++itr;
                        if (vol->getMaterialNumber( *itr ) <= MAT_DOMAIN7 ) {// only connect to an LC element
                          connectedVolumes[i] = *itr;
                        }
                    }
                }
                else {
                    throw std::runtime_error(format("Triangle {} is connected to {} tetrahedra, but expected 2, "
                                                    "in {}, {}", i, final.size(), __FILE__, __func__));
                }
        } // end exclude electrode only surfaces
    }
}

double det3D(Vec3 c[4]) {
  Vec3 v1 = c[1] - c[0];
  Vec3 v2 = c[2] - c[0];
  Vec3 v3 = c[3] - c[0];

  Vec3 cross = v2.cross(v3);
  return v1.dot(cross);
}

bool Mesh::containsCoordinate(idx elem, const Coordinates &coordinates, const Vec3 p) const {
  Vec3 c[4] = {coordinates.getPoint(getNode(elem, 0)), coordinates.getPoint(getNode(elem, 1)),
               coordinates.getPoint(getNode(elem, 2)), coordinates.getPoint(getNode(elem, 3))};

  const double eps = qlc3d_GLOBALS::GLOBAL_COORD_EPS;
  // first check if p is one of the corner nodes
  if (c[0].equals(p, eps) || c[1].equals(p, eps) || c[2].equals(p, eps) || c[3].equals(p, eps)) {
    return true;
  }

  // then check if p is an internal node
  // sum of 4 "sub-tet" determinants will equal totalDet if p is an internal node
  double totalDet = std::fabs(det3D(c));
  double sumSubDets = 0;
  for (int i = 0; i < 4; i++) {
    Vec3 cs[4] = {c[0], c[1], c[2], c[3]};
    cs[i] = p;
    sumSubDets += std::fabs(det3D(cs));
  }

  return std::fabs(totalDet - sumSubDets) <= eps;
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


void Mesh::calculateDeterminants3D(const Coordinates &coords) {
  if (Dimension != 3) {
    RUNTIME_ERROR("3D determinants can only be calculated for 3D meshes.");
  }
  if (getnElements() == 0) {
    RUNTIME_ERROR("Cannot calculate determinants for empty mesh.");
  }

  determinants.resize(getnElements(), 0.0);
  std::fill(determinants.begin(), determinants.end(), 0.0);

  TotalSize = 0.0; // reset total mesh volume
  Log::info("Calculating {} 3D determinants.", getnElements());

  for (idx i = 0 ; i < getnElements() ; i++) {
    unsigned int n[4] = {getNode(i,0), getNode(i,1), getNode(i,2), getNode(i,3)};
    Vec3 corners[4] = {coords.getPoint(n[0]), coords.getPoint(n[1]), coords.getPoint(n[2]), coords.getPoint(n[3])};
    determinants[i] = det3D(corners);

    if (determinants[i] < 0 ) { // should reorder nodes to ensure positive determinant
      Log::warn("Negative determinant for element {}, negating sign and continuing.", i);
      determinants[i] = -1.0 * determinants[i];
    }

    TotalSize += determinants[i] ; // add determinant of this element to total

    if (determinants[i] == 0) {
      RUNTIME_ERROR(format("Zero determinant in element {} with nodes: {} = ({}), {} = ({}), {} = ({}), {} = ({}).",
                           i, n[0], corners[0], n[1], corners[1], n[2], corners[2], n[3], corners[3]));
    }
  }
  TotalSize = TotalSize / 6.0; // scale tet determinant to volume
}

void Mesh::calculateSurfaceNormals(const Coordinates &coords, Mesh* tets) {
  // CALCULATES THE SURFACE NORMAL FOR A TRIANGULAR MESH
  if ((Dimension != 2) ||  (getnElements() < 1) || (getnNodes() < 1)) {
    RUNTIME_ERROR("Can not calculate surface normals.");
  }

  if (connectedVolumes.size() < getnElements()) {
    RUNTIME_ERROR("Can not calculate surface normals, connectedVolumes vector does not match number of elements.");
  }

  Log::info("Calculating {} surface normals and 2D determinants.", getnElements());

  surfaceNormals.resize(getnElements(), {0.0, 0.0, 0.0});
  std::fill(surfaceNormals.begin(), surfaceNormals.end(), Vec3(0.0, 0.0, 0.0));

  determinants.resize(getnElements(), 0.0);
  std::fill(determinants.begin(), determinants.end(), 0.0);

  // Loop over each element and calculate determinant and surface normal
  TotalSize = 0;
  for ( idx i = 0 ; i < getnElements() ; i ++) {
    unsigned int n[3] = {getNode(i,0), getNode(i,1), getNode(i,2)};

    Vec3 p0 = coords.getPoint(n[0]);
    Vec3 p1 = coords.getPoint(n[1]);
    Vec3 p2 = coords.getPoint(n[2]);

    Vec3 A = p1 - p0;
    Vec3 B = p2 - p0;

    Vec3 normalVec = A.cross(B);
    double det = normalVec.norm();
    if (det < 0 ) { // if det is negative make it not...
      Log::warn("Negative determinant for triangle {}, negating it and continuing.", i);
      det = -1.0 * det;
    }

    determinants[i] = det;
    TotalSize+= det / 2.0;
    normalVec.normalize();

    // CHECK THAT SURFACE NORMAL POINTS TOWARDS LC REGION. I.E. INWARDS
    // IF NOT REVERSE ITS ORIENTATION
    unsigned int t = getConnectedVolume(i); // index to neighbouring tet.

    if (t != NOT_AN_INDEX) { // index may only exist for LC tet connections.
      Vec3 triBary = (p0 + p1 + p2) / 3.;

      // calculate tet barycentre
      idx ntet[4] = {tets->getNode(t, 0), tets->getNode(t, 1), tets->getNode(t, 2), tets->getNode(t, 3)};

      Vec3 c0 = coords.getPoint(ntet[0]);
      Vec3 c1 = coords.getPoint(ntet[1]);
      Vec3 c2 = coords.getPoint(ntet[2]);
      Vec3 c3 = coords.getPoint(ntet[3]);

      Vec3 tetBary = (c0 + c1 + c2 + c3) / 4.0;

      // vector from tri barycentre to tet barycentre
      Vec3 v1 = tetBary - triBary;

      // dot product between surface normal and v1 determines orientation of surface normal
      double dot = normalVec.dot(v1);

      if (dot < 0.0) {
        // IF SURFACE NORMAL IS IN WRONG DIRECTION,
        // REVERSE NODE ORDER AND INVERT VECTOR DIRECTION
        idx indexA = i * getnNodes() + 1; // 2nd node of triangle
        idx indexB = i * getnNodes() + 2; // 3rd node of triangle
        iter_swap(nodes.begin() + indexA, nodes.begin() + indexB); // swap 2nd and 3rd nodes

        normalVec *= -1;
      }
    }

    setSurfaceNormal(i, normalVec);
  }
}

void Mesh::calcLocCoords(const idx elem, const Coordinates &coordinates, const Vec3 &targetPoint, double localCoordinates[4]) const {
  unsigned int elemNodes[4];
  loadNodes(elem, elemNodes);
  const double determinant = getDeterminant(elem) * 1e18;

  // node 0 - opposite to [1, 2, 3]
  Vec3 nodes[4] = {targetPoint,
                   coordinates.getPoint(elemNodes[1]),
                   coordinates.getPoint(elemNodes[2]),
                   coordinates.getPoint(elemNodes[3])};
  double det1 = std::abs(det3D(nodes));
  localCoordinates[0] = det1 / determinant;

  // node 1 - opposite to [0, 2, 3]
  nodes[1] = coordinates.getPoint(elemNodes[0]);
  nodes[2] = coordinates.getPoint(elemNodes[2]);
  nodes[3] = coordinates.getPoint(elemNodes[3]);
  double det2 = std::abs(det3D(nodes));
  localCoordinates[1] = det2 / determinant;

  // node 2 - opposite to [0, 1, 3]
  nodes[1] = coordinates.getPoint(elemNodes[0]);
  nodes[2] = coordinates.getPoint(elemNodes[1]);
  nodes[3] = coordinates.getPoint(elemNodes[3]);
  double det3 = std::abs(det3D(nodes));
  localCoordinates[2] = det3 / determinant;

  // node 3 - opposite to [0, 1, 2]
  nodes[1] = coordinates.getPoint(elemNodes[0]);
  nodes[2] = coordinates.getPoint(elemNodes[2]);
  nodes[3] = coordinates.getPoint(elemNodes[1]);
  double det4 = std::abs(det3D(nodes));
  localCoordinates[3] = det4 / determinant;

  // debug sanity check - sum of sub-tet determinants should equal total determinant
#ifndef NDEBUG
  double sum = det1 + det2 + det3 + det4;
  double error = std::fabs(sum - getDeterminant(elem) * 1e18);
  if (error > 1e-9) {
    RUNTIME_ERROR(format("Error in local coordinates calculation for element {}. Sum of sub-determinants ({}) does not equal total determinant ({}).",
                         elem, sum, getDeterminant(elem) * 1e18));
  }
#endif
}

void Mesh::loadNodes(idx elementIndex, idx *nodesOut) const {
  for (unsigned int i = 0; i < getnNodes(); i++) {
    nodesOut[i] = getNode(elementIndex, i);
  }
}

Vec3 Mesh::elementCentroid(unsigned int i, const Coordinates &coordinates) const {
  Vec3 centroid(0.0, 0.0, 0.0);
  for (unsigned int j = 0; j < getnNodes(); j++) {
    centroid += coordinates.getPoint(getNode(i, j));
  }
  centroid /= getnNodes();
  return centroid;
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

//
//void Mesh::CopySurfaceNormal(idx i, double* norm) const {
//#ifdef DEBUG
//    assert(i < nElements);
//    assert(SurfaceNormal);
//#endif
//
//    norm[0] = SurfaceNormal[i*3+0];
//    norm[1] = SurfaceNormal[i*3+1];
//    norm[2] = SurfaceNormal[i*3+2];
//}

void Mesh::removeElements(std::set<idx>& index) {
  if (index.empty()) {
    return;
  } else if (index.size() == getnElements()) {
    this->ClearMesh();
    return;
  }

  unsigned int maxIndex = *max_element(index.begin(), index.end());
  if (maxIndex >= getnElements()) {
    RUNTIME_ERROR(fmt::format("Element {} does not exist. Number of elements = {}.", maxIndex, getnElements()));
  }

  std::vector<idx> keepNodes;
  std::vector<idx> keepMaterials;
  std::vector<double> keepDeterminants;
  std::vector<Vec3> keepSurfaceNormals; // only for surface meshes
  std::vector<idx> keepConnectedVolumes; // only for surface meshes

  bool isSurfaceMesh = getDimension() == 2;

  unsigned int counter = 0;
  for (unsigned int i = 0; i < getnElements(); i++) {
    bool removeElement = index.count(i) > 0;

    if (removeElement) {
      continue; // skip this element
    } else {
      for (unsigned int j = 0; j < getnNodes(); j++) {
        keepNodes.push_back(getNode(i, j));
      }
      keepMaterials.push_back(getMaterialNumber(i));
      keepDeterminants.push_back(getDeterminant(i));

      if (isSurfaceMesh) {
        keepConnectedVolumes.push_back(getConnectedVolume(i));
        keepSurfaceNormals.push_back(getSurfaceNormal(i));
      }

      counter++;
    }
  }

  this->nodes = keepNodes;
  this->materials = keepMaterials;
  this->determinants = keepDeterminants;
  this->surfaceNormals = keepSurfaceNormals;
  this->connectedVolumes = keepConnectedVolumes;
}

void Mesh::ClearMesh(){
  nElements	= 0;
  nodes.clear();
  materials.clear();
  determinants.clear();
  surfaceNormals.clear();
  connectedVolumes.clear();
}

void Mesh::appendElements(const vector<idx> &nodeValues, const vector<idx> &materialValues) {
    // CHECK THAT EQUAL SIZE OF MATERIAL NUMBERS AND ELEMENT NODES
    idx nodesPerElem = (unsigned int) this->getnNodes();
    idx numNewElements = nodeValues.size() / nodesPerElem;

    if (numNewElements != materialValues.size() ) {
        throw std::runtime_error(format("Number of elements {} does not match number of materials {} in {}, {}",
                                        numNewElements, materialValues.size(), __FILE__, __func__ ));
    }

    for (auto n : nodeValues) {
      nodes.push_back(n);
    }

    for (auto m : materialValues) {
      materials.push_back(m);
    }

    this->nElements = materials.size();
}

void Mesh::CopyMesh(Mesh* rhs) {
    // require that number of dimensions and nodes per element match between this and rhs mesh
    if ( (getDimension() != rhs->getDimension() ) || (getnNodes() != rhs->getnNodes() ) ) {
      RUNTIME_ERROR(format("Dimension or nodes per element do not match in {}, {}."));
    }

    setnElements(rhs->getnElements() );		// set number of elements
    TotalSize = rhs->TotalSize;

    nodes = rhs->nodes;
    materials = rhs->materials;
    determinants = rhs->determinants;
    connectedVolumes = rhs->connectedVolumes;
    surfaceNormals = rhs->surfaceNormals;
}

void Mesh::ScaleDeterminants( const double& s) {
  for (idx i = 0; i < getnElements(); i ++) {
    determinants[i] = determinants[i] * s;
  }

  TotalSize = TotalSize * s;
}

void Mesh::gen_p_to_elem(vector<set<idx>> &p_to_elem) const {
    // find largest node number
    idx maxNodeNumber = 0;
    for (const auto &n : nodes) {
      maxNodeNumber = std::max(maxNodeNumber, n);
    }
    p_to_elem.clear();
    p_to_elem.reserve(maxNodeNumber + 1);

    // initialise vector to allow [index] access later
    for (idx i = 0; i < maxNodeNumber + 1; i++) {
        set <idx> empty;
        p_to_elem.push_back( empty );
    }

    for (idx i = 0; i <  this->getnElements(); i++) { // for every element
        for (idx j = 0; j < this->getnNodes(); j++) {// for every node
            idx n = this->getNode( i , j );
            p_to_elem[n].insert( i );
        }
    }
}
