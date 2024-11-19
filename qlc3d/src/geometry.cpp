#include <geometry.h>
#include <util/logging.h>
#include <util/exception.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>

const idx Geometry::NOT_AN_INDEX = std::numeric_limits<idx>::max();

Geometry::Geometry():
    regularGrid(nullptr) {
    npLC        = 0;
    t           = Mesh::tetMesh();
    e           = Mesh::triangleMesh();
    boundingBox = AABox();
    left_right_is_periodic = false;
    front_back_is_periodic = false;
    top_bottom_is_periodic = false;
}

Geometry::~Geometry() {
    delete regularGrid;
}

void Geometry::setTo(Geometry *geom) {
    this->ClearGeometry();
    npLC    = geom->getnpLC();                  // number of LC nodes
    boundingBox = geom->boundingBox;
    t->CopyMesh(geom->t.get());
    e->CopyMesh(geom->e.get());

    nodeNormals = geom->nodeNormals;
    coordinates_ = geom->getCoordinates().clone();

    left_right_is_periodic = geom->getleft_right_is_periodic();
    front_back_is_periodic = geom->getfront_back_is_periodic();
    top_bottom_is_periodic = geom->gettop_bottom_is_periodic();
    this->periNodes_.clear();
    periNodes_.insert(periNodes_.end(), geom->periNodes_.begin(), geom->periNodes_.end());
    if (this->regularGrid) delete regularGrid;
    if (geom->regularGrid)
        regularGrid = new RegularGrid(*geom->regularGrid);
}

void Geometry::ClearGeometry() {
    npLC = 0;
    nodeNormals.clear();
    boundingBox.setBoundsX(0, 0);
    boundingBox.setBoundsY(0, 0);
    boundingBox.setBoundsZ(0, 0);
    left_right_is_periodic = false;
    front_back_is_periodic = false;
    top_bottom_is_periodic = false;
    peri_equ_nodes.clear(); // this is never used anyways ??
    t->ClearMesh();
    e->ClearMesh();
    if (coordinates_ != nullptr) {
      coordinates_->clear();
    }
}

void Geometry::setCoordinates(const std::shared_ptr<Coordinates>& coordinates) {
  coordinates_ = coordinates;
  setnpLC(coordinates->size());

  boundingBox = coordinates->findBoundingBox();

  // reset node normals
  nodeNormals.clear();
  nodeNormals.resize(getnp());
}

void Geometry::setTetrahedra(const std::shared_ptr<Mesh> &tetrahedra) {
  t = tetrahedra;
}

void Geometry::setTriangles(const std::shared_ptr<Mesh> &triangles) {
  e = triangles;
}

unsigned int Geometry::getnp() const {
  return coordinates_ == nullptr ? 0 : coordinates_->size();
}

void Geometry::addCoordinates(const vector<double> &coords) {
  coordinates_->append(coords);
  setnpLC(coordinates_->size());
}

void Geometry::setMeshData(const std::shared_ptr<Coordinates> &coordinates,
                 std::vector<unsigned int> &&tetNodes, std::vector<unsigned int> &&tetMaterials,
                 std::vector<unsigned int> &&triNodes, std::vector<unsigned int> &&triMaterials) {
  setCoordinates(coordinates);
  t->setElementData(std::move(tetNodes), std::move(tetMaterials));
  e->setElementData(std::move(triNodes), std::move(triMaterials));
  ReorderDielectricNodes();
  e->setConnectedVolume(t.get()); // neighbour index tri -> tet
  t->calculateDeterminants3D(getCoordinates()); // calculate tetrahedral determinants
  t->ScaleDeterminants(1e-18); // scale to microns

  e->calculateSurfaceNormals(getCoordinates(), t.get()); // calculate triangle determinants and surface normal vectors
  e->ScaleDeterminants(1e-12); // scale to microns

  calculateNodeNormals();
  initialisePeriodicity(); // also makes periodic node indexes
}


unsigned int findMaxNodeNumber(const Mesh &m) {
  unsigned int maxNodeNumber = m.getNode(0, 0);
  for (unsigned int i = 0; i < m.getnElements(); i++) {
    for (unsigned int j = 0; j < m.getnNodes(); j++) {
      maxNodeNumber = std::max(maxNodeNumber, m.getNode(i, j));
    }
  }
  return maxNodeNumber;
}

void Geometry::updateMaxNodeNumbers() {
    unsigned int maxNodeNumberTets = findMaxNodeNumber(getTetrahedra());

    // for tets this should equal total number of coordinates
    if (maxNodeNumberTets != coordinates_->size() - 1) {
      RUNTIME_ERROR(fmt::format("Tet mesh max node number does not match number of coordinates: {} vs {}", maxNodeNumberTets, coordinates_->size()));
    }
}

void Geometry::setnpLC(int n) {
    npLC = n;
}

void Geometry::calculateNodeNormals() {
    assert(e != nullptr);
    assert(e->getnElements() > 0);

    Log::info("Calculating surface normals for {} alignment layer triangles", e->getnElements());

    nodeNormals.clear();
    nodeNormals.resize(3 * getnp());

    for (idx i = 0 ; i < e->getnElements() ; i ++) { // add neighbouring surface normals
        int m = e->getMaterialNumber(i);
        if (MATNUM_TO_FIXLC_NUMBER(m)) { // IF FIXLC SURFACE
            Vec3 elemNormal = e->getSurfaceNormal(i);
            for (idx j = 0; j < e->getnNodes() ; j++) {
              unsigned int nodeIdx = e->getNode(i, j);
              nodeNormals[nodeIdx] += elemNormal;
            }
        }
    }

    for (size_t i = 0; i < getnp(); i++) { // normalise length
      nodeNormals[i].normalize();
    }
}

void Geometry::ReorderDielectricNodes() {
    assert(t != nullptr);
    assert(t->getnElements() > 0);
    assert(coordinates_ != nullptr);
    assert(getnp() > 0);

    // check if dielectric elements exist
    bool DE_exist = false;
    const Mesh &tets = getTetrahedra();
    for (idx i = 0 ; i < t->getnElements() ; i++) {
        if (tets.getMaterialNumber(i) >= MAT_DIELECTRIC1) { // if material nuber > LC number
            DE_exist = true;
            break;
        }
    }
    setnpLC(getnp());
    if (!DE_exist) { // if no dielectric materials -> no need to reorder = exit
        this->updateMaxNodeNumbers();
        return;
    }
    //
    //GENERATE LIST OF MATERIAL NUMBERS FOR NODES. IN CASE OF DUAL VALUE LC PRECEDES
    //
    // mark all LC nodes as 1 and others as 0
    idx *lcde = (idx *) malloc(getnp() * sizeof(idx));
    memset(lcde , 0 , getnp() * sizeof(int));  // start with everything 0
    for (idx i = 0 ; i < t->getnElements() ; i ++) {
        if (t->getMaterialNumber(i) == MAT_DOMAIN1) {
            for (idx j = 0 ; j < t->getnNodes() ; j ++)  // loop over all nodes
                lcde[t->getNode(i , j)] = 1; // LC --> 1
        }
    }
    npLC = 0;
    vector <int> v_mat_index;
    for (size_t i = 0; i < getnp() ; i++) // first add all nodes marked as LC
        if (lcde[i] == 1) {
            v_mat_index.push_back(i);
            npLC++;
        }
    for (size_t i = 0; i < getnp() ; i++) // then add all non-LC nodes
        if (lcde[i] == 0)
            v_mat_index.push_back(i);

    free(lcde);
    //make inverse map
    vector <int> v_invmap;
    v_invmap.resize(v_mat_index.size() , -1);
    for (size_t i = 0 ; i < getnp() ; i ++)
        v_invmap[v_mat_index[i]] = i; //v_mat_index[i];//= i;

    // Reordered coordinates
    std::vector<Vec3> vecs;
    vecs.reserve(coordinates_->size());

    for (size_t i = 0; i < coordinates_->size(); i++) {
        vecs.push_back(coordinates_->getPoint(v_mat_index[i]));
    }
    coordinates_ = std::make_shared<Coordinates>(std::move(vecs));

    //REORDER TETRAHEDRA
    idx *newt = (idx *)malloc(t->getnElements() * t->getnNodes() * sizeof(idx));
    for (idx i = 0 ; i < t->getnElements() ; i++) {
        for (idx j = 0 ; j < t->getnNodes() ; j ++)  // loop over all nodes of element i
            newt[i * t->getnNodes() + j ] = v_invmap[ t->getNode(i , j) ];
    }
    t->setAllNodes(newt);  // copy node numbers
    free(newt);

    //REORDER TRIANGLES
    idx *newe = (idx *)malloc(e->getnElements() * e->getnNodes() * sizeof(idx));
    for (idx i = 0 ; i < e->getnElements() ; i++) {
        for (idx j = 0 ; j < e->getnNodes() ; j++)
            newe[i * e->getnNodes() + j ] = v_invmap[ e->getNode(i, j) ];
    }
    e->setAllNodes(newe);  // copy node numbers
    free(newe);
    this->updateMaxNodeNumbers();
}

void Geometry::setFacePeriNodes(list<size_t> &face0,
                                list<size_t> &face1,
                                const int &norm) {
  // SETS INDEX VALUES IN periNodes_ FOR FACES
  // norm is face normal, need coordinates to vectors parallel to it
  // norm = 0 -> COMPARE Y,Z
  // norm = 1 -> COMPARE X,Z
  // notm = 2 -> COMPARE X,Y
  int ind1 = (norm + 1) % 3;   // ind is a pre-calculated offeset to coordinate comparison in p.
  int ind2 = (norm + 2) % 3;   // Selects comparison of x,y or z coordinate. e.g. norm = 0 -> ind1 = 1, ind2 = 2
  double eps = 1e-5; // accuracy of coordinate comparison
  // SEARCH FOR NODE EQUIVALENCIES BY COMPARING COORDINATES
  // THAT ARE PERPENDICULAR TO FACE NORMAL
  list <size_t>:: iterator F0;    // FACE 0 NODES
  list <size_t>:: iterator F1;    // FACE 1 NODES
  int fc, bc; // debug counters
  for (F0 = face0.begin(), fc = 0; F0 != face0.end() ; ++F0, ++fc) { // LOOP OVER FACE 0
    bool found = false;

    int ind_n  = 0;     // index to neares (debug)
    double dist = 1000000;

    double f1 = coordinates_->getPoint(*F0).getDimension(ind1); // p[3 * (*F0) + ind1 ]; // coordinates of node F2 in face0
    double f2 = coordinates_->getPoint(*F0).getDimension(ind2); // p[3 * (*F0) + ind2 ];
    for (F1 = face1.begin() , bc = 0; F1 != face1.end() ; ++F1, ++bc) {
      double fa = coordinates_->getPoint(*F1).getDimension(ind1); // p[3 * (*F1) + ind1]; // coordinates of node F1 in face 1
      double fb = coordinates_->getPoint(*F1).getDimension(ind2); //p[3 * (*F1) + ind2];
      //compare coordinates
      double dist1 = fabs(f1 - fa); // distances in plane
      double dist2 = fabs(f2 - fb);
      double tdist = dist1 * dist1 + dist2 * dist2;

      if (tdist < dist) { // debug info only, keep track of nearest found node
        dist = tdist; // nearest distance
        ind_n = *F1;  // index to nearest distance
      }

      if (tdist < eps * eps) { // compare squared distances
        this->periNodes_[*F1] = *F0;
        found = true;
        break;
      }
    }// end for B
    if (!found) {
      Vec3 c0 = coordinates_->getPoint(*F0);
      Vec3 c1 = coordinates_->getPoint(ind_n);

      Log::error("normal={}.", norm);
      Log::error("coordinate {} FACE0: ({}).", *F0, c0);
      Log::error("nearest match in FACE1: index {}, distance = {}, at ({}).", ind_n, dist, c1);
      RUNTIME_ERROR("Periodic boundaries do not match");
    }
  }
}

void Geometry::setEdgePeriNodes(list<size_t> &edge0,
                                list<size_t> &edge1,
                                const int &dim) {
  const double BIGNUM = 1e99;
  list<size_t> :: iterator e0;
  list<size_t> :: iterator e1;
  double eps = 1e-5;
  // LOOP OVER EACH NODE IN EDGE 0 AND FIND ITS EQUIVALENT IN OTHER EDGES
  for (e0 = edge0.begin() ; e0 != edge0.end() ; ++e0) {
    // outer corner loop - corn0
    bool found = false;
    double dist = 0;

    double p0 = coordinates_->getPoint(*e0).getDimension(dim); // dim coordinate of node c0

    int iNearest = -1;
    double minDist = BIGNUM;
    // COMPARISON WITH C1
    for (e1 = edge1.begin() ; e1 != edge1.end() ; ++e1) {
      // inner corner loop - corn1
      double p1 = coordinates_->getPoint(*e1).getDimension(dim);

      dist = fabs(p0 - p1);
      if (dist <= minDist) {
        minDist = dist;
        iNearest = *e1;
      }
      if (dist <= eps) { // if same coordinate -> match!
        periNodes_[*e1] = *e0;
        found = true;
        break;
      }
    }
    if (!found) {
      Vec3 node = coordinates_->getPoint(*e0);
      Vec3 nearest = coordinates_->getPoint(iNearest);
      Log::error("Edge node 1 not found.");
      Log::error("Edge node 0 index = {}, coordinates = ({}).", *e0, node);
      Log::error("Nearest node index = {}, coordinates = ({}).", iNearest, nearest);
      RUNTIME_ERROR("Edge nodde 1 not found.");
    }
  }// end loop over all nodes in edge0
}

void Geometry::makePeriEquNodes() {
    // GENERATES INDEX OF PERIODIC NODE EQUIVALENCIES
    // LIST OF INDEXES TO ALL PERIODIC NODES
    vector <unsigned int> nodes;
    e->listNodesOfMaterial(nodes , MAT_PERIODIC);
    size_t numPeri = nodes.size();
    if (!numPeri) {
        periNodes_.clear();
        return;
    }
    periNodes_.clear();
    periNodes_.reserve(getnp());
    for (size_t i = 0 ; i < getnp() ; i++)
        periNodes_.push_back(i);
    // NEED TO CONSIDER 3 DIFFERENT CASES, DEPENDING ON TYPE OF PERIODICITY OF MODELLING WINDOW
    const double eps = 1e-5; // accuracy for coordinate comparisons
    const double xmin = boundingBox.getXMin();
    const double xmax = boundingBox.getXMax();
    const double ymin = boundingBox.getYMin();
    const double ymax = boundingBox.getYMax();
    const double zmin = boundingBox.getZMin();
    const double zmax = boundingBox.getZMin();
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
            !gettop_bottom_is_periodic()) {
        //SEPARATE NODES INTO TWO LISTS FOR FRONT AND BACK SURFACES
        list <size_t> front;
        list <size_t> back;
        for (size_t i = 0 ; i < nodes.size() ; i ++) {
            unsigned int n = nodes[i];
            if (FRONT) { // check if node i is on front surface
                front.push_back(n) ;
            } else if (BACK) { // check if node i is on back surface
                back.push_back(n);
            } else { // ERROR
                throw std::runtime_error(fmt::format("Expected node {} on front/back surface in {}, {}.",
                                                n, __FILE__, __func__));
            }
        }//end for i
        /// MAKE SURE EQUAL NUMBER OF NODES HAVE BEEN FOUND ON BOTH SURFACES
        if (front.size() != back.size()) {
            throw std::runtime_error(fmt::format("Periodic front surface contains {} nodes and back surface {} in {}, {}",
                                            front.size(), back.size(), __FILE__, __func__ ));
        }
        // SEARCH FOR NODE EQUIVALENCIES BY COMPARING X AND Z COORDINATES
        // BACK NODES MAP TO FRONT NODES
        this->setFacePeriNodes(front, back, 1);
    }
    // CASE 2 FRONT-BACK AND LEFT-RIGHT ARE PERIODIC
    else if (getfront_back_is_periodic() &&
             getleft_right_is_periodic() &&
             !gettop_bottom_is_periodic()) {
        // separate nodes into 8 lists, 4 x corners left/right and front/back planes
        list <size_t> edge0; //x = 0, y = 0
        list <size_t> edge1; //x = 0, y = max
        list <size_t> edge2; //x = max, y = max
        list <size_t> edge3; //x = max, y = 0
        list <size_t> front; //x = 0
        list <size_t> back;  //x = max
        list <size_t> right; //y = max
        list <size_t> left;  //y = 0
        for (size_t i = 0 ; i < nodes.size() ; i++) { // loop over all nodes and insert to correct list
            int n = nodes[i];
            if (LEFT && FRONT) { // edge0
                edge0.push_back(n);
            } else if (LEFT && BACK) { // edge1
                edge1.push_back(n);
            } else if (RIGHT  &&  BACK) { // edge2
                edge2.push_back(n);
            } else if (RIGHT && FRONT) { // edge3
                edge3.push_back(n);
            } else if (FRONT) { // front surface
                front.push_back(n);
            } else if (BACK) { // back surface
                back.push_back(n);
            } else if (LEFT) { // left surface
                left.push_back(n);
            } else if (RIGHT) { // right surface
                right.push_back(n);
            }
        }
        if ((edge0.size() != edge1.size()) ||
                (edge0.size() != edge2.size())  ||
                (edge0.size() != edge3.size())  ||
                (left.size() != right.size()) ||
                (front.size() != back.size())) {
            Log::error("Periodic node counts don't match.");
            Log::error("Corners 0, 1, 2, 3 contain [{}, {}, {}, {}] nodes.",
                       edge0.size(), edge1.size(), edge2.size(), edge3.size());
            Log::error("Front/back, left/right surfaces contain {}/{}, {}/{} nodes.",
                       front.size(), back.size(), left.size(), right.size());
            throw std::runtime_error(fmt::format("Periodic node counts don't match in {}, {}", __FILE__, __func__ ));
        }
        setEdgePeriNodes(edge0, edge1, 2);
        setEdgePeriNodes(edge0, edge2, 2);
        setEdgePeriNodes(edge0, edge3, 2);
        setFacePeriNodes(left  , right, 0); // COMPARE Y, Z
        setFacePeriNodes(front , back , 1);  // COMPARE X, Z
    }// END CASE 2
    // CASE 3 FRONT-BACK, LEFT-RIGHT AND TOP-BOTTOM ARE PERIODIC
    else if (getfront_back_is_periodic() &&
             getleft_right_is_periodic() &&
             gettop_bottom_is_periodic()) {
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
        int corner_nodes[8] = { -1, -1, -1, -1, -1, -1, -1, -1};
        for (size_t i = 0 ; i < nodes.size() ; i++) {
            int n = nodes[i];
            // CORNER NODES TAKE PRECEDENCE OVER OTHER NODES
            // FRONT LEFT BOTTOM
            if (FRONT && LEFT && BOTTOM)
                corner_nodes[0] = n;
            // FRONT RIGHT BOTTOM
            else if (FRONT && RIGHT && BOTTOM)
                corner_nodes[1] = n;
            else if (FRONT && LEFT && TOP)    /// FRONT LEFT TOP
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
            else if (LEFT && FRONT) { // edge0
                edge0.push_back(n);
            } else if (LEFT  && BACK) { // edge1
                edge1.push_back(n);
            } else if (RIGHT  && BACK) { // edge2
                edge2.push_back(n);
            } else if (RIGHT && FRONT) { // edge3
                edge3.push_back(n);
            }
            // 4 x Horizontal along X
            else if (FRONT && BOTTOM) {  // ymin and zmin
                edgea.push_back(n);
            } else if (BACK && BOTTOM) { // ymax and zmin
                edgeb.push_back(n);
            } else if (BACK && TOP) { // ymax and zmax
                edgec.push_back(n);
            } else if (FRONT && TOP) { // ymin and zmax
                edged.push_back(n);
            }
            // 4 x Horizontal along Y
            else if (LEFT && BOTTOM) {
                edgeA.push_back(n);   // xmin and zmin
            } else if (RIGHT && BOTTOM) {
                edgeB.push_back(n);   // xmax and zmin
            } else if (RIGHT && TOP) {
                edgeC.push_back(n);   // xmax and zmax
            } else if (LEFT && TOP) {
                edgeD.push_back(n);   // xmin and zmax
            }
            // FRONT/BACK, LEFT/RIGHT, TOP/BOTTOM FACES
            else if (FRONT) { // front surface
                front.push_back(n);
            } else if (BACK) { // back surface
                back.push_back(n);
            } else if (LEFT) { // left surface
                left.push_back(n);
            } else if (RIGHT) { // right surface
                right.push_back(n);
            } else if (BOTTOM) { // bottom surface
                bottom.push_back(n);
            } else if (TOP) { // top surface
                top.push_back(n);
            } else {
                throw std::runtime_error(fmt::format("Periodic node {} is not on an external surface in {}, {}.",
                                                n, __FILE__, __func__));
            }
        }// end for i, loop over all nodes
        // CHECK THAT OPPOSITE FACES HAVE EQUAL NUMBER OF NODES
        {
            // start dummy scope
            int mincorner = *min_element(corner_nodes, corner_nodes + 8);
            if (mincorner < 0) {
                throw std::runtime_error(fmt::format("Periodic corner nodes not found in {}, {}. "
                                                     "Indices are [{}. {}, {}, {}, {}, {}, {}, {}].",
                                                __FILE__, __func__, corner_nodes[0], corner_nodes[1], corner_nodes[2],
                                                corner_nodes[3], corner_nodes[4], corner_nodes[5], corner_nodes[6], corner_nodes[7]));
            }
            if (top.size() != bottom.size()) {
                throw std::runtime_error(fmt::format("Periodic top/bottom surfaces contain {}/{} nodes in {}, {}.",
                                                top.size(), bottom.size(), __FILE__, __func__));
            }
            if (left.size() != right.size()) {
                throw std::runtime_error(fmt::format("Periodic left/right surfaces contain {}/{} nodes in {}, {}.",
                                                left.size(), right.size(), __FILE__, __func__));
            }
            if (front.size() != back.size()) {
                throw std::runtime_error(fmt::format("Periodic front/back surfaces contain {}/{} nodes in {}, {}.",
                                                front.size(), back.size(), __FILE__, __func__));

            }
            // CHECK ALL CORNERS HAVE CORRECT NUMBER OF NODES
            size_t s0, s1, s2, s3;
            s0 = edge0.size(); s1 = edge1.size(); s2 = edge2.size() ; s3 = edge3.size();
            if ((s1 != s0) || (s2 != s0) || (s3 != s0)) {
                throw std::runtime_error(fmt::format("Periodic edges along z contain {}, {}, {}, {} nodes in {}, {}",
                                                s0, s1, s2, s3, __FILE__, __func__));
            }
            s0 = edgea.size(); s1 = edgeb.size(); s2 = edgec.size(); s3 = edged.size();
            if ((s1 != s0) || (s2 != s0) || (s3 != s0)) {
                throw std::runtime_error(fmt::format("Periodic edges along x contain {}, {}, {}, {} nodes in {}, {}",
                                                s0, s1, s2, s3, __FILE__, __func__));
            }
            s0 = edgeA.size(); s1 = edgeB.size(); s2 = edgeC.size(); s3 = edgeD.size();
            if ((s1 != s0) || (s2 != s0) || (s3 != s0)) {
                throw std::runtime_error(fmt::format("Periodic edges along y contain {}, {}, {}, {} nodes in {}, {}",
                                                s0, s1, s2, s3, __FILE__, __func__));
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
        setEdgePeriNodes(edge0, edge1 , 2);   // VERTICAL EDGES
        setEdgePeriNodes(edge0, edge2 , 2);
        setEdgePeriNodes(edge0, edge3 , 2);
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

void Geometry::initialisePeriodicity() {
    Log::info("Initialising peridioc surfaces");
    // CHECKS FOR TYPES OF PERIODICITY PRESENT IN CURRENT STRUCTURE.
    // POSSIBLE PERIODIC SURFACES ARE:
    //      LEFT/RIGHT
    //      FRONT/BACK
    //      TOP/BOTTOM
  const double EPS = 1e-7;
  for (idx i = 0 ; i < e->getnElements() ; i++) {
    if (e->getMaterialNumber(i) == MAT_PERIODIC) { // if surface is periodic
      // check to see which side surface is on by looking at the surface normal
      Vec3 snorm = e->getSurfaceNormal(i);

      if (fabs(fabs(snorm.x()) - 1.0) < EPS) {
        left_right_is_periodic = true;
      }
      else if (fabs(fabs(snorm.y()) - 1.0) < EPS) {
        front_back_is_periodic = true;
      }
      else if (fabs(fabs(snorm.z()) - 1.0) < EPS) {
        top_bottom_is_periodic = true;
      }
      else {
        RUNTIME_ERROR(fmt::format("Periodic surface element {} has invalid normal ({})", i, snorm));
      }
      // IF ALL SURFACES HAVE ALREADY BEEN IDENTIFIED AS PERIODIC
      // NO NEED TO CHECK FURTHER TRIANGLES
      if ((getleft_right_is_periodic()) &&
          (getfront_back_is_periodic()) &&
          (gettop_bottom_is_periodic()))
        break;
    }
  }
  // IF ANY PERIODIC TRIANGLES WERE DETECTED
  if (getleft_right_is_periodic()
      ||  getfront_back_is_periodic()
      ||  gettop_bottom_is_periodic()) {
    makePeriEquNodes();
  }
}

void Geometry::makeRegularGrid(const size_t &nx,
                               const size_t &ny,
                               const size_t &nz) {
    // CREATES REUGLAR GRID OBJECT USED FOR INTERPOLATING VALUES FROM
    // TETRAHEDRAL MESH ONTO REGULARLY SPACED GRID
    if (nx == 0 || ny == 0 || nz == 0) {
        return;
    }
    Log::info("Generating regular grid lookup with grid size nx={}, ny={}, nz={}", nx, ny, nz);

    if (regularGrid) {
        delete regularGrid;
    }
    regularGrid = new RegularGrid();

    regularGrid->createFromTetMesh(nx, ny, nz, this);
}

bool Geometry::brute_force_search(unsigned int &ind,             // return index
                                  const Vec3 &crd,              // search coordinate
                                  const bool &terminateOnError, // terminate if not found>
                                  const bool &requireLCEelement // only LC element index may be returned
                                 ) {
    // BRUTE FORCE DEBUG SEARCH FOR TETRAHEFRON THAN CONTAINS POINT WITH COORDINATES IN coord
    // coord IS ASSUMED TO BE OF LENGTH 3, FOR x, y, z
    // loop over each element
    for (idx i = 0 ; i < t->getnElements() ; i++) {
      bool found = t->containsCoordinate(i, getCoordinates(), crd);
      if (found) {    // If coord is in tet i
            if (requireLCEelement) { // WANT LC
                if (t->getMaterialNumber(i) <= MAT_DOMAIN7) {  // IF LC
                    ind = i;
                    return true;
                }
            } else { // DON'T CARE WHETHER LC OR DE
                ind = i ;
                return true; // exit function when found
            }
        }
    }// end for loop over all elems
    // IF COORDINATE WAS NOT FOUND APPLICATION MAY NEED TO BE TERMINATED
    if (terminateOnError) {
      RUNTIME_ERROR(fmt::format("Brute force search could not find coordinate at ({})", crd));
    }
    // SIGNAL A NON-FOUND COORDINATE BY RETURNING FALSE
    return false;
}

size_t Geometry::recursive_neighbour_search(
        const Vec3 &targetPoint,
        const vector<set<unsigned int> > &p_to_t,
        const size_t &currentTetIndex,
        std::set<size_t> &tetHistory,
        const bool &requireLCElement) {
  // Recursive search for tetrahedron that contains targetPoint, move to tet whose centroid is closest to targetPoint
  bool found = t->containsCoordinate(currentTetIndex, getCoordinates(), targetPoint);
  if (found) {
    if (!requireLCElement) { // if not worried about whether LC or DE element
      return currentTetIndex;
    } else { // LC element is required
      if (this->t->getMaterialNumber(currentTetIndex) <= MAT_DOMAIN7) {  // if LC element
        return currentTetIndex;
      }
    }
  }

  tetHistory.insert(currentTetIndex);   // history should be used to avoid visiting same element multiple times (this isn't implemented yet)
  // RECURSION LIMIT. UGLY SOLUTION TO STOP OUT OF MEMORY (?) CRASH FOR LARGE GRIDS AND/OR MESHES
  idx RECURSION_LIMIT = 1000;
  if ((idx) tetHistory.size() > RECURSION_LIMIT) return NOT_AN_INDEX;

  // CREATE LIST OF NEIGHBOURING ELEMENTS
  std::vector <unsigned int> neighs;
  for (idx i = 0 ; i < t->getnNodes() ; i++) {
    int n = t->getNode(currentTetIndex, i);
    neighs.insert(neighs.end(), p_to_t[n].begin(), p_to_t[n].end());
  }
  sort(neighs.begin(), neighs.end());
  auto itr = unique(neighs.begin(),neighs.end());
  neighs.erase(itr, neighs.end());
  // REMOVE REFERENCE TO SELF TOO....
  vector <double> dists;
  for (size_t i = 0 ; i < neighs.size() ; i++) {
    Vec3 centroid = getTetrahedra().elementCentroid(neighs[i], getCoordinates());
    dists.push_back(centroid.distanceSquared(targetPoint));
  }
  // FIND INDEX TO NEIGHBOUR THAT IS NEAREST
  size_t indn = min_element(dists.begin(), dists.end()) - dists.begin();

  while (dists[indn] < DBL_MAX) {
    size_t indt = neighs[indn]; // actual element number
    size_t indexFound = NOT_AN_INDEX;
    if (tetHistory.find(indt) == tetHistory.end()) {
      indexFound =  recursive_neighbour_search(targetPoint,
                                               p_to_t,
                                               indt,
                                               tetHistory,
                                               requireLCElement);
      if (indexFound != NOT_AN_INDEX) {
        return indexFound;
      }
    }

    // Mark this element as visited
    dists[indn] = DBL_MAX;
    indn = min_element(dists.begin() , dists.end()) - dists.begin();
  }
  // IF ALL NEIGHBOURS FAIL
  return NOT_AN_INDEX;
}

void Geometry::genIndToTetsByCoords(vector<unsigned int> &returnIndex,   // return index
                                    const Coordinates &targetCoordinates, // coordinates, but not necessarily from this same geometry
                                    const bool &terminateOnError,// whther to terminate app. if coordinate not found. default = true;
                                    const bool &requireLCElement) { // only LC element can be re returned
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
    if (t == nullptr || getTetrahedra().getnElements() == 0) {
      RUNTIME_ERROR("No tetrahedra elements defined");
    }

  returnIndex.clear();
  unsigned int nt = (unsigned int) this->t->getnElements();
  returnIndex.assign(targetCoordinates.size(), nt);   // assing with a value that is one too much initially
  vector <set<unsigned int>> p_to_t;
  t->gen_p_to_elem(p_to_t);

  // find most central tet, this is used as starting tet in other searches
  Vec3 structureCentroid = boundingBox.centre();

  std::set<size_t> tetHistory;
  unsigned int midTet = recursive_neighbour_search(structureCentroid, p_to_t,0,tetHistory);
  if (midTet == NOT_AN_INDEX) { // starting index at centre of structure not found (probably a hole or concave mesh)
    midTet = 0;
  }

  // for each target coordinate, find which tet contains it, starting each search from midTet
  for (unsigned int n = 0; n < targetCoordinates.size(); n++) { // for each coord
    std::set<size_t> searchHistory; // keeps track of tested tets to avoid repeating work
    Vec3 targetPoint = targetCoordinates.getPoint(n);

    // nearest neighbour search
    size_t t0 = recursive_neighbour_search(targetPoint, p_to_t,midTet,searchHistory, requireLCElement);
    if (t0 != NOT_AN_INDEX) {
      returnIndex[n] = t0;
    } else { // recursive search failed, try brute force
      unsigned int tetIndex = 0;
      if (brute_force_search(tetIndex, targetPoint, terminateOnError, requireLCElement)) {
        returnIndex[n] = tetIndex;
      } else { // BRUTE FORCE FAIL IS ALLOWED (NODE MAY BE OUTSIDE MESH)
        returnIndex[n] = Geometry::NOT_AN_INDEX; // MARK INDEX AS INVALID
        Log::info("Could not find regular grid point {} at ({}) in volume mesh. Assuming it is outside the mesh and continuing.",
                  n, targetPoint);
      }
    }
  }// end for each target coordinate
}

double Geometry::getAbsXDist(int i , double x) {
  double cx = getCoordinates().getPoint(i).x();
  return fabs(cx - x);
}

double Geometry::getAbsYDist(int i , double y) {
  double cy = getCoordinates().getPoint(i).y();
  return fabs(cy - y);
}

double Geometry::getAbsZDist(int i, double z) {
  double cz = getCoordinates().getPoint(i).z();
  return fabs(cz - z);
}

const Coordinates& Geometry::getCoordinates() const {
  if (!coordinates_) {
    RUNTIME_ERROR("Coordinates not set");
  }
  return *coordinates_;
}

bool Geometry::getleft_right_is_periodic() const {
    return left_right_is_periodic;
}

bool Geometry::getfront_back_is_periodic() const {
    return front_back_is_periodic;
}

bool Geometry::gettop_bottom_is_periodic() const {
    return top_bottom_is_periodic;
}

const std::vector<Vec3>& Geometry::getNodeNormals() const {
  return nodeNormals;
}

Vec3 Geometry::getNodeNormal(unsigned int i) const {
  return nodeNormals[i];
}

