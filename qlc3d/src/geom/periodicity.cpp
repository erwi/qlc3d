#include <geom/periodicity.h>
#include <geom/vec3.h>
#include <geom/coordinates.h>
#include <mesh.h>
#include <material_numbers.h>
#include <geom/aabox.h>
#include <util/logging.h>
#include <util/exception.h>

PeriodicityType::PeriodicityType(const Mesh &triangles) {

  if (!triangles.hasSurfaceNormals()) {
    RUNTIME_ERROR("Mesh does not have surface normals.");
  }

  left_right_is_periodic = false;
  front_back_is_periodic = false;
  top_bottom_is_periodic = false;

  const double EPS = 1e-7;
  for (unsigned int i = 0; i < triangles.getnElements(); i++) {
    unsigned int material = triangles.getMaterialNumber(i);
    if (material != MAT_PERIODIC) {
      continue;
    }

    auto normal = triangles.getSurfaceNormal(i);

    if (fabs(fabs(normal.x()) - 1.0) < EPS) {
      left_right_is_periodic = true;
    }
    else if (fabs(fabs(normal.y()) - 1.0) < EPS) {
      front_back_is_periodic = true;
    }
    else if (fabs(fabs(normal.z()) - 1.0) < EPS) {
      top_bottom_is_periodic = true;
    }
    else {
      RUNTIME_ERROR(fmt::format("Periodic surface element {} has invalid normal ({})", i, normal));
    }
    // IF ALL SURFACES HAVE ALREADY BEEN IDENTIFIED AS PERIODIC
    // NO NEED TO CHECK FURTHER TRIANGLES
    if (left_right_is_periodic && front_back_is_periodic && top_bottom_is_periodic) {
      break;
    }
  }
}


set<unsigned int> findPeriodicNodes(const Mesh &tris) {

  set<unsigned int> leftFace, rightFace, frontFace, backFace, topFace, bottomFace;


  auto equals = [](double a, double b) { return fabs(a - b) < 1e-7; };
  auto isLeftFace = [equals](const Vec3 &normal) { return equals(normal.x(), 1.); };
  auto isRightFace = [equals](const Vec3 &normal) { return equals(normal.x(), -1.); };
  auto isFrontFace = [equals](const Vec3 &normal) { return equals(normal.y(), 1.); };
  auto isBackFace = [equals](const Vec3 &normal) { return equals(normal.y(), -1.); };
  auto isTopFace = [equals](const Vec3 &normal) { return equals(normal.z(), 1.); };
  auto isBottomFace = [equals](const Vec3 &normal) { return equals(normal.z(), -1.); };

  const Vec3 rightNormal(-1, 0, 0);
  const Vec3 frontNormal(0, 1, 0);
  const Vec3 backNormal(0, -1, 0);
  const Vec3 topNormal(0, 0, 1);
  const Vec3 bottomNormal(0, 0, -1);

  for (unsigned int i = 0; i < tris.getnElements(); i++) {
    unsigned int material = tris.getMaterialNumber(i);
    if (material != MAT_PERIODIC) {
      continue;
    }

    auto normal = tris.getSurfaceNormal(i);

    if (isLeftFace(normal)) {
      leftFace.insert(i);
    } else if (isRightFace(normal)) {
      rightFace.insert(i);
    } else if (isFrontFace(normal)) {
      frontFace.insert(i);
    } else if (isBackFace(normal)) {
      backFace.insert(i);
    } else if (isTopFace(normal)) {
      topFace.insert(i);
    } else if (isBottomFace(normal)) {
      bottomFace.insert(i);
    }
  }
  return leftFace;
}


PeriodicNodesMapping::PeriodicNodesMapping(const Mesh &tris,
                                           const Coordinates &coords) : periodicityType(tris){

  initialisePeriodicNodes(tris, coords);

  if (!periodicityType.isAnyPeriodic()) {
    Log::warn("No periodic surfaces detected. Periodic nodes mapping will not be created.");
    return;
  }

  std::set<unsigned int> periodicNodes = findPeriodicNodes(tris);
}

void PeriodicNodesMapping::initialisePeriodicNodes(const Mesh &e, const Coordinates &coords) {
  Log::info("Initialising periodic surfaces");

  if (periodicityType.isAnyPeriodic()) {
    Log::info("Finding periodic nodes for periodic surfaces on left-right={}, front-back={}, top-bottom={}",
              periodicityType.isLeftRightPeriodic(), periodicityType.isFrontBackPeriodic(), periodicityType.isTopBottomPeriodic());
    makePeriEquNodes(periodicityType, e, coords);
  }
}

unsigned int PeriodicNodesMapping::getPeriodicNode(unsigned int node) const {
  if (periodicityType.isAnyPeriodic()) {
    return periNodes_[node];
  } else {
    return node;
  }
}

void PeriodicNodesMapping::setFacePeriNodes(std::list<unsigned int> &face0,
                                std::list<unsigned int> &face1,
                                const int &norm,
                                const Coordinates &coordinates) {
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
  //list <size_t>:: iterator F0;    // FACE 0 NODES
  //list <size_t>:: iterator F1;    // FACE 1 NODES
  unsigned int fc, bc = 0; // debug counters

  const AABox boundingBox = coordinates.findBoundingBox();
  for (auto F0 = face0.begin(); F0 != face0.end() ; ++F0, ++fc) { // LOOP OVER FACE 0
    bool found = false;

    int ind_n  = 0;     // index to neares (debug)
    double dist = 1000000;

    double f1 = coordinates.getPoint(*F0).getDimension(ind1); // p[3 * (*F0) + ind1 ]; // coordinates of node F2 in face0
    double f2 = coordinates.getPoint(*F0).getDimension(ind2); // p[3 * (*F0) + ind2 ];
    for (auto F1 = face1.begin(); F1 != face1.end() ; ++F1, ++bc) {
      double fa = coordinates.getPoint(*F1).getDimension(ind1); // p[3 * (*F1) + ind1]; // coordinates of node F1 in face 1
      double fb = coordinates.getPoint(*F1).getDimension(ind2); //p[3 * (*F1) + ind2];
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
      Vec3 c0 = coordinates.getPoint(*F0);
      Vec3 c1 = coordinates.getPoint(ind_n);

      Log::error("normal={}.", norm);
      Log::error("coordinate {} FACE0: ({}).", *F0, c0);
      Log::error("nearest match in FACE1: index {}, distance = {}, at ({}).", ind_n, dist, c1);
      RUNTIME_ERROR("Periodic boundaries do not match");
    }
  }
}


void PeriodicNodesMapping::setEdgePeriNodes(std::list<unsigned int> &edge0,
                                std::list<unsigned int> &edge1,
                                const int &dim,
                                const Coordinates &coordinates) {
  const double BIGNUM = 1e99;
  list<unsigned int> :: iterator e0;
  list<unsigned int> :: iterator e1;
  double eps = 1e-5;
  // LOOP OVER EACH NODE IN EDGE 0 AND FIND ITS EQUIVALENT IN OTHER EDGES
  for (e0 = edge0.begin() ; e0 != edge0.end() ; ++e0) {
    // outer corner loop - corn0
    bool found = false;
    double dist = 0;

    double p0 = coordinates.getPoint(*e0).getDimension(dim); // dim coordinate of node c0

    int iNearest = -1;
    double minDist = BIGNUM;
    // COMPARISON WITH C1
    for (e1 = edge1.begin() ; e1 != edge1.end() ; ++e1) {
      // inner corner loop - corn1
      double p1 = coordinates.getPoint(*e1).getDimension(dim);

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
      Vec3 node = coordinates.getPoint(*e0);
      Vec3 nearest = coordinates.getPoint(iNearest);
      Log::error("Edge node 1 not found.");
      Log::error("Edge node 0 index = {}, coordinates = ({}).", *e0, node);
      Log::error("Nearest node index = {}, coordinates = ({}).", iNearest, nearest);
      RUNTIME_ERROR("Edge nodde 1 not found.");
    }
  }// end loop over all nodes in edge0
}

void PeriodicNodesMapping::makePeriEquNodes(const PeriodicityType &periodicityType,
                                            const Mesh &e,
                                            const Coordinates& coordinates) {
  // GENERATES INDEX OF PERIODIC NODE EQUIVALENCIES
  // LIST OF INDEXES TO ALL PERIODIC NODES
  vector <unsigned int> nodes;
  e.listNodesOfMaterial(nodes , MAT_PERIODIC);
  size_t numPeri = nodes.size();
  if (!numPeri) {
    periNodes_.clear();
    return;
  }
  periNodes_.clear();
  const unsigned int np = coordinates.size();
  periNodes_.reserve(np);
  for (size_t i = 0 ; i < np ; i++)
    periNodes_.push_back(i);
  // NEED TO CONSIDER 3 DIFFERENT CASES, DEPENDING ON TYPE OF PERIODICITY OF MODELLING WINDOW
  const double eps = 1e-5; // accuracy for coordinate comparisons

  AABox boundingBox = coordinates.findBoundingBox();

  const double xmin = boundingBox.getXMin();
  const double xmax = boundingBox.getXMax();
  const double ymin = boundingBox.getYMin();
  const double ymax = boundingBox.getYMax();
  const double zmin = boundingBox.getZMin();
  const double zmax = boundingBox.getZMax();
  /// PROBABLY EVIL, BUT SO CONVENIENT...
#define LEFT    ( fabs(coordinates.getPoint(n).x() - xmin) <= eps )
#define RIGHT   ( fabs(coordinates.getPoint(n).x() - xmax) <= eps )
#define FRONT   ( fabs(coordinates.getPoint(n).y() - ymin) <= eps )
#define BACK    ( fabs(coordinates.getPoint(n).y() - ymax) <= eps )
#define BOTTOM  ( fabs(coordinates.getPoint(n).z() - zmin) <= eps )
#define TOP     ( fabs(coordinates.getPoint(n).z() - zmax) <= eps )
  //CASE 1 FRONT BACK IS PERIODIC ONLY
  if (periodicityType.isFrontBackPeriodic() &&
      !periodicityType.isLeftRightPeriodic() &&
      !periodicityType.isTopBottomPeriodic()) {
    //SEPARATE NODES INTO TWO LISTS FOR FRONT AND BACK SURFACES
    list <unsigned int> front;
    list <unsigned int> back;
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
    setFacePeriNodes(front, back, 1, coordinates);
  }
    // CASE 2 FRONT-BACK AND LEFT-RIGHT ARE PERIODIC
  else if (periodicityType.isFrontBackPeriodic() &&
           periodicityType.isLeftRightPeriodic() &&
           !periodicityType.isTopBottomPeriodic()) {
    // separate nodes into 8 lists, 4 x corners left/right and front/back planes
    list <unsigned int> edge0; //x = 0, y = 0
    list <unsigned int> edge1; //x = 0, y = max
    list <unsigned int> edge2; //x = max, y = max
    list <unsigned int> edge3; //x = max, y = 0
    list <unsigned int> front; //x = 0
    list <unsigned int> back;  //x = max
    list <unsigned int> right; //y = max
    list <unsigned int> left;  //y = 0
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
    setEdgePeriNodes(edge0, edge1, 2, coordinates);
    setEdgePeriNodes(edge0, edge2, 2, coordinates);
    setEdgePeriNodes(edge0, edge3, 2, coordinates);
    setFacePeriNodes(left  , right, 0, coordinates); // COMPARE Y, Z
    setFacePeriNodes(front , back , 1, coordinates);  // COMPARE X, Z
  }// END CASE 2
    // CASE 3 FRONT-BACK, LEFT-RIGHT AND TOP-BOTTOM ARE PERIODIC
  else if (periodicityType.isFrontBackPeriodic() &&
           periodicityType.isLeftRightPeriodic() &&
           periodicityType.isTopBottomPeriodic()) {
    // separate nodes into lists, 12 x edges left/right, front/back and top/bottom planes
    // Vertical corners along Z
    // Additionally, 7 corner nodes must point to bottom left (origin xmin,ymin,zmin) corner
    list <unsigned int> edge0; //x = 0, y = 0
    list <unsigned int> edge1; //x = 0, y = max
    list <unsigned int> edge2; //x = max, y = max
    list <unsigned int> edge3; //x = max, y = 0
    // Horizontal edges along X
    list <unsigned int> edgea; // y = 0, z = 0
    list <unsigned int> edgeb; // y = max, z = 0
    list <unsigned int> edgec; // y = max, z = max
    list <unsigned int> edged; // y = 0, z = max
    // Horizontal edges along Y
    list <unsigned int> edgeA; // x = 0, z = 0
    list <unsigned int> edgeB; // x = max, z = 0
    list <unsigned int> edgeC; // x = max, z = max
    list <unsigned int> edgeD; // x = 0, z = max
    list <unsigned int> front; //x = 0
    list <unsigned int> back;  //x = max
    list <unsigned int> right; //y = max
    list <unsigned int> left;  //y = 0
    list <unsigned int> top;   //z = max
    list <unsigned int> bottom;//z = 0;
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
        auto p = coordinates.getPoint(n);
        auto isOnTop = boundingBox.topFaceContains(p);
        throw std::runtime_error(fmt::format("Periodic node {} at {} is not on an external surface in {}, {}.",
                                             n, p, __FILE__, __func__));
      }
    }// end for i, loop over all nodes
    // CHECK THAT OPPOSITE FACES HAVE EQUAL NUMBER OF NODES
    {
      // start dummy scope
      int mincorner = *std::min_element(corner_nodes, corner_nodes + 8);
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
    setEdgePeriNodes(edge0, edge1 , 2, coordinates);   // VERTICAL EDGES
    setEdgePeriNodes(edge0, edge2 , 2, coordinates);
    setEdgePeriNodes(edge0, edge3 , 2, coordinates);
    setEdgePeriNodes(edgea, edgeb, 0, coordinates); // HORIZONTAL ALONG X
    setEdgePeriNodes(edgea, edgec, 0, coordinates);
    setEdgePeriNodes(edgea, edged, 0, coordinates);
    setEdgePeriNodes(edgeA, edgeB, 1, coordinates);  // HORIZONTAL ALONG Y
    setEdgePeriNodes(edgeA, edgeC, 1, coordinates);
    setEdgePeriNodes(edgeA, edgeD, 1, coordinates);
    setFacePeriNodes(left, right, 0, coordinates);   // COMPARE Y,Z
    setFacePeriNodes(front, back, 1, coordinates);   // COMPARE X,Z
    setFacePeriNodes(bottom, top, 2, coordinates);   // COMPARE X,Y
  }// end if 3 different periodicity cases
}

