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

PeriodicNodesMapping::PeriodicNodesMapping(const Mesh &tris,
                                           const Coordinates &coords) : periodicityType(tris) {
  Log::info("Initialising periodic surfaces");
  makePeriEquNodes(tris, coords);
}

/*
unsigned int PeriodicNodesMapping::getPeriodicNode(unsigned int node) const {
  if (periodicityType.isAnyPeriodic()) {
    return periNodes_[node];
  } else {
    return node;
  }
}
*/

/**
 * Set the periodic nodes for the given mesh and coordinates.
 * @param norm index to face normal axis, 0 = x, 1 = y, 2 = z. Anything else will cause error.
 * @param coordinates The coordinates.
 */
void PeriodicNodesMapping::setFacePeriNodes(const std::vector<unsigned int> &face0,
                                const std::vector<unsigned int> &face1,
                                int norm,
                                const Coordinates &coordinates) {
  if (norm < 0 || norm > 2) {
    RUNTIME_ERROR(fmt::format("Invalid norm value {}.", norm));
  }

  const double eps = 1e-5; // accuracy of coordinate comparison

  for (auto f0 : face0) {
    bool found = false;

    unsigned int indNearestSoFar  = std::numeric_limits<unsigned int>::max();
    double minDistanceSoFar = std::numeric_limits<double>::max();
    Vec3 p0 = coordinates.getPoint(f0);

    for (auto &f1 : face1) {
      Vec3 p1 = coordinates.getPoint(f1);

      // set the normal axis dimension to same value for both points. The, point-to-point distance is an in-plane calculation.
      if (norm == 0) {
        p1.set(p0.x(), p1.y(), p1.z());
      } else if (norm == 1) {
        p1.set(p1.x(), p0.y(), p1.z());
      } else { // norm == 2
        p1.set(p1.x(), p1.y(), p0.z());
      }

      double currentDistance = p0.distanceSquared(p1);

      if (currentDistance < minDistanceSoFar) { // debug info only, keep track of nearest found node
        minDistanceSoFar = currentDistance; // nearest distance
        indNearestSoFar = f1;  // index to nearest distance
      }

      if (currentDistance < eps * eps) { // compare squared distances
        this->periNodes_[f1] = f0;
        found = true;
        break;
      }
    }

    if (!found) {
      Vec3 c0 = coordinates.getPoint(f0);
      Vec3 c1 = coordinates.getPoint(indNearestSoFar);

      Log::error("normal={}.", norm);
      Log::error("coordinate {} FACE0: ({}).", f0, c0);
      Log::error("nearest match in FACE1: index {}, distance = {}, at ({}).", indNearestSoFar, minDistanceSoFar, c1);
      RUNTIME_ERROR("Periodic boundaries do not match");
    }
  }
}

/**
 * Maps nodes on edge1 to nodes on edge0.
 * @param edge0
 * @param edge1
 * @param dim index to dimension that coincides with edge direction, 0 = x, 1 = y, 2 = z. Anything else will cause error.
 * @param coordinates The coordinates
 */
void PeriodicNodesMapping::setEdgePeriNodes(const std::vector<unsigned int> &edge0,
                                            const std::vector<unsigned int> &edge1,
                                            int dim,
                                            const Coordinates &coordinates) {

  const double eps = 1e-5;
  for (auto &e0 : edge0) {
    bool found = false;
    unsigned int iNearest = std::numeric_limits<unsigned int>::max();
    double minDist = std::numeric_limits<double>::max();

    double p0 = coordinates.getPoint(e0).getDimension(dim); // dim coordinate of node c0

    for (auto &e1 : edge1) {
      // inner corner loop - corn1
      double p1 = coordinates.getPoint(e1).getDimension(dim);

      double dist = fabs(p0 - p1);
      if (dist <= minDist) {
        minDist = dist;
        iNearest = e1;
      }

      if (dist <= eps) { // if same coordinate -> match!
        periNodes_[e1] = e0;
        found = true;
        break;
      }
    }
    if (!found) {
      Vec3 node = coordinates.getPoint(e0);
      Vec3 nearest = coordinates.getPoint(iNearest);
      Log::error("Edge node 1 not found.");
      Log::error("Edge node 0 index = {}, coordinates = ({}).", e0, node);
      Log::error("Nearest node index = {}, coordinates = ({}).", iNearest, nearest);
      RUNTIME_ERROR("Edge node 1 not found.");
    }
  }
}

void PeriodicNodesMapping::makePeriEquNodes(const Mesh &e,
                                            const Coordinates& coordinates) {
  // GENERATES INDEX OF PERIODIC NODE EQUIVALENCIES
  // LIST OF INDEXES TO ALL PERIODIC NODES
  vector <unsigned int> nodes;
  e.listNodesOfMaterial(nodes , MAT_PERIODIC);
  periNodes_.clear();
  const unsigned int np = coordinates.size();
  periNodes_.reserve(np);
  for (size_t i = 0 ; i < np ; i++) {
    periNodes_.push_back(i);
  }

  if (!periodicityType.isAnyPeriodic()) {
    Log::info("No periodic surfaces detected.");
    // no periodic nodes, each node is equivalent to itself
    return;
  }

  Log::info("Finding periodic nodes for periodic surfaces on left-right={}, front-back={}, top-bottom={}",
            periodicityType.isLeftRightPeriodic(), periodicityType.isFrontBackPeriodic(), periodicityType.isTopBottomPeriodic());

  // NEED TO CONSIDER 3 DIFFERENT CASES, DEPENDING ON TYPE OF PERIODICITY OF MODELLING WINDOW
  const double eps = 1e-5; // accuracy for coordinate comparisons

  AABox boundingBox = coordinates.findBoundingBox();

  // Case 1 - Front-Back faces are periodic only
  if (periodicityType.isFrontBackPeriodic() &&
      !periodicityType.isLeftRightPeriodic() &&
      !periodicityType.isTopBottomPeriodic()) {
    std::vector<unsigned int> front;
    std::vector<unsigned int> back;
    for (auto &n : nodes) {
      Vec3 p = coordinates.getPoint(n);
      if (boundingBox.frontFaceContains(p, eps)) {
        front.push_back(n) ;
      }
      else if (boundingBox.backFaceContains(p, eps)) {
        back.push_back(n);
      } else {
        RUNTIME_ERROR(fmt::format("Expected node {} at {} on front/back surface.", n, p));
      }
    }
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
    std::vector<unsigned int> edge0; //x = 0, y = 0
    std::vector<unsigned int> edge1; //x = 0, y = max
    std::vector<unsigned int> edge2; //x = max, y = max
    std::vector<unsigned int> edge3; //x = max, y = 0
    std::vector<unsigned int> front; //x = 0
    std::vector<unsigned int> back;  //x = max
    std::vector<unsigned int> right; //y = max
    std::vector<unsigned int> left;  //y = 0

    for (auto &n : nodes) {
      const Vec3 &p = coordinates.getPoint(n);

      if (boundingBox.leftFaceContains(p, eps) && boundingBox.frontFaceContains(p, eps)) {
        edge0.push_back(n);
      } else if (boundingBox.leftFaceContains(p, eps) && boundingBox.backFaceContains(p, eps)) {
        edge1.push_back(n);
      } else if (boundingBox.rightFaceContains(p, eps)  && boundingBox.backFaceContains(p, eps)) {
        edge2.push_back(n);
      } else if (boundingBox.rightFaceContains(p, eps) && boundingBox.frontFaceContains(p, eps)) {
        edge3.push_back(n);
      } else if (boundingBox.frontFaceContains(p, eps)) {
        front.push_back(n);
      } else if (boundingBox.backFaceContains(p, eps)) {
        back.push_back(n);
      } else if (boundingBox.leftFaceContains(p, eps)) {
        left.push_back(n);
      } else if (boundingBox.rightFaceContains(p, eps)) {
        right.push_back(n);
      } else {
        RUNTIME_ERROR(fmt::format("Expected node {} at {} on front/back/left/right surface.", n, p));
      }
    }
    if ((edge0.size() != edge1.size()) ||
        (edge0.size() != edge2.size())  ||
        (edge0.size() != edge3.size())  ||
        (left.size() != right.size()) ||
        (front.size() != back.size())) {
      Log::error("Periodic node counts don't match.");
      Log::error("Corners 0, 1, 2, 3 contain [{}, {}, {}, {}] nodes.", edge0.size(), edge1.size(), edge2.size(), edge3.size());
      Log::error("Front/back, left/right surfaces contain {}/{}, {}/{} nodes.", front.size(), back.size(), left.size(), right.size());
      RUNTIME_ERROR("Periodic node counts don't match.");
    }
    setEdgePeriNodes(edge0, edge1, 2, coordinates);
    setEdgePeriNodes(edge0, edge2, 2, coordinates);
    setEdgePeriNodes(edge0, edge3, 2, coordinates);
    setFacePeriNodes(left  , right, 0, coordinates); // Compare y, z
    setFacePeriNodes(front , back , 1, coordinates); // Compare x, z
  }
  // CASE 3 FRONT-BACK, LEFT-RIGHT AND TOP-BOTTOM ARE PERIODIC
  else if (periodicityType.isFrontBackPeriodic() &&
           periodicityType.isLeftRightPeriodic() &&
           periodicityType.isTopBottomPeriodic()) {
    // separate nodes into lists, 12 x edges left/right, front/back and top/bottom planes
    // Vertical corners along Z
    // Additionally, 7 corner nodes must point to bottom left (origin xmin,ymin,zmin) corner
    std::vector<unsigned int> edge0; //x = 0, y = 0
    std::vector<unsigned int> edge1; //x = 0, y = max
    std::vector<unsigned int> edge2; //x = max, y = max
    std::vector<unsigned int> edge3; //x = max, y = 0
    // Horizontal edges along X
    std::vector<unsigned int> edgea; // y = 0, z = 0
    std::vector<unsigned int> edgeb; // y = max, z = 0
    std::vector<unsigned int> edgec; // y = max, z = max
    std::vector<unsigned int> edged; // y = 0, z = max
    // Horizontal edges along Y
    std::vector<unsigned int> edgeA; // x = 0, z = 0
    std::vector<unsigned int> edgeB; // x = max, z = 0
    std::vector<unsigned int> edgeC; // x = max, z = max
    std::vector<unsigned int> edgeD; // x = 0, z = max
    std::vector<unsigned int> front; //x = 0
    std::vector<unsigned int> back;  //x = max
    std::vector<unsigned int> right; //y = max
    std::vector<unsigned int> left;  //y = 0
    std::vector<unsigned int> top;   //z = max
    std::vector<unsigned int> bottom;//z = 0;
    // LOOP OVER ALL NODES AND INSERT TO CORRECT LIST
    unsigned int corner_nodes[8] = { std::numeric_limits<unsigned int>::max(),
                                     std::numeric_limits<unsigned int>::max(),
                                     std::numeric_limits<unsigned int>::max(),
                                     std::numeric_limits<unsigned int>::max(),
                                     std::numeric_limits<unsigned int>::max(),
                                     std::numeric_limits<unsigned int>::max(),
                                     std::numeric_limits<unsigned int>::max(),
                                     std::numeric_limits<unsigned int>::max()};

    for (auto n : nodes) {
      const auto &p = coordinates.getPoint(n);
      // CORNER NODES TAKE PRECEDENCE OVER OTHER NODES
      // front-left-bottom corner is independent
      if (boundingBox.frontFaceContains(p, eps) && boundingBox.leftFaceContains(p, eps) && boundingBox.bottomFaceContains(p, eps)) {
        corner_nodes[0] = n;
      } else if (boundingBox.frontFaceContains(p, eps) && boundingBox.rightFaceContains(p, eps) && boundingBox.bottomFaceContains(p, eps)) {
        corner_nodes[1] = n;
      } else if (boundingBox.frontFaceContains(p, eps) && boundingBox.leftFaceContains(p, eps) && boundingBox.topFaceContains(p, eps)) {
        corner_nodes[2] = n;
      } else if (boundingBox.frontFaceContains(p, eps) && boundingBox.rightFaceContains(p, eps) && boundingBox.topFaceContains(p, eps)) {
        corner_nodes[3] = n;
      } else  if (boundingBox.backFaceContains(p, eps) && boundingBox.leftFaceContains(p, eps) && boundingBox.bottomFaceContains(p, eps)) {
        corner_nodes[4] = n;
      } else  if (boundingBox.backFaceContains(p, eps) && boundingBox.rightFaceContains(p, eps) && boundingBox.bottomFaceContains(p, eps)) {
        corner_nodes[5] = n;
      } else  if (boundingBox.backFaceContains(p, eps) && boundingBox.leftFaceContains(p, eps) && boundingBox.topFaceContains(p, eps)) {
        corner_nodes[6] = n;
      } else  if (boundingBox.backFaceContains(p, eps) && boundingBox.rightFaceContains(p, eps) && boundingBox.topFaceContains(p, eps)) {
        corner_nodes[7] = n;
        // EDGE NODES
        // 4 x Vertical Corners
      } else if (boundingBox.leftFaceContains(p, eps) && boundingBox.frontFaceContains(p, eps)) {
        edge0.push_back(n);
      } else if (boundingBox.leftFaceContains(p, eps) && boundingBox.backFaceContains(p, eps)) {
        edge1.push_back(n);
      } else if (boundingBox.rightFaceContains(p, eps) && boundingBox.backFaceContains(p, eps)) {
        edge2.push_back(n);
      } else if (boundingBox.rightFaceContains(p, eps) && boundingBox.frontFaceContains(p, eps)) { // edge3
        edge3.push_back(n);
      }
        // 4 x Horizontal along X
      else if (boundingBox.frontFaceContains(p, eps) && boundingBox.bottomFaceContains(p, eps)) {
        edgea.push_back(n);
      } else if (boundingBox.backFaceContains(p, eps) && boundingBox.bottomFaceContains(p, eps)) {
        edgeb.push_back(n);
      } else if (boundingBox.backFaceContains(p, eps) && boundingBox.topFaceContains(p, eps)) {
        edgec.push_back(n);
      } else if (boundingBox.frontFaceContains(p, eps) && boundingBox.topFaceContains(p, eps)) {
        edged.push_back(n);
      }
        // 4 x Horizontal along Y
      else if (boundingBox.leftFaceContains(p, eps) && boundingBox.bottomFaceContains(p, eps)) {
        edgeA.push_back(n);   // xmin and zmin
      } else if (boundingBox.rightFaceContains(p, eps) && boundingBox.bottomFaceContains(p, eps)) {
        edgeB.push_back(n);
      } else if (boundingBox.rightFaceContains(p, eps) && boundingBox.topFaceContains(p, eps)) {
        edgeC.push_back(n);
      } else if (boundingBox.leftFaceContains(p, eps) && boundingBox.topFaceContains(p, eps)) {
        edgeD.push_back(n);
      }

      // FRONT/BACK, LEFT/RIGHT, TOP/BOTTOM FACES
      else if (boundingBox.frontFaceContains(p, eps)) {
        front.push_back(n);
      } else if (boundingBox.backFaceContains(p, eps)) {
        back.push_back(n);
      } else if (boundingBox.leftFaceContains(p, eps)) {
        left.push_back(n);
      } else if (boundingBox.rightFaceContains(p, eps)) {
        right.push_back(n);
      } else if (boundingBox.bottomFaceContains(p, eps)) {
        bottom.push_back(n);
      } else if (boundingBox.topFaceContains(p, eps)) {
        top.push_back(n);
      } else {
        RUNTIME_ERROR(fmt::format("Periodic node {} at {} is not on an external surface.", n, p));
      }
    }// end for i, loop over all nodes
    // CHECK THAT OPPOSITE FACES HAVE EQUAL NUMBER OF NODES
    {
      // start dummy scope
      unsigned int mincorner = *std::min_element(corner_nodes, corner_nodes + 8);
      if (mincorner == std::numeric_limits<unsigned int>::max()) {
        RUNTIME_ERROR(fmt::format("Periodic corner nodes not found. "
                                             "Indices are [{}. {}, {}, {}, {}, {}, {}, {}].",
                                             corner_nodes[0], corner_nodes[1], corner_nodes[2],
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

