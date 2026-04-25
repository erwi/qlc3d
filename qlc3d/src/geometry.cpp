#include <geometry.h>
#include <material_numbers.h>
#include <util/logging.h>
#include <util/exception.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <geom/periodicity.h>
#include <mesh/mesh.h>

const idx Geometry::NOT_AN_INDEX = std::numeric_limits<idx>::max();

Geometry::Geometry():
    regularGrid(nullptr) {
  npLC = 0;
  t = Mesh::tetMesh();
  e = Mesh::triangleMesh();
  boundingBox = AABox();
  periodicNodesMapping.reset();
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

    if (this->regularGrid) delete regularGrid;
    if (geom->regularGrid)
        regularGrid = new RegularGrid(*geom->regularGrid);

    periodicNodesMapping.reset();
}

void Geometry::ClearGeometry() {
    npLC = 0;
    nodeNormals.clear();
    boundingBox.setBoundsX(0, 0);
    boundingBox.setBoundsY(0, 0);
    boundingBox.setBoundsZ(0, 0);
    t->ClearMesh();
    e->ClearMesh();
    periodicNodesMapping.reset();
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

void Geometry::setMeshData(unsigned int elementOrder, const std::shared_ptr<Coordinates> &coordinates,
                 std::vector<unsigned int> &&tetNodes, std::vector<unsigned int> &&tetMaterials,
                 std::vector<unsigned int> &&triNodes, std::vector<unsigned int> &&triMaterials) {
  if (elementOrder != 1 && elementOrder != 2) {
    RUNTIME_ERROR(fmt::format("Invalid element order: {}", elementOrder));
  }

  setCoordinates(coordinates);
  t->setElementData(elementOrder == 1 ? ElementType::LINEAR_TETRAHEDRON : ElementType::QUADRATIC_TETRAHEDRON, std::move(tetNodes), std::move(tetMaterials));
  e->setElementData(elementOrder == 1 ? ElementType::LINEAR_TRIANGLE : ElementType::QUADRATIC_TRIANGLE, std::move(triNodes), std::move(triMaterials));
  ReorderDielectricNodes();
  e->setConnectedVolume(t.get()); // neighbour index tri -> tet
  t->calculateDeterminants3D(getCoordinates()); // calculate tetrahedral determinants
  t->ScaleDeterminants(qlc3d::units::CUBIC_MICROMETER_TO_CUBIC_METER);

  e->calculateSurfaceNormals(getCoordinates(), t.get()); // calculate triangle determinants and surface normal vectors
  e->ScaleDeterminants(qlc3d::units::SQUARE_MICROMETER_TO_SQUARE_METER);

  calculateNodeNormals();
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

    // Generate a new node order where LC nodes come first and dielectric-only nodes follow.
    std::vector<idx> isLcNode(getnp(), 0);
    for (idx i = 0 ; i < t->getnElements() ; i ++) {
        if (t->getMaterialNumber(i) == MAT_DOMAIN1) {
            for (idx j = 0 ; j < t->getnNodes() ; j ++) {
                isLcNode[t->getNode(i, j)] = 1;
            }
        }
    }

    npLC = 0;
    std::vector<idx> reorderedNodeIndices;
    reorderedNodeIndices.reserve(getnp());
    for (idx i = 0; i < getnp(); i++) {
        if (isLcNode[i] == 1) {
            reorderedNodeIndices.push_back(i);
            npLC++;
        }
    }
    for (idx i = 0; i < getnp(); i++) {
        if (isLcNode[i] == 0) {
            reorderedNodeIndices.push_back(i);
        }
    }

    if (reorderedNodeIndices.size() != getnp()) {
        RUNTIME_ERROR(fmt::format("Failed to reorder dielectric nodes: expected {} nodes, built {}.", getnp(), reorderedNodeIndices.size()));
    }

    std::vector<idx> inverseNodeMap(getnp(), NOT_AN_INDEX);
    for (idx i = 0; i < getnp(); i++) {
        inverseNodeMap[reorderedNodeIndices[i]] = i;
    }

    std::vector<Vec3> reorderedCoords;
    reorderedCoords.reserve(coordinates_->size());
    for (idx i = 0; i < getnp(); i++) {
        reorderedCoords.push_back(coordinates_->getPoint(reorderedNodeIndices[i]));
    }
    coordinates_ = std::make_shared<Coordinates>(std::move(reorderedCoords));

    std::vector<idx> reorderedTetNodes(t->getnElements() * t->getnNodes());
    std::vector<idx> tetMaterials(t->getnElements());
    std::vector<idx> oldTetNodes(t->getnNodes());
    for (idx i = 0 ; i < t->getnElements() ; i++) {
        t->loadNodes(i, oldTetNodes.data());
        tetMaterials[i] = t->getMaterialNumber(i);
        for (idx j = 0 ; j < t->getnNodes() ; j++) {
            const idx mappedNode = inverseNodeMap[oldTetNodes[j]];
            if (mappedNode == NOT_AN_INDEX) {
                RUNTIME_ERROR(fmt::format("Encountered unmapped tetrahedral node {} while reordering dielectric nodes.", oldTetNodes[j]));
            }
            reorderedTetNodes[i * t->getnNodes() + j] = mappedNode;
        }
    }
    t->setElementData(t->getElementType(), std::move(reorderedTetNodes), std::move(tetMaterials));

    std::vector<idx> reorderedTriNodes(e->getnElements() * e->getnNodes());
    std::vector<idx> triMaterials(e->getnElements());
    std::vector<idx> oldTriNodes(e->getnNodes());
    for (idx i = 0 ; i < e->getnElements() ; i++) {
        e->loadNodes(i, oldTriNodes.data());
        triMaterials[i] = e->getMaterialNumber(i);
        for (idx j = 0 ; j < e->getnNodes() ; j++) {
            const idx mappedNode = inverseNodeMap[oldTriNodes[j]];
            if (mappedNode == NOT_AN_INDEX) {
                RUNTIME_ERROR(fmt::format("Encountered unmapped triangle node {} while reordering dielectric nodes.", oldTriNodes[j]));
            }
            reorderedTriNodes[i * e->getnNodes() + j] = mappedNode;
        }
    }
    e->setElementData(e->getElementType(), std::move(reorderedTriNodes), std::move(triMaterials));

    this->updateMaxNodeNumbers();
}

void Geometry::makeRegularGrid(const size_t &nx,
                               const size_t &ny,
                               const size_t &nz) {
    // CREATES REUGLAR GRID OBJECT USED FOR INTERPOLATING VALUES FROM
    // TETRAHEDRAL MESH ONTO REGULARLY SPACED GRID
    if (nx == 0 || ny == 0 || nz == 0) {
      Log::info("Regular grid generation disabled, nx={}, ny={}, nz={}", nx, ny, nz);
      return;
    }
    Log::info("Generating regular grid lookup with grid size nx={}, ny={}, nz={}", nx, ny, nz);

    auto elementType = getTetrahedra().getElementType();
    if (elementType != ElementType::LINEAR_TETRAHEDRON) {
      throw NotYetImplementedException("Regular grid generation is only implemented for first order elements, got elementType=" +
                                               toString(elementType));
    }

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

const Coordinates& Geometry::getCoordinates() const {
  if (!coordinates_) {
    RUNTIME_ERROR("Coordinates not set");
  }
  return *coordinates_;
}

const std::vector<Vec3>& Geometry::getNodeNormals() const {
  return nodeNormals;
}

Vec3 Geometry::getNodeNormal(unsigned int i) const {
  return nodeNormals[i];
}

PeriodicNodesMapping& Geometry::createPeriodicNodesMapping() {
  if (periodicNodesMapping == nullptr) {
    periodicNodesMapping = std::make_unique<PeriodicNodesMapping>(getTriangles(), getCoordinates());
  }

  return *periodicNodesMapping.get();
}

void Geometry::clearPeriodicNodesMapping() {
  periodicNodesMapping.reset();
}