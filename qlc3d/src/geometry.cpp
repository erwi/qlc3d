#include <geometry.h>
#include <geom/tet-mesh-search.h>
#include <material_numbers.h>
#include <util/logging.h>
#include <util/exception.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <geom/periodicity.h>
#include <mesh/mesh.h>
#include <cassert>

const idx Geometry::NOT_AN_INDEX = std::numeric_limits<idx>::max();

Geometry::Geometry() {
  npLC = 0;
  t = Mesh::tetMesh();
  e = Mesh::triangleMesh();
  boundingBox = AABox();
  periodicNodesMapping.reset();
}

Geometry::~Geometry() = default;

void Geometry::setTo(Geometry *geom) {
    this->ClearGeometry();
    npLC    = geom->getnpLC();
    boundingBox = geom->boundingBox;
    t->CopyMesh(geom->t.get());
    e->CopyMesh(geom->e.get());

    nodeNormals = geom->nodeNormals;
    coordinates_ = geom->getCoordinates().clone();


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

void Geometry::genIndToTetsByCoords(vector<unsigned int> &returnIndex,
                                     const Coordinates &targetCoordinates,
                                     const bool &terminateOnError,
                                     const bool &requireLCElement) {
    TetMeshSearch search(getTetrahedra(), getCoordinates(), getBoundingBox());
    search.genIndToTetsByCoords(returnIndex, targetCoordinates, terminateOnError, requireLCElement);
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