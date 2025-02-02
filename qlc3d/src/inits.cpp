#include <iostream>
#include <memory>
#include <inits.h>
#include <simu.h>
#include <eventlist.h>
#include <refinement.h>
#include <io/meshreader.h>
#include <util/logging.h>
#include <util/exception.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <resultio.h>

#include <lc-representation.h>

/**
 * Goes through all triangle material numbers and tries to check that all is well.
 * @param triMaterials triangle material numbers
 * @param ne number of triangles
 */
void validateTriangleMaterials(const std::vector<idx> &triMaterials, const Electrodes &electrodes, const Alignment &alignment) {

    std::set<idx> uniqueMaterials(triMaterials.begin(), triMaterials.end());
    Log::info("Checking mesh surface materials. Found {} material numbers: {}",
              uniqueMaterials.size(), fmt::format("{}", fmt::join(uniqueMaterials, ", ")));

    for (idx i = 0; i < triMaterials.size(); ++i) {
        idx m = triMaterials[i];
        if ((m == MAT_PERIODIC) || (m == MAT_NEUMANN)) {
            continue;
        }
        size_t eNum = MATNUM_TO_ELECTRODE_NUMBER((size_t) m);
        size_t fNum = MATNUM_TO_FIXLC_NUMBER((size_t) m);
        if (m < MAT_ELECTRODE1) {
            RUNTIME_ERROR(fmt::format("Triangle {}, invalid material number {} < {}.", i, m, MAT_ELECTRODE1));
        }
        else if (eNum > 9) {
            RUNTIME_ERROR(fmt::format("Triangle {}, invalid Electrode number {} > 9", i, m));
        }
        else if (fNum > 9) {
            RUNTIME_ERROR(fmt::format("Triangle {}, invalid FixLC number {} > 9", i, m));
        }
        else if (eNum > electrodes.getnElectrodes()) {
          if (!electrodes.hasElectricField()) { // if no fixed electric field, the electrode must be properly defined
            RUNTIME_ERROR(
                    fmt::format("Triangle {} has electrode number {} but only {} electrode(s) have been defined.", i,
                                eNum, electrodes.getnElectrodes()));
          }
        }

        // If alignment surface material, check alignment has been defined for it
        size_t fixLCNum = MATNUM_TO_FIXLC_NUMBER(m);
        if (fixLCNum > 0 && !alignment.hasSurface(fixLCNum)) {
            RUNTIME_ERROR(fmt::format("Triangle {} has material number {} which is FixLC{} but no alignment surface has been defined for it.", i, m, fixLCNum));
        }
    }
}

/**
 * Goes through each material number for tetrahedra. Throws runtime_exception if invalid ones are found.
 */
void validateTetrahedralMaterials(const std::vector<idx> matt) {
    for (idx i = 0; i < matt.size(); ++i) {
        const idx m = matt[i];
        if ((m == MAT_DOMAIN1) ||
                (m==MAT_DIELECTRIC1) || (m==MAT_DIELECTRIC2) || (m==MAT_DIELECTRIC3) ||
                (m==MAT_DIELECTRIC4) || (m==MAT_DIELECTRIC5) || (m==MAT_DIELECTRIC6) ||
                (m==MAT_DIELECTRIC7) ) {
            continue;
        }
        else {
            RUNTIME_ERROR(fmt::format("Invalid material number {} found in mesh for tetrahedron {}. "
                                       "Valid material numbers are Domain1 and Dielectrics 1-7.", m, i));
        }
    }
}

/**
 * Find the index of the node in the nodes vector that is nearest to the target point.
 */
int findNearestNodeIndex(const Vec3 &target, const std::vector<Vec3> &nodes) {
  double minDist = std::numeric_limits<double>::max();
  int nearestIndex = -1;
  for (size_t i = 0; i < nodes.size(); i++) {
    double dist = target.distanceSquared(nodes[i]);
    if (dist < minDist) {
      minDist = dist;
      nearestIndex = i;
    }
  }
  return nearestIndex;
}

/**
 * Reorder the tetrahedral element node indices if the current ordering is the GMSH ordering.
 */
void reorderQuadraticTetNodeOrder(vector<idx> &tetNodes, const Coordinates &coords) {

  size_t numTets = tetNodes.size() / 10;
  Log::info("Checking quadratic tetrahedron element node oder for {} elements.", numTets);

  std::vector<Vec3> elemCoords;
  elemCoords.resize(10, Vec3());
  for (size_t i = 0; i < numTets; i++) {
    idx ind0 = i * 10;

    coords.loadCoordinates(&tetNodes[ind0], &tetNodes[ind0 + 10], &elemCoords[0]);

    // swap the two last node indices in tet element if current ordering is the GMSH ordering.
    // The expected mid-edge node positions are:
    // p8 = (p1 + p3) / 2
    // p9 = (p2 + p3) / 2
    Vec3 p8 = (elemCoords[1] + elemCoords[3]) * 0.5;
    Vec3 p9 = (elemCoords[2] + elemCoords[3]) * 0.5;

    // find indices to the actual nearest ones. This avoids comparison tolerance complications
    // but may not detect if the nodes are really far from expected positions.
    int indNearest8 = findNearestNodeIndex(p8, elemCoords);
    int indNearest9 = findNearestNodeIndex(p9, elemCoords);

    bool isExpectedOrder = indNearest8 == 8 && indNearest9 == 9;
    bool isGmsOrder = indNearest8 == 9 && indNearest9 == 8;
    if (isExpectedOrder) {
      continue;
    } else if (isGmsOrder) {
      std::swap(tetNodes[ind0 + 8], tetNodes[ind0 + 9]);
    } else {
      RUNTIME_ERROR(fmt::format("Tetrahedron {} has unexpected node ordering. Expected: 8, 9. Found: {}, {}.",
                                i, indNearest8, indNearest9));
    }
  }
}

void prepareGeometry(Geometry &geom,
                     RawMeshData &rawMeshData,
                     Electrodes &electrodes,
                     const Alignment &alignment,
                     const Vec3 &stretchVector,
                     unsigned int regularGridCountX,
                     unsigned int regularGridCountY,
                     unsigned int regularGridCountZ) {
  Log::info("mesh element order = {}, triangles count = {}, tetrahedra count = {}, nodes count = {}",
            rawMeshData.elementOrder, rawMeshData.triMaterials.size(), rawMeshData.tetMaterials.size(), rawMeshData.points.size());
  validateTriangleMaterials(rawMeshData.triMaterials, electrodes, alignment);
  validateTetrahedralMaterials(rawMeshData.tetMaterials);

  auto coordinates = std::make_shared<Coordinates>(std::move(rawMeshData.points));

  if (rawMeshData.elementOrder == 2) {
    reorderQuadraticTetNodeOrder(rawMeshData.tetNodes, *coordinates);
  }

  coordinates->scale(stretchVector);

  geom.setMeshData(rawMeshData.elementOrder, coordinates,
                   std::move(rawMeshData.tetNodes), std::move(rawMeshData.tetMaterials),
                   std::move(rawMeshData.triNodes), std::move(rawMeshData.triMaterials));

  geom.makeRegularGrid(regularGridCountX, regularGridCountY, regularGridCountZ);
}

void prepareGeometryWithDefaultBoundaries(Geometry& geom,
                                          const std::filesystem::path &meshFileName) {
  // read mesh data from file. Allocates the data arrays.
  RawMeshData rawMeshData = MeshReader::readMesh(meshFileName);

  // find electrodes from raw mesh data
  std::unordered_set<size_t> electrodeNumbers = findElectrodeNumbers(rawMeshData.triMaterials.begin(), rawMeshData.triMaterials.end());
  Electrodes electrodes = Electrodes::withInitialPotentials(std::vector<unsigned int>(electrodeNumbers.begin(), electrodeNumbers.end()), 0.0);

  // find alignment surfaces from raw mesh data
  std::unordered_set<size_t> fixLcNumbers = findFixlcNumbers(rawMeshData.triMaterials.begin(), rawMeshData.triMaterials.end());
  Alignment alignment = Alignment();
  for (auto fixLcNumber : fixLcNumbers) {
    alignment.addSurface(Surface::ofStrongAnchoring(fixLcNumber, 0, 0));
  }

  prepareGeometry(geom, rawMeshData, electrodes, alignment, {1, 1, 1}, 0, 0, 0);
}

void prepareGeometry(Geometry& geom,
                     const std::filesystem::path &meshFileName,
                     Electrodes& electrodes,
                     const Alignment& alignment,
                     const Vec3 &stretchVector,
                     unsigned int regularGridCountX,
                     unsigned int regularGridCountY,
                     unsigned int regularGridCountZ) {

    // read mesh data from file. Allocates the data arrays.
    RawMeshData rawMeshData = MeshReader::readMesh(meshFileName);
    prepareGeometry(geom, rawMeshData, electrodes, alignment, stretchVector, regularGridCountX, regularGridCountY, regularGridCountZ);
}

FILE* createOutputEnergyFile(Simu& simu) {
    FILE* fid = nullptr;
    if (simu.getOutputEnergy() == 1) {
      std::filesystem::path energyFilePath = simu.getSaveDirAbsolutePath() / "energy.m";
      fid = fopen( energyFilePath.string().c_str(), "w");
      if (fid == nullptr) {
        RUNTIME_ERROR(fmt::format("could not open output file for free energy by filename = {}.", energyFilePath));
      }
    }
    return fid;
}

void initialiseLcSolutionVector(SolutionVector &q,
                                const Simu &simu,
                                const LC &lc,
                                const InitialVolumeOrientation &boxes,
                                Alignment &alignment,
                                Geometry &geom) {
  const double S0 = lc.S0();
  boxes.setVolumeQ(q, S0, geom.getCoordinates());
  if (!simu.getLoadQ().empty()) {
    ResultIO::ReadResult(simu.getLoadQ(), q);
  }
  setSurfacesQ(q, alignment, S0, geom);
  q.initialiseLcBoundaries(geom, alignment);
}

void createMeshRefinementEvents(const MeshRefinement &meshRefinement,
                                EventList &eventListOut) {
    const unsigned int repRefIter = meshRefinement.getRepRefIter();
    const double repRefTime = meshRefinement.getRepRefTime();
    bool hasPeriodicRefinement = repRefIter > 0 || repRefTime > 0;
    Log::info("Creating periodic mesh refinement events at every {} iterations, {} seconds", repRefIter, repRefTime);
    eventListOut.setRepRefIter(repRefIter);
    eventListOut.setRepRefTime(repRefTime);

    if (hasPeriodicRefinement) {
        // create refinement event for periodically occurring mesh refinement.
        // These have no explicitly defined iterations or times
        for (auto &ref : meshRefinement.getRefinementConfig()) {
            if (!ref.occursPeriodically()) {
                continue;
            }
            RefInfo *info = RefInfo::ofPeriodicMeshRefinement(ref.type_, ref.values_, ref.x_, ref.y_, ref.z_);
            Event *event = Event::ofPeriodicMeshRefinement(info);
            eventListOut.addRepRefInfo(event);
        }

        // e.g. if user has defined RepRefIter > 0 or RepRefTime > 0, but no refinement objects with
        // empty explicit times/iterations lists
        if (eventListOut.getNumPeriodicRefinementObjects() == 0) {
            RUNTIME_ERROR("No refinement objects defined for periodic mesh refinement");
        }
    }

    // add refinement events occurring at explicitly defined iterations/times
    for (auto &ref : meshRefinement.getRefinementConfig()) {
        if (ref.occursPeriodically()) {
            continue;
        }

        // convert each explicitly defined refinement iteration to an *Event and add to event list
        for (unsigned int iter : ref.iterations_) {
            RefInfo * refInfo = RefInfo::make(ref.type_, iter, -1,
                                              ref.values_, ref.x_, ref.y_, ref.z_);
            auto *event = new Event(EVENT_REFINEMENT, iter, (void*) refInfo);
            eventListOut.insertIterEvent(event);
        }

        // convert each explicitly defined refinement time to an *Event adn add to event list
        for (double time : ref.times_) {
            RefInfo *refInfo = RefInfo::make(ref.type_, -1, time,
                                             ref.values_, ref.x_, ref.y_, ref.z_);
            auto *event = new Event(EVENT_REFINEMENT, time, (void*) refInfo);
            eventListOut.insertTimeEvent(event);
        }
    }
}

void createElectrodeSwitchingEvents(const Electrodes &electrodes,
                                    EventList &eventListOut) {
  // pass ownership of event to eventlist
  std::vector<Event*> events = electrodes.createSwitchingEvents();
  for (auto event : events) {
    eventListOut.insertTimeEvent(event);
  }
}