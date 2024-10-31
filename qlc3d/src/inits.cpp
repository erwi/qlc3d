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
            RUNTIME_ERROR(fmt::format("Triangle {} has electrode number {} but only {} electrode(s) have been defined.", i, eNum, electrodes.getnElectrodes()));
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

    Log::info("num points = {}", rawMeshData.points.size());

    // Throw exception if invalid element material numbers are detected.
    validateTriangleMaterials(rawMeshData.triMaterials, electrodes, alignment);
    validateTetrahedralMaterials(rawMeshData.tetMaterials);

    auto coordinates = std::make_shared<Coordinates>(std::move(rawMeshData.points));
    coordinates->scale(stretchVector);

    geom.setMeshData(coordinates,
                   std::move(rawMeshData.tetNodes), std::move(rawMeshData.tetMaterials),
                   std::move(rawMeshData.triNodes), std::move(rawMeshData.triMaterials));

    geom.makeRegularGrid(regularGridCountX, regularGridCountY, regularGridCountZ);
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

void initialiseLcSolutionVector(SolutionVector &q, const Simu &simu, const LC &lc, const Boxes &boxes, const Alignment &alignment, const Geometry &geom) {
  const double S0 = lc.S0();
  setVolumeQ(q, S0, boxes, geom.getCoordinates());
  if (!simu.getLoadQ().empty()) {
    ResultIO::ReadResult(simu.getLoadQ(), q);
  }
  setSurfacesQ(q, alignment, S0, geom);
  q.setFixedNodesQ(alignment, geom.getTriangles());  // set fixed surface anchoring
  q.setPeriodicEquNodes(geom);          // periodic nodes
  q.EnforceEquNodes(geom);                // makes sure values at periodic boundaries match
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