#include <io/orientation-loader.h>
#include <io/orientation-reader.h>
#include <io/orientation-assignment.h>
#include <solutionvector.h>
#include <geom/vec3.h>
#include <simu.h>
#include <util/logging.h>
#include <util/exception.h>

namespace qlc3d {
  namespace {
    void applyCurrentOrderParameterOverride(std::vector<OrientationSample> &samples, double s0,
                                           Simu::LoadInitialOrientationS0Mode mode) {
      if (mode != Simu::LoadInitialOrientationS0Mode::Current) {
        return;
      }

      for (auto &sample : samples) {
        const auto director = sample.tensor.toDirector();
        sample.tensor = TTensor::fromDirector(qlc3d::Director(director.nx(), director.ny(), director.nz(), s0));
      }
    }
  }

  void loadInitialOrientation(const Simu &simu, double s0, const Coordinates &meshCoordinates, SolutionVector &q) {
    const std::string &loadQ = simu.getLoadQ();
    const std::string &loadOrientation = simu.getLoadOrientation();

    if (!loadQ.empty() && !loadOrientation.empty()) {
      RUNTIME_ERROR("Both loadQ and loadOrientation are set in the settings file. "
                    "loadQ is deprecated - remove it and use only loadOrientation.");
    }

    std::string file;
    if (!loadQ.empty()) {
      Log::warn("loadQ is deprecated and will be removed in a future release. Use loadOrientation instead.");
      file = loadQ;
    } else if (!loadOrientation.empty()) {
      file = loadOrientation;
    } else {
      return; // nothing to load
    }

    if (std::filesystem::path f(file); !std::filesystem::exists(f)) {
      RUNTIME_ERROR("Orientation file " + file + " for initial LC configuration does not exist.");
    }
    Log::info("Reading initial orientation from file {}", file);

    auto reader = qlc3d::createOrientationReader(file);
    auto samples = reader->read(file, s0);
    applyCurrentOrderParameterOverride(samples, s0, simu.getLoadInitialOrientationS0Mode());
    if (reader->producesLocations()) {
      Vec3 stretch = simu.getStretchVector();
      for (auto &sample : samples) {
        const Vec3 &p = sample.location.value();
        sample.location = Vec3(p.x() * stretch.x(), p.y() * stretch.y(), p.z() * stretch.z());
      }
      Log::info("Assigning {} orientation samples to {} mesh nodes using nearest-neighbor assignment.",
                samples.size(), q.getnDoF());
      qlc3d::NearestNeighborAssignment{}.assign(samples, meshCoordinates, q);
    } else {
      Log::info("Assigning {} orientation samples to {} mesh nodes using exact-order assignment.",
                samples.size(), q.getnDoF());
      qlc3d::ExactOrderAssignment{}.assign(samples, meshCoordinates, q);
    }
  }
}
