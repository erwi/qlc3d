#include <algorithm>
#include <cctype>

#include <io/orientation-loader.h>
#include <io/orientation-reader.h>
#include <io/orientation-assignment.h>
#include <io/yarf.h>
#include <solutionvector.h>
#include <geom/vec3.h>
#include <simu.h>
#include <util/logging.h>
#include <util/exception.h>

namespace qlc3d {
  namespace {
    /**
     * Normalises the extension to lowercase so YARF detection remains case-insensitive.
     * @param path Candidate file path.
     * @return Lowercase file extension including the leading dot.
     */
    std::string lowercaseExtension(const std::filesystem::path &path) {
      std::string extension = path.extension().string();
      std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
      });
      return extension;
    }

    /**
     * Replaces the scalar order in a loaded sample set with the active material S0 when the current
     * mode is selected while preserving the sample director orientation.
     * @param samples Samples to rewrite.
     * @param s0 Requested scalar order.
     * @param mode Active S0 override mode.
     */
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

    /**
     * Replaces the scalar order in a loaded SolutionVector with the active material S0 when the current
     * mode is selected while preserving each node director orientation.
     * @param q Solution to rewrite.
     * @param s0 Requested scalar order.
     * @param mode Active S0 override mode.
     */
    void applyCurrentOrderParameterOverride(SolutionVector &q, double s0,
                                           Simu::LoadInitialOrientationS0Mode mode) {
      if (mode != Simu::LoadInitialOrientationS0Mode::Current) {
        return;
      }

      for (idx node = 0; node < q.getnDoF(); ++node) {
        const auto director = q.getDirector(node);
        q.setValue(node, TTensor::fromDirector(qlc3d::Director(director.nx(), director.ny(), director.nz(), s0)));
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

    const auto lowerExtension = lowercaseExtension(std::filesystem::path(file));

    if (lowerExtension == ".yarf") {
      const auto loaded = qlc3d::yarf::YarfFileReader{}.readQTensorSolution(std::filesystem::path(file));
      if (loaded.getnDoF() != q.getnDoF()) {
        RUNTIME_ERROR("YARF q_tensor node count (" + std::to_string(loaded.getnDoF()) +
                      ") does not match the target SolutionVector size (" + std::to_string(q.getnDoF()) + ").");
      }
      if (loaded.getnDimensions() != q.getnDimensions()) {
        RUNTIME_ERROR("YARF q_tensor dimension count (" + std::to_string(loaded.getnDimensions()) +
                      ") does not match the target SolutionVector dimensions (" + std::to_string(q.getnDimensions()) + ").");
      }
      q = loaded;
      applyCurrentOrderParameterOverride(q, s0, simu.getLoadInitialOrientationS0Mode());
      return;
    }

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
