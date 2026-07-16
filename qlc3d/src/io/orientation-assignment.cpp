#include <io/orientation-assignment.h>
#include <solutionvector.h>
#include <geom/coordinates.h>
#include <util/exception.h>
#include <globals.h>
#include <fmt/format.h>
#include <limits>

namespace qlc3d {
  void ExactOrderAssignment::assign(const std::vector<OrientationSample> &samples,
                                     const Coordinates &meshCoordinates, SolutionVector &q) const {
    if (samples.size() != q.getnDoF()) {
      RUNTIME_ERROR(fmt::format("The loaded result file size {} does not match the expected size {}",
                                 samples.size(), q.getnDoF()));
    }
    for (idx i = 0; i < q.getnDoF(); i++) {
      q.setValue(i, samples[i].tensor);
    }
  }

  void NearestNeighborAssignment::assign(const std::vector<OrientationSample> &samples,
                                          const Coordinates &meshCoordinates, SolutionVector &q) const {
    if (samples.empty()) {
      RUNTIME_ERROR("No orientation samples available for nearest-neighbor assignment");
    }
    for (const auto &sample : samples) {
      if (!sample.location.has_value()) {
        RUNTIME_ERROR("NearestNeighborAssignment requires all samples to have a location");
      }
    }

    for (idx i = 0; i < q.getnDoF(); i++) {
      const Vec3 &nodePoint = meshCoordinates.getPoint(i);
      double bestDistanceSquared = std::numeric_limits<double>::max();
      const OrientationSample *closest = nullptr;

      for (const auto &sample : samples) {
        double d = nodePoint.distanceSquared(*sample.location);
        if (d < bestDistanceSquared) {
          bestDistanceSquared = d;
          closest = &sample;
        }
      }
      q.setValue(i, closest->tensor);
    }
  }
}
