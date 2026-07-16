#ifndef PROJECT_QLC3D_ORIENTATION_ASSIGNMENT_H
#define PROJECT_QLC3D_ORIENTATION_ASSIGNMENT_H
#include <vector>
#include "orientation-sample.h"
class SolutionVector;
class Coordinates;

namespace qlc3d {
  /** Writes parsed OrientationSample values onto the LC nodes of a SolutionVector. */
  class OrientationAssignmentStrategy {
  public:
    virtual ~OrientationAssignmentStrategy() = default;
    virtual void assign(const std::vector<OrientationSample> &samples,
                         const Coordinates &meshCoordinates,
                         SolutionVector &q) const = 0;
  };

  /** For mesh-matched formats: samples carry no location and must be assigned in file order,
   *  1:1 onto q's LC nodes. Throws a runtime_error (RUNTIME_ERROR) on any count mismatch,
   *  matching legacy behavior/wording. */
  class ExactOrderAssignment : public OrientationAssignmentStrategy {
  public:
    void assign(const std::vector<OrientationSample> &samples,
                const Coordinates &meshCoordinates, SolutionVector &q) const override;
  };

  /** For point-cloud formats: samples carry locations that don't need to coincide with mesh nodes.
   *  For each LC mesh node, assigns the tensor value of the closest sample by Euclidean distance
   *  (brute-force search — acceptable for expected point-cloud sizes; the interface allows a smarter
   *  spatial index to replace this later without touching callers).
   *  IMPORTANT: sample locations must already be in the same (stretched) coordinate space as
   *  meshCoordinates before calling assign() -- this class does not itself apply StretchVector scaling;
   *  that is the caller's (orientation-loader's) responsibility, done once for the whole point cloud. */
  class NearestNeighborAssignment : public OrientationAssignmentStrategy {
  public:
    /** On exact ties in distance, the first sample encountered in list order wins (deterministic, not
     *  undefined). */
    void assign(const std::vector<OrientationSample> &samples,
                const Coordinates &meshCoordinates, SolutionVector &q) const override;
  };
}

#endif //PROJECT_QLC3D_ORIENTATION_ASSIGNMENT_H
