#ifndef PROJECT_QLC3D_ORIENTATION_LOADER_H
#define PROJECT_QLC3D_ORIENTATION_LOADER_H
class Simu;
class SolutionVector;
class Coordinates;

namespace qlc3d {
  /** Top-level entry point replacing the old ResultIO::ReadResult call site. Resolves loadQ/loadOrientation
   *  precedence (RUNTIME_ERROR if both are set), logs the loadQ deprecation warning if used, dispatches to the
   *  correct OrientationReader based on file content/extension, and assigns the parsed samples onto q's LC nodes
   *  using the assignment strategy appropriate for that format: exact order for mesh-matched formats (no
   *  locations, e.g. LCView), or nearest-neighbor (with sample locations scaled by simu's StretchVector before
   *  matching) for point-cloud formats (e.g. director CSV). Does nothing if neither loadQ nor loadOrientation is set.
   *  @param simu holds the (mutually exclusive) loadQ/loadOrientation settings.
   *  @param s0 default order parameter to use for formats/rows that don't specify one explicitly.
   *  @param meshCoordinates coordinates of the (already stretched) mesh nodes, used by the assignment strategy.
   *  @param q the solution vector whose LC nodes are populated.
   */
  void loadInitialOrientation(const Simu &simu, double s0, const Coordinates &meshCoordinates, SolutionVector &q);
}

#endif //PROJECT_QLC3D_ORIENTATION_LOADER_H
