#ifndef PROJECT_QLC3D_ORIENTATION_SAMPLE_H
#define PROJECT_QLC3D_ORIENTATION_SAMPLE_H
#include <optional>
#include <lc-representation.h>
#include <geom/vec3.h>

namespace qlc3d {
  /** One parsed orientation value: a Q-tensor (via TTensor) plus an optional 3D location.
   *  location is nullopt for "mesh-matched" formats (no coordinates in file, e.g. LCView),
   *  and set for "point-cloud" formats (e.g. director CSV). Whether a format produces samples
   *  with or without locations is intrinsic to that format (see OrientationReader::producesLocations()). */
  struct OrientationSample {
    TTensor tensor;
    std::optional<Vec3> location;
  };
}

#endif //PROJECT_QLC3D_ORIENTATION_SAMPLE_H
