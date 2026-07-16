#ifndef PROJECT_QLC3D_DIRECTOR_CSV_READER_H
#define PROJECT_QLC3D_DIRECTOR_CSV_READER_H
#include "orientation-reader.h"

namespace qlc3d {
  /** Reads a director CSV point cloud: header row (any order, case-insensitive) naming comma-separated
   *  columns x, y, z, nx, ny, nz, and an optional S column. Director vectors need not be pre-normalized
   *  (this reader normalizes them). If S is omitted, every sample's S is set to the s0 parameter passed to
   *  read(). Locations are returned in the file's own (unstretched) coordinate space -- StretchVector scaling
   *  is applied later by the caller (see NearestNeighborAssignment contract), not by this reader. */
  class DirectorCsvReader : public OrientationReader {
  public:
    [[nodiscard]] std::vector<OrientationSample> read(const std::string &fileName, double s0) const override;
    [[nodiscard]] bool producesLocations() const override { return true; }
  };
}

#endif //PROJECT_QLC3D_DIRECTOR_CSV_READER_H
