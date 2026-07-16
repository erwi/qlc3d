#ifndef PROJECT_QLC3D_ORIENTATION_READER_H
#define PROJECT_QLC3D_ORIENTATION_READER_H
#include <vector>
#include <string>
#include <memory>
#include "orientation-sample.h"

namespace qlc3d {
  /** Parses a file on disk into a list of OrientationSample. Implementations must not touch
   *  SolutionVector/mesh directly, so they can be unit tested with plain files. */
  class OrientationReader {
  public:
    virtual ~OrientationReader() = default;
    /** @param s0 default order parameter to use for formats/rows that don't specify one explicitly
     *           (ignored by formats that always carry an explicit S, e.g. LCView). */
    [[nodiscard]] virtual std::vector<OrientationSample> read(const std::string &fileName, double s0) const = 0;
    /** true if this format's samples always carry a location (point-cloud format);
     *  false if samples never carry a location (mesh-matched format). Never mixed within one format. */
    [[nodiscard]] virtual bool producesLocations() const = 0;
  };

  /** Reads the legacy LCView text result file format (no coordinates, mesh-matched). */
  class LcViewTextReader : public OrientationReader {
  public:
    [[nodiscard]] std::vector<OrientationSample> read(const std::string &fileName, double s0) const override;
    [[nodiscard]] bool producesLocations() const override { return false; }
  };

  /** Reads the legacy LCView binary result file format (no coordinates, mesh-matched). */
  class LcViewBinaryReader : public OrientationReader {
  public:
    [[nodiscard]] std::vector<OrientationSample> read(const std::string &fileName, double s0) const override;
    [[nodiscard]] bool producesLocations() const override { return false; }
  };

  /** Sniffs fileName (content, e.g. "RAW FLOAT TRI" marker - same as today's ResultIO::ReadResult) to decide
   *  between LcViewTextReader and LcViewBinaryReader. */
  [[nodiscard]] std::unique_ptr<OrientationReader> createLcViewReader(const std::string &fileName);

  /** Decides which OrientationReader to use for fileName: `.csv` extension (case-insensitive) routes to
   *  DirectorCsvReader; anything else falls back to the existing LCView text/binary content-sniffing
   *  (createLcViewReader). */
  [[nodiscard]] std::unique_ptr<OrientationReader> createOrientationReader(const std::string &fileName);
}

#endif //PROJECT_QLC3D_ORIENTATION_READER_H
