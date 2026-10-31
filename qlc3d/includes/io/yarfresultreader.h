
#ifndef PROJECT_QLC3D_YARFRESULTREADER_H
#define PROJECT_QLC3D_YARFRESULTREADER_H
#include <io/yarfresultwriter.h>
#include <io/orientation-reader.h>


/**
 * Reads mesh-matched q_tensor data from a YARF result snapshot. Ignores any other fields
 * The yarf file spec is given in #file: qlc3d/doc-impl/yarf-result-file-format.md
 */
class YarfResultReader : public qlc3d::OrientationReader {
public:
    /** Reads a YARF result file and returns the q_tensor samples as a vector of OrientationSample.
     *  @param fileName Path to the YARF result file.
     *  @param s0 Default order parameter to use for formats/rows that don't specify one explicitly
     *            (ignored by formats that always carry an explicit S, e.g. LCView).
     *  @return Vector of OrientationSample containing the q_tensor data from the YARF file.
     */
    [[nodiscard]] std::vector<qlc3d::OrientationSample> read(const std::string &fileName, double s0) const override;
    [[nodiscard]] bool producesLocations() const override { return false; }

    static Yarf::YarfHeader readYarfHeader(const std::string &fileName) ;
};
#endif //PROJECT_QLC3D_YARFRESULTREADER_H
