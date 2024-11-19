#ifndef PROJECT_QLC3D_LCVIEW_RESULT_OUTPUT_H
#define PROJECT_QLC3D_LCVIEW_RESULT_OUTPUT_H

#include <io/result-output.h>

class Mesh;
class Coordinates;

class LcViewResultFormatWriter : public ResultFormatWriter {
  [[nodiscard]] static std::string numberedMeshName(const SimulationState &simulationState, const std::string &meshName) ;

  static void writeMeshFile(const Coordinates &coordinates,
                            const Mesh &t,
                            const Mesh &e,
                            idx np,
                            const std::filesystem::path &fileName);
protected:
  double S0_;
  std::string meshName_;
  int lastMeshNumber_;
  /** path to the last written mesh file. Should be updated whenever a new mesh file is output */
  std::filesystem::path writtenMeshPath_;

  LcViewResultFormatWriter(const std::filesystem::path &outputDir, const std::string &meshName, double S0);// :
          //ResultFormatWriter(outputDir), meshName_(meshName), lastMeshNumber_(-1), S0_(S0) {};

  /** writes mesh file if required and sets writtenMeshPath_ */
  void writeMeshIfRequired(const Geometry &geom, const SimulationState &simulationState);
};

//=============================================
class LcViewBinaryResultFormatWriter : public LcViewResultFormatWriter {
  static void writeBinaryResultFile(const SolutionVector *v,
                             const SolutionVector *q,
                             double currentTime,
                             double S0,
                             const std::string &meshFileName,
                             const std::filesystem::path &filePath);

public:
  LcViewBinaryResultFormatWriter(const std::filesystem::path &outputDir, const std::string &meshName, double S0) :
          LcViewResultFormatWriter(outputDir, meshName, S0) {};

  [[nodiscard]] const std::string formatName() const override { return "LcView"; };
  void writeResult(const Geometry &geom, const SimulationState &simulationState) override;
};
//=============================================
class LcViewTxtResultFormatWriter : public LcViewResultFormatWriter {
  static constexpr char LCVIEW_TEXT_FORMAT_STRING[] = "%i %f %f %f %f %f %f\n";
  static void writeTextResultFile(const Coordinates &coordinates,
                           const Mesh *t,
                           const Mesh *e,
                           const SolutionVector &v,
                           const std::vector<qlc3d::Director> &dir,
                           double currentTime,
                           const std::string &meshFileName,
                           const std::filesystem::path &filePath);
public:
  LcViewTxtResultFormatWriter(const std::filesystem::path &outputDir, const std::string &meshName, double S0) :
          LcViewResultFormatWriter(outputDir, meshName, S0) {};

  [[nodiscard]] bool isDirectorRequired() const override { return true; };
  [[nodiscard]] const std::string formatName() const override { return "LcViewTxt"; };
  void writeResult(const Geometry &geom, const SimulationState &simulationState) override;
};

#endif //PROJECT_QLC3D_LCVIEW_RESULT_OUTPUT_H
