#ifndef PROJECT_QLC3D_RESULT_OUTPUT_H
#define PROJECT_QLC3D_RESULT_OUTPUT_H

#include <string>
#include <vector>
#include <set>
#include <simu.h>
#include <memory>
#include <filesystem>

// forward declarations
class SolutionVector;
class Geometry;
class SimulationState;

/**
 * Abstract base class for result output formats. Each output format should implement this interface.
 */
class ResultFormatWriter {
protected:
  const double *director = nullptr;
  const SolutionVector *qTensor = nullptr;
  const SolutionVector *potential = nullptr;
  const std::filesystem::path outputDirectory;

  std::string static iterationAsString(const SimulationState &simulationState);

public:

  explicit ResultFormatWriter(const std::filesystem::path &outputDirectory) : outputDirectory(outputDirectory) {};
  [[nodiscard]] virtual bool requiresDirector() const = 0;

  void setDirector(const double *dir) { this->director = dir; };
  void setQTensor(const SolutionVector &q) {this->qTensor = &q; };
  void setPotential(const SolutionVector &pot) {this->potential = &pot; }

  /** write the result to the output format. TODO: Geometry should be const */
  virtual void writeResult(const Geometry &geom, const SimulationState &simulationState) = 0;
  virtual const std::string formatName() const = 0;
  virtual ~ResultFormatWriter() = default;
};

class ResultOutput {
public:
  ResultOutput(const std::set<Simu::SaveFormats> &saveFormats, const std::string &meshName, double S0, const std::filesystem::path &outputDir);
  void writeResults(const Geometry &geom, const SolutionVector &potential, const SolutionVector &qtensor, const SimulationState &simulationState);

private:
  std::vector<std::shared_ptr<ResultFormatWriter>> outputFormatWriters_;
  /** output directory, relative to the application working directory */
  std::string outputDirectory_;

  [[nodiscard]] bool isDirectorRequired() const;
};

class RegularVTKFormatWriter : public ResultFormatWriter {
public:
  explicit RegularVTKFormatWriter(const std::filesystem::path &outputDir) : ResultFormatWriter(outputDir) {};
  [[nodiscard]] bool requiresDirector() const override { return true; };
  [[nodiscard]] const std::string formatName() const override { return "RegularVTK"; };
  void writeResult(const Geometry &geom, const SimulationState &simulationState) override;
  ~RegularVTKFormatWriter() override = default;
};

class RegularVecMatFormatWriter : public ResultFormatWriter {
public:
  explicit RegularVecMatFormatWriter(const std::filesystem::path &outputDir) : ResultFormatWriter(outputDir) {};
  [[nodiscard]] bool requiresDirector() const override { return true; };
  [[nodiscard]] const std::string formatName() const override { return "RegularVecMat"; };
  void writeResult(const Geometry &geom, const SimulationState &simulationState) override;
  ~RegularVecMatFormatWriter() override = default;
};

class DirStackZFormatWriter : public ResultFormatWriter {
public:
  explicit DirStackZFormatWriter(const std::filesystem::path &outputDir) : ResultFormatWriter(outputDir) {};
  [[nodiscard]] bool requiresDirector() const override { return true; };
  [[nodiscard]] const std::string formatName() const override { return "DirStackZ"; };
  void writeResult(const Geometry &geom, const SimulationState &simulationState) override;
  ~DirStackZFormatWriter() override = default;
};

class CsvUnstructuredFormatWriter : public ResultFormatWriter {
public:
  explicit CsvUnstructuredFormatWriter(const std::filesystem::path &outputDir) : ResultFormatWriter(outputDir) {};
  [[nodiscard]] bool requiresDirector() const override { return false; };
  [[nodiscard]] const std::string formatName() const override { return "CsvUnstructured"; };
  void writeResult(const Geometry &geom, const SimulationState &simulationState) override;
  ~CsvUnstructuredFormatWriter() override = default;
};

class VtkUnstructuredAsciiGridFormatWriter : public ResultFormatWriter {
public:
  explicit VtkUnstructuredAsciiGridFormatWriter(const std::filesystem::path &outputDir) : ResultFormatWriter(outputDir) {};
  [[nodiscard]] bool requiresDirector() const override { return false; }
  [[nodiscard]] const std::string formatName() const override { return "VtkUnstructuredAsciiGrid"; };
  void writeResult(const Geometry &geom, const SimulationState &simulationState) override;
  ~VtkUnstructuredAsciiGridFormatWriter() override = default;
};

#endif //PROJECT_QLC3D_RESULT_OUTPUT_H
