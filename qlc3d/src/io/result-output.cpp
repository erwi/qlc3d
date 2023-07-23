#include <io/result-output.h>
#include <io/lcview-result-output.h>
#include <geometry.h>
#include <solutionvector.h>
#include <filesystem>
#include <simulation-state.h>
#include <qlc3d.h>
#include <resultio.h>
#include <simu.h>
#include <util/logging.h>
#include "util/exception.h"

namespace fs = std::filesystem;

//<editor-fold desc="ResultOutput">
ResultOutput::ResultOutput(const std::set<Simu::SaveFormats> &saveFormats,
                           const std::string &meshName,
                           double S0,
                           const std::filesystem::path &outputDir) {
  outputFormatWriters_ = std::vector<std::shared_ptr<ResultFormatWriter>>();

  for (const auto &s : saveFormats) {
    if (s == Simu::SaveFormats::LCview) {
      outputFormatWriters_.push_back(std::shared_ptr<ResultFormatWriter>(new LcViewBinaryResultFormatWriter(outputDir, meshName, S0)));
    }
    else if (s == Simu::SaveFormats::LCviewTXT) {
      outputFormatWriters_.push_back(std::shared_ptr<ResultFormatWriter>(new LcViewTxtResultFormatWriter(outputDir, meshName, S0)));
    }
    else if (s == Simu::SaveFormats::CsvUnstructured) {
      outputFormatWriters_.push_back(std::shared_ptr<ResultFormatWriter>(new CsvUnstructuredFormatWriter(outputDir)));
    }
    else if (s == Simu::SaveFormats::DirStackZ) {
      outputFormatWriters_.push_back(std::shared_ptr<ResultFormatWriter>(new DirStackZFormatWriter(outputDir)));
    }
    else if (s == Simu::SaveFormats::RegularVTK) {
      outputFormatWriters_.push_back(std::shared_ptr<ResultFormatWriter>(new RegularVTKFormatWriter(outputDir)));
    }
    else if (s == Simu::SaveFormats::RegularVecMat) {
      outputFormatWriters_.push_back(std::shared_ptr<ResultFormatWriter>(new RegularVecMatFormatWriter(outputDir)));
    }
    else if (s == Simu::SaveFormats::VTKUnstructuredAsciiGrid) {
      outputFormatWriters_.push_back(std::shared_ptr<ResultFormatWriter>(new VtkUnstructuredAsciiGridFormatWriter(outputDir)));
    }
    else {
      RUNTIME_ERROR("Unknown save format: " + std::to_string(s));
    }
  }
}

void ResultOutput::writeResults(const Geometry &geom,
                                const SolutionVector &potential,
                                const SolutionVector &qtensor,
                                const SimulationState &simulationState) {
  // if any of current output format writers requires director, calculate director
  double *director = nullptr;
  if (isDirectorRequired()) {
    director = tensortovector(qtensor.Values, geom.getnpLC());
  }

  for (auto outputFormatWriter : outputFormatWriters_) {
    if (outputFormatWriter->isDirectorRequired()) {
      outputFormatWriter->setDirector(director);
    }
    outputFormatWriter->setQTensor(qtensor);
    outputFormatWriter->setPotential(potential);
    Log::info("Writing result with format {}", outputFormatWriter->formatName());
    outputFormatWriter->writeResult(geom, simulationState);
  }

  if (director != nullptr) {
    delete[] director;
  }
}

bool ResultOutput::isDirectorRequired() const {
  for (auto outputFormatWriter : outputFormatWriters_) {
    if (outputFormatWriter->isDirectorRequired()) {
      return true;
    }
  }
  return false;
}

bool ResultOutput::isRegularGridRequired() const {
  for (auto outputFormatWriter : outputFormatWriters_) {
    if (outputFormatWriter->isRegularGridRequired()) {
      return true;
    }
  }
  return false;
}
//</editor-fold>

std::string ResultFormatWriter::iterationAsString(const SimulationState &simulationState) {
  char numberChar[9];
  sprintf(numberChar, "%08d", simulationState.currentIteration());
  string fileNameNumber(numberChar);
  if (simulationState.state() == RunningState::COMPLETED) { // after completion, output filename with special counter value "final"
    fileNameNumber = "-final";
  }
  return fileNameNumber;
}

//<editor-fold desc="RegularVtkFormatWriter">
void RegularVTKFormatWriter::writeResult(const Geometry &geom, const SimulationState &simulationState)  {
  std::string filename = "regularvtk" + iterationAsString(simulationState) + ".vtk";
  std::filesystem::path filePath = outputDirectory / filename;

  RegularGrid* rGrid = geom.regularGrid;
  if (rGrid == nullptr) {
    RUNTIME_ERROR("RegularVTKFormatWriter requires regular grid");
  }

  rGrid->writeVTKGrid(filePath.c_str(),
                     potential->Values,
                     director,
                     geom.getnpLC());
}
//</editor-fold>

//<editor-fold desc="RegularVecMatFormatWriter">
void RegularVecMatFormatWriter::writeResult(const Geometry &geom, const SimulationState &simulationState)  {
  fs::path filePath = outputDirectory / ("regularvec" + iterationAsString(simulationState) + ".m");

  RegularGrid *rGrid = geom.regularGrid;
  if (rGrid == nullptr) {
    RUNTIME_ERROR("RegularVecMatFormatWriter requires regular grid");
  }
  rGrid->writeVecMat(filePath.c_str(),       // WRITE REGULAR GRID RESULT FILE
                    potential->Values,
                    director,
                    geom.getnpLC(),
                    simulationState.currentTime());
}
//</editor-fold>

//<editor-fold desc="DirStackZFormatWriter">
void DirStackZFormatWriter::writeResult(const Geometry &geom, const SimulationState &simulationState) {
  std::string filename = "dirstacksz" + iterationAsString(simulationState) + ".csv";

  RegularGrid *rGrid = geom.regularGrid;
  if (rGrid == nullptr) {
    RUNTIME_ERROR("DirStackZFormatWriter requires regular grid");
  }
  fs::path filePath = outputDirectory / filename;
  rGrid->writeDirStackZ(filePath.c_str(),
                       director,
                       geom.getnpLC(),
                       simulationState.currentTime());
}
//</editor-fold>

//<editor-fold desc="CsvUnstructuredFormatWriter">
void CsvUnstructuredFormatWriter::writeResult(const Geometry &geom, const SimulationState &simulationState) {
  fs::path filePath = outputDirectory / ("unstructured.csv." + std::to_string(simulationState.currentIteration()));
  ResultIO::writeCsvUnstructured(geom.getPtrTop(), *potential, *qTensor, filePath.string());
}
//</editor-fold>

//<editor-fold desc="VtkUnstructuredAsciiGridFormatWriter">
void VtkUnstructuredAsciiGridFormatWriter::writeResult(const Geometry &geom, const SimulationState &simulationState) {
  std::string fileName = "unstructured" + iterationAsString(simulationState) + ".vtk";
  std::filesystem::path filePath = outputDirectory / fileName;
  ResultIO::writeVtkUnstructuredAsciiGrid(
          geom.getPtrTop(), geom.getnp(), geom.getnpLC(), geom.getTetrahedra(), *potential, *qTensor, filePath.string());
}
//</editor-fold>
