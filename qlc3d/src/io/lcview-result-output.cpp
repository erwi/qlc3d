#include <io/lcview-result-output.h>
#include <io/meshwriter.h>
#include <util/logging.h>
#include <simulation-state.h>
#include <geometry.h>
#include <resultio.h>
#include <util/exception.h>
#include <lc-representation.h>

namespace fs = std::filesystem;

using namespace std;

//<editor-fold desc="LcViewResultFormatWriter">
LcViewResultFormatWriter::LcViewResultFormatWriter(const std::filesystem::path &outputDir, const std::string &meshName, double S0) :
  ResultFormatWriter(outputDir), S0_{ S0 }, lastMeshNumber_{ -1 } {

  // make sure mesh name does not contain any prefix path before file name
  meshName_ = fs::path(meshName).filename().string();
}

std::string LcViewResultFormatWriter::numberedMeshName(const SimulationState &simulationState, const std::string &meshName) {
  // calculate mesh name with mesh number appended. e.g. mesh.txt -> mesh0.txt
  std::string numberedMeshName = meshName;
  std::stringstream ss;
  ss << simulationState.meshNumber();
  std::string num;
  ss >> num;
  size_t pos = numberedMeshName.find_last_of("."); // position of separator point
  numberedMeshName.insert(pos, num);
  return numberedMeshName;
}

void LcViewResultFormatWriter::writeMeshIfRequired(const Geometry &geom, const SimulationState &simulationState) {
  std::string numberedMeshName = this->numberedMeshName(simulationState, meshName_);
  bool isFirstTime = lastMeshNumber_ == -1;
  bool isNewMesh = simulationState.meshNumber() > (unsigned int) lastMeshNumber_;
  if (isFirstTime || isNewMesh) {
    writtenMeshPath_ = outputDirectory / numberedMeshName;
    Log::info("Writing mesh file {}", writtenMeshPath_.string());
    writeMeshFile(geom.getCoordinates(), geom.getTetrahedra(), geom.getTriangles(), writtenMeshPath_);
    lastMeshNumber_ = simulationState.meshNumber();
  }
}


//</editor-fold>

//<editor-fold desc=LcViewBinaryFormatWriter>
void LcViewBinaryResultFormatWriter::writeResult(const Geometry &geom, const SimulationState &simulationState) {
  writeMeshIfRequired(geom, simulationState);

  fs::path outputFilePath = outputDirectory;
  if (simulationState.state() == RunningState::RUNNING || simulationState.state() == RunningState::INITIALISING) {
    std::stringstream ss;
    ss << std::setw(5) << std::setfill('0') << simulationState.currentIteration();
    outputFilePath = outputDirectory / ("result" + ss.str() + ".dat");
  } else if (simulationState.state() == RunningState::COMPLETED) {
    outputFilePath = outputDirectory / "result-final.dat";
  } else {
    RUNTIME_ERROR("Simulation state is not RUNNING or COMPLETED");
  }
  Log::info("Writing binary LcView result file {}", outputFilePath.string());
  writeBinaryResultFile(potential,
                        qTensor,
                        simulationState.currentTime().getTime(),
                        S0_,
                        writtenMeshPath_.filename().string(),
                        outputFilePath);
}

void LcViewBinaryResultFormatWriter::writeBinaryResultFile(const SolutionVector *v,
                                                           const SolutionVector *q,
                                                           double currentTime,
                                                           double S0,
                                                           const std::string &meshFileName,
                                                           const std::filesystem::path &filePath) {
  const int npLC = q->getnDoF();
  const int np = v->getnDoF();
  FILE *fid = fopen(filePath.string().c_str(), "wb");
  char time[20];
  sprintf(time, "%1.9f\n", currentTime);
  string text = "** Result Time :   ";
  text.append(time);
  fprintf(fid, "%s\n", text.c_str());
  fprintf(fid, "** z Compression Ratio :  1.00000\n");
  fprintf(fid, "%s\n", meshFileName.c_str());
  fprintf(fid, "RAW FLOAT TRI - S0, np, nsols\n");
  fprintf(fid, "%g %d %d\r\n", S0, np, 6);
  for (int i = 0;  i < np; i++) {
    if (i < npLC) {
      // WRITE LC REGIONS
      float value;
      value = (float)q->getValue(i, 0);
      fwrite((void *) &value, sizeof(float), 1, fid);
      value = (float)q->getValue(i, 1);
      fwrite((void *) &value, sizeof(float), 1, fid);
      value = (float)q->getValue(i, 2);
      fwrite((void *) &value, sizeof(float), 1, fid);
      value = (float)q->getValue(i, 4);
      fwrite((void *) &value, sizeof(float), 1, fid);
      value = (float)q->getValue(i, 3);
      fwrite((void *) &value, sizeof(float), 1, fid);
      value = v->getValue(i);
      fwrite((void *) &value , sizeof(float), 1 , fid);
    } else {
      // WRITE DIELECTRIC REGIONS
      float value = 0;
      for (int qc = 0; qc < 5; qc++) // q-component
        fwrite((void *) &value, sizeof(float), 1, fid);
      value = v->getValue(i);
      fwrite((void *) &value, sizeof(float), 1, fid);
    }
  }
  fclose(fid);
}
//</editor-fold>

//<editor-fold desc=LcViewTxtFormatWriter>
void LcViewTxtResultFormatWriter::writeResult(const Geometry &geom, const SimulationState &simulationState) {
  writeMeshIfRequired(geom, simulationState);

  fs::path filePath;
  if (simulationState.state() == RunningState::RUNNING || simulationState.state() == RunningState::INITIALISING) {
    std::stringstream ss;
    ss << std::setw(5) << std::setfill('0') << simulationState.currentIteration();
    filePath = outputDirectory / ("result-t-" + ss.str() + ".dat");
  } else if (simulationState.state() == RunningState::COMPLETED) {
    filePath = outputDirectory / "result-t-final.dat";
  } else {
    RUNTIME_ERROR("Simulation state is not RUNNING or COMPLETED");
  }

  writeTextResultFile(*potential,
                      *directors,
                      simulationState.currentTime().getTime(),
                      writtenMeshPath_.filename().string(),
                      filePath);
}

void LcViewTxtResultFormatWriter::writeTextResultFile(//const Coordinates& coordinates,
                                                      //const Mesh *t,
                                                      //const Mesh *e,
                                                      const SolutionVector &v,
                                                      const std::vector<qlc3d::Director> &dir,
                                                      double currentTime,
                                                      const std::string &meshFileName,
                                                      const std::filesystem::path &filePath) {

  FILE *fid = fopen(filePath.string().c_str(), "w");
  if (fid != nullptr) {
    const idx np = v.getnDoF();
    const idx npLC = dir.size();
    char str[15];
    sprintf(str, "%f", currentTime);
    std::string text = "** Result Time :    ";
    text.append(str);
    text.append("\n** z Compression Ratio :  1.00000\n");
    text.append(meshFileName + "\n");
    fprintf(fid, "%s", text.c_str()); //** Result Time :    0.00000000\n** z Compression Ratio :  1.00000\nmeshout.txt\n");
    for (unsigned int i = 0; i < np; i++) {
      if (i < npLC) {
        const qlc3d::Director &n = dir[i];
        fprintf(fid, LCVIEW_TEXT_FORMAT_STRING, i + 1, n.nx(), n.ny(), n.nz(), v.getValue(i), n.S(), n.S()); // NOTE: same S value is written twice. Should be the two different order parameters, or biaxiality P?
      } else {
        fprintf(fid, LCVIEW_TEXT_FORMAT_STRING, i + 1, 0., 0., 0., v.getValue(i), 0., 0.);
      }
    }
    fclose(fid);
  } else {
    RUNTIME_ERROR("Could not open result file: " + filePath.string());
  }
}
//</editor-fold>
