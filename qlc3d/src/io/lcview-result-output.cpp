#include <io/lcview-result-output.h>
#include <util/logging.h>
#include <simulation-state.h>
#include <geometry.h>
#include <resultio.h>
#include "util/exception.h"

namespace fs = std::filesystem;

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
  if (simulationState.meshNumber() > lastMeshNumber_|| lastMeshNumber_ == -1) {
    writtenMeshPath_ = outputDirectory / numberedMeshName;
    Log::info("Writing mesh file {}", writtenMeshPath_.string());
    //ResultIO::writeMesh(geom.getPtrTop(), geom.t, geom.e, geom.getnp(), writtenMeshPath_);
    writeMeshFile(geom.getPtrTop(), geom.t, geom.e, geom.getnp(), writtenMeshPath_);
    lastMeshNumber_ = simulationState.meshNumber();
  }
}

void LcViewResultFormatWriter::writeMeshFile(const double *p, Mesh *t, Mesh *e, idx np,
                                             const std::filesystem::path &fileName) {
  idx i;
  FILE *fid = fopen(fileName.string().c_str(), "w");
  if (fid != NULL) {
    fputs("MESH    dimension 3 ElemType Tetrahedra  Nnode 4\nCoordinates\n", fid);
    for (i = 0; i < np; i++) {
      fprintf(fid, "%i\t%f\t%f\t%f\n", i + 1, p[3 * i], p[3 * i + 1], p[3 * i + 2]);
    }
    fprintf(fid, "end coordinates\n\nElements\n");
    for (i = 0 ; i < t->getnElements() ; i++)
      fprintf(fid, "%i\t%i\t%i\t%i\t%i\t%i\n", i + 1, t->getNode(i, 0) + 1, t->getNode(i, 1) + 1, t->getNode(i, 2) + 1, t->getNode(i, 3) + 1, t->getMaterialNumber(i));
    fprintf(fid, "end elements\nMESH    dimension 3 ElemType Triangle  Nnode 3\nCoordinates\nend coordinates\n\nElements\n");
    for (i = 0 ; i < e->getnElements() ; i++)
      fprintf(fid, "%i\t%i\t%i\t%i\t%i\n", i + 1, e->getNode(i, 0) + 1, e->getNode(i, 1) + 1, e->getNode(i, 2) + 1, e->getMaterialNumber(i));
    fprintf(fid, "end elements\n");
    fclose(fid);
  } else {
    RUNTIME_ERROR(fmt::format("Could not open file for output mesh: {}", fileName));
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

  writeBinaryResultFile(geom.getPtrTop(),
                        geom.t,
                        geom.e,
                        potential,
                        qTensor,
                        simulationState.currentTime(),
                        S0_,
                        writtenMeshPath_.filename().string(),
                        outputFilePath);
}

void LcViewBinaryResultFormatWriter::writeBinaryResultFile(const double *p,
                                                           const Mesh *t,
                                                           const Mesh *e,
                                                           const SolutionVector *v,
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

  writeTextResultFile(geom.getPtrTop(),
                      geom.t,
                      geom.e,
                      potential,
                      director,
                      geom.getnpLC(),
                      simulationState.currentTime(),
                      writtenMeshPath_.filename().string(),
                      filePath);
}

void LcViewTxtResultFormatWriter::writeTextResultFile(const double *p,
                                                      const Mesh *t,
                                                      const Mesh *e,
                                                      const SolutionVector *v,
                                                      const double* n,
                                                      idx npLC,
                                                      double currentTime,
                                                      const std::string &meshFileName,
                                                      const std::filesystem::path &filePath) {
  idx np = v->getnDoF();
  char str[15];
  FILE *fid = fopen(filePath.string().c_str(), "w");

  if (fid != nullptr) {
    sprintf(str, "%f", currentTime);
    std::string text = "** Result Time :    ";
    text.append(str);
    text.append("\n** z Compression Ratio :  1.00000\n");
    text.append(meshFileName + "\n");
    fprintf(fid, "%s", text.c_str()); //** Result Time :    0.00000000\n** z Compression Ratio :  1.00000\nmeshout.txt\n");
    for (int i = 0; i < np; i++) {
      if (i < npLC)
        fprintf(fid, LCVIEW_TEXT_FORMAT_STRING, i + 1, n[i], n[i + npLC], n[i + 2 * npLC], v->Values[i], n[i + 3 * npLC], n[i + 4 * npLC]);
      else
        fprintf(fid, LCVIEW_TEXT_FORMAT_STRING, i + 1, 0., 0., 0., v->Values[i], 0., 0.);
    }
    fclose(fid);

  } else {
    RUNTIME_ERROR("Could not open result file: " + filePath.string());
  }
}
//</editor-fold>
