
#include <fstream>
#include <io/yarfresultwriter.h>
#include <io/meshwriter.h>
#include <util/logging.h>
#include <geometry.h>
#include <lc-representation.h>
#include <simulation-state.h>
#include <util/exception.h>
#include <solutionvector.h>
#include <algorithm>
#include <array>
#include <cstring>
#include <limits>

namespace {

/** inserts mesh number before the file extension */
static std::string makeNumberedMeshName(const std::string &meshName, const unsigned int meshNumber) {
  std::stringstream ss;
  ss << meshNumber;
  std::string num;
  ss >> num;
  size_t pos = meshName.find_last_of("."); // position of separator point
  if (pos == std::string::npos) {
    RUNTIME_ERROR("Mesh name does not contain a file extension: " + meshName);
  }
  std::string numberedMeshName = meshName;
  numberedMeshName.insert(pos, num);
  return numberedMeshName;
}

static std::filesystem::path createOutputFilePath(const std::filesystem::path &outputDirectory, const unsigned int iteration, const RunningState state) {
  if (state == RunningState::RUNNING || state == RunningState::INITIALISING) {
    std::stringstream ss;
    ss << std::setw(8) << std::setfill('0') << iteration;
    return outputDirectory / ("result" + ss.str() + ".yarf");
  } else if (state == RunningState::COMPLETED) {
    return outputDirectory / "result-final.yarf";
  } else {
    RUNTIME_ERROR("Simulation state is not RUNNING or COMPLETED");
  }
}


/** Calculate total payload size in bytes based on the last field's offset and size. */
static std::uint64_t calculatePayloadSizeBytes(const std::vector<Yarf::FieldDescriptor> &fields) {
  if (fields.empty()) {
    return 0;
  }

  std::uint64_t totalSize = 0;
  totalSize = fields.back().byte_offset;
  totalSize += fields.back().length * (fields.back().datatype == Yarf::DataType::Float64 ? sizeof(double) : sizeof(float));
  return totalSize;
}

/** Returns a non-owning pointer to a named field descriptor, or nullptr when absent. */
static const Yarf::FieldDescriptor *findField(const Yarf::YarfHeader &header, const std::string &name) {
  const auto it = std::find_if(header.fields.begin(), header.fields.end(), [&name](const Yarf::FieldDescriptor &field) {
    return field.name == name;
  });
  return it == header.fields.end() ? nullptr : &*it;
}

/** Buffers double values as little-endian IEEE-754 bytes before writing them. */
class LittleEndianDoubleWriter {
  std::ofstream &out_;
  std::array<char, 64 * 1024> buffer_{};
  std::size_t bufferedBytes = 0;

  void flushBuffer() {
    if (bufferedBytes == 0) {
      return;
    }
    out_.write(buffer_.data(), static_cast<std::streamsize>(bufferedBytes));
    if (!out_) {
      RUNTIME_ERROR("Failed to write YARF payload");
    }
    bufferedBytes = 0;
  }

public:
  explicit LittleEndianDoubleWriter(std::ofstream &out) : out_(out) {
    static_assert(sizeof(double) == sizeof(std::uint64_t), "YARF requires 64-bit doubles");
    static_assert(std::numeric_limits<double>::is_iec559, "YARF requires IEEE-754 doubles");
  }

  void write(const double value) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    for (unsigned int byte = 0; byte < sizeof(bits); ++byte) {
      buffer_[bufferedBytes++] = static_cast<char>((bits >> (byte * 8U)) & 0xffU);
      if (bufferedBytes == buffer_.size()) {
        flushBuffer();
      }
    }
  }

  void finish() {
    flushBuffer();
  }
};

}

std::vector<Yarf::FieldDescriptor> YarfResultWriter::createFieldDescriptors() const {
  std::vector<Yarf::FieldDescriptor> fields;
  std::uint64_t offsetNumBytes = 0;

  if (qTensor != nullptr) {
    const std::uint64_t numValues = qTensor->getnDoF() * 5;

    fields.push_back(Yarf::FieldDescriptor{"q_tensor",
      Yarf::FieldAssociation::LcNodes,
      Yarf::DataType::Float64,
      {"Qxx", "Qxy", "Qxz", "Qyy", "Qyz"},
      numValues,
      offsetNumBytes
    });
    offsetNumBytes += numValues * sizeof(double);
  }

  if (potential != nullptr) {
    const std::uint64_t numValues = potential->getnDoF();

    fields.push_back(Yarf::FieldDescriptor{"electric_potential",
      Yarf::FieldAssociation::AllNodes,
      Yarf::DataType::Float64,
      {"phi"},
      numValues,
      offsetNumBytes
    });
    offsetNumBytes += numValues * sizeof(double);
  }
  return fields;
}

Yarf::YarfHeader YarfResultWriter::writeYarfHeader(std::ofstream &out,
                                       const std::string &meshFileName,
                                       const SimulationState &simulationState,
                                       const Geometry &geom,
                                       const double S0) {
  if (!out.is_open() || !out.good()) {
    RUNTIME_ERROR("Output stream is not open or good for writing YARF header");
  }

  // Json header format is defined in #file:doc-impl/yarf-result-file-format.md
  Yarf::YarfHeader header;
  header.fields = createFieldDescriptors();

  header.simulation.mesh_file = meshFileName;
  header.simulation.S0 = S0;
  header.simulation.iteration = simulationState.currentIteration();
  header.simulation.time = simulationState.currentTime().getTime();
  header.simulation.npLC = geom.getnpLC();
  header.simulation.num_nodes = geom.getnp();
  header.simulation.num_elements = geom.getTetrahedra().getnElements();

  header.format.version = "1.0";
  header.format.endianness = "little";
  header.format.payload_size_bytes = calculatePayloadSizeBytes(header.fields);

  // convert header to json string
  nlohmann::json jsonHeader = header;
  const std::string jsonText = jsonHeader.dump();

  // write YARF magic number
  out.write("YAR1", 4);
  // JSON section length is a 4-byte little-endian unsigned integer
  if (jsonText.size() > std::numeric_limits<std::uint32_t>::max()) {
    RUNTIME_ERROR("YARF JSON header exceeds the maximum supported size");
  }
  const auto jsonLength = static_cast<std::uint32_t>(jsonText.size());
  char jsonLengthBytes[sizeof(jsonLength)];
  for (std::size_t i = 0; i < sizeof(jsonLength); ++i) {
    jsonLengthBytes[i] = static_cast<char>((jsonLength >> (i * 8U)) & 0xffU);
  }
  out.write(jsonLengthBytes, sizeof(jsonLengthBytes));
  out.write(jsonText.data(), static_cast<std::streamsize>(jsonText.size()));
  if (!out) {
    RUNTIME_ERROR("Failed to write YARF header");
  }
  return header;
}

void YarfResultWriter::writeQTensorPayload(std::ofstream &out,
                                           const Yarf::YarfHeader &header,
                                           const std::streampos payloadStart) const {
  const auto *field = findField(header, "q_tensor");
  if (field == nullptr) {
    RUNTIME_ERROR("YARF header is missing the q_tensor field descriptor");
  }
  if (qTensor == nullptr || field->datatype != Yarf::DataType::Float64 ||
      field->length != static_cast<std::uint64_t>(qTensor->getnDoF()) * 5U) {
    RUNTIME_ERROR("YARF q_tensor field descriptor does not match the solution vector");
  }
  if (field->byte_offset > static_cast<std::uint64_t>(std::numeric_limits<std::streamoff>::max())) {
    RUNTIME_ERROR("YARF q_tensor payload offset exceeds the supported stream range");
  }
  out.seekp(payloadStart + static_cast<std::streamoff>(field->byte_offset));
  if (!out) {
    RUNTIME_ERROR("Failed to seek to q_tensor payload");
  }

  LittleEndianDoubleWriter payloadWriter(out);
  for (idx node = 0; node < qTensor->getnDoF(); ++node) {
    const qlc3d::TTensor tracelessTensor(qTensor->getValue(node, 0),
                                        qTensor->getValue(node, 1),
                                        qTensor->getValue(node, 2),
                                        qTensor->getValue(node, 3),
                                        qTensor->getValue(node, 4));
    const qlc3d::QTensor tensor = tracelessTensor.toQTensor();
    payloadWriter.write(tensor.q1());
    payloadWriter.write(tensor.q3());
    payloadWriter.write(tensor.q5());
    payloadWriter.write(tensor.q2());
    payloadWriter.write(tensor.q4());
  }
  payloadWriter.finish();
}

void YarfResultWriter::writePotentialPayload(std::ofstream &out,
                                             const Yarf::YarfHeader &header,
                                             const std::streampos payloadStart) const {
  const auto *field = findField(header, "electric_potential");
  if (field == nullptr) {
    RUNTIME_ERROR("YARF header is missing the electric_potential field descriptor");
  }
  if (potential == nullptr || field->datatype != Yarf::DataType::Float64 ||
      field->length != static_cast<std::uint64_t>(potential->getnDoF())) {
    RUNTIME_ERROR("YARF electric_potential field descriptor does not match the solution vector");
  }
  if (field->byte_offset > static_cast<std::uint64_t>(std::numeric_limits<std::streamoff>::max())) {
    RUNTIME_ERROR("YARF electric_potential payload offset exceeds the supported stream range");
  }
  out.seekp(payloadStart + static_cast<std::streamoff>(field->byte_offset));
  if (!out) {
    RUNTIME_ERROR("Failed to seek to electric_potential payload");
  }

  LittleEndianDoubleWriter payloadWriter(out);
  for (idx node = 0; node < potential->getnDoF(); ++node) {
    payloadWriter.write(potential->getValue(node));
  }
  payloadWriter.finish();
}

void YarfResultWriter::writeMeshIfRequired(const Geometry &geom, const SimulationState &simulationState) {
  bool isFirstMesh = lastMeshNumber == YarfResultWriter::UNINITIALIZED;
  bool isNewMeshRequired = simulationState.meshNumber() > static_cast<unsigned int>(lastMeshNumber);
  if (isFirstMesh || isNewMeshRequired) {
    const std::string numberedMeshName = makeNumberedMeshName(meshName, simulationState.meshNumber());
    const std::filesystem::path writtenMeshPath = outputDirectory / numberedMeshName;
    Log::info("Writing mesh file {}", writtenMeshPath.string());
    writeMeshFile(geom.getCoordinates(), geom.getTetrahedra(), geom.getTriangles(), writtenMeshPath);
    lastMeshNumber = static_cast<int>(simulationState.meshNumber());
  }
}

void YarfResultWriter::writeResult(const Geometry &geom, const SimulationState &simulationState) {
  if (!std::filesystem::exists(outputDirectory)) {
    RUNTIME_ERROR("Output directory does not exist: " + outputDirectory.string());
  }
  if (qTensor != nullptr &&
      (qTensor->getnDimensions() != 5 || qTensor->getnDoF() != geom.getnpLC())) {
    RUNTIME_ERROR("YARF q_tensor solution dimensions do not match the geometry");
  }
  if (potential != nullptr &&
      (potential->getnDimensions() != 1 || potential->getnDoF() != geom.getnp())) {
    RUNTIME_ERROR("YARF electric_potential solution dimensions do not match the geometry");
  }
  writeMeshIfRequired(geom, simulationState);


  std::filesystem::path outputFilePath = createOutputFilePath(outputDirectory, simulationState.currentIteration(), simulationState.state());
  Log::info("Writing YARF result file {}", outputFilePath.string());

  std::ofstream out(outputFilePath, std::ios::out | std::ios::binary | std::ios::trunc);
  if (!out.is_open() || !out.good()) {
    RUNTIME_ERROR("Failed to open YARF result file: " + outputFilePath.string());
  }

  const auto header = writeYarfHeader(out, makeNumberedMeshName(meshName, simulationState.meshNumber()),
                                      simulationState, geom, S0);
  const auto payloadStart = out.tellp();
  if (payloadStart == std::streampos(-1)) {
    RUNTIME_ERROR("Failed to locate YARF payload start");
  }

  if (qTensor != nullptr) {
    writeQTensorPayload(out, header, payloadStart);
  }
  if (potential != nullptr) {
    writePotentialPayload(out, header, payloadStart);
  }

  out.close();
  if (!out) {
    RUNTIME_ERROR("Failed to finalize YARF result file: " + outputFilePath.string());
  }
}
