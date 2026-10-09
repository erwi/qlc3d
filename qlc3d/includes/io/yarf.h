#ifndef PROJECT_QLC3D_YARF_H
#define PROJECT_QLC3D_YARF_H

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iosfwd>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <geometry.h>
#include <geom/coordinates.h>
#include <simulation-state.h>
#include <solutionvector.h>
#include <json.hpp>

namespace qlc3d::yarf {

/**
 * Base exception for all YARF validation and I/O issues.
 * The message includes file and location context when available to make
 * malformed YARF files easier to diagnose.
 */
class YarfException : public std::runtime_error {
 public:
  explicit YarfException(std::string message,
                        const std::filesystem::path &filePath = {},
                        const std::string &location = {})
      : std::runtime_error(buildMessage(std::move(message), filePath, location)) {}

 private:
  static std::string buildMessage(const std::string &message,
                                 const std::filesystem::path &filePath,
                                 const std::string &location) {
    std::ostringstream text;
    text << message;
    if (!filePath.empty()) {
      text << " [file: " << filePath.string() << "]";
    }
    if (!location.empty()) {
      text << " [location: " << location << "]";
    }
    return text.str();
  }
};

/**
 * Raised when the YARF file cannot be opened, read, or written.
 */
class YarfFileException : public YarfException {
 public:
  using YarfException::YarfException;
};

/**
 * Raised when the JSON header cannot be parsed or decoded.
 */
class YarfJsonException : public YarfException {
 public:
  using YarfException::YarfException;
};

/**
 * Raised when the YARF schema or file-size validation fails.
 */
class YarfSchemaException : public YarfException {
 public:
  using YarfException::YarfException;
};

/**
 * Raised when a field association is unsupported by the current YARF subset.
 */
class YarfUnsupportedFieldException : public YarfException {
 public:
  using YarfException::YarfException;
};

/**
 * Raised when a datatype is unsupported by the current YARF subset.
 */
class YarfUnsupportedDatatypeException : public YarfException {
 public:
  using YarfException::YarfException;
};

/**
 * Supported YARF field associations.
 */
enum class FieldAssociation {
  LcNodes,
  AllNodes,
  LcElements,
  AllElements
};

/**
 * Supported binary scalar datatypes for YARF payload fields.
 */
enum class DataType {
  Float64,
  Float32
};

/**
 * Format-level metadata used to validate and decode a YARF file.
 */
struct FormatBlock {
  std::string version = "1.0";
  std::string endianness = "little";
  std::uint64_t payload_size_bytes = 0;
};

/**
 * Simulation metadata captured at the snapshot level.
 */
struct SimulationBlock {
  std::string mesh_file;
  std::int64_t iteration = 0;
  double time = 0.0;
  double S0 = 0.0;
  std::int64_t num_nodes = 0;
  std::int64_t num_elements = 0;
  std::int64_t npLC = 0;
};

/**
 * Describes a single field inside a YARF payload.
 *
 * @param name Logical field name such as q_tensor or electric_potential.
 * @param association Domain on which the field is defined.
 * @param datatype Underlying scalar datatype stored in the payload.
 * @param components Component labels in AoS order.
 * @param length Total scalar values in the field.
 * @param byte_offset Offset relative to the payload block start.
 */
struct FieldDescriptor {
  std::string name;
  FieldAssociation association = FieldAssociation::LcNodes;
  DataType datatype = DataType::Float64;
  std::vector<std::string> components;
  std::uint64_t length = 0;
  std::uint64_t byte_offset = 0;
};

/**
 * Full YARF JSON header plus parsed field descriptors.
 */
struct YarfHeader {
  FormatBlock format;
  SimulationBlock simulation;
  nlohmann::json metadata = nlohmann::json::object();
  std::vector<FieldDescriptor> fields;
};

/**
 * Adapters for the qlc3d Q-tensor buffer representation.
 */
class YarfQlc3dAdapter {
 public:
  static constexpr std::array<const char *, 5> kQTensorComponents{{"Qxx", "Qxy", "Qxz", "Qyy", "Qyz"}};
  static constexpr std::array<const char *, 1> kPotentialComponents{{"phi"}};

  /**
   * Exports a Q-tensor SolutionVector into a contiguous YARF AoS scalar buffer.
   * The qlc3d state is stored in the T-tensor basis, but the YARF payload uses the
   * canonical external order [Qxx, Qxy, Qxz, Qyy, Qyz] for each LC node.
   * @param q Q-tensor state for the LC nodes stored in the traceless T-tensor basis.
   * @return Flattened scalar array in canonical YARF component order per node.
   */
  static std::vector<double> exportQTensor(const SolutionVector &q);
  /**
   * Imports a Q-tensor scalar buffer back into a qlc3d SolutionVector.
   * The incoming YARF values must be encoded in the canonical order [Qxx, Qxy, Qxz, Qyy, Qyz].
   * @param values Flattened YARF scalar buffer in canonical component order per node.
   * @param q Target vector to populate. Q-tensor values are stored in the traceless T-tensor representation.
   */
  static void importQTensor(const std::vector<double> &values, SolutionVector &q);
  /**
   * Exports a nodal scalar potential SolutionVector into a contiguous YARF scalar buffer.
   * @param v Potential state for all nodes stored as a 1-component SolutionVector.
   * @return Flattened scalar array in nodal order.
   */
  static std::vector<double> exportPotential(const SolutionVector &v);
  /**
   * Imports a scalar potential buffer back into a qlc3d SolutionVector.
   * @param values Flattened nodal potential values.
   * @param v Target vector to populate.
   */
  static void importPotential(const std::vector<double> &values, SolutionVector &v);
  /**
   * Creates the supported field descriptor for a q_tensor payload.
   * @param npLC Number of LC nodes.
   * @param byteOffset Offset in bytes from the start of the payload.
   * @return Field descriptor matching the YARF contract for q_tensor.
   */
  static FieldDescriptor makeQTensorField(std::uint64_t npLC, std::uint64_t byteOffset = 0);
  /**
   * Creates the supported field descriptor for an electric_potential payload.
   * @param numNodes Number of nodes in the snapshot mesh.
   * @param byteOffset Offset in bytes from the start of the payload.
   * @return Field descriptor matching the YARF contract for electrostatic potential.
   */
  static FieldDescriptor makePotentialField(std::uint64_t numNodes, std::uint64_t byteOffset = 0);
  /**
   * Validates a field descriptor against the current q-tensor subset contract.
   * @param field Descriptor to validate.
   * @param simulation Simulation metadata used to check node counts.
   */
  static void validateQTensorField(const FieldDescriptor &field, const SimulationBlock &simulation);
  /**
   * Validates a field descriptor against the current electric-potential subset contract.
   * @param field Descriptor to validate.
   * @param simulation Simulation metadata used to check node counts.
   */
  static void validatePotentialField(const FieldDescriptor &field, const SimulationBlock &simulation);
};

/**
 * Converts a field association into its canonical YARF string.
 * @param association Association to stringify.
 * @return Canonical string identifier.
 */
std::string toString(FieldAssociation association);
/**
 * Converts a datatype into its canonical YARF string.
 * @param datatype Datatype to stringify.
 * @return Canonical string identifier.
 */
std::string toString(DataType datatype);
/**
 * Parses a field-association string.
 * @param value String value from JSON.
 * @param association Receives the parsed enum value.
 * @return True if the value is supported.
 */
bool tryParseFieldAssociation(const std::string &value, FieldAssociation &association);
/**
 * Parses a datatype string.
 * @param value String value from JSON.
 * @param datatype Receives the parsed enum value.
 * @return True if the value is supported.
 */
bool tryParseDataType(const std::string &value, DataType &datatype);

/**
 * Reader interface for reading a YARF file into typed metadata and field descriptors.
 */
class IYarfReader {
 public:
  virtual ~IYarfReader() = default;
  /**
   * Reads a YARF file from disk.
   * @param path Path to a YARF file.
   * @return Parsed metadata and field declarations.
   */
  [[nodiscard]] virtual YarfHeader readFile(const std::filesystem::path &path) const = 0;
};

/**
 * Writer interface for serializing a YARF header to disk.
 */
class IYarfWriter {
 public:
  virtual ~IYarfWriter() = default;
  /**
   * Writes a YARF file to disk.
   * @param path Output path for the YARF file.
   * @param header Header metadata to serialize.
   */
  virtual void writeFile(const std::filesystem::path &path, const YarfHeader &header) const = 0;
};

/**
 * Converts a YARF header between native struct form and JSON.
 */
class YarfJsonCodec {
 public:
  /**
   * Serializes a YARF header into a JSON object.
   * @param header Native YARF header.
   * @return JSON representation of the header.
   */
  static nlohmann::json encode(const YarfHeader &header);
  /**
   * Parses a YARF JSON object into native header data.
   * @param root JSON root object.
   * @return Native header representation.
   */
  static YarfHeader decode(const nlohmann::json &root);
  /**
   * Serializes to a JSON string.
   * @param header Header to encode.
   * @return Compact JSON string.
   */
  static std::string encodeString(const YarfHeader &header);
  /**
   * Parses a JSON string into a header.
   * @param jsonText JSON text to parse.
   * @return Native header representation.
   */
  static YarfHeader decodeString(const std::string &jsonText);
};

/**
 * Contains validation helpers for YARF file integrity and required schema checks.
 */
class YarfValidator {
 public:
  static constexpr std::array<char, 4> kMagicHeader{{'Y', 'A', 'R', '1'}};

  /**
   * Validates a single field descriptor against the supported subset and payload bounds.
   * @param field Descriptor to validate.
   * @param simulation Simulation metadata used for length checks.
   * @param payloadSizeBytes Declared payload size in bytes.
   * @param error Receives an error description on failure.
   * @return True when the descriptor is valid for the supported subset.
   */
  static bool validateFieldDescriptor(const nlohmann::json &field,
                                      const SimulationBlock &simulation,
                                      std::uint64_t payloadSizeBytes,
                                      std::string &error);
  /**
   * Checks whether the stream begins with the YAR1 magic header.
   * @param stream Input stream positioned at the start of a candidate file.
   * @return True when the stream starts with the expected magic.
   */
  static bool hasMagicHeader(std::ifstream &stream);
  /**
   * Reads the JSON-length prefix from the magic header.
   * @param stream Input stream positioned after the magic string.
   * @return Declared JSON length in bytes.
   */
  static std::uint32_t readJsonLength(std::ifstream &stream);
  /**
   * Validates the required top-level schema keys.
   * @param root Parsed JSON root object.
   * @param error Receives an error description on failure.
   * @return True when the schema passes the required validation steps.
   */
  static bool validateHeaderObject(const nlohmann::json &root, std::string &error);
  /**
   * Validates that the on-disk size matches the JSON and payload declarations.
   * @param actualFileSize File size, in bytes, on disk.
   * @param jsonLength Declared JSON length.
   * @param payloadSizeBytes Declared payload length.
   * @param error Receives an error description on mismatch.
   * @return True when the file length matches the declarations.
   */
  static bool validateFileSize(std::uint64_t actualFileSize,
                              std::uint32_t jsonLength,
                              std::uint64_t payloadSizeBytes,
                              std::string &error);
  /**
   * Validates a full file against its declared header and size.
   * @param path YARF file path.
   * @param root Parsed root object.
   * @param error Receives an error description on failure.
   * @return True when the file passes scope-limited WP1 validation.
   */
  static bool validateFile(const std::filesystem::path &path, const nlohmann::json &root, std::string &error);
};

/**
 * Reads a YARF file from disk using the public validation and codec pipeline.
 */
class YarfFileReader : public IYarfReader {
 public:
  /**
   * Reads and validates the YARF file.
   * @param path Input path.
   * @return Parsed header metadata.
   */
  [[nodiscard]] YarfHeader readFile(const std::filesystem::path &path) const override;
  /**
   * Reads the flattened scalar payload for a named field from disk.
   * @param path YARF file path.
   * @param fieldName Field name to read.
   * @return Scalar payload for the requested field in native little-endian form.
   */
  [[nodiscard]] std::vector<double> readFieldValues(const std::filesystem::path &path,
                                                   const std::string &fieldName) const;
  /**
   * Convenience reader for the Q-tensor field payload.
   * @param path YARF file path.
   * @return Flattened Q-tensor values in AoS order.
   */
  [[nodiscard]] std::vector<double> readQTensorValues(const std::filesystem::path &path) const;
  /**
   * Convenience reader that materializes a qlc3d Q-tensor vector.
   * @param path YARF file path.
   * @return Q-tensor state loaded into a SolutionVector.
   */
  [[nodiscard]] SolutionVector readQTensorSolution(const std::filesystem::path &path) const;
  /**
   * Convenience reader for the electric potential field payload.
   * @param path YARF file path.
   * @return Flattened potential values in nodal order.
   */
  [[nodiscard]] std::vector<double> readPotentialValues(const std::filesystem::path &path) const;
  /**
   * Convenience reader that materializes a qlc3d scalar potential vector.
   * @param path YARF file path.
   * @return Potential state loaded into a SolutionVector.
   */
  [[nodiscard]] SolutionVector readPotentialSolution(const std::filesystem::path &path) const;
};

/**
 * Writes a YARF snapshot with optional binary payload support.
 */
class YarfFileWriter : public IYarfWriter {
 public:
  /**
   * Writes the YARF header and payload sentinel to disk.
   * @param path Output path.
   * @param header Header metadata to write.
   */
  void writeFile(const std::filesystem::path &path, const YarfHeader &header) const override;
  /**
   * Writes a YARF file containing a q_tensor payload for a qlc3d SolutionVector.
   * @param path Output path.
   * @param header Header metadata to write.
   * @param q Q-tensor state to serialize.
   */
  void writeQTensorFile(const std::filesystem::path &path,
                        const YarfHeader &header,
                        const SolutionVector &q) const;
  /**
   * Writes a YARF file containing a q_tensor payload from a simulation block and metadata.
   * @param path Output path.
   * @param simulation Simulation metadata.
   * @param q Q-tensor state to serialize.
   * @param metadata Free-form JSON metadata object.
   */
  void writeQTensorFile(const std::filesystem::path &path,
                        const SimulationBlock &simulation,
                        const SolutionVector &q,
                        const nlohmann::json &metadata = nlohmann::json::object()) const;
  /**
   * Writes a YARF file containing an electric_potential payload for a qlc3d SolutionVector.
   * @param path Output path.
   * @param header Header metadata to write.
   * @param v Scalar potential state to serialize.
   */
  void writePotentialFile(const std::filesystem::path &path,
                         const YarfHeader &header,
                         const SolutionVector &v) const;
  /**
   * Writes a YARF file containing an electric_potential payload from simulation metadata.
   * @param path Output path.
   * @param simulation Simulation metadata.
   * @param v Scalar potential state to serialize.
   * @param metadata Free-form JSON metadata object.
   */
  void writePotentialFile(const std::filesystem::path &path,
                         const SimulationBlock &simulation,
                         const SolutionVector &v,
                         const nlohmann::json &metadata = nlohmann::json::object()) const;
  /**
   * Writes a YARF file containing both q_tensor and electric_potential payloads.
   * @param path Output path.
   * @param header Header metadata to write.
   * @param q Q-tensor state to serialize.
   * @param v Scalar potential state to serialize.
   */
  void writeQTensorAndPotentialFile(const std::filesystem::path &path,
                                   const YarfHeader &header,
                                   const SolutionVector &q,
                                   const SolutionVector &v) const;
  /**
   * Writes a YARF file containing both q_tensor and electric_potential payloads from simulation metadata.
   * @param path Output path.
   * @param simulation Simulation metadata.
   * @param q Q-tensor state to serialize.
   * @param v Scalar potential state to serialize.
   * @param metadata Free-form JSON metadata object.
   */
  void writeQTensorAndPotentialFile(const std::filesystem::path &path,
                                   const SimulationBlock &simulation,
                                   const SolutionVector &q,
                                   const SolutionVector &v,
                                   const nlohmann::json &metadata = nlohmann::json::object()) const;
  /**
   * Writes a generic binary payload described by the header.
   * @param path Output path.
   * @param header Header metadata to write.
   * @param payload Scalars written to the payload block in little-endian order.
   */
  void writeFieldPayload(const std::filesystem::path &path,
                         const YarfHeader &header,
                         const std::vector<double> &payload) const;
  /**
  * Writes a generic binary payload by delegating byte emission to a stream callback.
  * @param path Output path.
  * @param header Header metadata to write.
  * @param payloadWriter Callback that writes the payload bytes directly to the output stream.
  */
  void writeFieldPayload(const std::filesystem::path &path,
                         const YarfHeader &header,
                         const std::function<void(std::ostream &stream)> &payloadWriter) const;
};

}  // namespace qlc3d::yarf

namespace qlc3d {

/**
 * Materialised qlc3d state read back from a YARF snapshot file.
 *
 * The snapshot keeps the parsed header plus any supported q-tensor and electric-potential
 * SolutionVector payloads that were present in the file. The default runtime registry remains
 * untouched and read-side wiring stays explicit.
 */
struct YarfSnapshot {
 qlc3d::yarf::YarfHeader header;
 SolutionVector qTensor;
 SolutionVector electricPotential;
 bool hasQTensor = false;
 bool hasElectricPotential = false;
};

/**
 * Reads a YARF snapshot file using the qlc3d runtime naming convention for iteration-based files.
 *
 * A snapshot may be selected either by an explicit file path or by resolving the canonical naming
 * pattern `resultNNNNNNNN.yarf` inside the configured input directory.
 */
class YarfSnapshotReader {
 public:
 /**
  * Constructs a reader bound to an optional input directory.
  * @param inputDirectory Directory containing snapshot files, if not using an absolute path.
  */
 explicit YarfSnapshotReader(std::filesystem::path inputDirectory = {});

 /**
  * Resolves the canonical per-iteration filename used by the YARF writer seam.
  * @param simulationState Simulation state used to determine the iteration suffix.
  * @return Path to the snapshot under the configured directory.
  */
 [[nodiscard]] std::filesystem::path snapshotPath(const SimulationState &simulationState) const;

 /**
  * Resolves a snapshot path from either an explicit path or a state-based filename.
  * @param snapshotPath Optional user-specified path. When empty, the state-based name is used.
  * @param simulationState Simulation state used to name the file when no explicit path is supplied.
  * @return Resolved snapshot path.
  */
 [[nodiscard]] std::filesystem::path resolveSnapshotPath(const std::filesystem::path &snapshotPath,
                                                       const SimulationState &simulationState = {}) const;

 /**
  * Reads a complete YARF snapshot from a concrete file path.
  * @param snapshotPath File to parse.
  * @return Parsed snapshot header and materialised runtime vectors.
  */
 [[nodiscard]] YarfSnapshot readSnapshot(const std::filesystem::path &snapshotPath) const;

 /**
  * Reads a YARF snapshot using the qlc3d iteration naming convention.
  * @param simulationState Simulation state used to resolve the snapshot filename.
  * @return Parsed snapshot header and materialised runtime vectors.
  */
 [[nodiscard]] YarfSnapshot readSnapshot(const SimulationState &simulationState) const;

 private:
 std::filesystem::path inputDirectory_;
};

/**
 * Convenience helper for loading a YARF snapshot from a direct file path.
 * @param snapshotPath Snapshot file to open and parse.
 * @return Materialised snapshot data and parsed header.
 */
[[nodiscard]] YarfSnapshot readYarfSnapshot(const std::filesystem::path &snapshotPath);

/**
 * Convenience helper for loading the canonical snapshot file for a given simulation state.
 * @param inputDirectory Directory that contains the snapshot files.
 * @param simulationState Simulation state used to resolve the iteration-based filename.
 * @return Materialised snapshot data and parsed header.
 */
[[nodiscard]] YarfSnapshot readYarfSnapshot(const std::filesystem::path &inputDirectory,
                                         const SimulationState &simulationState);

/**
 * Writes a YARF snapshot from the qlc3d runtime containers without wiring YARF into the default
 * result-output registry.
 *
 * This seam is intentionally opt-in and stays distinct from the normal save-format plumbing so later
 * runtime exposure can be added without changing the default output paths.
 */
class YarfSnapshotWriter {
 public:
  /**
   * Constructs a writer bound to a specific output directory.
   * @param outputDirectory Directory that receives generated `.yarf` snapshots.
   */
  explicit YarfSnapshotWriter(std::filesystem::path outputDirectory = {});

  /**
   * Writes a combined YARF snapshot from runtime geometry and fields.
   * @param geometry Static geometry used to populate simulation counts.
   * @param potential Scalar potential to serialize.
   * @param qTensor Q-tensor state to serialize.
   * @param simulationState Simulation progress metadata and file naming state.
   * @param meshPath Optional mesh path stored in the snapshot metadata.
   * @param metadata Additional free-form metadata emitted into the YARF header.
   * @param writer Optional injected writer implementation used for testing failure translation.
   */
  void writeSnapshot(const Geometry &geometry,
                     const SolutionVector &potential,
                     const SolutionVector &qTensor,
                     const SimulationState &simulationState,
                     const std::filesystem::path &meshPath = {},
                     const nlohmann::json &metadata = nlohmann::json::object(),
                     const std::shared_ptr<qlc3d::yarf::IYarfWriter> &writer = {}) const;

  /**
   * Writes a snapshot through an injected writer callback for failure-path tests.
   * @param writeCallback Callback invoked in place of the default file writer.
   */
  void writeSnapshot(const Geometry &geometry,
                     const SolutionVector &potential,
                     const SolutionVector &qTensor,
                     const SimulationState &simulationState,
                     const std::function<void(const std::filesystem::path &path,
                                              const qlc3d::yarf::YarfHeader &header,
                                              const std::vector<double> &payload)> &writeCallback) const;

  /**
   * Produces the canonical YARF snapshot filename for a simulation state.
   * @param simulationState Simulation state used for iteration-based naming.
   * @return Path for the generated snapshot under the configured output directory.
   */
  [[nodiscard]] std::filesystem::path snapshotPath(const SimulationState &simulationState) const;

 private:
  std::filesystem::path outputDirectory_;
};

/**
 * Convenience function for a standalone qlc3d-to-YARF snapshot write.
 * The default runtime result output remains unchanged.
 *
 * @param outputDirectory Directory that receives the generated file.
 * @param geometry Geometry snapshot used for mesh counts.
 * @param potential Potential values to serialize.
 * @param qTensor Q-tensor values to serialize.
 * @param simulationState Simulation metadata and output naming.
 * @param meshPath Optional mesh path stored in the snapshot metadata.
 * @param metadata Free-form JSON metadata included in the YARF header.
 * @param writer Optional injected writer for failure-path tests.
 */
void writeYarfSnapshot(const std::filesystem::path &outputDirectory,
                       const Geometry &geometry,
                       const SolutionVector &potential,
                       const SolutionVector &qTensor,
                       const SimulationState &simulationState,
                       const std::filesystem::path &meshPath = {},
                       const nlohmann::json &metadata = nlohmann::json::object(),
                       const std::shared_ptr<qlc3d::yarf::IYarfWriter> &writer = {});

void writeYarfSnapshot(const std::filesystem::path &outputDirectory,
                       const Geometry &geometry,
                       const SolutionVector &potential,
                       const SolutionVector &qTensor,
                       const SimulationState &simulationState,
                       const std::function<void(const std::filesystem::path &path,
                                                const qlc3d::yarf::YarfHeader &header,
                                                const std::vector<double> &payload)> &writeCallback);

}  // namespace qlc3d

#endif  // PROJECT_QLC3D_YARF_H
