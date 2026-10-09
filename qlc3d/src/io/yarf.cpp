#include <io/yarf.h>

#include <array>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <lc-representation.h>
#include <system_error>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <util/exception.h>
#include <vector>

namespace qlc3d::yarf {
namespace {

/**
 * Reads the full contents of a YARF file into memory for header validation and payload extraction.
 * @param path Candidate YARF file.
 * @return File contents as a byte string.
 */
std::string readAllText(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw YarfFileException("Could not open YARF file for reading.", path, "read");
  }
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

/**
 * Decodes a 32-bit little-endian integer from a byte string.
 * @param bytes Source buffer.
 * @param offset Byte offset to start reading.
 * @param path File path used for diagnostics.
 * @return Little-endian uint32 value.
 */
std::uint32_t readUint32LittleEndian(const std::string &bytes, std::size_t offset, const std::filesystem::path &path = {}) {
  if (offset + 4 > bytes.size()) {
    throw YarfFileException("YARF file ended before the JSON length prefix could be read.", path, "JSON length prefix");
  }
  return static_cast<std::uint32_t>(static_cast<unsigned char>(bytes[offset]) |
                                   (static_cast<unsigned char>(bytes[offset + 1]) << 8) |
                                   (static_cast<unsigned char>(bytes[offset + 2]) << 16) |
                                   (static_cast<unsigned char>(bytes[offset + 3]) << 24));
}

/**
 * Encodes a 32-bit integer into a little-endian byte sequence for the file header.
 * @param stream Destination output stream.
 * @param value Integer to encode.
 */
void writeUint32LittleEndian(std::ostream &stream, std::uint32_t value) {
  const std::array<unsigned char, 4> encoded = {
      static_cast<unsigned char>(value & 0xFFu),
      static_cast<unsigned char>((value >> 8) & 0xFFu),
      static_cast<unsigned char>((value >> 16) & 0xFFu),
      static_cast<unsigned char>((value >> 24) & 0xFFu)};
  stream.write(reinterpret_cast<const char *>(encoded.data()), static_cast<std::streamsize>(encoded.size()));
}

/**
 * Validates that a JSON header fits within the 32-bit file-length prefix used by the format.
 * @param jsonText Serialized header text.
 * @return Encoded length in bytes.
 */
std::uint32_t jsonLengthFor(const std::string &jsonText) {
  if (jsonText.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    throw YarfSchemaException("YARF JSON header exceeds the 32-bit length limit.", {}, "JSON length");
  }
  return static_cast<std::uint32_t>(jsonText.size());
}

/**
 * Returns the scalar count for a descriptor when it is present.
 * @param field Descriptor being counted.
 * @return Number of scalar values represented by the field.
 */
std::uint64_t fieldScalarCount(const FieldDescriptor &field) {
  return field.length == 0 ? 0 : field.length;
}

/**
 * Writes a IEEE-754 double in little-endian order to a binary output stream.
 * @param stream Destination stream.
 * @param value Value to encode.
 */
void writeDoubleLittleEndian(std::ostream &stream, double value) {
  static_assert(sizeof(double) == 8, "double must be 8 bytes to match YARF float64 encoding");
  std::array<unsigned char, sizeof(double)> bytes = {{0, 0, 0, 0, 0, 0, 0, 0}};
  std::memcpy(bytes.data(), &value, sizeof(value));
  stream.write(reinterpret_cast<const char *>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
}

/**
 * Extracts a contiguous block of float64 values from the binary payload at a given offset.
 * @param payload Raw payload binary blob.
 * @param offset Offset in bytes from the payload start.
 * @param count Number of scalar values to decode.
 * @return Values decoded in native double precision.
 */
std::vector<double> readDoubleScalarPayload(const std::string &payload, std::uint64_t offset, std::uint64_t count) {
  std::vector<double> values;
  values.reserve(static_cast<std::size_t>(count));
  if (offset + count * sizeof(double) > payload.size()) {
    throw YarfSchemaException("YARF field payload exceeds the file payload bounds.", {}, "field payload bounds");
  }
  for (std::uint64_t i = 0; i < count; ++i) {
    const auto start = static_cast<std::size_t>(offset + i * sizeof(double));
    double value = 0.0;
    std::memcpy(&value, payload.data() + start, sizeof(value));
    values.push_back(value);
  }
  return values;
}

/**
 * Serialises a qlc3d solution vector into the canonical YARF Q-tensor scalar order.
 * @param stream Destination byte stream.
 * @param q Solution values stored in the internal T-tensor basis.
 */
void writeQTensorValues(std::ostream &stream, const SolutionVector &q) {
  const auto nodeCount = static_cast<std::size_t>(q.getnDoF());
  for (std::size_t node = 0; node < nodeCount; ++node) {
    const auto qTensor = TTensor(
        q.getValue(static_cast<idx>(node), 0),
        q.getValue(static_cast<idx>(node), 1),
        q.getValue(static_cast<idx>(node), 2),
        q.getValue(static_cast<idx>(node), 3),
        q.getValue(static_cast<idx>(node), 4))
        .toQTensor();

    writeDoubleLittleEndian(stream, qTensor.q1());
    writeDoubleLittleEndian(stream, qTensor.q3());
    writeDoubleLittleEndian(stream, qTensor.q5());
    writeDoubleLittleEndian(stream, qTensor.q2());
    writeDoubleLittleEndian(stream, qTensor.q4());
  }
}

/**
 * Serialises a qlc3d scalar potential solution into the YARF field payload order.
 * @param stream Destination byte stream.
 * @param v Potential values to write.
 */
void writePotentialValues(std::ostream &stream, const SolutionVector &v) {
  const auto nodeCount = static_cast<std::size_t>(v.getnDoF());
  for (std::size_t node = 0; node < nodeCount; ++node) {
    writeDoubleLittleEndian(stream, v.getValue(static_cast<idx>(node), 0));
  }
}

}  // namespace

std::vector<double> YarfQlc3dAdapter::exportQTensor(const SolutionVector &q) {
  const auto nodeCount = static_cast<std::size_t>(q.getnDoF());
  std::vector<double> values;
  values.reserve(nodeCount * 5U);
  for (std::size_t node = 0; node < nodeCount; ++node) {
    // Convert from qlc3d's internal T-tensor basis to the canonical external Q-tensor basis
    // used by YARF: [Qxx, Qxy, Qxz, Qyy, Qyz].
    const auto qTensor = TTensor(
        q.getValue(node, 0),
        q.getValue(node, 1),
        q.getValue(node, 2),
        q.getValue(node, 3),
        q.getValue(node, 4))
        .toQTensor();

    values.push_back(qTensor.q1());
    values.push_back(qTensor.q3());
    values.push_back(qTensor.q5());
    values.push_back(qTensor.q2());
    values.push_back(qTensor.q4());
  }
  return values;
}

void YarfQlc3dAdapter::importQTensor(const std::vector<double> &values, SolutionVector &q) {
  if (values.size() % 5U != 0U) {
    throw YarfSchemaException("YARF Q-tensor payload length is not divisible by 5.", {}, "q_tensor payload");
  }
  const auto nodeCount = static_cast<idx>(values.size() / 5U);
  SolutionVector next(nodeCount, 5);
  for (idx node = 0; node < nodeCount; ++node) {
    const auto qTensor = QTensor(
        values[static_cast<std::size_t>(node) * 5U + 0U],
        values[static_cast<std::size_t>(node) * 5U + 3U],
        values[static_cast<std::size_t>(node) * 5U + 1U],
        values[static_cast<std::size_t>(node) * 5U + 4U],
        values[static_cast<std::size_t>(node) * 5U + 2U]);

    next.setValue(node, TTensor::fromQTensor(qTensor));
  }
  q = std::move(next);
}

std::vector<double> YarfQlc3dAdapter::exportPotential(const SolutionVector &v) {
  const auto nodeCount = static_cast<std::size_t>(v.getnDoF());
  std::vector<double> values;
  values.reserve(nodeCount);
  for (std::size_t node = 0; node < nodeCount; ++node) {
    values.push_back(v.getValue(static_cast<idx>(node), 0));
  }
  return values;
}

void YarfQlc3dAdapter::importPotential(const std::vector<double> &values, SolutionVector &v) {
  const auto nodeCount = static_cast<idx>(values.size());
  SolutionVector next(nodeCount, 1);
  for (idx node = 0; node < nodeCount; ++node) {
    next.setValue(node, 0, values[static_cast<std::size_t>(node)]);
  }
  v = std::move(next);
}

FieldDescriptor YarfQlc3dAdapter::makeQTensorField(std::uint64_t npLC, std::uint64_t byteOffset) {
  FieldDescriptor field;
  field.name = "q_tensor";
  field.association = FieldAssociation::LcNodes;
  field.datatype = DataType::Float64;
  field.components = {"Qxx", "Qxy", "Qxz", "Qyy", "Qyz"};
  field.length = static_cast<std::uint64_t>(npLC) * 5U;
  field.byte_offset = byteOffset;
  return field;
}

FieldDescriptor YarfQlc3dAdapter::makePotentialField(std::uint64_t numNodes, std::uint64_t byteOffset) {
  FieldDescriptor field;
  field.name = "electric_potential";
  field.association = FieldAssociation::AllNodes;
  field.datatype = DataType::Float64;
  field.components = {"phi"};
  field.length = numNodes;
  field.byte_offset = byteOffset;
  return field;
}

void YarfQlc3dAdapter::validateQTensorField(const FieldDescriptor &field, const SimulationBlock &simulation) {
  if (field.name != "q_tensor") {
    throw YarfSchemaException("Q-tensor validation expected the q_tensor field name.", {}, "field name");
  }
  if (field.association != FieldAssociation::LcNodes) {
    throw YarfUnsupportedFieldException("YARF q_tensor must use lc_nodes association.", {}, "q_tensor association");
  }
  if (field.datatype != DataType::Float64) {
    throw YarfUnsupportedDatatypeException("YARF q_tensor requires float64 values.", {}, "q_tensor datatype");
  }
  if (field.components != std::vector<std::string>(kQTensorComponents.begin(), kQTensorComponents.end())) {
    throw YarfSchemaException("YARF q_tensor components must be [Qxx, Qxy, Qxz, Qyy, Qyz].", {}, "q_tensor components");
  }
  if (field.length != static_cast<std::uint64_t>(simulation.npLC) * 5U) {
    throw YarfSchemaException("YARF q_tensor length does not match npLC * 5.", {}, "q_tensor length");
  }
}

void YarfQlc3dAdapter::validatePotentialField(const FieldDescriptor &field, const SimulationBlock &simulation) {
  if (field.name != "electric_potential") {
    throw YarfSchemaException("Potential validation expected the electric_potential field name.", {}, "field name");
  }
  if (field.association != FieldAssociation::AllNodes) {
    throw YarfUnsupportedFieldException("YARF electric_potential must use all_nodes association.", {}, "electric_potential association");
  }
  if (field.datatype != DataType::Float64) {
    throw YarfUnsupportedDatatypeException("YARF electric_potential requires float64 values.", {}, "electric_potential datatype");
  }
  if (field.components != std::vector<std::string>(kPotentialComponents.begin(), kPotentialComponents.end())) {
    throw YarfSchemaException("YARF electric_potential components must be [phi].", {}, "electric_potential components");
  }
  if (field.length != static_cast<std::uint64_t>(simulation.num_nodes)) {
    throw YarfSchemaException("YARF electric_potential length does not match num_nodes.", {}, "electric_potential length");
  }
}

std::string toString(FieldAssociation association) {
  switch (association) {
    case FieldAssociation::LcNodes: return "lc_nodes";
    case FieldAssociation::AllNodes: return "all_nodes";
    case FieldAssociation::LcElements: return "lc_elements";
    case FieldAssociation::AllElements: return "all_elements";
  }
  throw YarfUnsupportedFieldException("Unsupported YARF field association.", {}, "FieldAssociation::toString");
}

std::string toString(DataType datatype) {
  switch (datatype) {
    case DataType::Float64: return "float64";
    case DataType::Float32: return "float32";
  }
  throw YarfUnsupportedDatatypeException("Unsupported YARF datatype.", {}, "DataType::toString");
}

bool tryParseFieldAssociation(const std::string &value, FieldAssociation &association) {
  if (value == "lc_nodes") {
    association = FieldAssociation::LcNodes;
    return true;
  }
  if (value == "all_nodes") {
    association = FieldAssociation::AllNodes;
    return true;
  }
  if (value == "lc_elements") {
    association = FieldAssociation::LcElements;
    return true;
  }
  if (value == "all_elements") {
    association = FieldAssociation::AllElements;
    return true;
  }
  return false;
}

bool tryParseDataType(const std::string &value, DataType &datatype) {
  if (value == "float64") {
    datatype = DataType::Float64;
    return true;
  }
  if (value == "float32") {
    datatype = DataType::Float32;
    return true;
  }
  return false;
}

nlohmann::json YarfJsonCodec::encode(const YarfHeader &header) {
  nlohmann::json root = nlohmann::json::object();
  root["format"] = {
      {"version", header.format.version},
      {"endianness", header.format.endianness},
      {"payload_size_bytes", header.format.payload_size_bytes}
  };

  root["simulation"] = {
      {"mesh_file", header.simulation.mesh_file},
      {"iteration", header.simulation.iteration},
      {"time", header.simulation.time},
      {"S0", header.simulation.S0},
      {"num_nodes", header.simulation.num_nodes},
      {"num_elements", header.simulation.num_elements},
      {"npLC", header.simulation.npLC}
  };

  root["metadata"] = header.metadata.empty() ? nlohmann::json::object() : header.metadata;

  nlohmann::json fields = nlohmann::json::array();
  for (const auto &field : header.fields) {
    fields.push_back({
        {"name", field.name},
        {"association", toString(field.association)},
        {"datatype", toString(field.datatype)},
        {"components", field.components},
        {"length", field.length},
        {"byte_offset", field.byte_offset}
    });
  }
  root["fields"] = fields;
  return root;
}

YarfHeader YarfJsonCodec::decode(const nlohmann::json &root) {
  if (!root.is_object()) {
    throw YarfSchemaException("YARF JSON root must be an object.", {}, "root");
  }

  std::string error;
  if (!YarfValidator::validateHeaderObject(root, error)) {
    throw YarfSchemaException(error, {}, "header validation");
  }

  YarfHeader header;
  const auto &format = root.at("format");
  const auto &simulation = root.at("simulation");

  header.format.version = format.at("version").get<std::string>();
  header.format.endianness = format.at("endianness").get<std::string>();
  header.format.payload_size_bytes = format.at("payload_size_bytes").get<std::uint64_t>();

  header.simulation.mesh_file = simulation.at("mesh_file").get<std::string>();
  header.simulation.iteration = simulation.at("iteration").get<std::int64_t>();
  header.simulation.time = simulation.at("time").get<double>();
  header.simulation.S0 = simulation.at("S0").get<double>();
  header.simulation.num_nodes = simulation.at("num_nodes").get<std::int64_t>();
  header.simulation.num_elements = simulation.at("num_elements").get<std::int64_t>();
  header.simulation.npLC = simulation.at("npLC").get<std::int64_t>();

  header.metadata = root.at("metadata");
  if (!header.metadata.is_object()) {
    throw YarfSchemaException("YARF metadata must be an object.", {}, "metadata");
  }

  for (const auto &entry : root.at("fields")) {
    std::string fieldError;
    if (!YarfValidator::validateFieldDescriptor(entry, header.simulation, header.format.payload_size_bytes, fieldError)) {
      throw YarfSchemaException(fieldError, {}, "field validation");
    }

    FieldDescriptor field;
    field.name = entry.at("name").get<std::string>();
    const auto associationValue = entry.at("association").get<std::string>();
    if (!tryParseFieldAssociation(associationValue, field.association)) {
      throw YarfUnsupportedFieldException("Unsupported YARF field association: " + associationValue, {}, "field association: " + field.name);
    }
    const auto datatypeValue = entry.at("datatype").get<std::string>();
    if (!tryParseDataType(datatypeValue, field.datatype)) {
      throw YarfUnsupportedDatatypeException("Unsupported YARF datatype: " + datatypeValue, {}, "field datatype: " + field.name);
    }
    if (entry.contains("components") && entry.at("components").is_array()) {
      for (const auto &component : entry.at("components")) {
        field.components.push_back(component.get<std::string>());
      }
    }
    field.length = entry.at("length").get<std::uint64_t>();
    field.byte_offset = entry.at("byte_offset").get<std::uint64_t>();
    header.fields.push_back(field);
  }

  return header;
}

std::string YarfJsonCodec::encodeString(const YarfHeader &header) {
  return YarfJsonCodec::encode(header).dump();
}

YarfHeader YarfJsonCodec::decodeString(const std::string &jsonText) {
  try {
    const auto root = nlohmann::json::parse(jsonText);
    return YarfJsonCodec::decode(root);
  } catch (const nlohmann::json::parse_error &exception) {
    throw YarfJsonException("Invalid YARF JSON at byte offset " + std::to_string(exception.byte) + ": " + exception.what(),
                           {},
                           "JSON parse");
  } catch (const nlohmann::json::exception &exception) {
    throw YarfJsonException("Invalid YARF JSON: " + std::string(exception.what()), {}, "JSON parse");
  }
}

bool YarfValidator::hasMagicHeader(std::ifstream &stream) {
  std::array<char, 4> magic = {{0, 0, 0, 0}};
  stream.read(magic.data(), static_cast<std::streamsize>(magic.size()));
  return stream && std::equal(magic.begin(), magic.end(), kMagicHeader.begin(), [](char actual, char expected) {
    return actual == expected;
  });
}

std::uint32_t YarfValidator::readJsonLength(std::ifstream &stream) {
  std::array<char, 4> bytes = {{0, 0, 0, 0}};
  stream.read(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  if (!stream) {
    throw YarfFileException("YARF file is missing the JSON length prefix.", {}, "JSON length prefix");
  }
  return static_cast<std::uint32_t>(
      static_cast<unsigned char>(bytes[0]) |
      (static_cast<unsigned char>(bytes[1]) << 8) |
      (static_cast<unsigned char>(bytes[2]) << 16) |
      (static_cast<unsigned char>(bytes[3]) << 24));
}

bool YarfValidator::validateFieldDescriptor(const nlohmann::json &field,
                                          const SimulationBlock &simulation,
                                          std::uint64_t payloadSizeBytes,
                                          std::string &error) {
  if (!field.is_object()) {
    error = "YARF field descriptor must be an object.";
    return false;
  }
  if (!field.contains("name") || !field.at("name").is_string()) {
    error = "YARF field name is required.";
    return false;
  }
  if (!field.contains("association") || !field.at("association").is_string()) {
    error = "YARF field association is required.";
    return false;
  }
  if (!field.contains("datatype") || !field.at("datatype").is_string()) {
    error = "YARF field datatype is required.";
    return false;
  }
  if (!field.contains("components") || !field.at("components").is_array()) {
    error = "YARF field components must be an array.";
    return false;
  }
  if (!field.contains("length") || !field.at("length").is_number_unsigned()) {
    error = "YARF field length must be an unsigned integer.";
    return false;
  }
  if (!field.contains("byte_offset") || !field.at("byte_offset").is_number_unsigned()) {
    error = "YARF field byte_offset must be an unsigned integer.";
    return false;
  }

  const std::string name = field.at("name").get<std::string>();
  if (name == "q_tensor") {
    const std::string association = field.at("association").get<std::string>();
    if (association != "lc_nodes") {
      error = "YARF q_tensor association must be lc_nodes.";
      return false;
    }
    if (field.at("datatype").get<std::string>() != "float64") {
      error = "YARF q_tensor datatype must be float64.";
      return false;
    }
    const auto components = field.at("components").get<std::vector<std::string>>();
    const std::vector<std::string> expected = {"Qxx", "Qxy", "Qxz", "Qyy", "Qyz"};
    if (components != expected) {
      error = "YARF q_tensor components must be [Qxx, Qxy, Qxz, Qyy, Qyz].";
      return false;
    }
    const auto length = field.at("length").get<std::uint64_t>();
    const std::uint64_t expectedLength = static_cast<std::uint64_t>(simulation.npLC) * 5U;
    if (length != expectedLength) {
      error = "YARF q_tensor length does not match npLC * 5.";
      return false;
    }
    const auto offset = field.at("byte_offset").get<std::uint64_t>();
    const std::uint64_t minimumEnd = offset + length * sizeof(double);
    if (minimumEnd > payloadSizeBytes) {
      error = "YARF q_tensor field payload exceeds the file payload bounds.";
      return false;
    }
  } else if (name == "electric_potential") {
    const std::string association = field.at("association").get<std::string>();
    if (association != "all_nodes") {
      error = "YARF electric_potential association must be all_nodes.";
      return false;
    }
    if (field.at("datatype").get<std::string>() != "float64") {
      error = "YARF electric_potential datatype must be float64.";
      return false;
    }
    const auto components = field.at("components").get<std::vector<std::string>>();
    const std::vector<std::string> expected = {"phi"};
    if (components != expected) {
      error = "YARF electric_potential components must be [phi].";
      return false;
    }
    const auto length = field.at("length").get<std::uint64_t>();
    const std::uint64_t expectedLength = static_cast<std::uint64_t>(simulation.num_nodes);
    if (length != expectedLength) {
      error = "YARF electric_potential length does not match num_nodes.";
      return false;
    }
    const auto offset = field.at("byte_offset").get<std::uint64_t>();
    const std::uint64_t minimumEnd = offset + length * sizeof(double);
    if (minimumEnd > payloadSizeBytes) {
      error = "YARF electric_potential field payload exceeds the file payload bounds.";
      return false;
    }
  }

  return true;
}

bool YarfValidator::validateHeaderObject(const nlohmann::json &root, std::string &error) {
  if (!root.is_object()) {
    error = "YARF JSON root must be an object.";
    return false;
  }

  const std::vector<std::string> requiredKeys = {"format", "simulation", "metadata", "fields"};
  for (const auto &key : requiredKeys) {
    if (!root.contains(key)) {
      error = "Missing required YARF top-level key: " + key;
      return false;
    }
  }

  if (!root.at("format").is_object()) {
    error = "YARF format block must be an object.";
    return false;
  }
  if (!root.at("simulation").is_object()) {
    error = "YARF simulation block must be an object.";
    return false;
  }
  if (!root.at("metadata").is_object()) {
    error = "YARF metadata block must be an object.";
    return false;
  }
  if (!root.at("fields").is_array()) {
    error = "YARF fields block must be an array.";
    return false;
  }

  const auto &format = root.at("format");
  if (!format.contains("version") || !format.at("version").is_string()) {
    error = "YARF format version is required.";
    return false;
  }
  if (format.at("version").get<std::string>() != "1.0") {
    error = "YARF format version must be '1.0'.";
    return false;
  }
  if (!format.contains("endianness") || !format.at("endianness").is_string()) {
    error = "YARF format endianness is required.";
    return false;
  }
  if (format.at("endianness").get<std::string>() != "little") {
    error = "YARF format endianness must be 'little'.";
    return false;
  }
  if (!format.contains("payload_size_bytes") || !format.at("payload_size_bytes").is_number_unsigned()) {
    error = "YARF payload size must be an unsigned integer.";
    return false;
  }

  return true;
}

bool YarfValidator::validateFileSize(std::uint64_t actualFileSize,
                                    std::uint32_t jsonLength,
                                    std::uint64_t payloadSizeBytes,
                                    std::string &error) {
  const std::uint64_t expectedFileSize = static_cast<std::uint64_t>(8) + jsonLength + payloadSizeBytes;
  if (actualFileSize != expectedFileSize) {
    error = "YARF actual file size does not match the declared JSON/payload lengths.";
    return false;
  }
  return true;
}

bool YarfValidator::validateFile(const std::filesystem::path &path, const nlohmann::json &root, std::string &error) {
  std::error_code ec;
  const auto fileSize = static_cast<std::uint64_t>(std::filesystem::file_size(path, ec));
  if (ec) {
    error = "Could not determine YARF file size: " + ec.message();
    return false;
  }

  if (!root.contains("format") || !root.at("format").contains("payload_size_bytes")) {
    error = "Missing payload_size_bytes in YARF format block.";
    return false;
  }

  const std::uint64_t declaredPayloadSize = root.at("format").at("payload_size_bytes").get<std::uint64_t>();
  std::uint32_t jsonLength = 0;
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    error = "Could not reopen YARF file for size validation.";
    return false;
  }
  input.seekg(4, std::ios::beg);
  std::array<char, 4> lengthBytes = {{0, 0, 0, 0}};
  input.read(lengthBytes.data(), static_cast<std::streamsize>(lengthBytes.size()));
  if (!input) {
    error = "YARF file ended before the JSON-length prefix.";
    return false;
  }
  jsonLength = static_cast<std::uint32_t>(
      static_cast<unsigned char>(lengthBytes[0]) |
      (static_cast<unsigned char>(lengthBytes[1]) << 8) |
      (static_cast<unsigned char>(lengthBytes[2]) << 16) |
      (static_cast<unsigned char>(lengthBytes[3]) << 24));

  return validateFileSize(fileSize, jsonLength, declaredPayloadSize, error);
}

YarfHeader YarfFileReader::readFile(const std::filesystem::path &path) const {
  const std::string bytes = readAllText(path);
  if (bytes.size() < 8) {
    throw YarfFileException("YARF file is too short to contain the magic header.", path, "file header");
  }

  if (bytes.substr(0, 4) != "YAR1") {
    throw YarfFileException("YARF file has an invalid magic header.", path, "file start");
  }

  const std::uint32_t jsonLength = readUint32LittleEndian(bytes, 4, path);
  if (bytes.size() < 8U + jsonLength) {
    throw YarfFileException("YARF file declares a JSON length longer than the actual file.", path, "JSON length");
  }

  const std::string jsonText = bytes.substr(8, jsonLength);
  try {
    const auto root = nlohmann::json::parse(jsonText);
    std::string error;
    if (!YarfValidator::validateHeaderObject(root, error)) {
      throw YarfSchemaException(error, path, "header validation");
    }
    if (!YarfValidator::validateFile(path, root, error)) {
      throw YarfSchemaException(error, path, "file validation");
    }
    const auto header = YarfJsonCodec::decode(root);
    if (!header.fields.empty()) {
      for (const auto &field : header.fields) {
        if (field.name == "q_tensor") {
          YarfQlc3dAdapter::validateQTensorField(field, header.simulation);
        } else if (field.name == "electric_potential") {
          YarfQlc3dAdapter::validatePotentialField(field, header.simulation);
        }
      }
    }
    return header;
  } catch (const nlohmann::json::parse_error &exception) {
    throw YarfJsonException("Invalid YARF JSON header at byte offset " + std::to_string(exception.byte) + ": " + exception.what(),
                            path,
                            "JSON payload");
  } catch (const nlohmann::json::exception &exception) {
    throw YarfJsonException("Invalid YARF JSON header: " + std::string(exception.what()), path, "JSON payload");
  }
}

std::vector<double> YarfFileReader::readFieldValues(const std::filesystem::path &path, const std::string &fieldName) const {
  const auto header = readFile(path);
  const auto bytes = readAllText(path);
  const auto jsonLength = readUint32LittleEndian(bytes, 4, path);
  const auto payloadOffset = static_cast<std::size_t>(8U + jsonLength);
  if (payloadOffset > bytes.size()) {
    throw YarfFileException("YARF file payload is truncated.", path, "payload read");
  }
  const std::string payload = bytes.substr(payloadOffset);

  for (const auto &field : header.fields) {
    if (field.name != fieldName) {
      continue;
    }
    if (field.datatype != DataType::Float64) {
      throw YarfUnsupportedDatatypeException("YARF field datatype is not supported in the current qlc3d adapter.", path, "field datatype: " + fieldName);
    }
    if (field.byte_offset + field.length * sizeof(double) > payload.size()) {
      throw YarfSchemaException("YARF field payload exceeds the file payload bounds.", path, "field payload bounds");
    }
    return readDoubleScalarPayload(payload, field.byte_offset, field.length);
  }

  throw YarfSchemaException("YARF field not found in file: " + fieldName, path, "field lookup");
}

std::vector<double> YarfFileReader::readQTensorValues(const std::filesystem::path &path) const {
  return readFieldValues(path, "q_tensor");
}

SolutionVector YarfFileReader::readQTensorSolution(const std::filesystem::path &path) const {
  const auto values = readQTensorValues(path);
  const auto header = readFile(path);
  if (header.simulation.npLC <= 0) {
    throw YarfSchemaException("YARF q_tensor simulation metadata is invalid for qlc3d import.", path, "simulation.npLC");
  }
  SolutionVector q(static_cast<idx>(header.simulation.npLC), 5);
  YarfQlc3dAdapter::importQTensor(values, q);
  return q;
}

std::vector<double> YarfFileReader::readPotentialValues(const std::filesystem::path &path) const {
  return readFieldValues(path, "electric_potential");
}

SolutionVector YarfFileReader::readPotentialSolution(const std::filesystem::path &path) const {
  const auto values = readPotentialValues(path);
  const auto header = readFile(path);
  if (header.simulation.num_nodes <= 0) {
    throw YarfSchemaException("YARF electric_potential simulation metadata is invalid for qlc3d import.", path, "simulation.num_nodes");
  }
  SolutionVector v(static_cast<idx>(header.simulation.num_nodes), 1);
  YarfQlc3dAdapter::importPotential(values, v);
  return v;
}

void YarfFileWriter::writeFile(const std::filesystem::path &path, const YarfHeader &header) const {
  const std::filesystem::path parent = path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }

  YarfHeader normalized = header;
  normalized.format.payload_size_bytes = 0;
  normalized.fields.clear();

  const std::string jsonText = YarfJsonCodec::encodeString(normalized);
  const std::uint32_t jsonLength = jsonLengthFor(jsonText);
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output) {
    throw YarfFileException("Could not open YARF file for writing.", path, "write");
  }

  output.write("YAR1", 4);
  writeUint32LittleEndian(output, jsonLength);
  output.write(jsonText.data(), static_cast<std::streamsize>(jsonText.size()));
  if (!output) {
    throw YarfFileException("Failed while writing YARF JSON header to disk.", path, "write");
  }
}

void YarfFileWriter::writeFieldPayload(const std::filesystem::path &path,
                                     const YarfHeader &header,
                                     const std::vector<double> &payload) const {
  YarfHeader normalized = header;
  if (normalized.format.payload_size_bytes != payload.size() * sizeof(double)) {
    normalized.format.payload_size_bytes = static_cast<std::uint64_t>(payload.size()) * sizeof(double);
  }

  writeFieldPayload(path, normalized, [&payload](std::ostream &output) {
    for (double value : payload) {
      writeDoubleLittleEndian(output, value);
    }
  });
}

void YarfFileWriter::writeFieldPayload(const std::filesystem::path &path,
                                     const YarfHeader &header,
                                     const std::function<void(std::ostream &stream)> &payloadWriter) const {
  const std::string jsonText = YarfJsonCodec::encodeString(header);
  const std::uint32_t jsonLength = jsonLengthFor(jsonText);
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output) {
    throw YarfFileException("Could not open YARF file for writing.", path, "write");
  }

  output.write("YAR1", 4);
  writeUint32LittleEndian(output, jsonLength);
  output.write(jsonText.data(), static_cast<std::streamsize>(jsonText.size()));
  if (!output) {
    throw YarfFileException("Failed while writing YARF JSON header to disk.", path, "write");
  }

  if (payloadWriter) {
    payloadWriter(output);
  }
  if (!output) {
    throw YarfFileException("Failed while writing YARF payload to disk.", path, "write");
  }
}

void YarfFileWriter::writeQTensorFile(const std::filesystem::path &path,
                                    const YarfHeader &header,
                                    const SolutionVector &q) const {
  YarfHeader normalized = header;
  normalized.format.version = "1.0";
  normalized.format.endianness = "little";
  normalized.fields.clear();
  normalized.fields.push_back(YarfQlc3dAdapter::makeQTensorField(static_cast<std::uint64_t>(q.getnDoF()), 0));
  normalized.format.payload_size_bytes = static_cast<std::uint64_t>(q.getnDoF()) * 5ULL * sizeof(double);
  writeFieldPayload(path, normalized, [&q](std::ostream &output) {
    writeQTensorValues(output, q);
  });
}

void YarfFileWriter::writeQTensorFile(const std::filesystem::path &path,
                                    const SimulationBlock &simulation,
                                    const SolutionVector &q,
                                    const nlohmann::json &metadata) const {
  YarfHeader header;
  header.simulation = simulation;
  header.metadata = metadata.empty() ? nlohmann::json::object() : metadata;
  header.format.version = "1.0";
  header.format.endianness = "little";
  header.fields.push_back(YarfQlc3dAdapter::makeQTensorField(static_cast<std::uint64_t>(simulation.npLC), 0));
  header.format.payload_size_bytes = static_cast<std::uint64_t>(simulation.npLC) * 5ULL * sizeof(double);
  writeFieldPayload(path, header, [&q](std::ostream &output) {
    writeQTensorValues(output, q);
  });
}

void YarfFileWriter::writePotentialFile(const std::filesystem::path &path,
                                     const YarfHeader &header,
                                     const SolutionVector &v) const {
  YarfHeader normalized = header;
  normalized.format.version = "1.0";
  normalized.format.endianness = "little";
  normalized.fields.clear();
  normalized.fields.push_back(YarfQlc3dAdapter::makePotentialField(static_cast<std::uint64_t>(v.getnDoF()), 0));
  normalized.format.payload_size_bytes = static_cast<std::uint64_t>(v.getnDoF()) * sizeof(double);
  writeFieldPayload(path, normalized, [&v](std::ostream &output) {
    writePotentialValues(output, v);
  });
}

void YarfFileWriter::writePotentialFile(const std::filesystem::path &path,
                                     const SimulationBlock &simulation,
                                     const SolutionVector &v,
                                     const nlohmann::json &metadata) const {
  YarfHeader header;
  header.simulation = simulation;
  header.metadata = metadata.empty() ? nlohmann::json::object() : metadata;
  header.format.version = "1.0";
  header.format.endianness = "little";
  header.fields.push_back(YarfQlc3dAdapter::makePotentialField(static_cast<std::uint64_t>(simulation.num_nodes), 0));
  header.format.payload_size_bytes = static_cast<std::uint64_t>(simulation.num_nodes) * sizeof(double);
  writeFieldPayload(path, header, [&v](std::ostream &output) {
    writePotentialValues(output, v);
  });
}

void YarfFileWriter::writeQTensorAndPotentialFile(const std::filesystem::path &path,
                                               const YarfHeader &header,
                                               const SolutionVector &q,
                                               const SolutionVector &v) const {
  YarfHeader normalized = header;
  normalized.format.version = "1.0";
  normalized.format.endianness = "little";
  normalized.fields.clear();
  const auto qScalarCount = static_cast<std::uint64_t>(q.getnDoF()) * 5ULL;
  const auto vScalarCount = static_cast<std::uint64_t>(v.getnDoF());
  normalized.fields.push_back(YarfQlc3dAdapter::makeQTensorField(static_cast<std::uint64_t>(q.getnDoF()), 0));
  normalized.fields.push_back(YarfQlc3dAdapter::makePotentialField(static_cast<std::uint64_t>(v.getnDoF()), qScalarCount * sizeof(double)));
  normalized.format.payload_size_bytes = (qScalarCount + vScalarCount) * sizeof(double);
  writeFieldPayload(path, normalized, [&q, &v](std::ostream &output) {
    writeQTensorValues(output, q);
    writePotentialValues(output, v);
  });
}

void YarfFileWriter::writeQTensorAndPotentialFile(const std::filesystem::path &path,
                                               const SimulationBlock &simulation,
                                               const SolutionVector &q,
                                               const SolutionVector &v,
                                               const nlohmann::json &metadata) const {
  YarfHeader header;
  header.simulation = simulation;
  header.metadata = metadata.empty() ? nlohmann::json::object() : metadata;
  header.format.version = "1.0";
  header.format.endianness = "little";
  const auto qScalarCount = static_cast<std::uint64_t>(simulation.npLC) * 5ULL;
  const auto vScalarCount = static_cast<std::uint64_t>(simulation.num_nodes);
  header.fields.push_back(YarfQlc3dAdapter::makeQTensorField(static_cast<std::uint64_t>(simulation.npLC), 0));
  header.fields.push_back(YarfQlc3dAdapter::makePotentialField(static_cast<std::uint64_t>(simulation.num_nodes), qScalarCount * sizeof(double)));
  header.format.payload_size_bytes = (qScalarCount + vScalarCount) * sizeof(double);
  writeFieldPayload(path, header, [&q, &v](std::ostream &output) {
    writeQTensorValues(output, q);
    writePotentialValues(output, v);
  });
}

}  // namespace qlc3d::yarf

namespace qlc3d {
namespace {

std::string yarfIterationSuffix(const SimulationState &simulationState) {
  char numberChar[9];
  std::snprintf(numberChar, sizeof(numberChar), "%08d", simulationState.currentIteration());
  return std::string(numberChar);
}

}  // namespace

YarfSnapshotReader::YarfSnapshotReader(std::filesystem::path inputDirectory)
    : inputDirectory_(std::move(inputDirectory)) {}

std::filesystem::path YarfSnapshotReader::snapshotPath(const SimulationState &simulationState) const {
  const auto fileName = std::filesystem::path("result" + yarfIterationSuffix(simulationState) + ".yarf");
  if (inputDirectory_.empty()) {
    return fileName;
  }
  const auto canonical = inputDirectory_ / fileName;
  return canonical;
}

std::filesystem::path YarfSnapshotReader::resolveSnapshotPath(const std::filesystem::path &snapshotPath,
                                                             const SimulationState &simulationState) const {
  if (snapshotPath.empty()) {
    return this->snapshotPath(simulationState);
  }

  const auto candidate = std::filesystem::path(snapshotPath);
  if (candidate.is_absolute() || candidate.has_parent_path() || inputDirectory_.empty()) {
    return candidate;
  }

  return inputDirectory_ / candidate;
}

YarfSnapshot YarfSnapshotReader::readSnapshot(const std::filesystem::path &snapshotPath) const {
  const auto resolvedPath = resolveSnapshotPath(snapshotPath);
  if (resolvedPath.empty()) {
    RUNTIME_ERROR("YARF snapshot path is empty.");
  }

  try {
    YarfSnapshot snapshot;
    qlc3d::yarf::YarfFileReader reader;
    snapshot.header = reader.readFile(resolvedPath);

    for (const auto &field : snapshot.header.fields) {
      if (field.name == "q_tensor") {
        snapshot.qTensor = reader.readQTensorSolution(resolvedPath);
        snapshot.hasQTensor = true;
      } else if (field.name == "electric_potential") {
        snapshot.electricPotential = reader.readPotentialSolution(resolvedPath);
        snapshot.hasElectricPotential = true;
      }
    }

    return snapshot;
  } catch (const qlc3d::yarf::YarfException &exception) {
    RUNTIME_ERROR(std::string("YARF snapshot read failed for qlc3d input: ") + exception.what());
  }
}

YarfSnapshot YarfSnapshotReader::readSnapshot(const SimulationState &simulationState) const {
  return readSnapshot(resolveSnapshotPath({}, simulationState));
}

YarfSnapshot readYarfSnapshot(const std::filesystem::path &snapshotPath) {
  return YarfSnapshotReader().readSnapshot(snapshotPath);
}

YarfSnapshot readYarfSnapshot(const std::filesystem::path &inputDirectory, const SimulationState &simulationState) {
  return YarfSnapshotReader(inputDirectory).readSnapshot(simulationState);
}

YarfSnapshotWriter::YarfSnapshotWriter(std::filesystem::path outputDirectory)
    : outputDirectory_(std::move(outputDirectory)) {}

std::filesystem::path YarfSnapshotWriter::snapshotPath(const SimulationState &simulationState) const {
  const auto fileName = std::filesystem::path("result" + yarfIterationSuffix(simulationState) + ".yarf");
  if (outputDirectory_.empty()) {
    return fileName;
  }
  return outputDirectory_ / fileName;
}

void YarfSnapshotWriter::writeSnapshot(const Geometry &geometry,
                                      const SolutionVector &potential,
                                      const SolutionVector &qTensor,
                                      const SimulationState &simulationState,
                                      const std::filesystem::path &meshPath,
                                      const nlohmann::json &metadata,
                                      const std::shared_ptr<qlc3d::yarf::IYarfWriter> &writer) const {
  if (!outputDirectory_.empty()) {
    std::error_code ec;
    std::filesystem::create_directories(outputDirectory_, ec);
    if (ec) {
      RUNTIME_ERROR("Failed to create YARF output directory: " + ec.message());
    }
  }

  const auto outputPath = snapshotPath(simulationState);
  if (qTensor.getnDimensions() != 5) {
    throw qlc3d::yarf::YarfSchemaException(
        "YARF q_tensor solution must have exactly 5 components.",
        outputPath,
        "q-tensor dimensions");
  }
  if (potential.getnDimensions() != 1) {
    throw qlc3d::yarf::YarfSchemaException(
        "YARF electric_potential solution must have exactly 1 component.",
        outputPath,
        "potential dimensions");
  }
  if (static_cast<std::uint64_t>(geometry.getnpLC()) != static_cast<std::uint64_t>(qTensor.getnDoF())) {
    throw qlc3d::yarf::YarfSchemaException(
        "YARF q_tensor node count does not match the geometry's npLC count.",
        outputPath,
        "geometry q-tensor validation");
  }
  if (static_cast<std::uint64_t>(geometry.getCoordinates().size()) != static_cast<std::uint64_t>(potential.getnDoF())) {
    throw qlc3d::yarf::YarfSchemaException(
        "YARF electric_potential node count does not match the geometry's coordinate count.",
        outputPath,
        "geometry potential validation");
  }

  if (writer) {
    try {
      qlc3d::yarf::YarfHeader header;
      header.format.version = "1.0";
      header.format.endianness = "little";
      header.simulation.mesh_file = meshPath.empty() ? "" : meshPath.string();
      header.simulation.iteration = simulationState.currentIteration();
      header.simulation.time = simulationState.currentTime().getTime();
      header.simulation.S0 = 0.0;
      header.simulation.num_nodes = static_cast<std::int64_t>(geometry.getCoordinates().size());
      header.simulation.num_elements = static_cast<std::int64_t>(geometry.getTetrahedra().getnElements());
      header.simulation.npLC = static_cast<std::int64_t>(geometry.getnpLC());
      header.metadata = metadata.empty() ? nlohmann::json::object() : metadata;
      writer->writeFile(outputPath, header);
      return;
    } catch (const qlc3d::yarf::YarfException &exception) {
      RUNTIME_ERROR(std::string("YARF snapshot write failed for qlc3d output: ") + exception.what());
    }
  }

  qlc3d::yarf::SimulationBlock simulation;
  simulation.mesh_file = meshPath.empty() ? "" : meshPath.string();
  simulation.iteration = simulationState.currentIteration();
  simulation.time = simulationState.currentTime().getTime();
  simulation.S0 = 0.0;
  simulation.num_nodes = static_cast<std::int64_t>(geometry.getCoordinates().size());
  simulation.num_elements = static_cast<std::int64_t>(geometry.getTetrahedra().getnElements());
  simulation.npLC = static_cast<std::int64_t>(geometry.getnpLC());

  try {
    qlc3d::yarf::YarfFileWriter fileWriter;
    fileWriter.writeQTensorAndPotentialFile(outputPath, simulation, qTensor, potential, metadata);
  } catch (const qlc3d::yarf::YarfException &exception) {
    RUNTIME_ERROR(std::string("YARF snapshot write failed for qlc3d output: ") + exception.what());
  }
}

void YarfSnapshotWriter::writeSnapshot(const Geometry &geometry,
                                      const SolutionVector &potential,
                                      const SolutionVector &qTensor,
                                      const SimulationState &simulationState,
                                      const std::function<void(const std::filesystem::path &path,
                                                               const qlc3d::yarf::YarfHeader &header,
                                                               const std::vector<double> &payload)> &writeCallback) const {
  if (!outputDirectory_.empty()) {
    std::error_code ec;
    std::filesystem::create_directories(outputDirectory_, ec);
    if (ec) {
      RUNTIME_ERROR("Failed to create YARF output directory: " + ec.message());
    }
  }

  const auto outputPath = snapshotPath(simulationState);
  const auto qValues = qlc3d::yarf::YarfQlc3dAdapter::exportQTensor(qTensor);
  const auto pValues = qlc3d::yarf::YarfQlc3dAdapter::exportPotential(potential);
  qlc3d::yarf::YarfHeader header;
  header.simulation.mesh_file = outputPath.string();
  header.simulation.iteration = simulationState.currentIteration();
  header.simulation.time = simulationState.currentTime().getTime();
  header.simulation.npLC = static_cast<std::int64_t>(qTensor.getnDoF());
  header.simulation.num_nodes = static_cast<std::int64_t>(potential.getnDoF());
  header.fields.push_back(qlc3d::yarf::YarfQlc3dAdapter::makeQTensorField(static_cast<std::uint64_t>(qTensor.getnDoF()), 0));
  header.fields.push_back(qlc3d::yarf::YarfQlc3dAdapter::makePotentialField(static_cast<std::uint64_t>(potential.getnDoF()), static_cast<std::uint64_t>(qValues.size()) * sizeof(double)));
  std::vector<double> payload;
  payload.reserve(qValues.size() + pValues.size());
  payload.insert(payload.end(), qValues.begin(), qValues.end());
  payload.insert(payload.end(), pValues.begin(), pValues.end());
  header.format.payload_size_bytes = static_cast<std::uint64_t>(payload.size()) * sizeof(double);
  try {
    writeCallback(outputPath, header, payload);
  } catch (const qlc3d::yarf::YarfException &exception) {
    RUNTIME_ERROR(std::string("YARF snapshot write failed for qlc3d output: ") + exception.what());
  }
}

void writeYarfSnapshot(const std::filesystem::path &outputDirectory,
                       const Geometry &geometry,
                       const SolutionVector &potential,
                       const SolutionVector &qTensor,
                       const SimulationState &simulationState,
                       const std::filesystem::path &meshPath,
                       const nlohmann::json &metadata,
                       const std::shared_ptr<qlc3d::yarf::IYarfWriter> &writer) {
  YarfSnapshotWriter writerFacade(outputDirectory);
  writerFacade.writeSnapshot(geometry, potential, qTensor, simulationState, meshPath, metadata, writer);
}

void writeYarfSnapshot(const std::filesystem::path &outputDirectory,
                       const Geometry &geometry,
                       const SolutionVector &potential,
                       const SolutionVector &qTensor,
                       const SimulationState &simulationState,
                       const std::function<void(const std::filesystem::path &path,
                                                const qlc3d::yarf::YarfHeader &header,
                                                const std::vector<double> &payload)> &writeCallback) {
  YarfSnapshotWriter writerFacade(outputDirectory);
  writerFacade.writeSnapshot(geometry, potential, qTensor, simulationState, writeCallback);
}

}  // namespace qlc3d
