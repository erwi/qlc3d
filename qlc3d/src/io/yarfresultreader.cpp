#include <io/yarfresultreader.h>
#include <json.hpp>
#include <util/exception.h>
#include <util/logging.h>
#include <array>
#include <filesystem>
#include <fstream>
#include <string>

Yarf::FieldDescriptor findFieldDescriptor(const std::string fieldName, std::vector<Yarf::FieldDescriptor> fields) {
  for (auto field : fields) {
    if (field.name == fieldName) {
      return field;
    }
  }
  RUNTIME_ERROR("Field not found: " + fieldName);
}

double readDouble(std::ifstream& inputStream) {
  double value;
  inputStream.read(reinterpret_cast<char*>(&value), sizeof(value));
  if (!inputStream) {
    RUNTIME_ERROR("Failed to read double from stream");
  }
  return value;
}

int32_t readInt(std::ifstream& inputStream) {
  int32_t value;
  inputStream.read(reinterpret_cast<char*>(&value), sizeof(value));
  if (!inputStream) {
    RUNTIME_ERROR("Failed to read int32 from stream");
  }
  return value;
}

std::vector<qlc3d::OrientationSample> readQTensorPayload(const Yarf::YarfHeader &header, std::ifstream& inputStream, const std::string &fileName) {
  auto qTensorField = findFieldDescriptor("q_tensor", header.fields);
  auto numQTensors = static_cast<unsigned int>(qTensorField.length / 5); // each q_tensor sample has 5 components
  std::vector<qlc3d::OrientationSample> result;
  result.reserve(numQTensors);

  inputStream.seekg(header.payloadStartOffset + qTensorField.byte_offset, std::ios::beg); // seek to the start of the q-tensor payload
  if (!inputStream) {
    RUNTIME_ERROR("Failed to seek to q_tensor payload in file: " + fileName);
  }

  for (unsigned int i = 0; i < numQTensors; ++i) {
    double Qxx = readDouble(inputStream);
    double Qxy = readDouble(inputStream);
    double Qxz = readDouble(inputStream);
    double Qyy = readDouble(inputStream);
    double Qyz = readDouble(inputStream);

    result.push_back(qlc3d::OrientationSample {
      .tensor = qlc3d::TTensor::fromQTensor(qlc3d::QTensor(Qxx, Qyy, Qxy, Qyz, Qxz)),
      .location = std::nullopt
    });
  }

  return result;
}

std::vector<qlc3d::OrientationSample> YarfResultReader::read(const std::string &fileName, double s0) const {

  Yarf::YarfHeader header = YarfResultReader::readYarfHeader(fileName);

  std::ifstream inputStream(fileName, std::ios::binary);

  return readQTensorPayload(header, inputStream, fileName);
}

Yarf::YarfHeader YarfResultReader::readYarfHeader(const std::string &fileName) {
  const std::filesystem::path filePath(fileName);
  if (!std::filesystem::exists(filePath)) {
    RUNTIME_ERROR("File does not exist " + fileName);
  }

  std::ifstream inputStream(filePath, std::ios::binary);
  if (!inputStream.is_open() || !inputStream.good()) {
    RUNTIME_ERROR("Failed to open file: " + fileName);
  }

  std::array<char, 4> magic{};
  inputStream.read(magic.data(), static_cast<std::streamsize>(magic.size()));
  if (!inputStream || std::string(magic.data(), magic.size()) != "YAR1") {
    RUNTIME_ERROR("Invalid YARF file: " + fileName);
  }

  // read the 4-byte little-endian JSON length prefix
  int32_t jsonLength = readInt(inputStream);

  if (jsonLength <= 0) {
    RUNTIME_ERROR("JSON header length must be greater than zero in file " + fileName);
  }

  std::string jsonText(jsonLength, '\0');
  inputStream.read(jsonText.data(), static_cast<std::streamsize>(jsonText.size()));
  if (!inputStream) {
    RUNTIME_ERROR("Failed to read JSON header from file: " + fileName);
  }

  const nlohmann::json jsonHeader = nlohmann::json::parse(jsonText);
  Yarf::YarfHeader header = jsonHeader.get<Yarf::YarfHeader>();

  constexpr std::int64_t PREFIX_SIZE = 4 + 4; // 4 bytes for magic number + 4 bytes for JSON length prefix
  header.payloadStartOffset = PREFIX_SIZE + jsonLength; // set the payload start offset to the position after the magic number and JSON length prefix

  return header;
}
