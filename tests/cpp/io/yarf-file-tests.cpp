#include <catch.h>

#include <array>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <vector>

#include <io/yarf.h>
#include <lc-representation.h>
#include <test-util.h>

namespace fs = std::filesystem;

namespace {

void writeLittleEndianDouble(std::ostream &stream, double value) {
  static_assert(sizeof(double) == 8, "double must be 8 bytes to match YARF float64 encoding");
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(value));
  std::array<unsigned char, sizeof(double)> bytes = {{0, 0, 0, 0, 0, 0, 0, 0}};
  for (std::size_t index = 0; index < bytes.size(); ++index) {
    bytes[index] = static_cast<unsigned char>((bits >> (8U * index)) & 0xFFu);
  }
  stream.write(reinterpret_cast<const char *>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
}

class FailingYarfWriter final : public qlc3d::yarf::IYarfWriter {
 public:
  void writeFile(const fs::path &path, const qlc3d::yarf::YarfHeader &header) const override {
    (void)header;
    throw qlc3d::yarf::YarfFileException("Injected YARF write failure.", path, "injected writer");
  }
};
}  // namespace

TEST_CASE("YARF header-only round trip succeeds") {
  // GIVEN a valid YARF header with simulation metadata, no payload fields, and a well-formed metadata map.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "header-only.yarf";

  qlc3d::yarf::YarfHeader header;
  header.format.version = "1.0";
  header.format.endianness = "little";
  header.format.payload_size_bytes = 0;
  header.simulation.mesh_file = "mesh.msh";
  header.simulation.iteration = 7;
  header.simulation.time = 3.5;
  header.simulation.S0 = 0.7;
  header.simulation.num_nodes = 12;
  header.simulation.num_elements = 9;
  header.simulation.npLC = 10;
  header.metadata = {
      {"description", "header-only snapshot"},
      {"solver_converged", false},
  };

  // WHEN the header is written to disk and read back.
  qlc3d::yarf::YarfFileWriter writer;
  writer.writeFile(path, header);

  qlc3d::yarf::YarfFileReader reader;
  const auto readHeader = reader.readFile(path);

  // THEN the round-tripped values match the original header and no field payload is present.
  REQUIRE(readHeader.format.version == "1.0");
  REQUIRE(readHeader.format.endianness == "little");
  REQUIRE(readHeader.format.payload_size_bytes == 0);
  REQUIRE(readHeader.simulation.mesh_file == "mesh.msh");
  REQUIRE(readHeader.metadata.at("description") == "header-only snapshot");
  REQUIRE(readHeader.fields.empty());
  REQUIRE(readHeader.simulation.iteration == header.simulation.iteration);
  REQUIRE(readHeader.simulation.time == header.simulation.time);
  REQUIRE(readHeader.simulation.S0 == header.simulation.S0);
  REQUIRE(readHeader.simulation.num_nodes == header.simulation.num_nodes);
  REQUIRE(readHeader.simulation.num_elements == header.simulation.num_elements);
  REQUIRE(readHeader.simulation.npLC == header.simulation.npLC);
  REQUIRE(readHeader.metadata == header.metadata);
}

TEST_CASE("YARF reader rejects invalid magic header") {
  // GIVEN a file whose first four bytes do not match the YARF magic header.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "bad-magic.yarf";

  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  output.write("NOPE", 4);
  output.write("\x10\x00\x00\x00", 4);
  output << "{}";
  output.close();

  // WHEN the reader tries to parse the file.
  qlc3d::yarf::YarfFileReader reader;

  // THEN it raises a YARF file exception that reports the invalid magic header.
  REQUIRE_THROWS_AS(reader.readFile(path), qlc3d::yarf::YarfFileException);
  REQUIRE_THROWS_WITH(reader.readFile(path), Catch::Contains(path.string()) && Catch::Contains("magic header"));
}

TEST_CASE("YARF reader rejects unsupported format version and endianness") {
  // GIVEN a JSON header that declares a non-canonical format version or byte order.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "bad-format.yarf";

  const nlohmann::json root = {
      {"format", {{"version", "2.0"}, {"endianness", "big"}, {"payload_size_bytes", 0}}},
      {"simulation", {{"mesh_file", "mesh.msh"}, {"iteration", 0}, {"time", 0.0}, {"S0", 0.0}, {"num_nodes", 1}, {"num_elements", 1}, {"npLC", 1}}},
      {"metadata", {{"description", "bad format"}}},
      {"fields", nlohmann::json::array()}
  };

  const std::string jsonText = root.dump();
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  output.write("YAR1", 4);
  const std::uint32_t jsonLength = static_cast<std::uint32_t>(jsonText.size());
  output.write(reinterpret_cast<const char *>(&jsonLength), sizeof(jsonLength));
  output.write(jsonText.data(), static_cast<std::streamsize>(jsonText.size()));
  output.close();

  // WHEN the header validation enforces the accepted contract.
  qlc3d::yarf::YarfFileReader reader;

  // THEN the reader rejects the unsupported format metadata before reading payload data.
  REQUIRE_THROWS_AS(reader.readFile(path), qlc3d::yarf::YarfSchemaException);
  REQUIRE_THROWS_WITH(reader.readFile(path), Catch::Contains("format"));
}

TEST_CASE("YARF reader rejects mismatched size versus payload size") {
  // GIVEN a YARF file whose declared JSON length does not match the actual file size.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "bad-size.yarf";

  nlohmann::json root = nlohmann::json::object({
      {"format", {
           {"version", "1.0"},
           {"endianness", "little"},
           {"payload_size_bytes", 0}
        }},
      {"simulation", {
           {"mesh_file", "mesh.msh"},
           {"iteration", 1},
           {"time", 0.0},
           {"S0", 0.0},
           {"num_nodes", 4},
           {"num_elements", 2},
           {"npLC", 4}
        }},
      {"metadata", {{"description", "bad size"}}},
      {"fields", nlohmann::json::array()}
  });

  const std::string jsonText = root.dump();
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  output.write("YAR1", 4);
  output.write("\x00\x00\x00\x00", 4);
  output.write(jsonText.data(), static_cast<std::streamsize>(jsonText.size()));
  output.write("\xFF\xFF\xFF\xFF", 4);

  const std::uint32_t declaredJsonLength = static_cast<std::uint32_t>(jsonText.size());
  output.seekp(4, std::ios::beg);
  output.write(reinterpret_cast<const char *>(&declaredJsonLength), sizeof(declaredJsonLength));
  output.close();

  // WHEN the header is validated before the payload is read.
  qlc3d::yarf::YarfFileReader reader;

  // THEN the reader throws a schema exception that describes the file-size mismatch.
  REQUIRE_THROWS_AS(reader.readFile(path), qlc3d::yarf::YarfSchemaException);
  REQUIRE_THROWS_WITH(reader.readFile(path), Catch::Contains(path.string()) && Catch::Contains("actual file size"));
}

TEST_CASE("YARF reader rejects missing required top-level keys") {
  // GIVEN a YARF payload that omits a required top-level key such as the fields array.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "missing-keys.yarf";

  const nlohmann::json root = {
      {"format", {{"version", "1.0"}, {"endianness", "little"}, {"payload_size_bytes", 0}}},
      {"simulation", {{"mesh_file", "mesh.msh"}, {"iteration", 0}, {"time", 0.0}, {"S0", 0.0}, {"num_nodes", 1}, {"num_elements", 1}, {"npLC", 1}}},
      {"metadata", {"description", "missing fields"}}
  };

  const std::string jsonText = root.dump();
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  output.write("YAR1", 4);
  const std::uint32_t jsonLength = static_cast<std::uint32_t>(jsonText.size());
  output.write(reinterpret_cast<const char *>(&jsonLength), sizeof(jsonLength));
  output.write(jsonText.data(), static_cast<std::streamsize>(jsonText.size()));
  output.close();

  // WHEN the schema validator checks the JSON structure.
  qlc3d::yarf::YarfFileReader reader;

  // THEN it throws a schema exception that explains the missing required key.
  REQUIRE_THROWS_AS(reader.readFile(path), qlc3d::yarf::YarfSchemaException);
  REQUIRE_THROWS_WITH(reader.readFile(path), Catch::Contains(path.string()) && Catch::Contains("Missing required YARF top-level key"));
}

TEST_CASE("YARF accepts unknown metadata keys") {
  // GIVEN a valid YARF header that includes future and legacy metadata entries.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "unknown-metadata.yarf";

  qlc3d::yarf::YarfHeader header;
  header.format.version = "1.0";
  header.format.endianness = "little";
  header.format.payload_size_bytes = 0;
  header.simulation.mesh_file = "mesh.msh";
  header.simulation.iteration = 5;
  header.simulation.time = 0.0;
  header.simulation.S0 = 0.0;
  header.simulation.num_nodes = 3;
  header.simulation.num_elements = 2;
  header.simulation.npLC = 3;
  header.metadata = {{"future_key", "future_value"}, {"legacy_flag", 64}};

  // WHEN the file is written and then read back.
  qlc3d::yarf::YarfFileWriter writer;
  writer.writeFile(path, header);

  qlc3d::yarf::YarfFileReader reader;
  const auto readHeader = reader.readFile(path);

  // THEN the metadata is preserved and the reader accepts keys it does not explicitly know.
  REQUIRE(readHeader.metadata.at("future_key") == "future_value");
  REQUIRE(readHeader.metadata.at("legacy_flag") == 64);
}

TEST_CASE("YARF q_tensor payload round trip preserves values") {
  // GIVEN a simulation with a q-tensor solution vector and a valid YARF field definition.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "qtensor.yarf";

  qlc3d::yarf::SimulationBlock simulation;
  simulation.mesh_file = "mesh.msh";
  simulation.iteration = 11;
  simulation.time = 1.25;
  simulation.S0 = 0.7;
  simulation.num_nodes = 9;
  simulation.num_elements = 6;
  simulation.npLC = 3;

  SolutionVector q(static_cast<idx>(simulation.npLC), 5);
  for (idx node = 0; node < simulation.npLC; ++node) {
    for (idx component = 0; component < 5; ++component) {
      const double value = static_cast<double>(node * 10 + component + 0.25);
      q.setValue(node, component, value);
    }
  }

  // WHEN the data is written to a YARF file and read back through the q-tensor helpers.
  qlc3d::yarf::YarfFileWriter writer;
  writer.writeQTensorFile(path, simulation, q, {{"description", "qtensor snapshot"}});

  qlc3d::yarf::YarfFileReader reader;
  const auto header = reader.readFile(path);

  // THEN the field metadata is correct and every q-tensor component matches the original value.
  REQUIRE(header.fields.size() == 1U);
  REQUIRE(header.fields.front().name == "q_tensor");
  REQUIRE(header.fields.front().length == static_cast<std::uint64_t>(simulation.npLC) * 5U);

  const auto expectedScalarValues = qlc3d::yarf::YarfQlc3dAdapter::exportQTensor(q);
  const auto scalarValues = reader.readQTensorValues(path);
  REQUIRE(scalarValues.size() == static_cast<std::size_t>(simulation.npLC) * 5U);
  REQUIRE(scalarValues == expectedScalarValues);

  const auto qRoundTrip = reader.readQTensorSolution(path);
  for (idx node = 0; node < simulation.npLC; ++node) {
    for (idx component = 0; component < 5; ++component) {
      REQUIRE(qRoundTrip.getValue(node, component) == Approx(q.getValue(node, component)).epsilon(1e-15));
    }
  }
}

TEST_CASE("YARF streaming payload writer accepts direct stream emission") {
  // GIVEN a YARF header with a scalar potential field and a direct stream writer callback.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "streamed-payload.yarf";

  qlc3d::yarf::YarfHeader header;
  header.simulation.num_nodes = 3;
  header.simulation.npLC = 2;
  header.format.version = "1.0";
  header.format.endianness = "little";
  header.fields.push_back(qlc3d::yarf::YarfQlc3dAdapter::makePotentialField(3, 0));
  header.format.payload_size_bytes = 3ULL * sizeof(double);

  const std::vector<double> values = {1.5, -2.0, 3.25};

  // WHEN the payload bytes are emitted through the stream callback.
  qlc3d::yarf::YarfFileWriter writer;
  writer.writeFieldPayload(path, header, [&values](std::ostream &stream) {
    for (double value : values) {
      writeLittleEndianDouble(stream, value);
    }
  });

  // THEN the values round trip without changing their order or magnitude.
  qlc3d::yarf::YarfFileReader reader;
  const auto readValues = reader.readPotentialValues(path);
  REQUIRE(readValues == values);
}

TEST_CASE("YARF q_tensor uses canonical component ordering across T-tensor conversion") {
  // GIVEN a q-tensor solution held in the internal T-tensor basis.
  SolutionVector q(1, 5);
  q.setValue(0, 0, 1.5);
  q.setValue(0, 1, -2.0);
  q.setValue(0, 2, 3.25);
  q.setValue(0, 3, 0.5);
  q.setValue(0, 4, -4.75);

  // WHEN the q-tensor values are exported using the canonical file ordering and then imported back.
  const auto qTensor = qlc3d::yarf::YarfQlc3dAdapter::exportQTensor(q);
  const auto directQ = qlc3d::TTensor(q.getValue(0, 0), q.getValue(0, 1), q.getValue(0, 2), q.getValue(0, 3), q.getValue(0, 4)).toQTensor();
  const std::vector<double> canonicalOrder = {directQ.q1(), directQ.q3(), directQ.q5(), directQ.q2(), directQ.q4()};
  REQUIRE(qTensor == canonicalOrder);

  SolutionVector roundTrip(1, 5);
  qlc3d::yarf::YarfQlc3dAdapter::importQTensor(qTensor, roundTrip);

  // THEN the exported order matches the canonical Q-tensor layout and the round trip preserves each component.
  for (idx component = 0; component < 5; ++component) {
   REQUIRE(roundTrip.getValue(0, component) == Approx(q.getValue(0, component)).epsilon(1e-16));
  }
}

TEST_CASE("YARF electric_potential payload round trip preserves values") {
  // GIVEN a potential solution for a simulation with one scalar value at each node.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "potential.yarf";

  qlc3d::yarf::SimulationBlock simulation;
  simulation.mesh_file = "mesh.msh";
  simulation.iteration = 12;
  simulation.time = 2.5;
  simulation.S0 = 0.9;
  simulation.num_nodes = 7;
  simulation.num_elements = 5;
  simulation.npLC = 3;

  SolutionVector v(static_cast<idx>(simulation.num_nodes), 1);
  for (idx node = 0; node < simulation.num_nodes; ++node) {
   v.setValue(node, 0, static_cast<double>(node) * 0.75 - 1.5);
  }

  // WHEN the potential field is written to disk and read back through the YARF helpers.
  qlc3d::yarf::YarfFileWriter writer;
  writer.writePotentialFile(path, simulation, v, {{"description", "potential snapshot"}});

  qlc3d::yarf::YarfFileReader reader;
  const auto header = reader.readFile(path);

  // THEN the field metadata is correct and every potential value matches the original solution.
  REQUIRE(header.fields.size() == 1U);
  REQUIRE(header.fields.front().name == "electric_potential");
  REQUIRE(header.fields.front().length == static_cast<std::uint64_t>(simulation.num_nodes));

  const auto expectedValues = qlc3d::yarf::YarfQlc3dAdapter::exportPotential(v);
  const auto scalarValues = reader.readPotentialValues(path);
  REQUIRE(scalarValues.size() == static_cast<std::size_t>(simulation.num_nodes));
  REQUIRE(scalarValues == expectedValues);

  const auto vRoundTrip = reader.readPotentialSolution(path);
  for (idx node = 0; node < simulation.num_nodes; ++node) {
   REQUIRE(vRoundTrip.getValue(node, 0) == Approx(v.getValue(node, 0)).epsilon(1e-15));
  }
}

TEST_CASE("YARF q_tensor plus electric_potential combined snapshot round trip preserves values") {
  // GIVEN a simulation snapshot containing both a q-tensor field and an electric potential field.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "combined.yarf";

  qlc3d::yarf::SimulationBlock simulation;
  simulation.mesh_file = "mesh.msh";
  simulation.iteration = 13;
  simulation.time = 3.0;
  simulation.S0 = 0.8;
  simulation.num_nodes = 5;
  simulation.num_elements = 3;
  simulation.npLC = 3;

  SolutionVector q(static_cast<idx>(simulation.npLC), 5);
  for (idx node = 0; node < simulation.npLC; ++node) {
   for (idx component = 0; component < 5; ++component) {
     q.setValue(node, component, static_cast<double>(node * 10 + component + 1.0));
   }
  }

  SolutionVector v(static_cast<idx>(simulation.num_nodes), 1);
  for (idx node = 0; node < simulation.num_nodes; ++node) {
   v.setValue(node, 0, static_cast<double>(node) - 0.5);
  }

  // WHEN the combined snapshot is written and read back.
  qlc3d::yarf::YarfFileWriter writer;
  writer.writeQTensorAndPotentialFile(path, simulation, q, v, {{"description", "combined snapshot"}});

  qlc3d::yarf::YarfFileReader reader;
  const auto header = reader.readFile(path);

  // THEN both field entries remain present and each solution component survives the round trip.
  REQUIRE(header.fields.size() == 2U);
  REQUIRE(header.fields.at(0).name == "q_tensor");
  REQUIRE(header.fields.at(1).name == "electric_potential");

  const auto qRoundTrip = reader.readQTensorSolution(path);
  for (idx node = 0; node < simulation.npLC; ++node) {
   for (idx component = 0; component < 5; ++component) {
     REQUIRE(qRoundTrip.getValue(node, component) == Approx(q.getValue(node, component)).epsilon(1e-15));
   }
  }

  const auto vRoundTrip = reader.readPotentialSolution(path);
  for (idx node = 0; node < simulation.num_nodes; ++node) {
   REQUIRE(vRoundTrip.getValue(node, 0) == Approx(v.getValue(node, 0)).epsilon(1e-15));
  }
}

TEST_CASE("YARF snapshot reader resolves state-based filenames and materialises runtime vectors") {
  // GIVEN a combined snapshot written with the canonical resultNNNNNNNN.yarf naming convention.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path directory = tmpDir.path();

  qlc3d::yarf::SimulationBlock simulation;
  simulation.mesh_file = "mesh.msh";
  simulation.iteration = 42;
  simulation.time = 1.75;
  simulation.S0 = 0.2;
  simulation.num_nodes = 4;
  simulation.num_elements = 2;
  simulation.npLC = 2;

  SolutionVector q(static_cast<idx>(simulation.npLC), 5);
  for (idx node = 0; node < simulation.npLC; ++node) {
    for (idx component = 0; component < 5; ++component) {
      q.setValue(node, component, static_cast<double>(node * 10 + component + 2.0));
    }
  }

  SolutionVector v(static_cast<idx>(simulation.num_nodes), 1);
  for (idx node = 0; node < simulation.num_nodes; ++node) {
    v.setValue(node, 0, static_cast<double>(node) * 0.5 - 1.0);
  }

  const auto snapshotPath = directory / "result00000042.yarf";
  qlc3d::yarf::YarfFileWriter writer;
  writer.writeQTensorAndPotentialFile(snapshotPath, simulation, q, v, {{"description", "state-based load"}});

  SimulationState state;
  state.currentIteration(42);

  // WHEN the read seam resolves the snapshot by simulation-state naming and reads it back.
  const auto snapshot = qlc3d::readYarfSnapshot(directory, state);
  const auto directSnapshot = qlc3d::YarfSnapshotReader(directory).readSnapshot(snapshotPath);

  // THEN the snapshot header and materialised vectors match the original data for both entry points.
  REQUIRE(snapshot.header.simulation.iteration == 42);
  REQUIRE(snapshot.hasQTensor);
  REQUIRE(snapshot.hasElectricPotential);
  REQUIRE(snapshot.qTensor.getnDoF() == simulation.npLC);
  REQUIRE(snapshot.electricPotential.getnDoF() == simulation.num_nodes);
  REQUIRE(snapshot.header.metadata.at("description") == "state-based load");
  for (idx node = 0; node < simulation.npLC; ++node) {
    for (idx component = 0; component < 5; ++component) {
      REQUIRE(snapshot.qTensor.getValue(node, component) == Approx(q.getValue(node, component)).epsilon(1e-15));
      REQUIRE(directSnapshot.qTensor.getValue(node, component) == Approx(q.getValue(node, component)).epsilon(1e-15));
    }
  }
  for (idx node = 0; node < simulation.num_nodes; ++node) {
    REQUIRE(snapshot.electricPotential.getValue(node, 0) == Approx(v.getValue(node, 0)).epsilon(1e-15));
    REQUIRE(directSnapshot.electricPotential.getValue(node, 0) == Approx(v.getValue(node, 0)).epsilon(1e-15));
  }
}

TEST_CASE("YARF electric_potential rejects mismatched field length") {
  // GIVEN a potential field whose declared length does not match the payload length.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "potential-bad-length.yarf";

  qlc3d::yarf::SimulationBlock simulation;
  simulation.num_nodes = 4;
  simulation.npLC = 3;
  qlc3d::yarf::YarfHeader header;
  header.simulation = simulation;
  header.format.version = "1.0";
  header.format.endianness = "little";
  header.fields.push_back(qlc3d::yarf::YarfQlc3dAdapter::makePotentialField(4, 0));
  header.fields.front().length = 3;
  header.format.payload_size_bytes = 8 * 3;

  // WHEN the malformed file is written and the reader validates it.
  qlc3d::yarf::YarfFileWriter writer;
  writer.writeFieldPayload(path, header, std::vector<double>(3, 1.0));

  qlc3d::yarf::YarfFileReader reader;

  // THEN it throws a schema exception because the field length is inconsistent.
  REQUIRE_THROWS_AS(reader.readFile(path), qlc3d::yarf::YarfSchemaException);
}

TEST_CASE("YARF q_tensor rejects reordered component names") {
  // GIVEN a serialized q-tensor field whose component names are in the wrong order for the canonical layout.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "qtensor-bad-components.yarf";

  nlohmann::json root = nlohmann::json::object({
      {"format", {{"version", "1.0"}, {"endianness", "little"}, {"payload_size_bytes", 120}}},
      {"simulation", {{"mesh_file", "mesh.msh"}, {"iteration", 1}, {"time", 0.0}, {"S0", 0.0}, {"num_nodes", 3}, {"num_elements", 2}, {"npLC", 3}}},
      {"metadata", {{"description", "bad q tensor"}}},
      {"fields", nlohmann::json::array({
          {
            {"name", "q_tensor"},
            {"association", "lc_nodes"},
            {"datatype", "float64"},
            {"components", nlohmann::json::array({"Qxx", "Qxz", "Qxy", "Qyy", "Qyz"})},
            {"length", 15},
            {"byte_offset", 0}
          }
      })}
  });

  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  output.write("YAR1", 4);
  const std::uint32_t jsonLength = static_cast<std::uint32_t>(root.dump().size());
  output.write(reinterpret_cast<const char *>(&jsonLength), sizeof(jsonLength));
  output << root.dump();
  output.write("\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00", 40);
  output.close();

  // WHEN the reader validates the field schema.
  qlc3d::yarf::YarfFileReader reader;

  // THEN it rejects the file with a schema exception before the payload is accepted.
  REQUIRE_THROWS_AS(reader.readFile(path), qlc3d::yarf::YarfSchemaException);
}

TEST_CASE("YARF q_tensor rejects mismatched field length") {
  // GIVEN a q-tensor field whose declared length differs from the payload length.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "qtensor-bad-length.yarf";

  qlc3d::yarf::SimulationBlock simulation;
  simulation.npLC = 3;
  qlc3d::yarf::YarfHeader header;
  header.simulation = simulation;
  header.format.version = "1.0";
  header.format.endianness = "little";
  header.format.payload_size_bytes = 8 * 14;
  header.fields.push_back(qlc3d::yarf::YarfQlc3dAdapter::makeQTensorField(3, 0));
  header.fields.front().length = 14;

  // WHEN the file is parsed.
  qlc3d::yarf::YarfFileWriter writer;
  writer.writeFieldPayload(path, header, std::vector<double>(14, 1.0));

  qlc3d::yarf::YarfFileReader reader;

  // THEN the reader rejects it as a schema violation.
  REQUIRE_THROWS_AS(reader.readFile(path), qlc3d::yarf::YarfSchemaException);
}

TEST_CASE("YARF q_tensor rejects payload range beyond declared bounds") {
  // GIVEN a q-tensor field whose payload begins beyond the valid byte range and extends past the declared bounds.
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path path = tmpDir.path() / "qtensor-bad-range.yarf";

  qlc3d::yarf::SimulationBlock simulation;
  simulation.npLC = 3;
  qlc3d::yarf::YarfHeader header;
  header.simulation = simulation;
  header.format.version = "1.0";
  header.format.endianness = "little";
  header.fields.push_back(qlc3d::yarf::YarfQlc3dAdapter::makeQTensorField(3, 0));
  header.fields.front().byte_offset = 8;
  header.format.payload_size_bytes = 8 * 15;

  // WHEN the reader reads the tensor values.
  qlc3d::yarf::YarfFileWriter writer;
  writer.writeFieldPayload(path, header, std::vector<double>(15, 1.0));

  qlc3d::yarf::YarfFileReader reader;

  // THEN it throws a schema exception because the range exceeds the valid payload window.
  REQUIRE_THROWS_AS(reader.readQTensorValues(path), qlc3d::yarf::YarfSchemaException);
}

TEST_CASE("YARF qlc3d writer seam uses the canonical resultNNNNNNNN.yarf naming convention") {
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path outputDir = tmpDir.path() / "yarf-output";
  const fs::path meshPath = "mesh.msh";
  const fs::path expectedPath = outputDir / "result00000013.yarf";

  Geometry geometry;
  geometry.setCoordinates(Coordinates::from({0.0, 0.0, 0.0,
                                            1.0, 0.0, 0.0,
                                            0.0, 1.0, 0.0,
                                            0.0, 0.0, 1.0}));
  geometry.setnpLC(2);

  SolutionVector q(2, 5);
  q.setValue(0, 0, 0.10);
  q.setValue(0, 1, -0.20);
  q.setValue(0, 2, 0.30);
  q.setValue(0, 3, 0.40);
  q.setValue(0, 4, -0.50);
  q.setValue(1, 0, 0.11);
  q.setValue(1, 1, -0.21);
  q.setValue(1, 2, 0.31);
  q.setValue(1, 3, 0.41);
  q.setValue(1, 4, -0.51);

  SolutionVector potential(4, 1);
  for (idx node = 0; node < 4; ++node) {
    potential.setValue(node, 0, static_cast<double>(node) * 0.5 - 1.0);
  }

  SimulationState simulationState;
  simulationState.currentIteration(13);
  simulationState.setCurrentTime(2.5);
  simulationState.state(RunningState::COMPLETED);

  qlc3d::writeYarfSnapshot(outputDir, geometry, potential, q, simulationState, meshPath, { {"description", "integration snapshot"} });

  REQUIRE(fs::exists(expectedPath));
  REQUIRE(!fs::exists(outputDir / "result-final.yarf"));
}

TEST_CASE("YARF qlc3d writer seam writes a combined snapshot from runtime containers") {
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path outputDir = tmpDir.path() / "yarf-output";
  const fs::path meshPath = "mesh.msh";
  const fs::path expectedPath = outputDir / "result00000013.yarf";

  Geometry geometry;
  geometry.setCoordinates(Coordinates::from({0.0, 0.0, 0.0,
                                            1.0, 0.0, 0.0,
                                            0.0, 1.0, 0.0,
                                            0.0, 0.0, 1.0}));
  geometry.setnpLC(2);

  SolutionVector q(2, 5);
  q.setValue(0, 0, 0.10);
  q.setValue(0, 1, -0.20);
  q.setValue(0, 2, 0.30);
  q.setValue(0, 3, 0.40);
  q.setValue(0, 4, -0.50);
  q.setValue(1, 0, 0.11);
  q.setValue(1, 1, -0.21);
  q.setValue(1, 2, 0.31);
  q.setValue(1, 3, 0.41);
  q.setValue(1, 4, -0.51);

  SolutionVector potential(4, 1);
  for (idx node = 0; node < 4; ++node) {
    potential.setValue(node, 0, static_cast<double>(node) * 0.5 - 1.0);
  }

  SimulationState simulationState;
  simulationState.currentIteration(13);
  simulationState.setCurrentTime(2.5);
  simulationState.state(RunningState::RUNNING);

  qlc3d::writeYarfSnapshot(outputDir, geometry, potential, q, simulationState, meshPath, { {"description", "integration snapshot"} });

  REQUIRE(fs::exists(expectedPath));
  qlc3d::yarf::YarfFileReader reader;
  const auto header = reader.readFile(expectedPath);
  REQUIRE(header.fields.size() == 2U);

  const auto qRoundTrip = reader.readQTensorSolution(expectedPath);
  for (idx node = 0; node < q.getnDoF(); ++node) {
    for (idx component = 0; component < 5; ++component) {
      REQUIRE(qRoundTrip.getValue(node, component) == Approx(q.getValue(node, component)).epsilon(1e-15));
    }
  }

  const auto potentialRoundTrip = reader.readPotentialSolution(expectedPath);
  for (idx node = 0; node < potential.getnDoF(); ++node) {
    REQUIRE(potentialRoundTrip.getValue(node, 0) == Approx(potential.getValue(node, 0)).epsilon(1e-15));
  }
}

TEST_CASE("YARF qlc3d writer seam translates injected YARF write failures to runtime errors") {
  TestUtil::TemporaryDirectory tmpDir;
  const fs::path outputDir = tmpDir.path() / "yarf-failure";

  Geometry geometry;
  geometry.setCoordinates(Coordinates::from({0.0, 0.0, 0.0,
                                            1.0, 0.0, 0.0,
                                            0.0, 1.0, 0.0,
                                            0.0, 0.0, 1.0}));
  geometry.setnpLC(1);

  SolutionVector q(1, 5);
  q.setValue(0, 0, 0.1);
  q.setValue(0, 1, 0.2);
  q.setValue(0, 2, 0.3);
  q.setValue(0, 3, 0.4);
  q.setValue(0, 4, 0.5);

  SolutionVector potential(4, 1);
  for (idx node = 0; node < 4; ++node) {
    potential.setValue(node, 0, static_cast<double>(node));
  }

  SimulationState simulationState;
  simulationState.currentIteration(7);
  simulationState.setCurrentTime(1.0);

  const auto failingWriter = std::make_shared<FailingYarfWriter>();
  REQUIRE_THROWS_WITH(
      qlc3d::writeYarfSnapshot(outputDir, geometry, potential, q, simulationState, "mesh.msh", {}, failingWriter),
      Catch::Contains("YARF snapshot write failed for qlc3d output") && Catch::Contains("Injected YARF write failure"));
}
