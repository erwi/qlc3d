#include <io/yarfresultwriter.h>
#include <catch.h>
#include <filesystem>
#include <fstream>
#include <geometry.h>
#include <inits.h>
#include <simulation-state.h>
#include <test-util.h>
#include <io/yarfresultreader.h>
#include <solutionvector.h>

namespace fs = std::filesystem;

TEST_CASE("The writer can be constructed") {
  YarfResultWriter resultWriter("output", "test-mesh.txt", 0.55);
}

TEST_CASE("The mesh writing with mesh number behavior is correct") {
  // GIVEN: Yarf result will be written to a temporary directory
  TestUtil::TemporaryDirectory tempDir;
  const fs::path& outputDir = tempDir.path();

  YarfResultWriter resultWriter(outputDir, "test-mesh.txt", 0.55);
  SimulationState simulationState;
  simulationState.currentIteration(0);

  // AND: Geometry loaded from a small mesh test resource
  Geometry geom;
  prepareGeometryWithDefaultBoundaries(geom, TestUtil::RESOURCE_SMALL_CUBE_GMSH_MESH);

  // WHEN write it first time
  resultWriter.writeResult(geom, simulationState);

  // THEN: The mesh file should be written
  const fs::path mesh0Path = outputDir / "test-mesh0.txt";
  REQUIRE(fs::exists(mesh0Path));

  // delete the mesh file so we can detect if it is written again
  fs::remove(mesh0Path);
  REQUIRE_FALSE(fs::exists(mesh0Path));

  // WHEN: Write it again with the same mesh number
  resultWriter.writeResult(geom, simulationState);

  // THEN: The mesh file should not be written again
  REQUIRE_FALSE(fs::exists(mesh0Path));

  // WHEN: The mesh number is incremented and new result is written
  simulationState.incrementMeshNumber();
  resultWriter.writeResult(geom, simulationState);

  // THEN: A mesh with new number should be written
  REQUIRE(fs::exists(outputDir / "test-mesh1.txt"));
}

TEST_CASE("YarfResultWriter writes result file with correct header") {
  TestUtil::TemporaryDirectory tempDir;
  const fs::path& outputDir = tempDir.path();
  const double S0 = 0.2;
  YarfResultWriter resultWriter(outputDir, "test-mesh.txt", S0);
  SimulationState simulationState;
  simulationState.currentIteration(10);
  simulationState.currentTime().setTime(20e-3);

  Geometry geom;
  prepareGeometryWithDefaultBoundaries(geom, TestUtil::RESOURCE_UNIT_CUBE_DIELECTRIC_NEUMAN_GMSH_MESH);

  // WHEN: Write result
  resultWriter.writeResult(geom, simulationState);

  // THEN: The result file should be written to the output directory
  const fs::path resultFilePath = outputDir / "result00000010.yarf"; // 10 iterations, padded to 8 digits
  REQUIRE(fs::exists(resultFilePath));
  REQUIRE(fs::file_size(resultFilePath) > 0);

  // AND: The header of the output file should contain the expected values
  YarfResultReader resultReader;
  Yarf::YarfHeader yarfHeader = YarfResultReader::readYarfHeader(resultFilePath.string());

  REQUIRE(yarfHeader.simulation.mesh_file == "test-mesh0.txt");
  REQUIRE(yarfHeader.simulation.S0 == S0);
  REQUIRE(yarfHeader.simulation.iteration == simulationState.currentIteration());
  REQUIRE(yarfHeader.simulation.time == simulationState.currentTime().getTime());
  REQUIRE(yarfHeader.simulation.npLC == geom.getnpLC());
  REQUIRE(yarfHeader.simulation.num_nodes == geom.getnp());
  REQUIRE(yarfHeader.simulation.num_elements == geom.getTetrahedra().getnElements());

  REQUIRE(yarfHeader.format.version == "1.0");
  REQUIRE(yarfHeader.format.endianness == "little");
  REQUIRE(yarfHeader.format.payload_size_bytes == 0); // No payload written in this test

  // file size in number of bytes should be equal to the sum of the magic number, JSON length prefix, and JSON header length
  const std::uint64_t expectedFileSize = yarfHeader.payloadStartOffset; // no payload written, so file size is just the offset to the start of the payload
  const std::uint64_t actualFileSize = fs::file_size(resultFilePath);
  REQUIRE(actualFileSize == expectedFileSize);
}

TEST_CASE("YarfResultWriter writes result file with payload") {
  // GIVEN a geometry, simulation state, and populated Q-tensor and potential vectors.
  TestUtil::TemporaryDirectory tempDir;
  const fs::path& outputDir = tempDir.path();
  const double S0 = 0.3;
  YarfResultWriter resultWriter(outputDir, "test-mesh.txt", S0);
  SimulationState simulationState;
  simulationState.currentIteration(5);
  simulationState.currentTime().setTime(15e-3);

  Geometry geom;
  prepareGeometryWithDefaultBoundaries(geom, TestUtil::RESOURCE_UNIT_CUBE_DIELECTRIC_NEUMAN_GMSH_MESH);

  // Create dummy solution vectors for q-tensor and potential
  SolutionVector qTensor(geom.getnpLC(), 5);
  SolutionVector potential(geom.getnp(), 1);
  for (idx i = 0; i < qTensor.getnDoF(); ++i) {
    qTensor.setValue(i, {i + 0.1, i + 0.2, i + 0.3, i + 0.4, i + 0.5});
  }
  for (idx i = 0; i < potential.getnDoF(); ++i) {
    potential.setValue(i, 0, i);
  }

  resultWriter.setQTensor(qTensor);
  resultWriter.setPotential(potential);

  // WHEN the result is written with both solution vectors assigned.
  resultWriter.writeResult(geom, simulationState);

  // THEN the header describes both fields and their payload ranges.
  const fs::path resultFilePath = outputDir / "result00000005.yarf"; // iteration padded to 8 digits
  REQUIRE(fs::exists(resultFilePath));
  REQUIRE(fs::file_size(resultFilePath) > 0);

  // AND: The header of the output file should contain the expected values
  YarfResultReader resultReader;
  Yarf::YarfHeader yarfHeader = YarfResultReader::readYarfHeader(resultFilePath.string());

  REQUIRE(yarfHeader.fields.size() == 2); // q-tensor and potential fields
  const auto &qTensorField = yarfHeader.fields.at(0);
  REQUIRE(qTensorField.name == "q_tensor");
  REQUIRE(qTensorField.association == Yarf::FieldAssociation::LcNodes);
  REQUIRE(qTensorField.datatype == Yarf::DataType::Float64);
  const std::vector<std::string> qTensorComponents = {"Qxx", "Qxy", "Qxz", "Qyy", "Qyz"};
  REQUIRE(qTensorField.components == qTensorComponents);
  REQUIRE(qTensorField.length == static_cast<std::uint64_t>(qTensor.getnDoF()) * 5U);
  REQUIRE(qTensorField.byte_offset == 0);

  const auto &potentialField = yarfHeader.fields.at(1);
  REQUIRE(potentialField.name == "electric_potential");
  REQUIRE(potentialField.association == Yarf::FieldAssociation::AllNodes);
  REQUIRE(potentialField.datatype == Yarf::DataType::Float64);
  REQUIRE(potentialField.components == std::vector<std::string>{"phi"});
  REQUIRE(potentialField.length == static_cast<std::uint64_t>(potential.getnDoF()));
  REQUIRE(potentialField.byte_offset == qTensorField.length * sizeof(double));
  REQUIRE(yarfHeader.format.payload_size_bytes ==
          qTensorField.length * sizeof(double) + potentialField.length * sizeof(double));

  // check that the file size in bytes equals the offset to the start of the payload plus the payload size
  const std::uint64_t expectedFileSize = yarfHeader.payloadStartOffset + yarfHeader.format.payload_size_bytes;
  const std::uint64_t actualFileSize = fs::file_size(resultFilePath);
  REQUIRE(actualFileSize == expectedFileSize);

  auto samples = resultReader.read(resultFilePath.string(), 0);
  REQUIRE(samples.size() == static_cast<std::size_t>(qTensor.getnDoF()));

  // THEN the written Q-tensors round-trip to the original traceless representation.
  for (idx node = 0; node < qTensor.getnDoF(); ++node) {
    REQUIRE(samples[node].tensor.t1() == Approx(qTensor.getValue(node, 0)).margin(1e-12));
    REQUIRE(samples[node].tensor.t2() == Approx(qTensor.getValue(node, 1)).margin(1e-12));
    REQUIRE(samples[node].tensor.t3() == Approx(qTensor.getValue(node, 2)).margin(1e-12));
    REQUIRE(samples[node].tensor.t4() == Approx(qTensor.getValue(node, 3)).margin(1e-12));
    REQUIRE(samples[node].tensor.t5() == Approx(qTensor.getValue(node, 4)).margin(1e-12));
  }

  // AND the potential payload contains each node's scalar value in mesh order.
  // Note that the reader does not return the potential values, so we read them directly from the binary file.
  std::ifstream payloadStream(resultFilePath, std::ios::binary);
  REQUIRE(payloadStream.is_open());
  payloadStream.seekg(static_cast<std::streamoff>(yarfHeader.payloadStartOffset + potentialField.byte_offset),
                      std::ios::beg);
  REQUIRE(payloadStream.good());
  for (idx node = 0; node < potential.getnDoF(); ++node) {
    double writtenPotential = 0.0;
    payloadStream.read(reinterpret_cast<char *>(&writtenPotential), sizeof(writtenPotential));
    REQUIRE(payloadStream.good());
    REQUIRE(writtenPotential == Approx(potential.getValue(node)).margin(1e-12));
  }
}