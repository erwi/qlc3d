#ifndef PROJECT_QLC3D_YARFRESULTWRITER_H
#define PROJECT_QLC3D_YARFRESULTWRITER_H
#include <fstream>
#include <json.hpp>

#include "result-output.h"

# pragma region "YARF Metadata"
namespace Yarf {
  /**
   * Supported YARF field associations.
   */
  enum class FieldAssociation {
    LcNodes,
    AllNodes,
    LcElements,
    AllElements
  };
  NLOHMANN_JSON_SERIALIZE_ENUM(FieldAssociation, {
    {FieldAssociation::LcNodes, "lc_nodes"},
    {FieldAssociation::AllNodes, "all_nodes"},
    {FieldAssociation::LcElements, "lc_elements"},
    {FieldAssociation::AllElements, "all_elements"}
  })
  /**
   * Supported binary scalar datatypes for YARF payload fields.
   */
  enum class DataType {
    Float64,
    Float32
  };
  NLOHMANN_JSON_SERIALIZE_ENUM(DataType, {
    {DataType::Float64, "float64"},
    {DataType::Float32, "float32"}
  })

  /**
   * Format-level metadata used to validate and decode a YARF file.
   */
  struct FormatBlock {
    std::string version = "1.0";
    std::string endianness = "little";
    std::uint64_t payload_size_bytes = 0;
  };
  NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FormatBlock, version, endianness, payload_size_bytes)

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
  NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SimulationBlock, mesh_file, iteration, time, S0, num_nodes, num_elements, npLC)

  /**
   * Describes a single field inside a YARF payload.
   *
   * @param name Logical field name such as q_tensor or electric_potential.
   * @param association Domain on which the field is defined.
   * @param datatype Underlying scalar datatype stored in the payload.
   * @param components Component labels in AoS order.
   * @param length Total scalar values count in the field (not number of bytes).
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
  NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FieldDescriptor, name, association, datatype, components, length, byte_offset)

  /**
   * Full YARF JSON header plus parsed field descriptors.
   *
   * @param payloadStartOffset Absolute byte offset from the start of the file to the start of the
   * payload (magic number + length prefix + JSON length). Populated only by the reader and not
   * serialized to JSON.
   */
  struct YarfHeader {
    FormatBlock format;
    SimulationBlock simulation;
    nlohmann::json metadata = nlohmann::json::object();
    std::vector<FieldDescriptor> fields;
    std::int64_t payloadStartOffset = 0;
  };
  NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(YarfHeader, format, simulation, metadata, fields)
}
# pragma endregion



/**
 * Writes simulation snapshots in the YARF header-and-binary-payload format.
 */
class YarfResultWriter : public ResultFormatWriter {
  static constexpr int UNINITIALIZED = -1;

protected:
  double S0;
  std::string meshName;
  int lastMeshNumber;
public:
  /**
   * Creates a YARF result writer.
   * @param outputDir Directory in which result and mesh files are written.
   * @param meshName Base name of the mesh file associated with results.
   * @param S0 Equilibrium order parameter recorded in the result header.
   */
  YarfResultWriter(const std::filesystem::path &outputDir, const std::string &meshName, double S0) :
          ResultFormatWriter(outputDir), S0(S0), meshName(meshName), lastMeshNumber(UNINITIALIZED) {};

  /** @return The display name of this result format. */
  [[nodiscard]] const std::string formatName() const override { return "Yarf"; };

  /**
   * Writes a YARF snapshot and any mesh file required for the current mesh number.
   * @param geom Geometry and mesh represented by the snapshot.
   * @param simulationState Iteration, time, and mesh-number metadata for the snapshot.
   */
  void writeResult(const Geometry &geom, const SimulationState &simulationState) override;

private:
  /**
   * Writes the mesh file when the current mesh number has not yet been written.
   * @param geom Geometry to write.
   * @param simulationState State containing the current mesh number.
   */
  void writeMeshIfRequired(const Geometry &geom, const SimulationState &simulationState);

  /**
   * Writes the YARF header to the output stream.
   * @param out Output stream to write to.
   * @param meshFileName Name of the mesh file associated with the simulation.
   * @param simulationState State metadata stored in the header.
   * @param geom Geometry metadata stored in the header.
   * @param S0 Scalar order parameter value for the simulation.
   * @return Header containing the field descriptors written to the stream.
   */
  Yarf::YarfHeader writeYarfHeader(std::ofstream &out,
    const std::string &meshFileName,
    const SimulationState &simulationState,
    const Geometry &geom,
    const double S0);

  /**
   * Writes the Q-tensor payload described by the header.
   * @param out Output stream positioned at the YARF file.
   * @param header Header containing the field descriptor and payload offset.
   * @param payloadStart Start of the binary payload in the output stream.
   */
  void writeQTensorPayload(std::ofstream &out, const Yarf::YarfHeader &header, std::streampos payloadStart) const;

  /**
   * Writes the electric-potential payload described by the header.
   * @param out Output stream positioned at the YARF file.
   * @param header Header containing the field descriptor and payload offset.
   * @param payloadStart Start of the binary payload in the output stream.
   */
  void writePotentialPayload(std::ofstream &out, const Yarf::YarfHeader &header, std::streampos payloadStart) const;

  /**
   * Creates descriptors for the solution vectors currently assigned to this writer.
   * @return Field descriptors with contiguous payload offsets.
   */
  std::vector<Yarf::FieldDescriptor> createFieldDescriptors() const;

};


#endif //PROJECT_QLC3D_YARFRESULTWRITER_H
