# YARF implementation notes

This document describes the qlc3d YARF implementation through the WP5 read-seam milestone, covering the Q-tensor and electric-potential field contracts, the qlc3d `SolutionVector` adapters, and the explicit read/write seams that remain opt-in for runtime integration.

## Scope

The implementation provides a small, testable YARF I/O package under `qlc3d/includes/io/` and `qlc3d/src/io/`:

- `YarfHeader` and supporting plain-data types
- `YarfJsonCodec` for serializing and parsing the JSON header
- `YarfValidator` for format integrity and descriptor validation
- `IYarfReader` / `IYarfWriter` interfaces
- `YarfFileReader` / `YarfFileWriter` concrete file implementations
- `YarfQlc3dAdapter` for qlc3d `SolutionVector` conversion and field validation

This remains a standalone library rather than a runtime save/load path, keeping the format logic isolated and testable.

## Supported subset

The implementation supports:

- `format.version` with the canonical value `"1.0"`
- `format.endianness` with the canonical value `"little"`
- `simulation` metadata as a required object
- `metadata` as an arbitrary JSON object that is preserved verbatim when read back
- header-only files with zero payload, and binary payload files containing the supported Q-tensor and electric-potential fields

The supported field contracts are:

- `q_tensor`
  - `association = "lc_nodes"`
  - `datatype = "float64"`
  - `components = ["Qxx", "Qxy", "Qxz", "Qyy", "Qyz"]`
  - `length = npLC * 5`
  - binary payload stored as contiguous AoS values in the component order above
- `electric_potential`
  - `association = "all_nodes"`
  - `datatype = "float64"`
  - `components = ["phi"]`
  - `length = num_nodes`
  - binary payload stored as one scalar per node in mesh order

The accepted field associations are still the general YARF subset:

- `lc_nodes`
- `all_nodes`
- `lc_elements`
- `all_elements`

The accepted datatypes are:

- `float64`
- `float32`

## Validation rules

The validator implements the header-first checks required by the format specification, plus the Q-tensor-specific invariants added in WP2:

- magic header is exactly `YAR1`
- JSON length is read as a little-endian `uint32`
- required top-level keys are present: `format`, `simulation`, `metadata`, `fields`
- actual file size must equal `8 + json_length + payload_size_bytes`
- a `q_tensor` field must use `lc_nodes`, `float64`, and the exact five component labels above
- an `electric_potential` field must use `all_nodes`, `float64`, and the exact scalar label `phi`
- q-tensor `length` must equal `npLC * 5`
- potential `length` must equal `num_nodes`
- payload ranges must stay within the file's declared payload size

## File layout

The writer emits:

1. 4-byte magic (`YAR1`)
2. 4-byte little-endian JSON length
3. serialized UTF-8 JSON header
4. optional raw binary payload, written as little-endian IEEE-754 `float64` values in AoS order

The qlc3d adapter flattens a `SolutionVector` of shape `(npLC, 5)` into the contiguous vector `[Qxx, Qxy, Qxz, Qyy, Qyz]` for each node, and a scalar potential `SolutionVector` of shape `(num_nodes, 1)` into the contiguous vector `[phi]` for each node.

## Integration seams

The supported write-side seam is:

- `qlc3d::YarfSnapshotWriter`
- `qlc3d::writeYarfSnapshot(...)`

The write helper accepts the runtime `Geometry`, potential `SolutionVector`, Q-tensor `SolutionVector`, and `SimulationState`, then emits a file named `result` + 8-digit iteration + `.yarf` under the requested output directory. This uses the canonical YARF naming convention, e.g. `result00000013.yarf`, and does not silently inject the regular output writer registry.

The supported read-side seam is:

- `qlc3d::YarfSnapshotReader`
- `qlc3d::readYarfSnapshot(...)`

The read helper resolves either an explicit snapshot path or the canonical state-based name under the configured directory, then parses the header and materialises `q_tensor` and `electric_potential` `SolutionVector` values when those fields are present. This is intentionally a qlc3d-facing adapter layer and not yet wired into the main runtime file-discovery or result-output registry.

All YARF-layer exceptions raised during the integration boundary are translated to `std::runtime_error` with the original YARF message preserved, so qlc3d callers can handle them consistently without importing the YARF-specific exception hierarchy.

The runtime wiring remains deferred by design. The save-side and read-side seams are intentionally explicit and opt-in so the default runtime output registry and restart flow stay stable while the YARF format continues to be exercised in isolated tests and staged integration work.

## Future extension

Later work packages can add additional supported fields and qlc3d integration seams, but the current implementation is intentionally explicit and strict around the supported q-tensor subset.
