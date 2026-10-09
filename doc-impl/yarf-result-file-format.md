# YARF (Yet Another Result Format) Specification

**Extension:** `.yarf`  
**Purpose:** A hybrid text-binary file format designed for storing finite element simulation snapshots, combining the extreme extensibility of JSON for metadata with the high-performance memory-mapping capabilities of raw binary payloads.

---

## 1. High-Level File Layout

A `.yarf` file consists of three contiguous blocks of data:

| Section | Size | Description |
| :--- | :--- | :--- |
| **Magic Header** | 8 bytes | Format identifier and JSON block size. |
| **JSON Header** | Variable | UTF-8 encoded JSON text containing all metadata and payload schemas. |
| **Binary Payload**| Variable | Contiguous raw binary arrays described by the JSON header. |

The total file size must exactly equal: `8 + (JSON Length) + (Payload Size)`.

---

## 2. Magic Header (8 Bytes)

The file begins with an 8-byte header to identify the format and dictate the size of the subsequent text block.

*   **Bytes 0–3:** Magic string `"YAR1"` (ASCII encoded). Hexadecimal: `0x59 0x41 0x52 0x31`.
*   **Bytes 4–7:** A 32-bit unsigned integer (`uint32`) in **little-endian** byte order. This integer represents the exact length of the JSON Header block in bytes.

---

## 3. JSON Header

Immediately following the 8-byte magic header is the JSON header.
*   It must be exactly the length specified in Bytes 4–7 of the Magic Header.
*   It must be valid UTF-8 encoded text.
*   It is **not** null-terminated.

### 3.1. Top-Level Structure

The JSON root must be an object containing the following keys:

```json
{
  "format": {},
  "simulation": {},
  "metadata": {},
  "fields": []
}
```

### 3.2. `format` Object
Describes the core format specification and physical layout of the binary payload.

*   `version` *(string)*: The version of the YARF specification this file adheres to (e.g., `"1.0"`).
*   `endianness` *(string)*: The byte order of the binary payload. Typical values are `"little"` or `"big"`. All binary data in the payload must adhere to this.
*   `payload_size_bytes` *(integer)*: The total size of the Binary Payload block in bytes. Used by readers to verify file integrity and detect incomplete writes.

### 3.3. `simulation` Object
Core invariant metadata for the simulation snapshot.

*   `mesh_file` *(string)*: Relative or absolute path to the external mesh geometry file (e.g., GiD mesh).
*   `iteration` *(integer)*: The simulation step or iteration number.
*   `time` *(number)*: The physical simulation time at this snapshot.
*   `S0` *(number)*: The equilibrium order parameter.
*   `num_nodes` *(integer)*: Total number of nodes in the referenced mesh.
*   `num_elements` *(integer)*: Total number of elements in the referenced mesh.
*   `npLC` *(integer)*: The number of nodes where the Liquid Crystal (LC) material is defined. **Domain Invariant:** LC nodes are strictly ordered such that they occupy indices `0` through `npLC - 1` in the external mesh.

### 3.4. `metadata` Object
An open-ended object for arbitrary, application-specific data. Readers should ignore unrecognized keys in this object. Examples include debug strings, solver convergence flags, or custom tolerances.

### 3.5. `fields` Array
An array of objects, where each object describes a single continuous array of data located in the Binary Payload section.

*   `name` *(string)*: The identifier of the field (e.g., `"q_tensor"`, `"electric_potential"`).
*   `association` *(string)*: The domain over which the data is defined. Standard values:
    *   `"lc_nodes"`: Data is defined only on the first `npLC` nodes.
    *   `"all_nodes"`: Data is defined on all `num_nodes` nodes.
    *   `"all_elements"`: Data is defined on all `num_elements` elements.
*   `datatype` *(string)*: The numeric type of the binary data. Standard values: `"float64"`, `"float32"`.
*   `components` *(array of strings)*: The names of the dimensions for each entry. (e.g., `["Qxx", "Qxy", "Qxz", "Qyy", "Qyz"]` for a 5-component tensor, or `["phi"]` for a scalar).
*   `length` *(integer)*: The total number of scalar values in the flattened 1D array. This is strictly equal to the number of entities defined by the `association` multiplied by the number of `components`. This allows readers to allocate memory immediately without inferring sizes.
*   `byte_offset` *(integer)*: The starting byte position of this field's data **relative to the start of the Binary Payload block**. (e.g., an offset of `0` means the data starts immediately after the JSON block).

---

## 4. Binary Payload

The binary payload begins exactly at byte index `8 + (JSON Length)` of the file.

*   Data is stored as raw, tightly packed binary arrays as defined by the `fields` array.
*   There is no enforced padding or alignment between fields, though writers may optionally inject padding for memory-alignment purposes by adjusting the `byte_offset` of the subsequent field.
*   Multidimensional data (like the 5-component Q-tensor) is interleaved/AoS (Array of Structures). For example, node 0's 5 components are followed by node 1's 5 components in contiguous memory.

---

## 5. File Integrity Validation

To ensure a `.yarf` file is uncorrupted and fully written, a reader should perform the following checks upon opening:
1.  Read the first 4 bytes. Assert they equal `"YAR1"`.
2.  Read the next 4 bytes as a little-endian `uint32` to get `json_length`.
3.  Calculate the absolute file size on disk (`actual_file_size`).
4.  Parse the JSON header and extract `format.payload_size_bytes`.
5.  Assert that `actual_file_size == 8 + json_length + payload_size_bytes`.

---

## Appendix A: Example JSON Header

```json
{
  "format": {
    "version": "1.0",
    "endianness": "little",
    "payload_size_bytes": 46000000
  },
  "simulation": {
    "mesh_file": "geometry/chamber_v2.msh",
    "iteration": 1500,
    "time": 0.05,
    "S0": 0.53,
    "num_nodes": 1000000,
    "npLC": 750000,
    "num_elements": 980000
  },
  "metadata": {
    "description": "Mid-point snapshot, transitioning to steady-state.",
    "solver_converged": false
  },
  "fields": [
    {
      "name": "q_tensor",
      "association": "lc_nodes",
      "datatype": "float64",
      "components": ["Qxx", "Qxy", "Qxz", "Qyy", "Qyz"],
      "length": 3750000,
      "byte_offset": 0
    },
    {
      "name": "electric_potential",
      "association": "all_nodes",
      "datatype": "float64",
      "components": ["phi"],
      "length": 1000000,
      "byte_offset": 30000000
    },
    {
      "name": "element_debug_scalar",
      "association": "all_elements",
      "datatype": "float32",
      "components": ["error_estimate"],
      "length": 980000,
      "byte_offset": 38000000
    }
  ]
}
```