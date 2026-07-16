# Initial LC Orientation Loading

## Overview

The initial LC orientation (Q-tensor) assigned to every LC mesh node before the simulation starts is built up in
three ordered stages inside `initialiseLcSolutionVector` (`qlc3d/src/inits.cpp`):

1. **Box initialisation** — `InitialVolumeOrientation::setVolumeQ` sets analytic orientation within configured
   cuboid regions (`Normal`, `Random`, `Hedgehog` types).
2. **File-based loading** (this document) — if the legacy `loadQ` settings key is set, the Q-tensor is overwritten
   from a file on disk.
3. **Surface/alignment conditions** — `setSurfacesQ` and `q.initialiseLcBoundaries` apply boundary conditions last.

This ordering must not change: file-based loading only overwrites node values already initialised by the box
stage, and surface conditions are always applied last.

This document describes the file-based loading abstraction (stage 2). Both the new `loadOrientation` settings key
and the legacy `loadQ` settings key route through this abstraction. Two file formats are supported: LCView
text/binary (mesh-matched: no coordinates, exact node-count/order match required) and director CSV (point-cloud:
coordinates + director, nearest-neighbor assignment onto mesh nodes).

---

## Source Files

| File | Responsibility |
|------|---------------|
| `qlc3d/includes/io/orientation-sample.h` | `qlc3d::OrientationSample` — common intermediate representation: a `TTensor` plus an optional 3D `location`. |
| `qlc3d/includes/io/orientation-reader.h` | `qlc3d::OrientationReader` interface, `LcViewTextReader`/`LcViewBinaryReader` implementations, `createLcViewReader` LCView content-sniffing dispatch, and `createOrientationReader` top-level extension-based dispatch. |
| `qlc3d/src/io/orientation-reader.cpp` | Implementation of the LCView text/binary readers (migrated from the former `ResultIO::readTextLcViewResultFile`/`readBinaryLcViewResultFile`) and the binary-marker-based format sniff (migrated from `ResultIO::ReadResult`); `createOrientationReader`'s `.csv` extension check. |
| `qlc3d/includes/io/director-csv-reader.h` | `qlc3d::DirectorCsvReader` — `OrientationReader` implementation for the director CSV point-cloud format. |
| `qlc3d/src/io/director-csv-reader.cpp` | CSV header parsing (case-insensitive, any column order), required/optional column validation, per-row parsing/normalization, and `RUNTIME_ERROR` on any malformed input. |
| `qlc3d/includes/io/orientation-assignment.h` | `qlc3d::OrientationAssignmentStrategy` interface, `ExactOrderAssignment`, and `NearestNeighborAssignment` implementations. |
| `qlc3d/src/io/orientation-assignment.cpp` | `ExactOrderAssignment::assign` — validates sample count against `SolutionVector::getnDoF()` and copies tensors onto LC nodes in order. `NearestNeighborAssignment::assign` — brute-force nearest-sample lookup per mesh node. |
| `qlc3d/includes/io/orientation-loader.h` | `qlc3d::loadInitialOrientation` — top-level orchestration entry point. |
| `qlc3d/src/io/orientation-loader.cpp` | Resolves `loadQ`/`loadOrientation` precedence (fatal error if both set, deprecation warning if `loadQ` used), dispatches to a reader via `createOrientationReader`, applies `StretchVector` scaling to point-cloud sample locations, and assigns samples onto the `SolutionVector` using the assignment strategy appropriate for the format (`producesLocations()`). |
| `qlc3d/src/inits.cpp` (`initialiseLcSolutionVector`) | Orchestration: calls `qlc3d::loadInitialOrientation(simu, S0, geom.getCoordinates(), q)`. |

## Architecture

```
file on disk --[OrientationReader::read]--> vector<OrientationSample> --[OrientationAssignmentStrategy::assign]--> SolutionVector
```

- **`OrientationSample`**: one parsed orientation value — a `qlc3d::TTensor` plus `std::optional<Vec3> location`.
  `location` is `std::nullopt` for mesh-matched formats (no coordinates in the file, e.g. LCView text/binary).
  Whether a format's samples carry locations is intrinsic to the format, exposed via
  `OrientationReader::producesLocations()`.

- **`OrientationReader`**: parses a file on disk into `vector<OrientationSample>`, without touching
  `SolutionVector`/mesh directly — this keeps readers independently unit-testable with plain files.
  - `LcViewTextReader`: parses the legacy LCView text result format (3 header lines, then
    `id nx ny nz v S S` rows). Stops at the first all-zero director row (marks "end of LC region"), matching
    legacy behavior. Does **not** validate row count against any mesh size — that is the assignment strategy's
    job.
  - `LcViewBinaryReader`: parses the legacy LCView binary result format (5 header lines, then
    `S0 np nsol` header, then `np` records of `nsol` floats each, of which the first 5 are the Q-tensor
    components `q1 q2 q3 q5 q4` in that on-disk order). Reads exactly `np` samples (the file's own declared
    count), not a target mesh size.
  - `createLcViewReader(fileName)`: checks the file exists, then sniffs the first 5 lines for the
    `"RAW FLOAT TRI"` marker to decide between `LcViewBinaryReader` and `LcViewTextReader` — same sniffing logic
    as the former `ResultIO::ReadResult`.
  - `DirectorCsvReader`: parses the director CSV point-cloud format. First non-empty line is the header: split on
    `,`, trimmed, matched case-insensitively against the recognized column names `x, y, z, nx, ny, nz` (required)
    and `s` (optional). Any header column not in this set, or any duplicate/missing required column, is a
    `RUNTIME_ERROR` (fail loud — this also surfaces typos in column names, e.g. a misspelled `nx`, which would
    otherwise silently become both a missing-required-column error and a silently-ignored extra column). Each data
    row is split on `,` and parsed positionally using the header-derived column indices; a zero-length `(nx, ny,
    nz)` vector is a `RUNTIME_ERROR` (unlike LCView text, there is no "all-zero row marks end of data" convention
    for CSV). Non-zero director vectors are normalized. `S` uses the row's own value if the column is present,
    otherwise the `s0` parameter passed to `read()`. An empty or header-only file is a `RUNTIME_ERROR`.
  - `createOrientationReader(fileName)`: the top-level dispatcher used by `loadInitialOrientation`. Routes `.csv`
    (case-insensitive extension match) to `DirectorCsvReader`; everything else falls back to `createLcViewReader`'s
    content-sniffing.

- **`OrientationAssignmentStrategy`**: writes parsed samples onto the LC nodes of a `SolutionVector`.
  - `ExactOrderAssignment`: for mesh-matched formats (samples carry no location). Throws
    (`RUNTIME_ERROR`, i.e. `std::runtime_error`) if `samples.size() != q.getnDoF()`, with the message
    `"The loaded result file size {actual} does not match the expected size {expected}"` — same wording as the
    legacy behavior. On success, copies `samples[i].tensor` onto `q` node `i` for all LC nodes, in file order.
  - `NearestNeighborAssignment`: for point-cloud formats (samples carry a `location`). For each LC mesh node,
    brute-force scans all samples and assigns the tensor of the closest one by Euclidean distance
    (`Vec3::distanceSquared`). On exact ties, the first sample encountered in list order wins (deterministic).
    Throws `RUNTIME_ERROR` if `samples` is empty, or if any sample lacks a `location` (programmer error — this
    strategy must only ever be invoked with located samples). Callers are responsible for pre-scaling sample
    locations by `StretchVector` into the same coordinate space as `meshCoordinates` before calling `assign()` —
    this class does not itself apply any scaling.

- **Orchestration** (`qlc3d::loadInitialOrientation`, `qlc3d/src/io/orientation-loader.cpp`, called from
  `initialiseLcSolutionVector` in `qlc3d/src/inits.cpp`): resolves which of `simu.getLoadQ()` /
  `simu.getLoadOrientation()` is set:
  - If both are set: `RUNTIME_ERROR` (fatal) — ambiguous configuration, since `loadQ` is deprecated in favor of
    `loadOrientation`.
  - If only `loadQ` is set: logs a `Log::warn` deprecation message suggesting `loadOrientation`, then proceeds
    using that file.
  - If only `loadOrientation` is set: proceeds using that file.
  - If neither is set: does nothing.
  - Once a file is resolved, a reader is selected via `createOrientationReader`, and samples are parsed via
    `reader->read(file, s0)`. If `reader->producesLocations()` is true (point-cloud format, e.g. director CSV),
    each sample's location is first scaled component-wise by `simu.getStretchVector()` (mapping the file's own
    unstretched coordinate space into the same stretched space as `meshCoordinates`), then
    `NearestNeighborAssignment` writes the samples onto `q`. Otherwise (mesh-matched format, e.g. LCView),
    `ExactOrderAssignment` writes them directly. This runs strictly between box initialisation and
    surface/alignment condition application, per the ordering requirement above.

## Director CSV Format Summary

- Header row, comma-delimited, case-insensitive column names, any column order.
- Required columns: `x, y, z, nx, ny, nz`.
- Optional column: `s` — if omitted, every row uses the `s0` (equilibrium order parameter) passed by the caller.
- Director vectors need not be pre-normalized; the reader normalizes them.
- A single data row is valid and, combined with `NearestNeighborAssignment`, applies one uniform orientation to
  every mesh node — the mechanism behind "single point applies everywhere".
- All malformed input (missing/unrecognized/duplicate columns, wrong field count, unparsable numbers, zero-length
  director, empty/header-only file) is fatal (`RUNTIME_ERROR`) — there are no silent fallbacks.

## Tests

- `tests/cpp/io/orientation-reader-tests.cpp` — `LcViewTextReader`/`LcViewBinaryReader`/`createLcViewReader`,
  `ExactOrderAssignment` unit tests.
- `tests/cpp/io/orientation-assignment-tests.cpp` — `NearestNeighborAssignment` unit tests (synthetic samples).
- `tests/cpp/io/director-csv-reader-tests.cpp` — `DirectorCsvReader` parsing/validation unit tests, and
  `createOrientationReader` extension-based dispatch tests.
- `tests/cpp/io/orientation-loader-tests.cpp` — end-to-end `loadInitialOrientation` tests covering
  `loadQ`/`loadOrientation` precedence, LCView regression, and director CSV end-to-end scenarios (single point
  applies everywhere, multiple points with nearest-neighbor selection, and `StretchVector` scaling of CSV
  locations before matching).
