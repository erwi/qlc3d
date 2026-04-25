# Quadratic Element Support — State of the Codebase

> Last updated: April 2026

This document assesses the current state of support for quadratic finite elements (10-node tetrahedra `TET10`, 6-node triangles `TRI6`) compared to the fully-supported linear baseline (4-node tetrahedra `TET4`, 3-node triangles `TRI3`).

---

## Table of Contents

1. [Background and Design Intent](#1-background-and-design-intent)
2. [Linear Element Baseline](#2-linear-element-baseline)
3. [What Is Implemented for Quadratic Elements](#3-what-is-implemented-for-quadratic-elements)
4. [Incomplete or Blocked Work](#4-incomplete-or-blocked-work)
5. [Test Coverage](#5-test-coverage)
6. [Technical Debt](#6-technical-debt)
7. [Modularity and Coupling Concerns](#7-modularity-and-coupling-concerns)
8. [Recommended Next Steps](#8-recommended-next-steps)

---

## 1. Background and Design Intent

The codebase was originally written for linear tetrahedral finite elements. Work to introduce quadratic (serendipity/Lagrange) TET10 and TRI6 elements has been partially completed. The core shape-function machinery and mesh conversion utilities are done; the principal remaining work is removing explicit `NotYetImplementedException` guards in the solver entry points and fixing a dependency-chain block in the regular-grid / spatial-lookup subsystem.

A critical architectural decision was made in `prepareGeometry()` in `inits.cpp`: **the mesh is always converted to quadratic (TET10/TRI6) in memory at geometry preparation time**, regardless of what was loaded from disk. A linear input mesh is automatically promoted via `convertLinearMeshDataToQuadratic()`. This means all downstream code only ever sees TET10/TRI6, but the solver entry points that sit between `prepareGeometry()` and the already-quadratic-aware assembly routines still throw before those routines are reached.

---

## 2. Linear Element Baseline

Linear elements (`ElementType::LINEAR_TETRAHEDRON`, `ElementType::LINEAR_TRIANGLE`) are fully supported end-to-end.

| Area | Status |
|---|---|
| `ElementType` enum, `getNodesPerElement()`, `getElementOrder()` | ✅ Complete |
| `Mesh` data model, node storage, determinants, surface normals | ✅ Complete |
| GiD mesh reader (4-node tet, 3-node tri) | ✅ Complete |
| Gmsh reader — element types 2 (TRI3), 4 (TET4) | ✅ Complete |
| `TetShapeFunction(1)`, `TriShapeFunction(1)` shape functions | ✅ Complete |
| Gauss quadrature points (`Keast0`, `Keast4`, `Keast8`, `Tri4thOrder`) | ✅ Complete |
| Q-tensor (LC) solver assembly | ✅ Complete |
| Potential solver assembly | ✅ Complete |
| Adaptive mesh refinement (red-green tet splitting) | ✅ Complete |
| Regular spatial grid (`makeRegularGrid`) | ✅ Complete |
| Q-tensor interpolation (`interpolateQTensor`) | ✅ Complete |
| Result output — VTK, LcView, result files | ✅ Complete |

---

## 3. What Is Implemented for Quadratic Elements

### 3.1 Data Model (`mesh/mesh.h`)

- `ElementType::QUADRATIC_TRIANGLE = 3` (6 nodes) and `ElementType::QUADRATIC_TETRAHEDRON = 4` (10 nodes) are defined.
- `getNodesPerElement()` returns 6 / 10 correctly.
- `getElementOrder()` returns 2 for both types.
- `Mesh::setElementData()` validates node-count consistency for quadratic types.
- Per-element loops throughout `Mesh` use `getnNodes()` dynamically, so they correctly include mid-edge nodes for TET10/TRI6.
- `calculateDeterminants3D()` and `calculateSurfaceNormals()` correctly use only the 4 (or 3) corner nodes, which is exact for straight-edged elements.
- `gen_p_to_elem()`, `listFixLCSurfaceNodes()`, `findElectrodeSurfaceNodes()`, `findNodesWhere()` all iterate `getnNodes()` and therefore include mid-edge nodes.

### 3.2 Shape Functions (`fe/gaussian-quadrature.h`)

Both quadratic shape-function families are fully implemented:

- **`TetShapeFunction(2)`** — `initialiseQuadraticTet()`: standard TET10 serendipity basis; computes all 10 `sh`, `shR`, `shS`, `shT` values; `initialiseElement()` maps natural → global coordinates via the full 10-node Jacobian.
- **`TriShapeFunction(2)`** — `setQuadraticTrianglePoints()`: standard TRI6 basis; computes all 6 `sh`, `shR`, `shS` values.
- Public accessors `N(i)`, `Nx(i)`, `Ny(i)`, `Nz(i)`, `sample()`, `sampleQ()`, `sampleQX/Y/Z()` etc. loop over `nodesPerElement` and work for both 4 and 10 nodes.

### 3.3 Gmsh Mesh I/O

- The Gmsh reader (`io/gmsh-read.cpp`) reads element types `9` (TRI6) and `11` (TET10).
- `RawMeshData` carries `elementOrder_` (1 or 2) from the file.
- Node count arithmetic is computed from `elementOrder_` throughout the reader.

### 3.4 Geometry Preparation and Mesh Conversion (`inits.cpp`, `mesh/element-split-convert.h`)

The full conversion pipeline is implemented:

| Function | Purpose |
|---|---|
| `convertLinearMeshDataToQuadratic()` | TET4→TET10 and TRI3→TRI6; mid-edge nodes at edge midpoints; respects shared edges between elements. Uses sparse IRC matrix for edge→node-index map. |
| `reorderQuadraticTetNodeOrder()` | Fixes Gmsh's TET10 node ordering to match qlc3d's internal convention (Gmsh swaps nodes 8 and 9). Detected at load time by proximity of actual coordinates to expected mid-edge positions. |
| `recombineLinearisedMeshToQuadratic()` | Heuristically re-constitutes TET10 from a mesh that was previously split 1 TET10 → 8 TET4. |
| `splitQuadraticTetrahedronToLinear()` | TET10 → 8 × TET4 (used for LcView output). |
| `splitQuadraticTriangleToLinear()` | TRI6 → 4 × TRI3. |
| `recombineLinearTetsToQuadratic()` / `recombineLinearTrianglesToQuadratic()` | Inverse of the two split functions above. |

`prepareGeometry()` always converts the in-memory mesh to quadratic before returning, so all downstream simulation code only ever sees TET10/TRI6.

### 3.5 Boundary-Tet Node Reordering (`fe/fe-util.h`)

`reorderQuadraticBoundaryTetNodes()` reorders the 10 nodes of a TET10 so that the triangle face nodes (3 corners + 3 mid-edge) come first — required for Neumann surface integrals. This exists alongside the linear-only `reorderBoundaryTetNodes()`.

### 3.6 LC Solver (Q-Tensor) Assembly (`lc/lc-solver.cpp`)

**Functionally ready for quadratic elements.** The assembly is fully parameterised by `npe = shapes.getNumPointsPerElement()`:

- `assembleLocalVolumeMatrix()` works for both 4 and 10 nodes.
- `assembleMatrixSystemVolumeTerms()` detects element order, creates `TetShapeFunction(elementOrder)`, and uses `Keast8` (8th-order integration, appropriate for quadratic basis products).
- `assembleMatrixSystemWeakAnchoring()` detects element order, creates `TriShapeFunction(elementOrder)` with `Tri4thOrder`.
- A `#ifndef NDEBUG` block in `assembleLocalVolumeMatrix()` verifies TET10 mid-edge node positions at runtime.

### 3.7 Potential Solver Assembly (`potential/potential-solver.cpp`)

The assembly routines are quadratic-aware:

- `assembleVolume()` creates `TetShapeFunction(elementOrder)` with `Keast4`.
- `assembleNeumann()` branches on `elementType` to call `reorderQuadraticBoundaryTetNodes` vs `reorderBoundaryTetNodes`.

However, the public entry point `solvePotential()` throws a `NotYetImplementedException` **before** reaching this code. See §4.

### 3.8 VTK Output (`io/vtkiofun.cpp`)

The VTK writer correctly outputs VTK cell type 24 (quadratic tet, 10 nodes) vs cell type 10 (linear tet, 4 nodes). It is gated upstream by the `ResultOutput::writeResults()` check (see §4).

### 3.9 LcView Output (`io/lcview-result-output.cpp`)

`writeQuadraticTetrahedra()` is implemented; it splits each TET10 into 8 TET4 for the LcView binary format (which does not support TET10 natively). Gated upstream by the same check described in §4.

---

## 4. Incomplete or Blocked Work

### 4.1 `NotYetImplementedException` Guards

Four explicit runtime-throw sites block use of quadratic elements in the full simulation pipeline:

| File | Entry point | What it prevents |
|---|---|---|
| `potential/potential-solver.cpp` | `PotentialSolver::solvePotential()` | All electrode potential solves. The assembly code **directly below this guard** already supports quadratic. |
| `io/result-output.cpp` | `ResultOutput::writeResults()` | All result output — VTK, LcView, and result files — even though both underlying writers already support TET10. |
| `geometry.cpp` | `Geometry::makeRegularGrid()` | Regular spatial grid construction. This cascades (see §4.2). |
| `autorefinement.cpp` | `autoref()` | All adaptive mesh refinement. |

Removing the first two guards (potential solver and result output) would immediately enable quadratic end-to-end runs on non-adaptive problems without electrodes.

### 4.2 Regular Grid and Spatial Lookup Chain

`Geometry::makeRegularGrid()` is blocked, which cascades to:

- `Geometry::genIndToTetsByCoords()` — unavailable
- `interpolateQTensor()` — cannot run (needs `genIndToTetsByCoords` for nearest-tet lookup)
- `autoref()` — cannot complete even if its own guard is removed

The fix for `makeRegularGrid()` itself is likely small: the regular grid only needs the 4 corner nodes of each TET10 to determine bounding boxes and containment (which is exact for straight-edged elements). The mid-edge nodes are not needed for spatial indexing.

### 4.3 Adaptive Mesh Refinement

The entire refinement module (`src/refinement/`, `tet-splitter`, `tet-classifier`, `midpoint-node-lookup`) operates on linear (4-node) elements. Red-green splitting produces 4-node child tets. After refinement there is no re-promotion to quadratic.

One practical approach for adaptive quadratic refinement: refine using the linear machinery, then call `convertLinearMeshDataToQuadratic()` on the result. This would need to be wired up after the regular-grid blocker is resolved (since the Q-tensor interpolation needed to map old solution values onto the new mesh also depends on the regular grid).

### 4.4 Q-Tensor Interpolation (`refinement/q-tensor-interpolator.cpp`)

`interpolateQTensor()` uses only corner nodes (indices 0–3) for both barycentric coordinate calculation and Q-tensor value sampling. For straight-edged TET10 this is **geometrically exact** (the same geometry as TET4), but it discards the higher-order field representation — interpolated Q values will be linear within each element even though the original field is quadratic. This is not documented in the source.

### 4.5 GiD Mesh Reader

`ReadGiDMesh3D()` reads only 4-node tets and 3-node tris. GiD meshes are always loaded as linear, then converted by `prepareGeometry()`. No path exists to load a native quadratic GiD mesh.

### 4.6 `reorderBoundaryTetNodes()` Incomplete (Linear Version)

The linear variant `reorderBoundaryTetNodes()` in `fe/fe-util.h` explicitly throws `NotYetImplementedException("Reordering only 4-node tetrahedra is supported")` for any input with more than 4 nodes. A `TODO` comment notes the generalised ordering for higher-order elements. The quadratic counterpart `reorderQuadraticBoundaryTetNodes()` exists and must be called explicitly; `assembleNeumann()` already branches correctly.

---

## 5. Test Coverage

### 5.1 Tests That Exist and Pass

| Test file | What is tested |
|---|---|
| `tests/cpp/fe/gaussian-integration-tests.cpp` | `TetShapeFunction(2)` — volume integrals (identity Jacobian, x⁴ polynomial, Q-tensor sampling, permittivity sampling); TET10 boundary integral with a real mesh; `TriShapeFunction(2)` — weight sums and area integrals |
| `tests/cpp/mesh/element-split-convert-tests.cpp` | TET10 ↔ 8×TET4 split/recombine; TRI6 ↔ 4×TRI3 split/recombine; `convertLinearMeshDataToQuadratic` for single tet, two-tet shared-face, and full Gmsh mesh |
| `tests/cpp/io/meshio-tests.cpp` | Loading a quadratic Gmsh mesh yields `ElementType::QUADRATIC_TETRAHEDRON` |

There is an explicit `TODO` comment in `gaussian-integration-tests.cpp` (near line 13): *"repeat this with quadratic element"* for the linear tet integration test case.

### 5.2 Missing Tests

The following areas have no tests or are actively blocked from being tested:

| Untested area | Notes |
|---|---|
| `PotentialSolver` with quadratic mesh | A test resource `RESOURCE_THIN_QUADRATIC_GMSH_MESH` is referenced in the potential solver test file, but the test hits the `NotYetImplementedException` guard immediately. The test exists but cannot pass. |
| LC solver (Q-tensor) end-to-end with quadratic mesh | No integration test runs the full LC solve loop with TET10. The unit tests for assembly are adequate, but there is no solve-to-convergence smoke test. |
| `reorderQuadraticTetNodeOrder()` | The Gmsh node-order fix in `inits.cpp` has no dedicated test. |
| `recombineLinearisedMeshToQuadratic()` round-trip | No test verifies the split→recombine cycle preserves connectivity. |
| Result output with quadratic mesh | No test exercises VTK or LcView output paths with a TET10 mesh. |
| `interpolateQTensor()` with quadratic mesh | Not tested. Corner-node-only interpolation should still be tested explicitly to document the approximation order. |
| Geometry preparation pipeline (`prepareGeometry`) | No standalone test; covered implicitly only via integration tests that load a Gmsh mesh. |

### 5.3 Suggested New Tests

1. **Potential solver smoke test** — once the `NotYetImplementedException` is removed, use `RESOURCE_THIN_QUADRATIC_GMSH_MESH` and verify that the solve produces the same voltage distribution (to within numerical tolerance) as the linear version of the same mesh.
2. **LC solve convergence test** — run the LC solver to convergence on a small TET10 mesh; check that the Q-tensor norms are finite and energy is decreasing.
3. **`reorderQuadraticTetNodeOrder()`** — load a Gmsh-exported TET10 mesh and verify that mid-edge node positions match expected midpoints after reordering.
4. **`recombineLinearisedMeshToQuadratic()` round-trip** — split a TET10 mesh to linear, run the recombine function, and verify the result matches the original connectivity.
5. **VTK output round-trip** — write a TET10 mesh to VTK and verify cell type 24 is emitted.
6. **`interpolateQTensor()` quadratic vs linear** — interpolate a known quadratic field; document/assert the expected (linear) accuracy.

---

## 6. Technical Debt

### 6.1 Dual Shape-Function API Coexistence

The old `GaussianQuadratureTet<NGP>` and `GaussianQuadratureTri<NGP>` template classes (linear-only, hardcoded to 4 and 3 nodes respectively) coexist with the new `TetShapeFunction` / `TriShapeFunction` classes. Several methods on the old templates are annotated `/** Deprecated */`. They are still instantiated in some code paths and would silently give wrong results if handed quadratic element data. These should be removed or completely replaced.

### 6.2 Per-Gauss-Point Jacobian Computation

`TetShapeFunction::initialiseElement()` recomputes the element inverse Jacobian at **every** Gauss point via the standard isoparametric summation over all `nodesPerElement` nodes.

Although the Jacobian is mathematically constant for straight-sided elements, the quadratic basis derivatives in the sum evaluate to slightly different floating-point values at different Gauss-point locations (relative variation ≈ 2×10⁻¹⁵). When a single cached value is used for all Gauss points, the consistent bias accumulates across ~10 000 elements and a ≈15 000-DOF linear system to produce a ~1.4×10⁻³ error in the solved potential — larger than the `3×10⁻⁴` test tolerance. Recomputing at every Gauss point avoids this systematic bias; the uncorrelated floating-point errors across Gauss points partially cancel in the Gauss integration, preserving the original accuracy.

This behaviour is documented in `plans/quadratic-jacobian-optimization.md`, which records the investigation and explains why the caching optimization was not retained.

### 6.3 `Mesh::Dimension` Redundancy

The `Mesh` constructor comment notes that the `dimension` field is redundant since it can be deduced from `elementType`. The field and its setter/getter remain, creating an inconsistency risk.

### 6.4 `Mesh::nElements` Double-Accounting

`nElements` is stored as both a field and implicitly via `nodes.size() / getnNodes()`. `getnElements()` computes from the latter. The two can become inconsistent. `setnElements()` exists solely to maintain the field. This should be unified.

### 6.5 `setAllNodes()` Deprecated Call Sites

`setAllNodes()` is marked deprecated but is still called in old code paths.

### 6.6 LcView Implicit Mesh Expansion

`writeQuadraticTetrahedra()` silently expands 1 TET10 → 8 TET4 in the output file. Result files for quadratic meshes therefore contain 8× more connectivity entries than the mesh in memory. This is undocumented in comments or user-facing documentation and could cause confusion during postprocessing.

### 6.7 `recombineLinearisedMeshToQuadratic()` Fragility

Uses heuristics (multiples-of-8/4 element counts, matching materials) to detect a pre-split quadratic mesh. On failure it silently falls back to `convertLinearMeshDataToQuadratic()`, which creates a different set of nodes. This silent fallback can corrupt a workflow that depends on stable node indices across a save/reload cycle.

### 6.8 Magic Scaling Constants

`1e-18` (tet volume, µm³→m³), `1e-12` (tri area, µm²→m²), and `1e-6` (per-coordinate scale) appear in multiple assembly locations without named constants. This makes unit analysis fragile and the constants are easy to miss during a review.

### 6.9 Stateful `setIntegrationPoints()` Lazy Initialiser

`TetShapeFunction::setIntegrationPoints()` returns silently if already called, making the object stateful in a non-obvious way. The current workaround in parallel assembly loops is to make `firstprivate` copies for each OpenMP thread. A cleaner design would be to make this purely functional (immutable after construction, or factory pattern) and avoid mutable lazy-init state.

---

## 7. Modularity and Coupling Concerns

### 7.1 `prepareGeometry()` Is the Implicit "Always Quadratic" Gateway

The decision to always convert to TET10/TRI6 is buried inside a general geometry-preparation function with no configuration switch. There is currently no way to run in strictly linear mode for debugging or performance comparison without modifying `prepareGeometry()` itself. Converting this to an explicit policy parameter would make the intent visible and testable.

### 7.2 Entry-Point Guards Disconnected from Assembly Readiness

`solvePotential()` and `ResultOutput::writeResults()` throw at their public entry points, even though the assembly and output code immediately below them is quadratic-aware. This creates the false impression that quadratic is unsupported throughout those subsystems. The guards are coarse-grained; any format- or solver-specific limitations should live inside the individual format writers or solver paths, not at the top level.

### 7.3 `makeRegularGrid` Block Cascades Through Unrelated Features

The cascade `makeRegularGrid` → `genIndToTetsByCoords` → `interpolateQTensor` → `autoref` means that one unresolved guard blocks spatially unrelated features (interpolation, refinement). These dependencies should be documented explicitly, and ideally the spatial lookup interface should be injectable so that individual features can be tested in isolation without requiring the full regular-grid infrastructure.

### 7.4 Result Output Format Limitations Should Live in Format Writers

The linear-only guard in `ResultOutput::writeResults()` prevents both `VtkResultFormatWriter` (which supports TET10 natively via VTK cell type 24) and `LcViewResultFormatWriter` (which supports TET10 by splitting) from being exercised. Moving the guard into the individual writers would allow the VTK path to work immediately, while the LcView path continues to document its implicit-expansion behaviour.

### 7.5 `interpolateQTensor()` Implicit Linear-Only Contract

`interpolateQTensor` uses only corner nodes and linear barycentric coordinates. For straight-edged TET10 this is exact for geometry but linear-order for field values. This implicit contract is not expressed in the function signature or documentation, making it a hidden coupling between the interpolator and the assumption that elements are geometrically linear (straight edges, no curved isoparametric mapping).

### 7.6 GiD Reader Always Forces Conversion

Every GiD mesh loaded goes through `convertLinearMeshDataToQuadratic()` unconditionally. The original linear mesh is discarded, the conversion cost is always paid, and there is no signalling to the caller that this happened. A cleaner design would have the GiD reader produce `RawMeshData` with `elementOrder_ == 1` and document that `prepareGeometry()` will promote it, or alternatively expose a GiD mesh upgrade utility.

---

## 8. Recommended Next Steps

In rough priority order:

1. **Remove the `solvePotential()` guard** — this immediately enables quadratic potential solves. The assembly code is ready. Add a smoke test using the existing quadratic test mesh resource.
2. **Remove the `ResultOutput::writeResults()` guard** — enables VTK output immediately. Move any format-specific guards into the individual writers.
3. **Fix `makeRegularGrid()` for TET10** — use only the 4 corner nodes for bounding-box / containment checks (correct for straight-edge TET10). This unblocks `interpolateQTensor` and the refinement pipeline.
4. **Wire up adaptive refinement for quadratic meshes** — after the grid blocker is resolved; refine using linear machinery then call `convertLinearMeshDataToQuadratic()` on the result and re-map the Q-tensor solution.
5. **Remove the old `GaussianQuadratureTet<NGP>` / `GaussianQuadratureTri<NGP>` deprecated classes** — reduce API surface and eliminate silent-wrong-answer risk.
6. **Jacobian caching** — investigated and found numerically unsafe with the current isoparametric assembly; see `plans/quadratic-jacobian-optimization.md` for details and possible future approaches.
7. **Add tests for all the missing areas** listed in §5.2.
8. **Document the LcView implicit expansion** explicitly in code and user documentation.
9. **Replace magic scaling constants** with named `constexpr` values.
10. **Redesign `setIntegrationPoints()` lazy-init** to be stateless / construction-time initialised.

