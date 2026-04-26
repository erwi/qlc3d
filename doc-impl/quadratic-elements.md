# Quadratic Element Support — Current Snapshot

> Last updated: April 2026

This document describes the current state of quadratic finite-element support in qlc3d for 10-node tetrahedra (`TET10`) and 6-node triangles (`TRI6`). It is a snapshot of what works now, what is still limited, and what work would unlock additional functionality.

---

## Table of Contents

1. [Current Scope](#1-current-scope)
2. [Implemented Quadratic Support](#2-implemented-quadratic-support)
3. [Current Limitations](#3-current-limitations)
4. [Test Coverage Snapshot](#4-test-coverage-snapshot)
5. [Current Caveats and Technical Debt](#5-current-caveats-and-technical-debt)
6. [Recommended Milestones](#6-recommended-milestones)

---

## 1. Current Scope

Quadratic support is active across the main geometry, assembly, search, and output paths.

The key project-wide behavior is in `prepareGeometry()` in `inits.cpp`:

- linear meshes are promoted in memory to quadratic meshes through `convertLinearMeshDataToQuadratic()`;
- native quadratic Gmsh meshes are accepted directly;
- downstream simulation code normally operates on quadratic tetrahedra and triangles.

For straight-edged elements, the codebase treats quadratic geometry as the same physical volume and surface geometry defined by the corner nodes, with mid-edge nodes carrying the higher-order field representation.

---

## 2. Implemented Quadratic Support

### 2.1 Mesh model and element metadata

- `ElementType::QUADRATIC_TRIANGLE` and `ElementType::QUADRATIC_TETRAHEDRON` are defined.
- `getNodesPerElement()` and `getElementOrder()` return the correct values for quadratic elements.
- `Mesh::setElementData()` validates node-count consistency for quadratic connectivity.
- Mesh traversal code that uses `getnNodes()` includes mid-edge nodes.
- `calculateDeterminants3D()` and `calculateSurfaceNormals()` use corner nodes, which is exact for straight-edged quadratic elements.

### 2.2 Shape functions and quadrature

- `TetShapeFunction(2)` implements the full 10-node basis and derivatives.
- `TriShapeFunction(2)` implements the full 6-node basis and derivatives.
- Sampling helpers such as `N(i)`, `Nx(i)`, `Ny(i)`, `Nz(i)`, `sample()`, `sampleQ()`, and derivative samplers work with quadratic element sizes.
- Quadratic LC assembly uses `Keast8` for volume terms and `Tri4thOrder` for surface terms.
- Quadratic potential assembly uses `Keast4` for volume terms and `Tri4thOrder` for Neumann surface terms.

### 2.3 Mesh I/O and geometry preparation

- The Gmsh reader accepts element type `9` (`TRI6`) and `11` (`TET10`).
- `RawMeshData` carries element order through the read/prepare pipeline.
- `prepareGeometry()`:
  - recombines previously linearised quadratic meshes when the split pattern is detected,
  - promotes first-order meshes to second order when needed,
  - reorders Gmsh TET10 mid-edge nodes into qlc3d's internal convention,
  - validates and snaps quadratic mid-edge nodes when they are within tolerance of the expected midpoint.

### 2.4 Conversion and split/recombine utilities

The following utilities exist and are in active use:

- `convertLinearMeshDataToQuadratic()`
- `reorderQuadraticTetNodeOrder()`
- `recombineLinearisedMeshToQuadratic()`
- `splitQuadraticTetrahedronToLinear()`
- `splitQuadraticTriangleToLinear()`
- `recombineLinearTetsToQuadratic()`
- `recombineLinearTrianglesToQuadratic()`

### 2.5 Solvers and boundary handling

- LC volume assembly and weak-anchoring assembly accept quadratic elements.
- `PotentialSolver::solvePotential()` accepts quadratic meshes and reaches quadratic-aware assembly code.
- Neumann assembly in the potential solver branches correctly between `reorderQuadraticBoundaryTetNodes()` and `reorderBoundaryTetNodes()` based on element type.

### 2.6 Spatial lookup and result output

- `RegularGrid` can be built for quadratic tetrahedral meshes.
- Spatial lookup for straight-edged TET10 remains geometrically exact because containment and barycentric coordinates use the corner-node geometry.
- The VTK writer emits VTK cell type `24` for quadratic tetrahedra.
- LcView mesh output supports quadratic input meshes by splitting each quadratic element into linear sub-elements during export.

---

## 3. Current Limitations

### 3.1 Adaptive refinement is still linear-only

`autoref()` in `autorefinement.cpp` rejects any tetrahedral mesh whose element type is not `ElementType::LINEAR_TETRAHEDRON`.

Current consequence:

- adaptive refinement cannot be run on the quadratic meshes produced by `prepareGeometry()`;
- a full quadratic simulation can solve on a prepared mesh, but it cannot continue through the refinement path.

### 3.2 Q-tensor transfer after mesh changes is only first-order on TET10

`interpolateQTensor()` uses only the four corner nodes of the containing tetrahedron:

- geometry lookup is still exact for straight-edged TET10;
- field transfer is linear within each element even when the source solution is quadratic.

Current consequence:

- any workflow that depends on carrying a quadratic solution between meshes loses the higher-order field representation.

### 3.3 Native quadratic GiD input is not supported

`ReadGiDMesh3D()` reads only 4-node tetrahedra and 3-node triangles.

Current consequence:

- GiD input always enters the system as linear data and is then promoted by `prepareGeometry()`;
- there is no direct path for loading native quadratic GiD connectivity.

### 3.4 Generic boundary-node reordering is not higher-order-safe

`reorderBoundaryTetNodes()` throws for anything other than a 4-node tetrahedron. Quadratic Neumann assembly works today because it explicitly calls `reorderQuadraticBoundaryTetNodes()`, but the generic helper is still linear-only.

Current consequence:

- higher-order callers must know to bypass the generic helper and call the quadratic-specific function.

---

## 4. Test Coverage Snapshot

### 4.1 Existing quadratic coverage

| Test file | Current coverage |
|---|---|
| `tests/cpp/fe/gaussian-integration-tests.cpp` | TET10 and TRI6 shape functions, volume/surface integration, and quadratic sampling helpers |
| `tests/cpp/mesh/element-split-convert-tests.cpp` | TET10/TRI6 split-recombine utilities and linear→quadratic conversion |
| `tests/cpp/io/meshio-tests.cpp` | Quadratic Gmsh mesh loading |
| `tests/cpp/init-geometry-tests.cpp` | `prepareGeometry()`, quadratic node-order correction, midpoint snapping/validation, and standalone regular-grid creation |
| `tests/cpp/potential/potential-solver-test.cpp` | Quadratic potential solves, including Neumann-boundary cases |
| `tests/cpp/lc/lc-solver-tests.cpp` | Steady-state LC solve on quadratic meshes |
| `tests/cpp/geom/geom-tests.cpp` | TET10 containment, barycentric coordinates, centroid behavior, and regular-grid interaction |
| `tests/cpp/io/result-output-test.cpp` | Quadratic LCView result I/O and quadratic VTK cell-type output |
| `tests/cpp/refinement/q-tensor-interpolator-tests.cpp` | Q-tensor interpolation size/copy behavior and exact preservation of linear fields through refinement |

### 4.2 Coverage gaps that still matter

- No test currently exercises adaptive refinement on an in-memory quadratic mesh, because `autoref()` rejects it.
- `interpolateQTensor()` is tested for linear fields, but there is no regression test that documents the expected loss of quadratic field order on TET10.
- LcView export is tested through result read/write behavior, but there is no explicit regression test that documents the 1 `TET10` → 8 `TET4` mesh expansion.
- `recombineLinearisedMeshToQuadratic()` does not have a regression test focused on failure behavior and node-index stability across save/reload-style workflows.

---

## 5. Current Caveats and Technical Debt

### 5.1 `prepareGeometry()` is implicitly "always quadratic"

`prepareGeometry()` performs the mesh-order promotion policy internally. There is no explicit switch for keeping the simulation in first-order mode.

Current consequence:

- linear-vs-quadratic comparisons, debugging, and performance studies require code changes or alternate setup paths instead of a documented runtime choice.

### 5.2 `recombineLinearisedMeshToQuadratic()` is heuristic and non-authoritative

`recombineLinearisedMeshToQuadratic()` infers whether a first-order mesh is really a split quadratic mesh from element-count divisibility and connectivity/material consistency checks.

Current consequence:

- a workflow that depends on stable quadratic node identity across write/read cycles is sensitive to whether the heuristic succeeds;
- if recombination does not trigger, `prepareGeometry()` falls back to fresh linear→quadratic promotion, which produces a different quadratic node set.

### 5.3 LcView export changes mesh topology representation

`writeQuadraticTetrahedra()` and `writeQuadraticTriangles()` export quadratic meshes by splitting them into linear sub-elements because the target format does not store quadratic cells.

Current consequence:

- exported mesh connectivity does not match the in-memory quadratic element count;
- post-processing workflows must treat LcView mesh files as linearised representations of the quadratic mesh.

### 5.4 Per-Gauss-point Jacobian recomputation is required for current numerical accuracy

`TetShapeFunction::initialiseElement()` recomputes the inverse Jacobian at every Gauss point. For the current straight-sided quadratic implementation, caching a single Jacobian for the whole element introduces a measurable bias in solved results.

Current consequence:

- Jacobian caching is not currently a safe optimization path without further numerical work.

### 5.5 `setIntegrationPoints()` remains stateful

`TetShapeFunction::setIntegrationPoints()` mutates the shape-function object and silently returns if the points are already initialised.

Current consequence:

- parallel assembly relies on thread-local copies rather than a clearly immutable quadrature object design.

---

## 6. Recommended Milestones

### Milestone 1 — Quadratic adaptive refinement

**Work required**

- remove the linear-only guard in `autoref()`;
- define the refinement path for quadratic meshes, including how refined meshes are re-promoted to TET10/TRI6;
- preserve boundary conditions, regular-grid rebuilds, and solution-vector remapping through that path.

**Unlocked functionality**

- adaptive LC and electrostatic simulations on quadratic meshes;
- end-to-end quadratic workflows that include mesh refinement instead of stopping at the initial solve.

### Milestone 2 — Higher-order Q-tensor transfer between meshes

**Work required**

- extend `interpolateQTensor()` so that field transfer uses the quadratic representation instead of only the four corner nodes;
- add regression tests that distinguish exact linear preservation from expected quadratic-field behavior.

**Unlocked functionality**

- refinement and remeshing workflows that preserve quadratic solution quality;
- more accurate restart and mesh-adaptation cycles for quadratic LC simulations.

### Milestone 3 — Explicit mesh-order input policy

**Work required**

- add a documented policy for mesh-order handling at input time;
- either support native quadratic GiD meshes or make first-order GiD ingestion plus promotion an explicit, inspectable step;
- expose whether a mesh was loaded as quadratic, recombined to quadratic, or promoted from linear input.

**Unlocked functionality**

- direct quadratic GiD workflows, if native support is added;
- reproducible debugging and benchmarking across linear and quadratic input paths.

### Milestone 4 — Stable and explicit quadratic export/save-reload semantics

**Work required**

- document and test the linearisation semantics of LcView export;
- harden the split/recombine workflow around `recombineLinearisedMeshToQuadratic()` so failures are explicit and node-identity-sensitive workflows are protected.

**Unlocked functionality**

- predictable quadratic post-processing;
- safer save/reload pipelines for workflows that depend on stable quadratic connectivity and node numbering.

### Milestone 5 — Higher-order-safe helper APIs

**Work required**

- remove or generalise linear-only helper assumptions such as `reorderBoundaryTetNodes()`;
- reduce the number of call sites that must know special-case quadratic helper functions by name.

**Unlocked functionality**

- easier extension of boundary-integral code paths to additional higher-order element workflows;
- less element-order-specific coupling in future solver and post-processing code.
