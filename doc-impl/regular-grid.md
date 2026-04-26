# Regular Grid Architecture

This document describes the regular-grid subsystem after the Part-2 decoupling refactor.

---

## Overview

The regular grid is a uniformly-spaced interpolation grid that overlays the
tetrahedral mesh. For each grid point the system precomputes which tetrahedron
contains it and the barycentric weights; later field values (potential, Q-tensor
director/order parameter) can be interpolated onto the grid in O(n) time.

---

## Class Responsibilities

### `TetMeshSearch`  (`geom/tet-mesh-search.h`, `geom/tet-mesh-search.cpp`)

Self-contained spatial search over a tetrahedral mesh.  Accepts `const Mesh&`,
`const Coordinates&` and `const AABox&` in its constructor; all operations are
read-only and thread-safe for the mesh data.

Key methods:
- `genIndToTetsByCoords()` — batch search: for every point in a `Coordinates`
  object find the index of the containing tetrahedron.
- `bruteForceSearch()` — single-point linear scan over all elements (fallback).
- `recursiveNeighbourSearch()` (private) — greedy centroid walk used as the
  primary search strategy.

Works for both `LINEAR_TETRAHEDRON` (TET4) and `QUADRATIC_TETRAHEDRON` (TET10)
element types because all geometry operations use only the four corner-node
indices (0–3) which are always the first four entries in the node list.

### `RegularGrid`  (`regulargrid.h`, `regulargrid.cpp`)

Owns the lookup table and provides interpolation methods.

**No longer depends on `Geometry`.**  The creation entry point is:

```cpp
bool createFromTetMesh(unsigned nx, unsigned ny, unsigned nz,
                       const Mesh& tets,
                       const Coordinates& coords,
                       const AABox& bounds);
```

Internally `generateLookupList()` creates a `TetMeshSearch` and calls
`genIndToTetsByCoords()`.

### `buildRegularGrid()`  (`regulargrid-factory.h`, `regulargrid-factory.cpp`)

Free function factory.  Accepts a `Geometry` reference and grid dimensions,
returns `unique_ptr<RegularGrid>` (or `nullptr` when any dimension is zero).
This is the preferred call site for all simulation code.

### `Geometry`  (`geometry.h`, `geometry.cpp`)

**No longer owns a `RegularGrid`.**  `makeRegularGrid()` and `getRegularGrid()`
have been removed.  `genIndToTetsByCoords()` and `brute_force_search()` are
thin wrappers around `TetMeshSearch`.

---

## Ownership and Lifetime

The canonical owner of the `RegularGrid` is `SimulationContainer` via its
`regGrid` member (`unique_ptr<RegularGrid>`).

`buildRegularGrid()` is called:
- After `prepareGeometry()` in `SimulationContainer::initialise()`.
- After mesh refinement inside `autoref()` in `autorefinement.cpp`.

The grid pointer is threaded as `RegularGrid*` through:

```
SimulationContainer::regGrid
  ↓  regGrid.get()
handleInitialEvents(…, regGrid)
  └─ handlePreRefinement(…, regGrid)   ← may update regGrid
handleEvents(…, regGrid)
  └─ handleMeshRefinement(…, regGrid) ← may update regGrid
       └─ autoref(…, regGrid)         ← rebuilds regGrid after refinement
ResultOutput::writeResults(geom, v, q, regGrid.get(), state)
  └─ ResultFormatWriter::setRegularGrid(regGrid)
```

After any mesh-refinement event, `autoref()` rebuilds the grid from the refined
geometry and writes it back through the `unique_ptr<RegularGrid>&` out-parameter.

---

## Result Output

Format writers that need the regular grid (`RegularVTKFormatWriter`,
`RegularVecMatFormatWriter`, `DirStackZFormatWriter`) retrieve it via
`regularGrid_` (a `RegularGrid*` set by `setRegularGrid()` in
`ResultOutput::writeResults()`).  Writers that do not use the grid
(`CsvUnstructuredFormatWriter`, `VtkUnstructuredAsciiGridFormatWriter`)
safely ignore it.

---

## TET10 Support

`TetMeshSearch` and `RegularGrid::generateLookupList()` work correctly for
quadratic (TET10) element meshes.  The geometry computations (`containsCoordinate`,
`calcLocCoords`, `elementCentroid`) operate only on the four corner-node
coordinates (indices 0–3), which are exact for straight-edged quadratic elements.
Field interpolation is therefore geometrically exact but first-order accurate in
field values — mid-edge Q-tensor values at indices 4–9 are not used.

