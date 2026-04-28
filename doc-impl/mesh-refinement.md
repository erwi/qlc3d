# Mesh Refinement

## Overview

qlc3d supports adaptive mesh refinement (AMR) of tetrahedral meshes during simulation. The refinement is geometry-only (no de-refinement/coarsening). After refinement, the Q-tensor and potential solution vectors are interpolated onto the new mesh, and the simulation continues from that state.

The algorithm is a variant of **red-green refinement** for tetrahedral meshes. Selected tetrahedra ("red" tets) are split into 8 children. Neighbouring tetrahedra that share bisected edges but are not themselves fully split are split in a compatible but simpler way ("green" tets), preserving mesh conformity.

---

## Source Files

| File | Responsibility |
|------|---------------|
| `qlc3d/includes/meshrefinement.h` | Configuration data structures: `RefinementConfig`, `MeshRefinement` |
| `qlc3d/src/meshrefinement.cpp` | Implementation of `RefinementConfig::validate()` and `toSpecs()`, `MeshRefinement` |
| `qlc3d/includes/refinement/refinement-spec.h` | `RefinementSpec` – unified config/runtime descriptor |
| `qlc3d/src/refinement/refinement-spec.cpp` | Implementation of `RefinementSpec` factory methods and helpers |
| `qlc3d/src/findrefelems.cpp` | Element selection: which tetrahedra need refinement |
| `qlc3d/src/refinement.cpp` | Core refinement pipeline: `Refine`, `modify_geometry` |
| `qlc3d/src/refinement2.cpp` | Thin new-element dispatcher and geometry updates |
| `qlc3d/includes/refinement/midpoint-node-lookup.h` | Abstract edge→midpoint lookup interface used by the splitters |
| `qlc3d/src/refinement/midpoint-node-lookup.cpp` | Map-backed midpoint lookup implementation and edge hashing helpers |
| `qlc3d/includes/refinement/tet-splitter.h` | Pure splitter API for tet and triangle child connectivity templates |
| `qlc3d/src/refinement/tet-splitter.cpp` | Pure, testable tet/triangle splitting templates |
| `qlc3d/includes/refinement/tet-classifier.h` | Public API for `ClassificationResult`, `classifyRefinement`, and the four pure helper functions |
| `qlc3d/src/refinement/tet-classifier.cpp` | Implementation of tet/triangle classification and the iterative expansion loop |
| `qlc3d/includes/refinement/periodic-edge-expander.h` | Public API for `PeriodicEdgeExpander`: catalogue construction and periodic edge expansion |
| `qlc3d/src/refinement/periodic-edge-expander.cpp` | Implementation of `PeriodicEdgeExpander`: catalogue, `findTranslation`, and `expand` |
| `qlc3d/includes/refinement/element-selector.h` | Public API for `ElementSelector`: abstract interface and factory methods |
| `qlc3d/src/refinement/element-selector.cpp` | Concrete implementations: `ChangeSelector`, `SphereSelector`, `BoxSelector` |
| `qlc3d/includes/refinement/q-tensor-interpolator.h` | Public API for `interpolateScalar` and `interpolateQTensor` |
| `qlc3d/src/refinement/q-tensor-interpolator.cpp` | Barycentric Q-tensor interpolation from an old mesh onto a new mesh |
| `qlc3d/src/autorefinement.cpp` | High-level `autoref()`: refinement loop, `interpolateQTensor`, geometry swap |
| `qlc3d/src/refinementhandler.cpp` | Event-system integration: `handleMeshRefinement`, `handlePreRefinement` |
| `qlc3d/includes/refinement.h` | Refinement declarations and refinement type constants |

### Refinement module layout

The pipeline is split into focused modules:

- `element-selector.h/.cpp` provides the `ElementSelector` abstract interface with three concrete implementations (`ChangeSelector`, `SphereSelector`, `BoxSelector`) plus static factory methods `makeChange`, `makeSphere`, `makeBox`, and `fromSpec`. Each selector is independently testable.
- `midpoint-node-lookup.h/.cpp` provides the abstract edge-to-midpoint lookup interface and its map-backed implementation.
- `tet-splitter.h/.cpp` provides pure splitter functions for tet and triangle child connectivity, plus the node-resolution helpers for green-2 and green-3 cases.
- `tet-classifier.h/.cpp` provides the complete tet/triangle classification pipeline, implemented as four pure helper functions plus the high-level `classifyRefinement` composer. `classifyRefinement` uses `PeriodicEdgeExpander` to handle periodic boundary mirroring.
- `periodic-edge-expander.h/.cpp` provides `PeriodicEdgeExpander`: catalogues boundary edges once at construction time and exposes `expand()` to add periodic mirror edges to the bisected-edge list.
- `q-tensor-interpolator.h/.cpp` provides `interpolateScalar` and `interpolateQTensor` for barycentric interpolation of Q-tensor fields from an old mesh onto a new mesh.
- `refinement2.cpp` is a thin adapter that resolves element node ordering, performs midpoint lookups, and dispatches to the pure splitters.
- `refinement.cpp` contains only `Refine` and `modify_geometry`.

Unit coverage lives under `tests/cpp/refinement/`:

- `midpoint-node-lookup-tests.cpp`
- `tet-splitter-tests.cpp`
- `tet-classifier-tests.cpp`
- `periodic-edge-expander-tests.cpp`
- `q-tensor-interpolator-tests.cpp`
- `refinement-spec-tests.cpp`
- `refinement-event-tests.cpp`
- `element-selector-tests.cpp`

The broader end-to-end regression for the full refinement pipeline is in `tests/cpp/meshrefinement-test.cpp`.

---

## Data Model

### `RefinementSpec`
Unified descriptor for a mesh refinement region. Created via the static factory methods:

- `RefinementSpec::makeExplicit(type, iteration, time, values, x, y, z)` — for an event that fires at a specific iteration or simulation time (exactly one of `iteration`/`time` must be positive).
- `RefinementSpec::makePeriodic(type, values, x, y, z)` — for an event that recurs on every period tick (no explicit schedule).

Key accessors:
- `getType()` — `Type::Change`, `Type::Sphere`, or `Type::Box`
- `isPeriodic()` — true iff constructed via `makePeriodic`
- `getIteration()` / `getTime()` — when the event fires (0 if not applicable)
- `getRefIter()` — number of successive refinement sub-passes: `values.size()` for Change/Sphere, `x.size()/2` for Box
- `getValue(i)`, `getX()`, `getY()`, `getZ()` — spatial/threshold parameters
- `clone()` — creates an independent deep copy (used when scheduling periodic repetitions)

### `RefinementConfig`
Plain data read from the settings file for a single `[REFINEMENT]` block. Contains:
- `type_` — `"change"`, `"sphere"`, or `"box"`
- `iterations_` / `times_` — when to trigger (both empty = periodic)
- `values_` — thresholds or radii per refinement iteration
- `x_`, `y_`, `z_` — spatial bounds (for `sphere`/`box`)

Validated at construction time via `validate()`. Converted to a list of `RefinementSpec` objects via `toSpecs()`:
- If periodic, produces a single `makePeriodic` spec.
- If explicit, produces one `makeExplicit` spec per iteration, plus one per time.

### `MeshRefinement`
Container stored on the `Simu` object. Holds a list of `RefinementConfig`s and the periodic repetition period `repRefIter_`.

---

## Triggering

Refinement is driven by the event system:

1. During setup, `createMeshRefinementEvents()` in `inits.cpp` iterates `RefinementConfig` objects, calls `toSpecs()`, and creates `Event` objects using `Event::makeRefinement(std::move(spec), iteration)` or `Event::makeRefinement(std::move(spec), time)`. Periodic specs are stored as `repRefinements_` in the `EventList`.
2. At each simulation step, `EventList::manageReoccurringEvents()` clones (`RefinementSpec::clone()`) each periodic spec and inserts a new timed event.
3. `eventhandler` checks which events are due and forwards them to `handleMeshRefinement()` (or `handlePreRefinement()` for initial refinement).

`handleMeshRefinement()` (in `refinementhandler.cpp`) collects `getRefinementSpec()` pointers from all events and passes them as `std::vector<const RefinementSpec*>` to `autoref()`. Events are deleted after refinement.

The `Event` class stores the `RefinementSpec` as a `std::unique_ptr<RefinementSpec>`, ensuring automatic lifetime management. The `void* eventData_` path is retained only for `EVENT_SWITCHING` events (electrode switching).

---

## High-Level Flow (`autoref`)

```
autoref(geom_orig, geom, q, v, specs, ...)
  │
  ├─ Determine maximum number of refinement passes (maxRefIter)
  │
  ├─ [Quadratic branch] If geom_orig contains TET10/TRI6 elements:
  │    ├─ splitQuadraticGeometryToLinear(geom_orig) → raw_linear (TET4/TRI3 RawMeshData)
  │    └─ Construct linear Geometry geom_orig_linear from raw_linear
  │    (geom_orig itself is left unchanged; the linear copy is used as the loop base)
  │
  ├─ Clone geom_orig [or geom_orig_linear for quadratic input] → geom_temp (working copy)
  │
  ├─ for refiter = 0 .. maxRefIter-1:
  │    ├─ get_index_to_tred(...)   ← select which tets to split
  │    │    ├─ interpolate q onto geom_temp if needed (Change type)
  │    │    └─ findTets_Change / findTets_Sphere / findTets_Box
  │    │
  │    └─ Refine(geom_temp, i_tet)  ← split selected tets (always TET4 inside loop)
  │
  ├─ [Quadratic branch] Re-promote refined linear mesh to TET10/TRI6:
  │    ├─ Extract RawMeshData raw_refined from geom_temp
  │    ├─ convertLinearMeshDataToQuadratic(raw_refined)   (adds fresh mid-edge nodes)
  │    └─ Reconstruct quadratic Geometry geom_temp from raw_refined
  │
  ├─ Post-refinement cleanup (runs on the quadratic geom_temp for TET10 input):
  │    ├─ calculateNodeNormals
  │    ├─ rebuild regular grid via `buildRegularGrid(...)`
  │    ├─ Resize & re-initialise potential SolutionVector v
  │    ├─ interpolate Q-tensor q onto new mesh via `interpolateQTensor`
  │    │    (uses 10-node quadratic basis when source mesh is TET10)
  │    └─ re-apply strong surface BC (setStrongSurfacesQ)
  │
  └─ geom ← geom_temp  (replace current geometry)
```

### Quadratic split–refine–recombine strategy

When `geom_orig` holds TET10/TRI6 elements the high-level flow uses a three-phase approach:

1. **Split** — `splitQuadraticGeometryToLinear` decomposes every TET10 into 8 TET4 children and every TRI6 into 4 TRI3 children, producing a linear `RawMeshData`. Coordinates are passed through unchanged; mid-edge nodes are discarded (they are recomputed after recombination).
2. **Refine** — the complete existing red-green refinement pipeline (`classifyRefinement`, `Refine`, `modify_geometry`) runs on the TET4 mesh without any modification.
3. **Recombine** — `convertLinearMeshDataToQuadratic` adds fresh mid-edge nodes at exact edge midpoints to produce a valid TET10/TRI6 mesh. Mid-edge node indices are not required to be stable across refinement cycles.

The linear refinement path is entirely unchanged: when `geom_orig` is TET4 the split/recombine phases are skipped.

---

## Element Selection (`findrefelems.cpp` + `element-selector.h/.cpp`)

Three strategies mark tetrahedra as "red" (requiring full 1→8 split). Each strategy is encapsulated as a concrete `ElementSelector` subclass created via `ElementSelector::fromSpec(spec)` or the individual factory methods. The legacy free functions in `findrefelems.cpp` delegate to the same logic.

| Strategy | Function | Selection criterion |
|----------|----------|---------------------|
| `Change` | `findTets_Change` | Max absolute change in any Q-tensor component across the tet's nodes exceeds threshold |
| `Sphere` | `findTets_Sphere` | Any LC node within a sphere of given centre and radius |
| `Box`    | `findTets_Box` | Any LC node within an axis-aligned bounding box |

All three take a `const RefinementSpec&` and populate `i_tet[]` with `RED_TET` (value 6) for selected elements.

---

## Core Refinement Algorithm (`Refine`, `refinement.cpp`)

The `Refine(geom, i_tet)` function performs one refinement pass on a geometry in-place:

### Step 1 – Classify (`tet-classifier.cpp`)

`classifyRefinement(geom, i_tet)` handles all classification and region-expansion in one call:

1. **`collectRedTetEdges`** — collect all 6 edges of every red tet into a deduplicated, sorted `vector<Line>`.
2. **`PeriodicEdgeExpander::expand`** — add periodic mirror edges for meshes with periodic boundaries (constructed once before the loop, called each iteration).
3. **`countEdgesPerElement`** — for each bisected edge, find all tets and triangles that share it; build per-element counts and the element→edge index `t_to_l`.
4. **`resolveGreen3RedAmbiguity`** — tets with 3 bisected edges are green-3 if those 3 edges form a closed face loop (3 unique nodes), or red if they span 4 unique nodes.
5. **`assignRefinementTypes`** — map edge counts to type codes; counts ≥ 4 are promoted to `RED_TET` (6).
6. Steps 1–5 repeat until no new red tets appear (iterative expansion to maintain mesh conformity).
7. After the loop, surface triangles are classified with one final `countEdgesPerElement` call.

Returns a `ClassificationResult` containing:

| Field | Contents |
|-------|----------|
| `i_tet` | Refinement type per tet (0 / 1 / 2 / 3 / 6) |
| `i_tri` | Refinement type per triangle |
| `lines` | Deduplicated sorted bisected edge list |
| `t_to_l` | Tet → bisected-edge index sets |
| `e_to_l` | Triangle → bisected-edge index sets |

Refinement type codes:

| Value | Meaning | Children created |
|-------|---------|-----------------|
| 0 | No refinement | — |
| 1 (`GREEN1_TET`) | Green-1: one bisected edge | 2 tets |
| 2 (`GREEN2_TET`) | Green-2: two bisected edges | 3 or 4 tets |
| 3 (`GREEN3_TET`) | Green-3: three edges forming a face loop | 4 tets |
| 6 (`RED_TET`) | Red: all 6 edges bisected | 8 tets |

### Step 2 – Create new elements (`refinement2.cpp` + `tet-splitter`)
`create_new_elements`:
1. Builds a midpoint-node lookup mapping `(nodeA, nodeB) → newMidpointNodeIndex` for all bisected lines.
2. `create_new_coordinates` computes midpoint positions.
3. Thin dispatcher in `refinement2.cpp` resolves node ordering, looks up midpoint nodes, and calls the pure splitter templates in `tet-splitter`.

### Step 3 – Modify geometry
`modify_geometry` updates the `Geometry` object:
- Appends new coordinate data
- Removes old split elements, appends new ones
- Reorders dielectric nodes to the end
- Recalculates determinants, surface normals, etc.

---

## Q-Tensor Interpolation (`q-tensor-interpolator.h/.cpp`)

The Q-tensor interpolation logic lives in the dedicated module
`qlc3d/includes/refinement/q-tensor-interpolator.h` / `qlc3d/src/refinement/q-tensor-interpolator.cpp`.

`interpolateQTensor(geomNew, geomOld, qOld) → SolutionVector`:
1. For each new LC node, finds the old tetrahedron containing it (`genIndToTetsByCoords`).
2. Computes barycentric coordinates within that old tet (`calcLocCoords`).
3. If the largest barycentric coordinate is ≈1 (i.e. the node coincides with an existing corner), copies the old value directly.
4. Otherwise, interpolates each of the 5 Q-tensor components:
   - **TET4 source**: linearly interpolates using 4 barycentric weights and corner values (`interpolateScalar`).
   - **TET10 source**: evaluates all 10 quadratic shape functions (`TetShapeFunction(2)`) at the reference point `(ξ, η, ζ) = (λ1, λ2, λ3)` derived from the 4 barycentric weights, then performs a 10-node weighted sum. This preserves quadratic field accuracy across mesh changes.

The lower-level helper `interpolateScalar(loc, S)` computes Σ loc[i]·S[i] for four barycentric weights and corner values (used by the TET4 path only).

Unit tests are in `tests/cpp/refinement/q-tensor-interpolator-tests.cpp`, including `interpolateQTensor_quadraticFieldOnTET10_isExact` which verifies exact interpolation of a quadratic field on a TET10 mesh.

---

## Periodicity Handling

When the mesh has periodic boundaries, bisecting an edge on one periodic face must be mirrored on the opposite face to maintain periodicity. `PeriodicEdgeExpander` encapsulates this logic:

- **Construction** scans all `MAT_PERIODIC` surface triangles and classifies their edges into per-face buckets (`front`, `back`, `left`, `right`, `top`, `bottom`) and per-corner buckets (12 corner-line vectors for the vertical and horizontal edge lines at box corners).
- **`expand(lines)`** iterates every bisected edge, detects which periodic face or corner it belongs to, finds the translational mirror in the opposite bucket via `findTranslation`, and appends it. Duplicates are removed before returning.
- **`findTranslation`** is a static predicate: given an edge and a candidate list, it returns the first candidate that is a translational copy (same shape, shifted in the coordinate direction(s) indicated by a 3-element `dir[]` mask).

**Known limitation**: top/bottom periodicity (`isTopBottomPeriodic`) always throws `RUNTIME_ERROR("Top-bottom periodicity not implemented in mesh refinement yet.")`.

---

## Known Issues and Technical Debt

1. **The higher-level dispatcher (`refinement2.cpp`) lacks direct unit tests** — the midpoint lookup, splitter, and classifier layers all have dedicated unit tests, but the dispatcher itself is validated only through end-to-end integration tests.

2. **Top/bottom periodicity not implemented** — `PeriodicEdgeExpander::expand` always throws `RUNTIME_ERROR("Top-bottom periodicity not implemented in mesh refinement yet.")` when the mesh has top/bottom periodic boundaries.
