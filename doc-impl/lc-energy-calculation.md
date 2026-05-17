# LC Free Energy Calculation

This document describes how the LC free energy is configured, computed, and written to output in qlc3d.

## Overview

The free energy calculation is an optional diagnostic that integrates the LC free energy density over the LC domain at every simulation iteration and writes the per-component results to a CSV file. It is entirely separate from the energy minimisation done by the solver; it only reads the current solution state and produces output.

The total free energy is:

$$
F = \int_\Omega (f_D + f_B - f_E)\, d\Omega + \int_\Gamma f_S\, d\Gamma
$$

Decomposed into four components:

| Component | Symbol | Description |
|-----------|--------|-------------|
| Elastic | Felastic | Total elastic distortion energy (L1..L6 combined) |
| Thermotropic | Fth | Landau-de Gennes bulk ordering energy |
| Electric | Fe | Dielectric + flexoelectric coupling to the field |
| Surface | Fs | Surface anchoring energy |

---

## Configuration

Energy output is controlled by a single boolean setting in the simulation settings file (`.qfg`):

```
outputEnergy = 1    # 0 = disabled (default), 1 = enabled
```

The setting key is defined in `qlc3d/includes/settings_file_keys.h` as `SFK_OUTPUT_ENERGY = "outputEnergy"` and is read in `qlc3d/src/settings-reader.cpp`. It is stored in the `Simu` class (`outputEnergy_` field, accessed via `Simu::getOutputEnergy()`). The default value is `0` (disabled).

---

## Architecture

The energy calculation is implemented as a set of independent, testable modules:

| File | Responsibility |
|------|---------------|
| `energy/energy-result.h` | `EnergyResult` struct with per-component fields and `total()` |
| `energy/lc-energy-density.h/cpp` | Pure Gauss-point energy density functions (no integration), including `surfaceEnergyDensity` |
| `energy/lc-energy-integrator.h/cpp` | Integrates all energy contributions; provides `integrateVolumeEnergy`, `integrateSurfaceEnergy`, and the combined `integrateEnergy` |
| `energy/lc-energy-calculator.h/cpp` | Top-level calculator; delegates to `integrateEnergy` and returns `EnergyResult` |
| `io/energy-csv-writer.h/cpp` | Writes `EnergyResult` rows to a CSV file |

`LcEnergyCalculator` is a stateless object injected into `SimulationContainer` via its constructor. CSV writing is handled separately by an `EnergyCsvWriter` that `SimulationContainer` creates internally during `initialise()` when `outputEnergy == 1`.

The volume and surface integrators (`integrateVolumeEnergy`, `integrateSurfaceEnergy`) are defined in the same `lc-energy-integrator.h` but remain individually callable for independent unit testing. All energy density functions live in `lc-energy-density.h/cpp`.

> **Note on absolute energy values:** No ground-state offsets are subtracted from any energy component. All reported energies (thermotropic, surface) are raw values and will generally be negative at the equilibrium state. The intended use is to monitor that the total energy decreases over simulation iterations, not to interpret the absolute numbers.

---

## Output File Format

The output file is `<saveDir>/energy.csv`. It is created at the start of the simulation. A row is appended at the start of each iteration (capturing the state from the previous iteration) and one final row is appended in `postSimulationTasks()` after all iterations have completed, capturing the final converged state.

CSV columns:

```
time_s,iteration,elastic_J,thermotropic_J,electric_J,surface_J,total_J
```

- `time_s` — simulation time in seconds.
- `iteration` — simulation iteration number.
- `elastic_J`, `thermotropic_J`, `electric_J`, `surface_J` — per-component energies in Joules.
- `total_J` — sum of all four components.

---

## Numerical Integration

### Quadrature

Integration uses the `TetShapeFunction` and `TriShapeFunction` classes from `<fe/gaussian-quadrature.h>`, which are the same classes used by the LC and potential solvers.

| Element | Shape function class | Integration points |
|---------|--------------------|--------------------|
| Linear tet (TET4) | `TetShapeFunction(1)` | `Keast8` (45 pts, 8th order) |
| Quadratic tet (TET10) | `TetShapeFunction(2)` | `Keast8` (45 pts, 8th order) |
| Linear tri (TRI3) | `TriShapeFunction(1)` | `Tri4thOrder` (7 pts, 4th order) |
| Quadratic tri (TRI6) | `TriShapeFunction(2)` | `Tri4thOrder` (7 pts, 4th order) |

The element order is read at runtime from the mesh via `tets.getElementType()` → `getElementOrder(elementType)`, so the same integrator code handles both linear and quadratic mesh types.

### Volume integral

`integrateVolumeEnergy(lc, geom, v, q)` (in `lc-energy-integrator.cpp`):

1. Iterates over all tetrahedra with `materialNumber <= MAT_DOMAIN7` (LC domain elements only).
2. For each element, calls `shapes.initialiseElement(nodes, det)` and loops over Gauss points via `shapes.hasNextPoint()` / `shapes.nextPoint()`.
3. At each Gauss point, fills a `GaussPointData` struct from shape function interpolation of `q` and `v` nodal values (including spatial derivatives via `shapes.Nx/Ny/Nz`).
4. Calls `elasticEnergyDensity`, `thermotropicEnergyDensity`, and `electricEnergyDensity` and accumulates `shapes.getWeight() * det * density`.
5. Returns a partially populated `EnergyResult` (`surface` left at 0.0).

### Surface integral

`integrateSurfaceEnergy(alignment, geom, q, S0)` (in `lc-energy-integrator.cpp`):

1. Retrieves all weak-anchoring surfaces from `alignment.getWeakSurfacesByFixLcNumber()`.
2. Iterates over all triangle elements; skips those whose `fixLcNumber` has no weak-anchoring entry.
3. For homeotropic surfaces (`usesSurfaceNormal() == true`), uses the per-node surface normal; otherwise, uses the fixed easy-axis vectors `v1`, `v2` from the `Surface` object.
4. Evaluates the raw Rapini-Papoular surface energy density (via `surfaceEnergyDensity`) at each Gauss point and accumulates the integral. No ground-state offset is applied.

---

## Energy Density Expressions

All density functions are implemented in `lc-energy-density.cpp` as pure functions taking `GaussPointData` and `EnergyMaterialParams`.

### Elastic energy

`elasticEnergyDensity` evaluates the total elastic distortion using the G1/G2/G4/G6 invariants:

| Invariant | Expression | Corresponds to |
|-----------|-----------|----------------|
| G1 | Q_ij,k Q_ij,k (in T-basis) | Isotropic elastic gradient |
| G2 | Q_ij,j Q_ik,k (in T-basis) | Divergence-squared |
| G4 | ε_ijk Q_il Q_jl,k | Chiral twist (proportional to n·curl n) |
| G6 | Q_lk Q_ij,l Q_ij,k | Non-linear elastic |

The L1..L6 Landau-de Gennes coefficients are mapped from the Frank constants K11, K22, K33 via the `LC` class.

### Thermotropic energy

```
f_B = A/2 * R  +  B/3 * tr(Q³)  +  C/4 * R²
```

No ground-state offset is applied. The value is negative at the equilibrium state.

### Electric energy

The dielectric term uses the anisotropic permittivity tensor:

```
ε_ij = (ε_⊥ + Δε/3) δ_ij  +  (2Δε / (3 S₀)) Q_ij
f_E = ε₀/2 * ε_ij * E_i * E_j
```

Flexoelectric polarisation:

```
P_i = ξ_a Q_ij,j  +  ξ_b Q_ij Q_jk,k
```

where `ξ_a = 2/(9*S0) * (e1 + 2*e3)` and `ξ_b = 4/(9*S0²) * (e1 - e3)`.

---

## Material Parameters

| Parameter | Source | Role |
|-----------|--------|------|
| S0 | `lc.S0()` | Equilibrium scalar order parameter |
| K11, K22, K33 | `lc.K11()`, etc. | Frank elastic constants |
| A, B, C | `lc.A()`, `lc.B()`, `lc.C()` | Landau thermotropic coefficients |
| ε_∥, ε_⊥ | `lc.eps_par()`, `lc.eps_per()` | Dielectric permittivities |
| e1, e3 | `lc.e1()`, `lc.e3()` | Flexoelectric coefficients |
| p0 | `lc.p0()` | Cholesteric pitch (0 = non-chiral), q0 = 2π/p0 |

---

## Element Selection

Only volume elements whose material number satisfies `materialNumber <= MAT_DOMAIN7` (LC domain) contribute to the volume integral. Surface elements are filtered by `fixLcNumber` against the set of weak-anchoring surfaces registered in the `Alignment` object.

---

## Call Sequence Summary

```
main-app-qlc3d.cpp
  └─ creates LcEnergyCalculator (stateless, injected into SimulationContainer)

SimulationContainer::initialise()
  └─ if outputEnergy == 1:
       energyCsvWriter_.emplace(savePath/"energy.csv")
            └─ opens energy.csv and writes CSV header

SimulationContainer::runIteration()
  └─ if energyCsvWriter_.has_value():
       result = energyCalculator_.calculate(lc, geom, v, q, alignment)
         └─ integrateEnergy(...)
              ├─ integrateVolumeEnergy(...)   → elastic, thermotropic, electric
              └─ integrateSurfaceEnergy(...)  → surface
       energyCsvWriter_->write(time, iteration, result)
         └─ appends one CSV row (state from the previous iteration)

SimulationContainer::postSimulationTasks()
  └─ if energyCsvWriter_.has_value():
       result = energyCalculator_.calculate(lc, geom, v, q, alignment)
       energyCsvWriter_->write(time, iteration, result)
         └─ appends final CSV row for the converged end state
```

---

## Symbolic Energy Expressions

The file `equations/lc-energy-derivation.py` contains symbolic SymPy derivations of all four energy densities (f_D, f_B, f_E, f_S) directly in the T-tensor basis (`t1..t5` and their spatial derivatives). The C++ variable mapping is:

| C++ variable | SymPy symbol | Note |
|-------------|-------------|------|
| `q1`..`q5` | `t1`..`t5` | T-tensor DOF values at Gauss point |
| `q1x`..`q5z` | `t1x`..`t5z` | Physical spatial gradients |
| `Vx`, `Vy`, `Vz` | `E_x`, `E_y`, `E_z` | Electric field = −∇φ |
| `S0` | `S_0` | Equilibrium order parameter |
