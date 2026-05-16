# LC Free Energy Calculation

This document describes how the LC free energy is configured, computed, and written to output in qlc3d.

## Overview

The free energy calculation is an optional diagnostic that integrates the LC free energy density over the LC domain at every simulation iteration and writes the per-component results to a MATLAB/Octave-compatible output file. It is entirely separate from the energy minimisation done by the solver; it only reads the current solution state and produces output.

The total free energy is:

$$
F = \int_\Omega (f_D + f_B - f_E)\, d\Omega + \int_\Gamma f_S\, d\Gamma
$$

The energy output currently covers the bulk volume integral only (no surface anchoring term), decomposed into:

| Component | Symbol | Energy contribution |
|-----------|--------|---------------------|
| Splay | F11 | K11 splay distortion |
| Twist | F22 | K22 twist distortion (including chirality) |
| Bend | F33 | K33 bend distortion |
| Thermotropic | Fth | Landau-de Gennes bulk ordering |
| Electric | Fe | Dielectric + flexoelectric coupling to the field |

---

## Configuration

Energy output is controlled by a single boolean setting in the simulation settings file (`.qfg`):

```
outputEnergy = 1    # 0 = disabled (default), 1 = enabled
```

The setting key is defined in `qlc3d/includes/settings_file_keys.h` as `SFK_OUTPUT_ENERGY = "outputEnergy"` and is read in `qlc3d/src/settings-reader.cpp`. It is stored in the `Simu` class (`outputEnergy_` field, accessed via `Simu::getOutputEnergy()`). The default value is `0` (disabled).

---

## File Lifecycle

### Opening the output file

At the end of `SimulationContainer::initialise()` (`simulation-container.cpp`), `createOutputEnergyFile()` is called (defined in `inits.cpp`):

```cpp
FILE* createOutputEnergyFile(Simu& simu) {
    // opens <saveDir>/energy.m for writing if outputEnergy == 1
}
```

The file is opened at `<saveDir>/energy.m` (within the simulation's save directory). If the file cannot be opened, a `RUNTIME_ERROR` is raised. The `FILE*` handle (`Energy_fid`) is stored as a member of `SimulationContainer`.

### Writing per-iteration data

At the start of each iteration in `SimulationContainer::runIteration()`:

```cpp
if (simu->getOutputEnergy()) {
    CalculateFreeEnergy(Energy_fid, currentIteration, currentTime, *lc, &geom1, &v, &q);
}
```

Energy is calculated **before** the Q-tensor and potential are updated for that iteration.

### Closing the file

`closeEnergyFile(fid, simu)` (in `energy.cpp`) is called from the simulation teardown path. It appends the MATLAB array terminator `];` and closes the file handle.

---

## Output File Format

The output file `energy.m` uses MATLAB/Octave matrix syntax:

```matlab
% columns are:
% time[s],splay,twist,bend,thermotropic,dielectric
F = [ ...
1.000000e-09    <F11>    <F22>    <F33>    <Fth>    <Fe>;
...
];
```

- First column: simulation time in seconds.
- Remaining columns: integrated energy components in Joules.
- The header comment is written only on the first iteration (`currentIteration == 1`).

---

## Numerical Integration

The integral over the LC domain uses **11-point 3D Gauss-Legendre quadrature** on linear tetrahedral elements (P1 elements). The quadrature is set up inside the `Energy` namespace in `energy.cpp`.

### Gauss points and weights

Constants and arrays in the `Energy` namespace:

| Symbol | Value | Description |
|--------|-------|-------------|
| `ngp` | 11 | Number of Gauss points |
| `gp[11][4]` | — | Barycentric coordinates of each Gauss point |
| `w[11]` | w11, w12, w13 | Corresponding quadrature weights |
| `w11` | −74/5625 | Weight for centre point |
| `w12` | 343/45000 | Weights for 4 vertex-adjacent points |
| `w13` | 56/2250 | Weights for 6 edge-midpoint points |

### Shape functions

`init_shape()` pre-computes the P1 shape functions and their parametric derivatives at each Gauss point:

| Array | Content |
|-------|---------|
| `sh1[ngp][4]` | Shape function values N_i(ξ_gp) |
| `sh1r[ngp][4]` | ∂N_i/∂r (r-derivative, constant ±1 for P1) |
| `sh1s[ngp][4]` | ∂N_i/∂s |
| `sh1t[ngp][4]` | ∂N_i/∂t |

### Jacobian and physical derivatives

For each element and Gauss point:

1. The 3×3 Jacobian is assembled from the node coordinates (converted from micrometres to metres via `MICROMETER_TO_METER`).
2. The determinant `Jdet` is read from the precomputed value on the `Mesh` object.
3. The inverse Jacobian `Jinv[3][3]` is computed analytically.
4. Physical shape function gradients `dSh[4][3]` are obtained by the chain rule: `dSh[i][k] = Σ_r sh1r[igp][i] * Jinv[r][k]`.

### Interpolation at Gauss points

At each Gauss point, the five T-tensor components and their spatial gradients are interpolated from the nodal DOFs stored in `SolutionVector *q`:

```
q_n     = Σ_i N_i * q->getValue(node_i, n)      (component value)
q_n_x   = Σ_i dSh[i][0] * q->getValue(node_i, n) (x-derivative)
q_n_y   = Σ_i dSh[i][1] * q->getValue(node_i, n)
q_n_z   = Σ_i dSh[i][2] * q->getValue(node_i, n)
```

Electric field components are similarly interpolated from `SolutionVector *v` (the potential), giving `(Vx, Vy, Vz) = -∇φ`.

---

## Energy Density Expressions

All expressions are in terms of the five T-tensor components `q1..q5` and their spatial derivatives. The T-basis is defined as in `equations.md`:

```
T1 = (3 ê_z⊗ê_z − I) / √6   (axial)
T2 = (ê_x⊗ê_x − ê_y⊗ê_y) / √2   (biaxial)
T3 = (ê_x⊗ê_y + ê_y⊗ê_x) / √2   (xy shear)
T4 = (ê_y⊗ê_z + ê_z⊗ê_y) / √2   (yz shear)
T5 = (ê_x⊗ê_z + ê_z⊗ê_x) / √2   (xz shear)
```

### Pre-computed scalar invariants

The code evaluates several intermediate scalar invariants before summing the energy contributions:

| Variable | Expression | Physical meaning |
|----------|-----------|------------------|
| `R` | q1²+q2²+q3²+q4²+q5² | tr(Q²) = Q_ij Q_ij |
| `G1` | L1-type gradient invariant | Isotropic elastic gradient squared: Q_ij,k Q_ij,k |
| `G2` | L2-type divergence invariant | Square of divergence of Q: Q_ij,j Q_ik,k |
| `G3` | L3-weighted K24 saddle-splay term | Q_ik,j Q_ij,k |
| `G4` | L4-type chiral twist invariant | ε_ijk Q_il Q_jl,k (proportional to n·curl n) |
| `G6` | L6-type non-linear elastic invariant | Q_lk Q_ij,l Q_ij,k |

### Elastic energy components

The Frank elastic constants K11, K22, K33 are mapped to Landau-de Gennes coefficients L1, L2, L3, L4, L6 (see [Material Parameters](#material-parameters) below). The energy calculation uses the inverse mapping to express K11, K22, K33 contributions directly:

```
F_splay = 4*G2/(9*S0²) - 2*G1/(27*S0²) - 4*G6/(27*S0³)
F_bend  = 2*G1/(27*S0²) + 4*G6/(27*S0³)
F_twist = ((G4 - (3/2)*R*q0) / (9*S0²) * 4)²   [squared, chiral offset included]

F11 += 0.5 * K11 * mul * F_splay
F22 += 0.5 * K22 * mul * F_twist
F33 += 0.5 * K33 * mul * F_bend
```

where `mul = w[igp] * Jdet` is the quadrature weight times Jacobian determinant, and `q0 = 2π/p0` is the chirality wavenumber (zero if `p0 == 0`).

### Thermotropic (bulk) energy

```
f_B = A/2 * R  +  B/3 * tr(Q³)  +  C/4 * R²
```

The `B/3 * tr(Q³)` cubic term is expanded explicitly in `q1..q5`. The reference ground-state value `f0` is subtracted so that the reported energy is zero at the equilibrium scalar order:

```
f0 = (3A/4)*S0² + (B/4)*S0³ + (9C/16)*S0⁴
Fth += mul * (f_B_elem - f0)
```

### Electric field energy

The dielectric contribution uses the anisotropic permittivity tensor:

```
ε_ij = (ε_⊥ + Δε/3) δ_ij  +  (2Δε / (3 S₀)) Q_ij
```

The electric energy density is:

```
f_E = ε₀/2 * ε_ij * E_i * E_j
```

which is implemented as:

```
Fel_elem = e0 * (-Vx² - Vy² - Vz²) * epsav * 0.5
         + e0 * deleps * (... terms in V×V×q ...)
```

where `epsav = ε_⊥/S0` and `deleps = (ε_∥ - ε_⊥)/S0`.

### Flexoelectric energy

The flexoelectric polarisation is:

```
P_i = ξ_a Q_ij,j  +  ξ_b Q_ij Q_jk,k
```

The two coefficients are computed from the flexoelectric parameters `e1` and `e3`:

```
ξ_a = efe  = 2/(9*S0) * (e1 + 2*e3)
ξ_b = efe2 = 4/(9*S0²) * (e1 - e3)
```

Each term is added to `Fflx` only if its coefficient is non-zero. The flexoelectric energy contribution is included in the `Fe` output column together with the dielectric energy.

---

## Material Parameters

### Frank elastic constants → L-coefficients

The `LC` class (in `lc.h` / `lc.cpp`) computes the Landau-de Gennes L-coefficients from the Frank constants and the equilibrium order parameter S0. S0 is found from the thermotropic parameters:

```
S0 = (-B + √(B² - 24AC)) / (6C)
```

The Frank-to-Landau conversions are:

| L-coefficient | Formula |
|---------------|---------|
| L1 | 2(K33 - K11 + 3 K22) / (27 S0²) |
| L2 | 4(K11 - K22) / (9 S0²) |
| L3 | 4 K24 / (9 S0²) |
| L4 | 8 q0 K22 / (9 S0²) — chirality wavenumber q0 = 2π/p0 |
| L6 | 4(K33 - K11) / (27 S0³) |

### Parameters used in the energy calculation

| Parameter | Source in code | Role |
|-----------|---------------|------|
| S0 | `lc.S0()` | Equilibrium scalar order parameter |
| K11, K22, K33 | `lc.K11()`, `lc.K22()`, `lc.K33()` | Frank elastic constants |
| A, B, C | `lc.A()`, `lc.B()`, `lc.C()` | Landau thermotropic coefficients |
| ε_∥, ε_⊥ | `lc.eps_par()`, `lc.eps_per()` | Dielectric permittivities (parallel/perpendicular) |
| e1, e3 | `lc.e1()`, `lc.e3()` | Flexoelectric coefficients |
| p0 | `lc.p0()` | Cholesteric pitch (0 = non-chiral) |

---

## Element Selection

Only elements whose material number satisfies `materialNumber <= MAT_DOMAIN7` (i.e., LC domain elements) contribute to the integral. Fixed-surface, electrode, or non-LC elements are excluded.

---

## Call Sequence Summary

```
SimulationContainer::initialise()
  └─ createOutputEnergyFile(*simu)        // opens <saveDir>/energy.m

SimulationContainer::runIteration()       // called every iteration
  └─ if (simu->getOutputEnergy())
       └─ CalculateFreeEnergy(fid, iter, time, lc, geom, v, q)
            ├─ init_shape()               // Gauss point shape functions
            ├─ for each LC tetrahedron:
            │    for each Gauss point:
            │      interpolate q1..q5, gradients, E-field
            │      evaluate G1, G2, G3, G4, G6, R
            │      accumulate F11, F22, F33, Fth, Fflx, Fe
            └─ fprintf(fid, "time F11 F22 F33 Fth Fe;\n")

SimulationContainer::postSimulationTasks() (or teardown)
  └─ closeEnergyFile(fid, simu)           // appends "];" and fclose
```

---

## Symbolic Energy Expressions

The file `equations/lc-energy-derivation.py` contains symbolic SymPy derivations of all four energy densities (f_D, f_B, f_E, f_S) directly in the T-tensor basis (`t1..t5` and their spatial derivatives):

| Function | Returns |
|----------|---------|
| `thermotropic_energy_density()` | f_B expanded in t1..t5 |
| `elastic_energy_density()` | f_D expanded in t1..t5 and their first derivatives, parameterised by L1,L2,L3,L4,L6 |
| `electric_energy_density()` | f_E with ε_ij and P_flexo as separate outputs |
| `surface_anchoring_energy_density()` | f_S with easy-axis vectors v1, v2 |

The script also includes **sanity checks** that verify the SymPy-derived Q matrix and ε tensor exactly match the runtime C++ implementations in `lc-representation.cpp` (via `TTensor::toQTensor()` and `DielectricPermittivity::fromTTensor()`). These symbolic expressions can serve as a foundation for **C++ code generation** to replace or validate the hand-written formulas in `energy.cpp`.

The key **C++ variable mapping** is:

| C++ variable | SymPy symbol | Note |
|-------------|-------------|------|
| `q1`..`q5` | `t1`..`t5` | T-tensor DOF values at Gauss point |
| `q1x`..`q5z` | `t1x`..`t5z` | Physical spatial gradients |
| `Vx`, `Vy`, `Vz` | `E_x`, `E_y`, `E_z` | Electric field = −∇φ |
| `S0` | `S_0` | Equilibrium order parameter |
| `epsav` | `ε_⊥/S0` | Scaled perpendicular permittivity |
| `deleps` | `Δε/S0` | Scaled permittivity anisotropy |
| `efe` | `ξ_a` | Linear flexoelectric coefficient |
| `efe2` | `ξ_b` | Quadratic flexoelectric coefficient |

