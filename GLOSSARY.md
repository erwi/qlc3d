# qlc3d Glossary

Terms used across the project documentation (`doc-impl` and `qlc3d/doc`), grouped by topic. The goal is a single, project-wide ubiquitous language. Terms are alphabetically ordered within each section.

## Simulation control and settings

- **dt** : Initial time step duration in seconds. When set to `0`, the steady-state Newton-Raphson solver is used.
- **dtFunction** : Four-element array controlling adaptive time-step scaling as a piecewise-linear function of `R = dQ/TargetdQ`.
- **dtLimits** : Two-element vector `[mindt, maxdt]` clamping the adaptive time step.
- **EndCriterion** : Settings key selecting how a simulation ends: `Iterations`, `Time`, or `Change`.
- **EndValue** : Numerical threshold paired with `EndCriterion`.
- **Event** : Scheduled action in the simulation loop (electrode switching, mesh refinement, or result save). Managed by `EventList`.
- **MaxError** : Accuracy parameter for the nonlinear time-stepping Newton iterations.
- **numAssemblyThreads** : Number of OpenMP threads used in FEM matrix assembly.
- **NumMatrixSolverThreads** : Number of threads used by the linear solvers.
- **Settings file** : Text file (`.qfg`) of `key = value` pairs controlling the simulation.
- **Simu** : Class holding simulation mode, time-stepping, output, and refinement settings.
- **TargetdQ** : Target maximum Q-tensor change per time step used by adaptive time stepping; also used as Newton damping when `dt = 0`.

## Liquid-crystal physics

- **Director** : The unit vector **n** describing the average orientation of LC molecules; related to the Q-tensor by `Q_ij = S/2 (3 n_i n_j − δ_ij)`.
- **Easy axis / Easy direction** : Preferred surface orientation specified by `FIXLCn.Easy` as tilt, twist, and rotation angles.
- **Elastic energy** : Distortion energy density `f_D` computed from the G1/G2/G4/G6 invariants (L1..L6 Landau-de Gennes coefficients).
- **Electric energy** : Field-induced energy density `f_E` from dielectric anisotropy and flexoelectric polarisation.
- **Flexoelectric polarisation** : `P_i = ξ_a Q_ij,j + ξ_b Q_ij Q_jk,k`; contributes to the potential RHS via bound charge density `ρ = P_i,i`.
- **Frank constants** : Elastic constants `K11` (splay), `K22` (twist), `K33` (bend); mapped internally to Landau-de Gennes coefficients L1..L6.
- **Free energy** : Total LC energy `F = ∫(f_D + f_B − f_E)dΩ + ∫f_S dΓ`; optionally written to `energy.csv` when `outputEnergy = 1`.
- **LC** : Liquid crystal; also the `LC` class that holds material constants.
- **Order parameter** : Scalar `S` (or tensor `Q`) measuring the degree of LC alignment.
- **Q-tensor** : Traceless symmetric tensor describing LC order, expanded in the five-dimensional T-basis as `Q = Σ q_i T_i`.
- **Surface anchoring energy** : Rapini-Papoular style energy density `f_S` at weak-anchoring surfaces.
- **T-basis / T-tensor basis** : Five orthonormal basis tensors `T_1`..`T_5` used to represent the Q-tensor; C++ variables `q1`..`q5` correspond to SymPy symbols `t1`..`t5`.
- **Thermotropic energy** : Bulk ordering energy density `f_B` from Landau-de Gennes coefficients A, B, C.
- **Tilt** : Polar angle component of the director or easy-axis orientation.
- **Twist** : Azimuthal angle component of the director or easy-axis orientation.

## Materials and boundary conditions

- **Anchoring** : Surface condition that constrains the LC director at a boundary. Implemented via `FIXLCn` settings and classified as *strong*, *weak*, *homeotropic*, *degenerate*, *freeze*, *polymerise*, or *manualNodes*.
- **Degenerate anchoring** : Anchoring type where the LC is free to rotate in the local surface plane; a negative `Strength` turns it into weak homeotropic anchoring.
- **Dielectric** : Non-LC volume region (`Dielectric1`..`Dielectric7`) with a relative permittivity set via `eps_dielectric`.
- **Domain** : Volume material region. `Domain1` is the LC region; `MAT_DOMAIN7` is used internally as the upper bound for LC domain material numbers.
- **Electrode** : Surface material (`Electrode1`..`Electrode9`) with a time-dependent potential defined by `En.Time` and `En.Pot`.
- **FixLC** : Surface material group (`FixLC1`..`FixLC99`) used to define LC anchoring conditions.
- **Homeotropic anchoring** : Anchoring type that aligns the LC along the local surface normal.
- **ManualNodes anchoring** : Anchoring type that uses `FIXLCn.Params` as explicit mesh node indices to fix the LC orientation.
- **MAT_DOMAIN7** : Internal upper bound on material numbers that identify LC domain volume elements; used to filter elements in volume integrals.
- **MAT_PERIODIC** : Internal material number for periodic boundary surface triangles.
- **Neumann** : Boundary condition for the potential solver and "free surface" for the LC director. *Note: `mesh.md` spells this "Neuman" in one place.*
- **Periodic boundary** : Translational symmetry boundary; supported in x, y, and/or z directions for axis-aligned cuboid meshes.
- **Strong anchoring** : Anchoring type that fixes the LC orientation to the easy direction.
- **Weak anchoring** : Anchoring type with finite strength `FIXLCn.Strength` allowing deviation from the easy direction.

## Mesh, geometry, and refinement

- **Adaptive mesh refinement (AMR)** : The geometry-only tetrahedral refinement performed during a simulation; no de-refinement/coarsening is supported. See also *red-green refinement*.
- **Autoref** : High-level adaptive refinement entry point (`autoref()` in `autorefinement.cpp`) that orchestrates element selection, splitting, geometry rebuild, and Q-tensor interpolation.
- **Barycentric coordinates** : Local coordinates (`loc`) within a tetrahedron used for interpolation and containment tests; for TET4 there are four weights, for TET10 the reference point `(ξ, η, ζ)` is derived from them.
- **Box (refinement)** : `REFINEMENTi.Type = Box` selects tetrahedra whose nodes lie inside an axis-aligned cuboid region.
- **Change (refinement)** : `REFINEMENTi.Type = Change` selects tetrahedra where the maximum absolute change in any Q-tensor component exceeds a threshold.
- **Geometry** : Object holding `Coordinates`, tetrahedral and triangle `Mesh` objects, node normals, and (via `SimulationContainer`) an optional `RegularGrid`.
- **Gmsh** : Primary supported mesh format (ASCII 4.1); element types include TET4 (`4`), TET10 (`11`), TRI3 (`2`), and TRI6 (`9`).
- **Green refinement** : Conforming refinement of tetrahedra that share bisected edges but are not fully split. Sub-types are *green-1*, *green-2*, and *green-3*.
- **MeshElementOrder** : Settings key controlling element order policy: `native`, `quadratic`, or `linear`.
- **MeshRefinement** : Container on the `Simu` object holding `RefinementConfig`s and the periodic repetition period `RepRefIter`.
- **Red refinement** : Full 1→8 split of a tetrahedron when all six edges are bisected.
- **Red-green refinement** : Tetrahedral refinement algorithm combining full red splits with conforming green splits.
- **RefinementConfig** : Plain data read from one `[REFINEMENT]` settings block.
- **RefinementSpec** : Unified runtime descriptor for a refinement region, created from `RefinementConfig`.
- **Regular grid** : Uniformly spaced Cartesian interpolation grid built over the tetrahedral mesh for fast field output.
- **RepRefIter** : Settings key defining the iteration period for repeating mesh refinement.
- **Sphere (refinement)** : `REFINEMENTi.Type = Sphere` selects tetrahedra whose nodes lie inside a sphere.
- **TET4 / TET10** : Linear and quadratic tetrahedral volume elements (4 and 10 nodes).
- **TetMeshSearch** : Self-contained spatial search over a tetrahedral mesh.
- **TRI3 / TRI6** : Linear and quadratic triangle surface elements (3 and 6 nodes).

## Solvers and numerics

- **Assembly** : Element-by-element construction of the sparse FEM matrix `K` and right-hand side `L`.
- **dQ** : Largest change in any Q-tensor component between successive iterations; used as a convergence measure and to drive adaptive time stepping.
- **Gauss point** : Integration point inside an element where shape functions, densities, and Jacobians are evaluated.
- **GMRES** : Generalised Minimal Residual linear solver used when the system matrix is non-symmetric.
- **ILCSolver** : Interface for LC solvers; concrete implementations are `SteadyStateLCSolver` and `TimeSteppingLCSolver`.
- **ILU** : Incomplete LU preconditioner.
- **Jacobi** : Diagonal preconditioner.
- **Newton-Raphson** : Steady-state solver mode used when `dt = 0`.
- **PCG** : Preconditioned Conjugate Gradient linear solver used when the system matrix is symmetric.
- **PotentialSolver** : Class that assembles and solves the Poisson equation for the electric potential.
- **SolutionVector** : Node-indexed array storing either Q (5 DOF/node) or V (1 DOF/node), with a DOF map distinguishing free from fixed nodes.
- **SteadyStateLCSolver** : LC solver that performs a single Newton step per call.
- **TimeSteppingLCSolver** : LC solver using an implicit predictor-corrector time-stepping scheme.

## Initial conditions and boxes

- **Hedgehog** : `BOXn.Type = Hedgehog` creates a +1 hedgehog defect at the box centre.
- **InitialVolumeOrientation** : Settings subsystem for defining initial LC orientation via `BOXn` structures (the term appears in the user docs alongside `BOX`).
- **LoadQ** : Settings key for loading a previously saved Q-tensor result file as the initial condition.
- **Normal (box type)** : `BOXn.Type = Normal` sets a uniform or spatially varying director orientation.
- **Random (box type)** : `BOXn.Type = Random` sets a randomized director orientation within the box.

## Output and I/O

- **DirStackZ** : CSV output format saving director components on a regular grid, stacked along Z.
- **LcEnergyCalculator** : Stateless top-level object that computes the optional LC free-energy diagnostic.
- **LcView** : External viewer program and the default binary result format (`LCview`).
- **RegularVecMat** : Output format saving regular-grid data as MATLAB-loadable text files.
- **RegularVTK** : VTK output format with Q-tensor and potential interpolated onto a regular grid.
- **SaveDir** : Settings key naming the results subdirectory.
- **SaveFormat** : Settings key selecting result output formats (`LCview`, `LCviewTXT`, `RegularVTK`, `RegularVecMat`, `DirStackZ`, `VTKUnstructuredAsciiGrid`).
- **SaveIter** : Iteration frequency for writing intermediate result files.
- **SaveTime** : Time-based frequency for writing intermediate result files.
- **VTK** : Output format family, including `VTKUnstructuredAsciiGrid` for ParaView.

## Noted inconsistencies and imperfections

1. **Energy output filename**: `qlc3d/doc/README.md` says `energy.m`, while `doc-impl/lc-energy-calculation.md` says `<saveDir>/energy.csv`.
2. **Energy setting key case**: `OutputEnergy` in `README.md` vs `outputEnergy` in `doc-impl/lc-energy-calculation.md`.
3. **LCView / LCview / LcView**: The viewer/format is spelled `LCView` (class), `LCview` (settings value), and `LcView` (docs) in different places.
4. **Neumann spelling**: `qlc3d/doc/mesh.md` spells the boundary condition "Neuman" (missing final *n*).
5. **InitialVolumeOrientation vs BOX**: The user docs use `InitialVolumeOrientation` as a subsystem name but the actual settings keys are `BOXn.*`; this may be a leftover from a previous naming scheme.
6. **Time-stepping scheme name**: `README.md` calls it "nonlinear Crank-Nicholson", while `doc-impl/algorithm-flow.md` calls it "implicit Adams–Bashforth-like".
7. **VTK format name**: `README.md` uses `VTKUnstructuredAsciiGrid`, while the implementation writer is `VtkUnstructuredAsciiGridFormatWriter` (camel-case difference).
8. **DirStackZ / DirstackZ**: `README.md` uses both capitalisations in the same section.
9. **dtFunction / DtFunction**: `README.md` uses both cases for the adaptive time-stepping function key.
10. **SaveDir / saveDir**: `README.md` uses `SaveDir`; `doc-impl/lc-energy-calculation.md` uses `saveDir`.
11. **Default values marked as TODO**: `README.md` lists `EndCriterion` default and `EndValue` default as `TODO: check this`.
12. **Flexoelectric sign convention**: `README.md` notes "(which sign convention?)" for `e1`/`e3`; this should be resolved and documented.
13. **Polymerise anchoring**: Described as a "Special undocumented secret feature" in `README.md`; if it is supported it should be properly documented.
14. **Domain numbering**: `mesh.md` documents only `Domain1` as the LC region, while the code/docs use `MAT_DOMAIN7` as the inclusive upper bound for LC domains (`Domain1`..`Domain7`).
15. **Periodic material naming**: User docs use `Periodic`; internal code uses `MAT_PERIODIC`.
16. **Q-tensor variable naming**: C++ uses `q1`..`q5` for T-basis coefficients, but these are also the names of the old Q-tensor components in some contexts; the SymPy derivation uses `t1`..`t5` to avoid ambiguity.
