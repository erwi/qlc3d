# Index
This directory contains documentation about implementation details split into separate files organised by topic:

- algorithm-flow.md: High-level algorithm flow of the qlc3d simulation from settings parsing through initialisation, main loop, solvers, events, and output.
- configuration-system.md: How the settings (.qfg/.txt) file is discovered, parsed by the Reader/SettingsReader classes, converted into Simu/LC/Electrodes/Alignment/Box/MeshRefinement/SolverSettings objects, and consumed by the simulation.
- known-bugs.md: Concrete, verified defects found, with file:line citations.
- equations.md: Mathematical notation and equations for the electric potential Poisson problem and LC free-energy terms (elastic, thermotropic, electric, surface).
- lc-energy-calculation.md: Configuration, architecture, numerical integration, and output format of the optional LC free-energy diagnostic.
- mesh-refinement.md: Adaptive red-green mesh refinement pipeline, element selection, classification, splitting, and Q-tensor interpolation.
- quadratic-elements.md: Current state of quadratic finite-element support for 10-node tetrahedra (TET10) and 6-node triangles (TRI6).
- regular-grid.md: Regular-grid subsystem architecture, spatial search, interpolation, and ownership after the Part-2 decoupling refactor.
