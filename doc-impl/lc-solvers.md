# LC solver architecture and extension points

The active LC solver stack in qlc3d is deliberately split into three layers:

- `ILCSolver` is the public strategy interface for the simulation loop.
- `ImplicitLCSolver` holds the common FEM assembly and linear-solver plumbing.
- Concrete solver classes implement the nonlinear algorithm itself.

This design keeps the solver selection at the simulation boundary while reusing the same local matrix assembly code for multiple algorithms.

## 1. The public interface

`ILCSolver` in `qlc3d/includes/lc/lc-solver.h` defines the single required entry point:

```cpp
virtual LCSolverResult solve(
    SolutionVector &q,
    const SolutionVector &v,
    const Geometry &geom,
    SimulationState &simulationState) = 0;
```

Every concrete LC solver returns a `LCSolverResult` with:

- `solverType` (`STEADY_STATE` or `TIME_STEPPING`)
- `iterations`
- `dq` (largest Q-tensor increment)
- `converged`
- `maxIterationsReached`
- `elapsedTimes` (assembly and solve timing)

The simulation does not know or care which concrete algorithm is used; it only calls the `ILCSolver` interface.

## 2. Shared implementation and matrix solve helper

`ImplicitLCSolver` contains the common state and assembly helpers used by both concrete implementations:

- `K` – global Jacobian matrix
- `L` – Euler-Lagrange right-hand side vector
- `X` – solution increment vector
- `assembleMatrixSystemVolumeTerms()`
- `assembleMatrixSystemWeakAnchoring()`
- `assembleMatrixSystem()`
- `solveMatrixSystem()`

The most important detail is the linear-solver policy inside `solveMatrixSystem()` (`lc-solver.cpp`):

- if `isSymmetricMatrix` is true, use `PCG` with a diagonal preconditioner
- otherwise, use `GMRES` with an LU-based preconditioner

The code sets:

```cpp
isSymmetricMatrix(lc.p0() == 0.0);
```

so the actual linear-solver choice is inferred from the material parameters, not from the `Q_Solver` settings key. In other words, the matrix symmetry is the runtime decision point, while the settings file is only a place where the values are stored.

## 3. Available concrete solvers today

The active solver classes are created in `qlc3d/src/main-app-qlc3d.cpp`:

```cpp
unique_ptr<ILCSolver> lcSolver =
    simu->simulationMode() == SteadyState ?
        unique_ptr<ILCSolver>(new SteadyStateLCSolver(...)) :
        unique_ptr<ILCSolver>(new TimeSteppingLCSolver(...));
```

This means the solver choice is made at the application level by the simulation mode, not by a polymorphic factory inside the solver layer.

### SteadyStateLCSolver

`SteadyStateLCSolver` performs a single Newton step per call:

1. assemble `K` and `L` from the current Q-tensor field
2. solve `K * ΔQ = L`
3. update `Q <- Q - ΔQ`

This is a "one-shot" nonlinear solve and is selected when the simulation is in steady-state mode (`dt = 0` in practice).

### TimeSteppingLCSolver

`TimeSteppingLCSolver` is the implicit time-stepping algorithm. It adds a mass-matrix term and performs a Newton loop within each time step:

1. predict the next state using the previous `dQ/dt`
2. assemble `K` and `L` for the current guess
3. modify the system with the implicit time-stepping terms
4. solve the modified system for the Newton increment
5. update `Q <- Q - ΔQ`
6. repeat until `ΔQ` falls below tolerance or the Newton iteration cap is reached

This solver also stores the previous RHS and previous Q-tensor state to compute the next time-step predictor.

## 4. Common FEM assembly pattern

Both solver classes share the same local-element assembly flow defined in `ImplicitLCSolver`:

- process each tetrahedral element in the LC region
- build the local `lK` and `lL` for that element
- accumulate contributions from thermotropic, elastic, chiral, dielectric, and weak-anchoring terms
- insert the local matrix into the global sparse matrix using `addToGlobalMatrix()`

This is the reason the shared base class is useful: the nonlinear algorithm changes, but the element-level physics and FEM bookkeeping remain the same.

## 5. How to add another solver in the future

The extension pattern is straightforward:

1. Add a new concrete class implementing `ILCSolver`.
2. If it reuses the standard local assembly and matrix solve logic, inherit from `ImplicitLCSolver`.
3. Implement `solve(...)` to do the nonlinear algorithm and return a `LCSolverResult`.
4. Instantiate it in `main-app-qlc3d.cpp` behind the same `ILCSolver` interface.
5. Keep the physics assembly in the shared base rather than duplicating it in each solver.

A minimal pattern looks like this:

```cpp
class NewLCSolver final : public ILCSolver, protected ImplicitLCSolver {
public:
  LCSolverResult solve(SolutionVector &q,
                       const SolutionVector &v,
                       const Geometry &geom,
                       SimulationState &simulationState) override {
    // assemble K/L, call solveMatrixSystem(), update q, and return LCSolverResult
  }
};
```

The extension point is intentionally narrow: the simulation loop only sees `ILCSolver`, while the numerical details live in the concrete class and the shared `ImplicitLCSolver` helpers.

## 6. Current state of the solver documentation

The implementation is reasonably well described in `doc-impl/algorithm-flow.md`, but the solver architecture itself is spread across a few files:

- `lc-solver.h` defines the interface and concrete classes
- `lc-solver.cpp` contains the actual nonlinear logic and matrix assembly
- `algorithm-flow.md` explains the high-level flow of the steady-state and time-stepping paths
- `known-bugs.md` records the places where the documented settings API does not match the live implementation

The key point for maintainers is that the real runtime solver choice in qlc3d is:

- choose the nonlinear algorithm by simulation mode (`SteadyState` vs time stepping)
- choose the linear backend internally by matrix symmetry (`PCG` or `GMRES`)

and not by the `Q_Solver` / `V_Solver` settings that are still described in older documentation.
