# Configuration System

This document describes how qlc3d reads its settings (`.qfg`/`.txt`) file, how the
values are validated and converted into in-memory objects, and how those objects are
consumed by the rest of the simulation. It covers the full path from the command line
to the solvers.

## 1. Entry point and file discovery

`main()` in `qlc3d/src/main-app-qlc3d.cpp` creates a `Configuration` object and calls
`parseArgs(argc, args, configuration)`:

- `args[1]` (optional) is the path to the settings file. If relative, it is resolved
  against the current working directory. The file must exist or the program aborts
  with `RUNTIME_ERROR`. (`main-app-qlc3d.cpp:30-43`)
- `args[2]` (optional) is a working directory. If given, it must be an absolute,
  existing directory; the process `chdir`s into it before anything else happens.
  (`main-app-qlc3d.cpp:20-27`)
- If no settings file argument is given, `Configuration`'s default path
  `./meshes/test.txt` is used (`configuration.cpp:9`).

`runSimulation()` then calls `configuration.readSettings()`, which does all the actual
parsing (`main-app-qlc3d.cpp:48`).

## 2. `Configuration` and `SettingsReader`

`Configuration` (`qlc3d/includes/configuration.h`, `qlc3d/src/configuration.cpp`) is a
thin container that owns:
- `Simu`, `LC`, `MeshRefinement`, `Electrodes`, `SolverSettings`, `Alignment`,
  `InitialVolumeOrientation`

`Configuration::readSettings()` constructs a `SettingsReader` for
`configuration.settingsFile()` and pulls each parsed object out of it
(`configuration.cpp:15-24`). The `SettingsReader` constructor immediately calls
`read()`, which:

1. Opens the file and asserts it exists (`settings-reader.cpp:19-22`).
2. Constructs a `Reader` (see below) configured as case-insensitive, with string
   values lower-cased, and with `${ENV_VAR}` substitution enabled
   (`settings-reader.cpp:24-30`).
3. Calls, in this order: `readSimu`, `readLC`, `readAlignment`, `readRefinement`,
   `readElectrodes`, `readSolverSettings`, then constructs an empty
   `InitialVolumeOrientation` and calls `readInitialVolumeOrientation`
   (`settings-reader.cpp:32-41`).

Each `read*` method builds one settings object using a `*Builder` class (where one
exists) and stores it in a `unique_ptr` member. Callers (`Configuration`) retrieve the
constructed objects via `simu()`, `lc()`, `refinement()`, `electrodes()`,
`solverSettings()`, `alignment()`, `initialVolumeOrientation()`; each of these
`std::move`s the object out and can only be called once per parse
(`settings-reader.h:45-51`).

## 3. The `Reader` key/value parser

`Reader` (`qlc3d/includes/reader.h`, header-only) implements the low-level text format:

- One `key = value` assignment per line; `#` starts a comment that runs to end of
  line; blank lines are skipped (`reader.h:237-349`).
- Exactly one `=` per line is required; both key and value must be non-empty after
  trimming whitespace (`splitByChar`, `reader.h:244-269`).
- Keys may not contain a `"` character (`validateKey`, `reader.h:703-708`).
- Duplicate key definitions in the same file are rejected with a `ReaderError`
  reporting the original line number (`reader.h:331-339`).
- If `setCaseSensitivity(false)` is used (as `SettingsReader` does), both keys and
  values are lower-cased while reading. If `setLowerCaseStringValues(true)` is also
  set, values are lower-cased even when case sensitivity is otherwise on
  (`reader.h:316-325`).
- If `setEnvironmentVariableSubstitution(true)` is used, any `${VAR}` sequence in a
  line's value is textually replaced with `getenv("VAR")`; an undefined environment
  variable throws (`reader.h:189-211`).
- Values can be scalars, quoted strings (`"..."`, no embedded whitespace check
  applies once quoted), or bracketed arrays `[a, b, c]` of scalars/strings
  (`reader.h:544-604`). Arrays must have exactly one `[` and one `]`, with `[` first
  and `]` last, and must not be empty (`reader.h:568-604`).
- `Reader::get<T>(key, default)` returns a default when the key is absent;
  `Reader::getOptional<T>(key)` returns `std::optional<T>`;
  `Reader::getValueByKey<T>(key)` throws `ReaderError` if the key is missing or the
  value cannot be converted to `T` (`reader.h:606-654`).
- `Reader::containsKeyWithPrefix(prefix)` is used to detect whether any
  `FIXLC*`/`REFINEMENT*` block is present at all before iterating numbered entries
  (`reader.h:484-491`).
- A separate, unrelated `Reader::readValidKeysFile()` / `isValidKey()` mechanism
  exists to restrict the settings file to a fixed list of valid keys (with `*` as a
  wildcard converted to a regex `.*`), but `SettingsReader` never calls
  `readValidKeysFile()`, so this whitelist feature is currently inactive in
  production use (`reader.h:351-412`).

## 4. Settings file keys

All settings-file key strings are declared as `const std::string SFK_*` constants in
`qlc3d/includes/settings_file_keys.h`. Numbered/wildcarded keys (e.g. `E1.Time`,
`FIXLC3.Anchoring`, `BOX2.Tilt`) use the literal `*` character as a placeholder and the
helper `wildcardToNum(base, n)` to substitute the actual number
(`settings_file_keys.h:92-104`). Note: the `SFK_E_TIME`, `SFK_E_POTS`, `SFK_E_FIELD`,
and all `SFK_FIXLC_*` constants are declared but **not actually used** by
`SettingsReader` — `readElectrodes` and `readAlignment` build the equivalent key
strings from raw literals instead (e.g. `"e" + std::to_string(i) + ".time"`,
`"FIXLC" + to_string(i)`), so these declared constants are effectively dead code
(`settings-reader.cpp:188-189, 286`).

## 5. `Simu` — general simulation settings

Read by `SettingsReader::readSimu` into a `SimuBuilder`
(`settings-reader.cpp:99-156`), then built into an immutable `Simu`
(`qlc3d/includes/simu.h`, `qlc3d/src/simu.cpp`).

| Key | Default | Notes |
|---|---|---|
| `MeshName` | *(required)* | Mesh file, resolved relative to `configuration.currentDirectory()` (`simulation-container.cpp:113`). |
| `dt` | `1e-9` | Initial time step. `simulationMode()` is derived, not stored: `dt > 0` ⇒ `TimeStepping`, otherwise ⇒ `SteadyState` (`simu.h:148`). |
| `QMatrixSolver` | `auto` | One of `auto`/`pcg`/`gmres`. See "Solver settings" below — this setting only affects the (currently dead/commented-out) code in `solve_pcg.cpp`; it has **no effect** on the active Q-tensor solver. |
| `TargetdQ` | `1e-3` | Must be `> 0`. |
| `MaxError` | `1e-3` | Must be `> 0`; Newton iteration convergence criterion, passed to `TimeSteppingLCSolver`. |
| `EndCriterion` | `Time` | One of `iterations`/`time`/`change`. |
| `EndValue` | `1e-3` | Must be `>= 0`. |
| `loadQ` | `""` | Deprecated; mutually exclusive with `loadOrientation` (throws `RUNTIME_ERROR` if both set) — see `orientation-loader.cpp:10-17`. |
| `loadOrientation` | `""` | Preferred replacement for `loadQ`. |
| `LoadInitialOrientationS0Mode` | `file` | Controls whether file-based initial orientation loads keep the file’s scalar order (`file`, legacy default) or replace it with the active material equilibrium `S0` (`current`, opt-in). |
| `saveDir` | `"res"` | Resolved to an absolute path at build time as `workingDir / saveDir` (`simu.cpp:83`), where `workingDir` is the current path *at the time `SimuBuilder` was constructed* (not necessarily the final working directory — see Known Bugs). |
| `dtLimits` | `[1e-9, 1e-4]` | Array of exactly 2 values, both `> 0`, high `>` low is only asserted as high `> 0` (label says should be `> min`, but code does not check this — see Known Bugs). |
| `dtFunction` | `[0.5, 0.8, 1.2, 10]` | Array of exactly 4 values; no validation performed (`SimuBuilder::dtFunction`, TODO comment in code). |
| `stretchVector` | `[1, 1, 1]` | Array of exactly 3 values, each `> 0`. |
| `RegularGridSize` | `[0, 0, 0]` | Array of exactly 3 values; either all zero or all `> 0`. |
| `SaveFormat` | `{}` (empty set) | Array of strings, each matched case-insensitively against `Simu::VALID_SAVE_FORMATS = {lcview, regularvtk, regularvecmat, dirstackz, lcviewtxt, csvunstructured, vtkunstructuredasciigrid}`. |
| `MeshElementOrder` | `native` | One of `native`/`quadratic`/`linear`. |
| `outputEnergy` | `0` | Boolean-as-int; if truthy, `SimulationContainer` opens `energy.csv` in the save directory. |
| `outputFormat` | `0` (binary) | Stored and exposed via `getOutputFormat()`, but **no code in the repository reads this value** — see Known Bugs. |
| `SaveIter` | `0` | Must be `>= 0`. |
| `SaveTime` | `0` | Must be `>= 0`. |
| `NumAssemblyThreads` | `0` | Also copied into `SolverSettings::nThreads` (see below) — the same key is read twice, once into `Simu` and once into `SolverSettings`. |
| `NumMatrixSolverThreads` | `0` | Stored in `Simu` as `numMatrixSolverThreads_`/`getMatrixSolverThreadCount()`, but **no code in the repository ever calls `getMatrixSolverThreadCount()`** — see Known Bugs. |

`Simu::getSaveFormatStrings()` converts the parsed `set<SaveFormats>` bitfield back to
strings for logging (`simu.cpp:53-59`). `PotentialConsistency` (`simu.h:19`) and the
`QMatrixSolvers` enum's `Auto` value are declared but `Auto` is never distinguished from
any other case in the code that reads it (again see Known Bugs, since the only
consumer is dead code).

## 6. `LC` — liquid crystal material parameters

Read by `SettingsReader::readLC` into an `LCBuilder`
(`settings-reader.cpp:158-175`), all fields optional with defaults tuned for the 5CB
liquid crystal (`qlc3d/includes/lc.h:75-87`):

| Key | Default |
|---|---|
| `K11`, `K22`, `K33` | `10e-12` |
| `K24` | `0` |
| `p0` | `0` (achiral) |
| `A` | `-1.2e5` |
| `B` | `-2.1333e6` |
| `C` | `1.7333e6` |
| `eps_par` | `18.5` |
| `eps_per` | `7.0` |
| `e1`, `e3` | `0` |
| `gamma1` | `0.0777` |

`LC`'s constructor immediately derives several quantities used by the solvers, all as a
function of `A`, `B`, `C` via `S0` (`qlc3d/src/lc.cpp:10-53`):

- `S0 = (-B + sqrt(B² - 24·A·C)) / (6·C)`. This must lie in `[0, 1]`, otherwise the
  constructor throws `std::invalid_argument` — this is the *only* validation
  performed on `A`/`B`/`C`; a negative discriminant produces a `NaN` value, which is
  caught by the same `isnan` check but with a generic error message (`lc.cpp:44-53`).
- `L1 = 2(K33 - K11 + 3·K22) / (27·S0²)`
- `L2 = 4(K11 - K22) / (9·S0²)`
- `L3 = 4·K24 / (9·S0²)`
- `L4 = 8·q0·K22 / (9·S0²)` where `q0 = 2π / p0`, or `0` if `p0 == 0`
- `L5 = 0` (unconditionally, no settings key controls this)
- `L6 = 4(K33 - K11) / (27·S0³)`
- `u1 = 2·gamma1 / (9·S0²)`
- `deleps() = (eps_par - eps_per) / S0` (computed on demand, not cached)

These `L1..L6`, `u1`, `S0` values (not the raw `K11`/`A`/`B`/`C`/... values) are what the
LC free-energy and Q-tensor assembly code actually consumes.

## 7. `Electrodes` — potentials and electric field

Read by `SettingsReader::readElectrodes` (`settings-reader.cpp:177-209`):

- `EField = [x, y, z]` (optional, exactly 3 components): if present, `Electrodes` is
  constructed via `withConstantElectricField`, and `hasElectricField()` becomes true.
  In this mode `getCurrentPotentials()` always returns an empty map and
  `isPotentialCalculationRequired()` returns `false` — the potential solver is
  skipped entirely and a fixed, spatially-uniform E-field is used instead
  (`electrodes.cpp:122-157`).
- Otherwise, for `i` in `1..99`, keys `E{i}.Time` (array of times) and `E{i}.Pot`
  (array of potentials, same length) define an `Electrode`. An electrode number `i`
  with no `Time`/`Pot` keys is simply skipped (not an error). `EField` and electrode
  potentials are mutually exclusive by construction (only one `Electrodes` factory is
  called).
- `Electrode`'s constructor requires `times.size() == potentials.size()` and that
  `times` is sorted ascending, or it throws (`electrodes.cpp:13-35`). If the first
  time is `> 0`, a `(0, 0)` pair is logged as a warning and prepended — see Known Bugs
  for a related edge-case in this logic.
- `Electrode::getPotentialAtTime(t)` performs a step/zero-order-hold lookup: it
  returns the potential of the latest switching time `<= t` (not interpolated)
  (`electrodes.cpp:37-53`). Querying a negative time, or a time before the first
  defined switching time, throws via `RUNTIME_ERROR`.
- `Electrode::createSwitchingEvents()` creates one `Event` per defined switching
  time; these are consumed by `EventList`/`SimulationContainer` to actually change
  the applied potential during the simulation (`electrodes.cpp:55-63`,
  `simulation-container.cpp:88`).
- `eps_dielectric` (relative permittivity of dielectric regions) has a getter/setter
  (`getDielectricPermittivity`/`setDielectricPermittivities`) and a hard-coded
  default of `{1.0}` set in the `Electrodes()` default constructor
  (`electrodes.cpp:71-73`), but **`SettingsReader::readElectrodes` never reads any
  `eps_dielectric*` key from the settings file**, so this value can currently only be
  changed by code that calls `setDielectricPermittivities()` directly, not via the
  settings file, despite `SFK_EPS_DIELECTRIC` being declared in
  `settings_file_keys.h` — see Known Bugs.

## 8. `Alignment` / `Surface` — FIXLC anchoring surfaces

Read by `SettingsReader::readAlignment` (`settings-reader.cpp:281-358`). For each
`FIXLC{i}` (`i` in `1..99`) block present, `FIXLC{i}.Anchoring` selects one of the
following (case-insensitive after lower-casing by the `Reader`), each mapped to a
`Surface` factory method in `qlc3d/src/alignment.cpp`:

| `Anchoring` value | Required keys | Behavior |
|---|---|---|
| `strong` | `Easy` (array) | If all elements of `Easy` parse as strings that are not purely numeric (see the `isValueArrayOfNumbers` bug below — in practice this branch is effectively always taken for numeric literals too), `Easy` is treated as `[tiltExpr, twistExpr]` string expressions evaluated per-point via `CartesianExpression` in **absolute** mesh coordinates. Otherwise (this branch is currently dead — see Known Bugs) `Easy[0]`/`Easy[1]` are read directly as fixed tilt/twist angles in degrees. In both cases, **only the first two elements of `Easy` are used**; a third element (documented in example settings files as a "rotation" angle) is silently ignored (`ofStrongAnchoring` always passes `rotDegrees = 0`, `alignment.cpp:220-232`). |
| `homeotropic` | none | Strong anchoring with the surface normal as the easy direction; `overrideVolume` is ignored (always `true`) since "it doesn't make sense in the strong homeotropic case" (`settings-reader.cpp:320-322`). |
| `weak` | `Easy` (2 or 3 elements), `Strength`, `K1`, `K2` | Fixed tilt/twist angles only (`Easy[0]`, `Easy[1]`); a 3rd `Easy` element is accepted by the size check (`== 2 || == 3`) but then never read (`settings-reader.cpp:324-332`, `ofWeakAnchoring`, `alignment.cpp:299-311`). |
| `degenerate` | `Strength` | If `Strength >= 0`, planar-degenerate anchoring (surface normal repels the director); if `Strength < 0`, silently reinterpreted as `ofWeakHomeotropic(i, -Strength, ...)` instead (`settings-reader.cpp:334-344`). |
| `weakhomeotropic` | `Strength` | Weak homeotropic anchoring. |
| `freeze` | none | "Freezes" whatever LC orientation the nodes already have at initialisation time; `overrideVolume` is forced to `false` (`ofFreeze`, `alignment.cpp:313-330`). |
| *(any other value)* | — | Throws `ReaderError` ("Invalid anchoring type ... This may be a typo ... or it has not yet been implemented"). |

Common optional key: `FIXLC{i}.overrideVolume` (bool) — whether the surface's easy
orientation overwrites the initial LC orientation at solver startup, versus leaving
whatever the box/loaded orientation already assigned. Defaults are type specific (see
`ofStrongAnchoring`/`ofWeakAnchoring`/... default parameter values); `homeotropic` and
`freeze` ignore this key entirely.

`AnchoringType::Polymerise` is defined in the `AnchoringType` enum
(`alignment.h:21-23`) but has no corresponding settings-file value, no factory method,
and no handling in `Surface`'s per-type logic anywhere else in the codebase — it is
unreachable from the settings file today.

Surface tilt/twist expressions (used only by `strong` anchoring) are evaluated with
`(x, y, z)` in **absolute mesh coordinates** (in the internal, SI unit system —
meters), unlike `Box` expressions below which use box-local normalized coordinates.

## 9. `InitialVolumeOrientation` / `Box` — initial bulk LC orientation

Read by `SettingsReader::readInitialVolumeOrientation` (`settings-reader.cpp:360-415`).
For each `BOX{i}` (`i` in `1..99`) whose `.Type` key is present, a `Box` is built via
`BoxBuilder`:

| Key | Default | Notes |
|---|---|---|
| `Box{i}.Type` | `Normal` | One of `Normal`/`Random`/`Hedgehog` (`Box::VALID_TYPES`, `box.cpp:17`). |
| `Box{i}.Params` | `[]` | Free-form parameter array; `Params[0]` (if present) is used as an exponent (power) in the tilt/twist formula below, default `1.0` (`Box::getParam`, `box.cpp:119-125,164-186`). |
| `Box{i}.X`, `.Y`, `.Z` | `[0, 0]` | Each must be exactly length 2 (min, max); define the box's `AABox` bounding region. |
| `Box{i}.Tilt`, `.Twist` | `[0, 0]` (as numbers) or `"0"` (as an expression) | If the value is a numeric array (checked via `Reader::isValueArray`, not the buggy `isValueArrayOfNumbers`), it must be length 2: `[base, gradient]`. If it is a string, it is treated as a `CartesianExpression`. |

For a numeric `Tilt = [base, gradient]` (and similarly `Twist`), the tilt/twist angle
at a point `p` inside the box is:

```
zNormalised = (p.z - boxZMin) / (boxZMax - boxZMin)
tiltDegrees = base + pow(zNormalised * gradient, Params[0] or 1.0)
```

i.e. only the z-coordinate is used for the numeric (non-expression) form; x/y
normalized coordinates are computed but unused in that path (`box.cpp:164-186`).
Expression-form `Tilt`/`Twist` are evaluated with `(x, y, z)` **normalized to `[0, 1]`
within the box's own bounding box** — the opposite convention from `Surface`
expressions, which use absolute coordinates (see section 8).

`InitialVolumeOrientation::setVolumeQ` iterates all mesh nodes and, for each box, sets
the Q-tensor at nodes whose coordinates fall inside that box's `AABox`
(`box.cpp:145-162`). Boxes are applied in the order they were added (i.e. increasing
`BOX{i}` number); a node inside multiple overlapping boxes ends up with whichever box
processed it last.

## 10. `MeshRefinement` — adaptive mesh refinement configuration

Read by `SettingsReader::readRefinement` (`settings-reader.cpp:211-252`):

- `RepRefIter` (default `0`): iteration period for periodically-recurring refinement.
- `RepRefTime`: explicitly **not supported** — if present and non-zero, throws
  `std::invalid_argument` telling the user to use `RepRefIter` instead
  (`settings-reader.cpp:219-224`).
- `REFINEMENT{i}` (`i` in `1..99`, only scanned if any `REFINEMENT*` key exists,
  and stopping at the first missing `i`): each defines a `RefinementConfig` with
  `.Type` (required: `Change`/`Sphere`/`Box`, case-insensitive), and optional
  `.X`/`.Y`/`.Z`/`.Iterations`/`.Times`/`.Values` arrays. Validation rules per type are
  enforced in `RefinementConfig::validate()` (`meshrefinement.cpp:9-36`):
  - `Change`: requires non-empty `Values`.
  - `Sphere`: requires non-empty `Values`, `X`, `Y`, `Z`.
  - `Box`: requires non-empty, equal-length `X`/`Y`/`Z`, and an even count (pairs
    per box region).
- A `RefinementConfig` with empty `Iterations` and empty `Times` is considered
  "periodic" (`occursPeriodically()`, `meshrefinement.h:42`) and is only actually
  scheduled if `RepRefIter > 0`; `createMeshRefinementEvents`
  (`qlc3d/src/inits.cpp:240-` ff.) throws `RUNTIME_ERROR` if `RepRefIter > 0` but no
  periodic `REFINEMENT*` object exists (`inits.cpp:261-265`).
- Explicit (non-periodic) configs are expanded into one `RefinementSpec` per listed
  iteration and per listed time (`RefinementConfig::toSpecs()`,
  `meshrefinement.cpp:42-55`), each becoming a distinct simulation `Event`.

## 11. `SolverSettings` — linear solver tuning

Read by `SettingsReader::readSolverSettings` (`settings-reader.cpp:254-279`); defaults
set in the `SolverSettings()` constructor (`qlc3d/src/solver-settings.cpp:4-28`):

| Key | Default | Actually consumed by |
|---|---|---|
| `NumAssemblyThreads` | `1` | `main-app-qlc3d.cpp` calls `omp_set_num_threads(solverSettings->getnThreads())` (`main-app-qlc3d.cpp:62-63`). Note this is the *same settings key* also stored separately on `Simu` (see section 5). |
| `Q_Newton_Panic_Iter` | `10` | Passed to `TimeSteppingLCSolver`'s constructor (`main-app-qlc3d.cpp:60`). |
| `Q_Newton_Panic_Coeff` | `0.1` | Stored, but no call site outside `solver-settings.cpp` reads `getQ_Newton_Panic_Coeff()`. |
| `Q_GMRES_Maxiter` | `2000` | Used by the active Q-tensor GMRES solve in `lc-solver.cpp:40` (`0` is treated specially: it falls back to the matrix's column count instead of literally 0 iterations). |
| `Q_GMRES_Restart` | `100` | Used in `lc-solver.cpp:41`. |
| `Q_GMRES_Toler` | `1e-7` | Used in `lc-solver.cpp:42`. |
| `Q_GMRES_Preconditioner` | `LUinc` (`2`) | Stored; not read outside `solver-settings.cpp` in the current codebase's active (non-commented) Q-solver path. |
| `Q_Solver` | `Q_Solver_GMRES` (`1`) | Stored via `getQ_Solver()`/`setQ_Solver()`, but **no active code path reads `getQ_Solver()`**. The actual Q-tensor solver choice (PCG vs GMRES) is made in `lc-solver.cpp` based on whether the Q-tensor matrix is symmetric (which follows from whether `p0 == 0`, i.e. achiral vs chiral), not from this setting. See Known Bugs. |
| `Q_PCG_Preconditioner`, `Q_PCG_Maxiter`, `Q_PCG_Toler` | `Diagonal`(`0`), `2000`, `1e-7` | Stored; not read by any active (non-commented-out) code. |
| `V_GMRES_Maxiter` | `2000` | Used in `potential-solver.cpp:450`. |
| `V_GMRES_Restart` | `100` | Used in `potential-solver.cpp:451`. |
| `V_GMRES_Toler` | `1e-7` | Used in `potential-solver.cpp:454`. |
| `V_GMRES_Preconditioner` | `LUinc` (`2`) | Stored; not read by any active code in the current codebase. |
| `V_Solver` | `V_Solver_GMRES` (`1`) | Stored via `getV_Solver()`/`setV_Solver()`, but **no active code path reads `getV_Solver()`** — the potential solver always uses GMRES (`potential-solver.cpp:448-457`). |
| `V_PCG_Preconditioner`, `V_PCG_Maxiter`, `V_PCG_Toler` | `Cholinc`(`1`), `2000`, `1e-7` | Stored; not read by any active code — the potential solver never uses PCG. |

In short: **only the `Q_GMRES_*` and `V_GMRES_*` settings currently have any effect**
on the running simulation. `Q_Solver`, `V_Solver`, all `*_PCG_*` settings, and both
`*_GMRES_Preconditioner` settings are parsed and validated but not consulted by any
active solver code path. See `known-bugs.md` for details.

## 12. Where settings end up being used

`SimulationContainer::initialise()` (`qlc3d/src/simulation-container.cpp`) is the main
consumer of the parsed configuration objects once `main-app-qlc3d.cpp` has built the
solver/event/state objects:

- `simu`, `lc`, `electrodes` are re-fetched from `Configuration` at the start of
  `initialise()` (`simulation-container.cpp:50-52`).
- Working directory is switched to `configuration.currentDirectory()` again
  (`simulation-container.cpp:54-58`).
- The results directory (`simu->getSaveDirAbsolutePath()`) is created if missing, and
  the settings file that was used is copied into it as `settings.qfg` for
  reproducibility (`simulation-container.cpp:61-107`).
- `createMeshRefinementEvents(*configuration.refinement(), eventList)` and
  `createElectrodeSwitchingEvents(*electrodes, eventList)` convert the parsed
  refinement/electrode configuration into scheduled simulation `Event`s
  (`simulation-container.cpp:87-88`).
- The mesh is loaded and geometry prepared using `simu->meshName()`,
  `simu->getStretchVector()`, and `simu->getMeshElementOrder()`
  (`simulation-container.cpp:113-120`).
- A regular grid is (optionally) built from `simu->getRegularGridXCount/YCount/ZCount()`
  (`simulation-container.cpp:121-124`); if any requested output format requires a
  regular grid but the grid size was never configured, a `std::runtime_error` is
  thrown (`simulation-container.cpp:80-82`).
- `configuration.getInitialVolumeOrientation()` (the parsed `Box`es) and
  `configuration.getAlignment()` (the parsed `FIXLC` `Surface`s) are used together to
  set the initial Q-tensor, first from the boxes, then optionally overwritten by a
  loaded orientation file (`simu->getLoadQ()`/`getLoadOrientation()`), then by the
  surfaces (`inits.cpp`, function `initialiseLcSolutionVector`, called from
  `simulation-container.cpp:149`).
