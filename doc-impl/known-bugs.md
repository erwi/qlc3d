# Known Bugs in the Configuration System

This file lists concrete defects found while investigating qlc3d (see `configuration-system.md` for the full description of the
intended/actual behavior). Each entry is a verified fact from reading the code, with
file:line citations. 

## 1. `SolverSettings::setnThreads` validates the wrong variable

`qlc3d/src/solver-settings.cpp:30-35`:

```cpp
void SolverSettings::setnThreads(int num) {
  if (nThreads < 0) {                 // checks the OLD value, not `num`
    throw std::runtime_error("Number of threads must be 0 or positive");
  }
  nThreads = num;
}
```

The guard checks the current (old) value of `nThreads`, which starts at the default
`1` and can never be negative through this setter, instead of checking the incoming
parameter `num`. As a result, setting `NumAssemblyThreads` to a negative value in the
settings file is silently accepted here, and the error is only raised later, with a
confusing "Number of threads must be 0 or positive" message, the next time
`getnThreads()` is called (`solver-settings.cpp:43-48`), rather than at settings-parse
time.

## 2. `Q_Solver` and `V_Solver` settings are read but never consulted

- `Q_Solver` is parsed in `settings-reader.cpp:259` into
  `SolverSettings::setQ_Solver`. `SolverSettings::getQ_Solver()`
  (`solver-settings.h:57`) has no callers anywhere in `qlc3d/src` or
  `qlc3d/includes` outside of `solver-settings.cpp` itself.
- `V_Solver` is parsed in `settings-reader.cpp:270` into
  `SolverSettings::setV_Solver`. `SolverSettings::getV_Solver()`
  (`solver-settings.h:58`) similarly has no callers anywhere else.

The actual choice of Q-tensor solver (PCG vs GMRES) is made in `lc-solver.cpp` purely
based on whether the assembled matrix is symmetric, which in turn follows from whether
the LC is chiral (`p0 != 0`) — not from the `Q_Solver` setting. The potential solver
(`potential-solver.cpp:448-457`) unconditionally uses GMRES and never checks
`V_Solver`. Users who set `Q_Solver = PCG` or `V_Solver = PCG` in their settings file
will see no effect and get no warning that the setting was ignored.

## 3. All `*_PCG_*` solver settings are dead configuration

`Q_PCG_Preconditioner`, `Q_PCG_Maxiter`, `Q_PCG_Toler`, `V_PCG_Preconditioner`,
`V_PCG_Maxiter`, `V_PCG_Toler` are all parsed in `settings-reader.cpp:259-278` and
stored in `SolverSettings`, but none of their corresponding getters
(`getQ_PCG_Maxiter`, `getQ_PCG_Preconditioner`, `getQ_PCG_Toler`,
`getV_PCG_Maxiter`, `getV_PCG_Preconditioner`, `getV_PCG_Toler`) are called from any
file other than `solver-settings.cpp` itself. The only place in the repository that
references Q-tensor PCG solving with these settings
(`qlc3d/src/solve_pcg.cpp`) is entirely commented out
(the whole function body, `solve_pcg.cpp:9-37`, is inside a `/* ... */` block), so
this file currently compiles to nothing. `Simu::QMatrixSolvers` (`Auto`/`PCG`/`GMRES`,
read from the separate `QMatrixSolver` key) is also referenced only inside this same
dead code (`solve_pcg.cpp:22,29`), so that setting is dead as well.

Also related: `Q_GMRES_Preconditioner` and `V_GMRES_Preconditioner` are parsed and
stored, but neither `getQ_GMRES_Preconditioner()` nor `getV_GMRES_Preconditioner()` is
called from `lc-solver.cpp` or `potential-solver.cpp` (the two places that actually run
GMRES solves), so the configured preconditioner choice has no effect on the active
GMRES solves either — only `Maxiter`, `Restart`, and `Toler` are actually read there.

## 4. `NumMatrixSolverThreads` setting has no effect

`SFK_NUM_MATRIX_SOLVER_THREADS` ("NumMatrixSolverThreads") is parsed in
`settings-reader.cpp:121` into `Simu`'s `numMatrixSolverThreads_` field, exposed via
`Simu::getMatrixSolverThreadCount()` (`simu.h:154`). No code outside `simu.h`/`simu.cpp`
calls `getMatrixSolverThreadCount()`. Only `Simu::getAssemblyThreadCount()` /
`SolverSettings::getnThreads()` (populated from the separate `NumAssemblyThreads` key)
are actually used to call `omp_set_num_threads()` (`main-app-qlc3d.cpp:62-63`).

## 5. `outputFormat` setting has no effect

`SFK_OUTPUT_FORMAT` ("outputFormat") is parsed in `settings-reader.cpp:115` into
`Simu::outputFormat_`, exposed via `Simu::getOutputFormat()` (`simu.h:177`). No code
anywhere in the repository calls `getOutputFormat()`. The field is already flagged in
its own declaration comment as "TODO should be part of list of save formats? Looks
like not used anywhere" (`simu.h:112`), confirming this is a known, still-unresolved
dead setting.

## 6. `eps_dielectric` setting is declared but never read from the settings file

`SFK_EPS_DIELECTRIC` ("eps_dielectric") is declared in `settings_file_keys.h:54`, and
`Electrodes` has a working `getDielectricPermittivity()`/`setDielectricPermittivities()`
API plus a hard-coded single-element default of `{1.0}` set in the default constructor
(`electrodes.cpp:71-73`). However, `SettingsReader::readElectrodes`
(`settings-reader.cpp:177-209`) never reads any `eps_dielectric` (or numbered
`Dielectric{i}.eps` style) key from the settings file, and never calls
`setDielectricPermittivities()`. There is currently no way for a user to configure a
non-default dielectric permittivity through the settings file, even though the data
model and constant for the key both exist.

## 7. `Reader::isValueArrayOfNumbers` checks the key string instead of the value string

`qlc3d/includes/reader.h:445-457`:

```cpp
bool Reader::isValueArrayOfNumbers(const std::string &key) const {
  if (!isValidNumber(key)) {     // BUG: validates `key` (e.g. "FIXLC1.Easy"), not its value
    return false;
  }
  try {
    auto array = getValueByKey<std::vector<double>>(key);
  } catch (...) {
    return false;
  }
  return true;
}
```

`isValidNumber()` is meant to check whether a *value string* looks like a number. Here
it is called on `key`, which is a settings-file key name such as `"FIXLC1.Easy"`. Since
key names contain letters and dots, `isValidNumber(key)` is always `false` for any
realistic key, so `isValueArrayOfNumbers()` always returns `false`, regardless of
whether the value is actually a numeric array.

The only call site is `SettingsReader::readAlignment` for `FIXLC{i}.Anchoring =
strong` surfaces (`settings-reader.cpp:296-316`):

```cpp
if (reader.isValueArrayOfNumbers(key)) {
  // numeric tilt/twist angle branch — currently unreachable
  ...
} else if (reader.isValueArrayOfStrings(key)) {
  // expression branch — always taken instead, even for plain numeric values
  ...
}
```

Because `isValueArrayOfNumbers` always returns `false`, the "plain numeric tilt/twist"
branch for `FIXLC*.Anchoring = Strong` is unreachable code, and every `FIXLC*.Easy`
value — even a purely numeric array like `[5.0, 90.0, 0.0]`, as used in
`examples/switching-dynamics-1d/settings.txt:33` — is instead parsed as a pair of
string *expressions* and evaluated through `CartesianExpression`/`tinyexpr`. This
happens to produce the same numeric result for constant expressions like `"5.0"`, so
the bug is currently benign for typical settings files, but it means the intended fast
path never executes and any settings file relying on the documented numeric-array
behavior is silently and invisibly routed through the expression evaluator instead.

## 8. Third `FIXLC*.Easy` value ("rotation") is accepted but always discarded

For `weak` anchoring, `SettingsReader::readAlignment` explicitly validates that
`FIXLC{i}.Easy` has 2 or 3 elements (`settings-reader.cpp:325`: `assertTrue(easyAngles.size()
== 2 || easyAngles.size() == 3, ...)`), implying a third ("rotation") value is a
supported input. However, only `easyAngles[0]` and `easyAngles[1]` are ever read
(`settings-reader.cpp:326-332`), and `Surface::ofWeakAnchoring` hard-codes the third
`easyAnglesDegrees` component to `0` (`alignment.cpp:299-311`). The same is true for
`strong` anchoring's numeric branch (`ofStrongAnchoring(unsigned int, double, double,
bool)`, `alignment.cpp:220-232`, which only takes `tiltDegrees`/`twistDegrees`, no
rotation parameter). Example settings files (e.g.
`examples/switching-dynamics-1d/settings.txt:32-33`) supply a third value in the
comment ("easy tilt, twist, rotation angles"), but that value has no effect for any
currently supported anchoring type.

## 9. `Electrode` constructor can double-count the initial time-0 sample

`qlc3d/src/electrodes.cpp:13-35`:

```cpp
Electrode::Electrode(unsigned int electrodeNumber, const std::vector<double> &times, const std::vector<double> &potentials) {
  ...
  if (times.empty()) {
    this->times_ = {0};
    this->potentials_ = {0};
  } else if (times[0] > 0) {
    Log::warn(...);
    this->times_ = {0};
    this->potentials_ = {0};
  }
  ...
  for (auto t : times) { this->times_.push_back(t); }
  for (auto p : potentials) { this->potentials_.push_back(p); }
}
```

When `times` is non-empty and `times[0] > 0`, the code first sets `times_`/`potentials_`
to `{0}`/`{0}`, then unconditionally appends *all* of the original `times`/`potentials`
afterwards — this is the intended "insert an implicit t=0, V=0 sample before the
user's data" behavior and works correctly in that case. However, if `times` is
non-empty and `times[0] == 0` (the normal/expected case where the user already
provides a t=0 sample), neither branch reassigns `times_`/`potentials_`, so they stay
at their default-constructed empty state, and the subsequent unconditional
`push_back` loop appends the user's data as-is — this case behaves correctly too. The
logic is fragile/non-obvious (three implicit paths achieving the same append pattern
via different means) but was not found to produce incorrect output for any of the
three cases (`empty`, `times[0] > 0`, `times[0] == 0`) during review; it is listed here
because the structure is easy to misread as double-inserting the `(0,0)` pair and would
become a real bug if the trailing `for` loops were ever changed without noticing the
implicit dependency on the preceding `if`/`else if`.

## 10. Inconsistent coordinate normalization between `Box` and `Surface` expressions

`Box::getTiltAt`/`getTwistAt` (`box.cpp:164-186`) evaluate `Box{i}.Tilt`/`Box{i}.Twist`
string expressions using coordinates normalized to `[0, 1]` within that box's own
bounding box. `Surface::getEasyTiltAngleAt`/`getEasyTwistAngleAt`
(`alignment.cpp:145-159`) evaluate `FIXLC{i}.Easy` string expressions using absolute
mesh coordinates (in meters), with no normalization. This is not a crash-causing bug,
but it is an undocumented inconsistency: an expression written for a `Box` cannot be
reused verbatim for a `Surface`, and vice versa, without accounting for the different
coordinate conventions.

## 11. `AnchoringType::Polymerise` is defined but unreachable

`AnchoringType::Polymerise` is declared in the enum (`alignment.h:22`), and a `TODO`
comment for a future `Surface::ofPolymerise()` factory exists (`alignment.h:108`), but:
- `SettingsReader::readAlignment` has no `else if (type == "polymerise")` branch
  (`settings-reader.cpp:281-358`), so this value can never be produced from a settings
  file — any settings file specifying `FIXLC{i}.Anchoring = Polymerise` hits the final
  `else` branch and throws a `ReaderError` ("Invalid anchoring type").
- `ManualNodes` is likewise defined in the enum but also has no reachable path in
  `SettingsReader::readAlignment` today.

## 12. `Reader::readValidKeysFile` / key whitelist mechanism is unused in production

`Reader` supports restricting the set of allowed settings keys via
`readValidKeysFile()`/`isValidKey()` (`reader.h:351-412`), with wildcard (`*`) support
converted to regex. `SettingsReader::read()` never calls `readValidKeysFile()`
(`settings-reader.cpp:19-42`), so `_validKeys` is always empty and the check at
`reader.h:328` (`if ((_validKeys.size() > 0) && !isValidKey(key))`) never fires in the
real application — any key name at all is currently accepted (and simply ignored if
unrecognized by `SettingsReader`), with no "unknown key" validation performed against
the actual, hard-coded set of keys `SettingsReader` understands.

