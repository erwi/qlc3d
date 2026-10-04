# Known Bugs in the Configuration System

This file lists concrete defects found while investigating qlc3d (see `configuration-system.md` for the full description of the
intended/actual behavior). Each entry is a verified fact from reading the code, with
file:line citations. 

## 1. `eps_dielectric` setting is declared but never read from the settings file

`SFK_EPS_DIELECTRIC` ("eps_dielectric") is declared in `settings_file_keys.h:54`, and
`Electrodes` has a working `getDielectricPermittivity()`/`setDielectricPermittivities()`
API plus a hard-coded single-element default of `{1.0}` set in the default constructor
(`electrodes.cpp:71-73`). However, `SettingsReader::readElectrodes`
(`settings-reader.cpp:177-209`) never reads any `eps_dielectric` (or numbered
`Dielectric{i}.eps` style) key from the settings file, and never calls
`setDielectricPermittivities()`. There is currently no way for a user to configure a
non-default dielectric permittivity through the settings file, even though the data
model and constant for the key both exist.

## 2. Third `FIXLC*.Easy` value ("rotation") is accepted but always discarded

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

## 3. `Electrode` constructor can double-count the initial time-0 sample

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

## 4. Inconsistent coordinate normalization between `Box` and `Surface` expressions

`Box::getTiltAt`/`getTwistAt` (`box.cpp:164-186`) evaluate `Box{i}.Tilt`/`Box{i}.Twist`
string expressions using coordinates normalized to `[0, 1]` within that box's own
bounding box. `Surface::getEasyTiltAngleAt`/`getEasyTwistAngleAt`
(`alignment.cpp:145-159`) evaluate `FIXLC{i}.Easy` string expressions using absolute
mesh coordinates (in meters), with no normalization. This is not a crash-causing bug,
but it is an undocumented inconsistency: an expression written for a `Box` cannot be
reused verbatim for a `Surface`, and vice versa, without accounting for the different
coordinate conventions.

## 5. `AnchoringType::Polymerise` is defined but unreachable

`AnchoringType::Polymerise` is declared in the enum (`alignment.h:22`), and a `TODO`
comment for a future `Surface::ofPolymerise()` factory exists (`alignment.h:108`), but:
- `SettingsReader::readAlignment` has no `else if (type == "polymerise")` branch
  (`settings-reader.cpp:281-358`), so this value can never be produced from a settings
  file — any settings file specifying `FIXLC{i}.Anchoring = Polymerise` hits the final
  `else` branch and throws a `ReaderError` ("Invalid anchoring type").
- `ManualNodes` is likewise defined in the enum but also has no reachable path in
  `SettingsReader::readAlignment` today.


