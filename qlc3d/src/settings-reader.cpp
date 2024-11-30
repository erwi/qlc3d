#include <settings-reader.h>
#include <settings_file_keys.h>
#include <meshrefinement.h>
#include <reader.h> // reads our custom key/value based settings file format
#include <utility>
#include <cassert>
#include <set>
#include <electrodes.h>
#include <geom/vec3.h>
#include "solver-settings.h"


SettingsReader::SettingsReader(const std::filesystem::path &fileName):
fileName_(fileName),
simu_(nullptr), lc_(nullptr), electrodes_(nullptr), meshRefinement_(nullptr), solverSettings_(nullptr), alignment_(nullptr) {
    read();
}

void SettingsReader::read() {
    // check that settings file exists
    std::ifstream f(fileName_);
    assertTrue(f.good(), "Settings file does not exist: " + fileName_.string());

    Reader reader;
    reader.setCaseSensitivity(false);
    reader.setLowerCaseStringValues(true); // all returned string values are lower case.

    // substitute environment variables into the settings file when encountering special formatting string ${ENV_VAR}
    reader.setEnvironmentVariableSubstitution(true);
    reader.readSettingsFile(fileName_);

    readSimu(reader);
    readLC(reader);

    readAlignment(reader);
    readRefinement(reader);
    readElectrodes( reader);
    readSolverSettings(reader);

    initialVolumeOrientation_ = std::make_unique<InitialVolumeOrientation>();
    readInitialVolumeOrientation(reader);
}

std::unique_ptr<Simu> SettingsReader::simu() {
    assert(simu_.get() != nullptr); // read not called
    return std::move(simu_);
}

std::unique_ptr<LC> SettingsReader::lc() {
    assert(lc_ != nullptr);
    return std::move(lc_);
}

std::unique_ptr<MeshRefinement> SettingsReader::refinement() {
    assert(meshRefinement_ != nullptr);
    return std::move(meshRefinement_);
}

std::unique_ptr<Electrodes> SettingsReader::electrodes() {
    assert(electrodes_ != nullptr);
    return std::move(electrodes_);
}

std::unique_ptr<SolverSettings> SettingsReader::solverSettings() {
    assert(solverSettings_ != nullptr);
    return std::move(solverSettings_);
}

std::unique_ptr<Alignment> SettingsReader::alignment() {
    assert(alignment_ != nullptr);
    return std::move(alignment_);
}

std::unique_ptr<InitialVolumeOrientation> SettingsReader::initialVolumeOrientation() {
    assert(initialVolumeOrientation_ != nullptr);
    return std::move(initialVolumeOrientation_);
}

// <editor-fold desc="Private Methods">
void SettingsReader::assertTrue(bool condition, const std::string &errorMsg) {
    if (!condition) {
        throw ReaderError(errorMsg, fileName_.string());
    }
}

/** reader reads optional value by given key. if a value is present, the consumer accepts it. */
void readDouble(Reader &r, const std::string& key, std::function<void(double)> consumer) {
    if (auto v = r.getOptional<double>(key)) {
        consumer(v.value());
    }
}

void readInt(Reader &r, const std::string& key, std::function<void(int)> consumer) {
    if (auto v = r.getOptional<int>(key)) {
        consumer(v.value());
    }
}

void SettingsReader::readSimu(Reader &reader) {
    using namespace std;
    SimuBuilder builder;
    // meshName is mandatory, all other settings in Simu are optional as sensible defaults exist.
    builder.meshFileName(reader.getValueByKey<string>(SFK_MESH_NAME));

    readDouble(reader, SFK_DT, [&](double v) { builder.initialTimeStep( v); });
    //if (auto t = reader.getOptional<double>(SFK_DT)) { builder.initialTimeStep(t.value()); }
    if (auto s = reader.getOptional<string>(SFK_Q_MATRIX_SOLVER)) { builder.qMatrixSolver(s.value()); }
    if (auto v = reader.getOptional<double>(SFK_TARGET_DQ)) { builder.targetDQ(v.value()); }
    if (auto v = reader.getOptional<double>(SFK_MAX_ERROR)) { builder.maxError(v.value()); }
    if (auto v = reader.getOptional<string>(SFK_END_CRITERION)) { builder.endCriterion(v.value()); }
    if (auto v = reader.getOptional<string>(SFK_LOAD_Q)) { builder.loadQ(v.value()); }
    if (auto v = reader.getOptional<string>(SFK_SAVE_DIR)) { builder.saveDir(v.value()); }
    if (auto v = reader.getOptional<double>(SFK_END_VALUE)) { builder.endValue(v.value()); }
    if (auto v = reader.getOptional<int>(SFK_OUTPUT_FORMAT)) { builder.outputFormat(v.value()); }
    if (auto v = reader.getOptional<int>(SFK_OUTPUT_ENERGY)) { builder.outputEnergy(v.value()); }
    //if (auto v = reader.getOptional<int>(SFK_SAVE_FORMAT)) { builder.saveFormat(v.value()); } // TODO list of enum
    if (auto v = reader.getOptional<int>(SFK_SAVE_ITER)) { builder.saveIter(v.value()); }
    if (auto v = reader.getOptional<double>(SFK_SAVE_TIME)) { builder.saveTime(v.value()); }
    if (auto v = reader.getOptional<int>(SFK_NUM_ASSEMBLY_THREADS)) { builder.numAssemblyThreads(v.value()); }
    if (auto v = reader.getOptional<int>(SFK_NUM_MATRIX_SOLVER_THREADS)) { builder.numMatrixSolverThreads(v.value()); }

    if (auto v = reader.getOptional<vector<double>>(SFK_DT_LIMITS)) {
        assertTrue(v.value().size() == 2, SFK_DT_LIMITS + " length should be 2");
        builder.dtLimits(v.value()[0], v.value()[1]);
    }

    if (auto v = reader.getOptional<vector<double>>(SFK_DT_FUNCTION)) {
        auto dtf = v.value();
        assertTrue(dtf.size() == 4, SFK_DT_FUNCTION + " length should be 4");
        builder.dtFunction(dtf[0], dtf[1], dtf[2], dtf[3]);
    }

    if (auto v = reader.getOptional<vector<double>>(SFK_STRETCH_VECTOR)) {
        auto sv = v.value();
        assertTrue(sv.size() == 3, SFK_STRETCH_VECTOR + " length should be 3");
        builder.stretchVector(sv[0], sv[1], sv[2]);
    }

    if (auto v = reader.getOptional<vector<double>>(SFK_REGULAR_GRID_SIZE)) {
        auto rgs = v.value();
        assertTrue(rgs.size() == 3, SFK_REGULAR_GRID_SIZE + " length should be 3");
        builder.regularGridSize(rgs[0], rgs[1], rgs[2]);
    }

    if (auto v = reader.getOptional<vector<string>>(SFK_SAVE_FORMAT)) {
        std::set<string> sf(v.value().begin(), v.value().end());
        builder.saveFormat(sf);
    }

    simu_.reset(builder.build());
}

void SettingsReader::readLC(Reader &reader) {
    LCBuilder builder;
    readDouble(reader, SFK_K11, [&](double k11) { builder.K11(k11); });
    readDouble(reader, SFK_K22, [&](double k22) { builder.K22(k22); });
    readDouble(reader, SFK_K24, [&](double k24) { builder.K24(k24); });
    readDouble(reader, SFK_K33, [&](double k33) { builder.K33(k33); });
    readDouble(reader, SFK_p0, [&](double p0) { builder.p0(p0); });
    readDouble(reader, SFK_A, [&](double a) { builder.A(a); });
    readDouble(reader, SFK_B, [&](double b) { builder.B(b); });
    readDouble(reader, SFK_C, [&](double c) { builder.C(c); });
    readDouble( reader, SFK_EPS_PAR, [&](double e) { builder.eps_par(e); });
    readDouble( reader, SFK_EPS_PER, [&](double e) { builder.eps_per(e); });
    readDouble( reader, SFK_E1, [&](double e) { builder.e1(e); });
    readDouble(reader, SFK_E3, [&](double e) { builder.e3(e); });
    readDouble( reader, SFK_GAMMA1, [&](double g) { builder.gamma1(g); });

    lc_.reset(builder.build());
}

void SettingsReader::readElectrodes(Reader &reader) {

  std::vector<double> eField = reader.getOptional<std::vector<double>>("efield").value_or(vector<double>());
  assertTrue(eField.size() == 0 || eField.size() == 3, "Unexpected number of components in EField. Expected 3.");
  std::shared_ptr<Vec3> eFieldVec = nullptr;
  if (eField.size() == 3) {
    eFieldVec = std::make_shared<Vec3>(eField[0], eField[1], eField[2]);
  }

  std::vector<std::shared_ptr<Electrode>> electrodesVector;
  for (unsigned int i = 1; i <= 99; i++) {
    std::string switchingTimesKey = "e" + std::to_string(i) + ".time";
    std::string switchingValuesKey = "e" + std::to_string(i) + ".pot";

    auto maybeTimes = reader.getOptional<vector<double>>(switchingTimesKey);
    auto maybePots = reader.getOptional<vector<double>>(switchingValuesKey);

    auto times = maybeTimes.value_or(vector<double>());
    auto pots = maybePots.value_or(vector<double>());
    assertTrue(times.size() == pots.size(), "Electrode " + std::to_string(i) + " has different number of switching times and potentials");

    if (!times.empty()) {
      std::shared_ptr<Electrode> electrodePtr = std::make_shared<Electrode>(i, times, pots);
      electrodesVector.push_back(electrodePtr);
    }
  }

  if (eFieldVec != nullptr) {
    electrodes_ = std::make_unique<Electrodes>(Electrodes::withConstantElectricField(*eFieldVec));
  } else {
    electrodes_ = std::make_unique<Electrodes>(Electrodes::withElectrodePotentials(electrodesVector));
  }
}

void SettingsReader::readRefinement(Reader &reader) {
    using namespace std;

    meshRefinement_ = std::make_unique<MeshRefinement>();

    // read periodically repeating refinement times/iterations
    meshRefinement_->setRepRefIter(reader.getOptional<unsigned int>("RepRefIter").value_or(0));
    meshRefinement_->setRepRefTime(reader.getOptional<double>("RepRefTime").value_or(0));

    // read all refinement objects
    vector<RefinementConfig> refinement;
    const vector<double> emptyDoubles;
    const vector<int> emptyInts;

    if (reader.containsKeyWithPrefix("REFINEMENT")) {
        for (int i = 1; i <= 99; i++) {
            string keyBase = "REFINEMENT" + to_string(i);

            if (!reader.containsKeyWithPrefix(keyBase)) {
                break;
            }

            auto type = reader.getValueByKey<string>(keyBase + ".Type"); // required, the other may be optional, depending on type
            auto x = reader.getOptional<vector<double>>(keyBase + ".X").value_or(emptyDoubles);
            auto y = reader.getOptional<vector<double>>(keyBase + ".Y").value_or(emptyDoubles);
            auto z = reader.getOptional<vector<double>>(keyBase + ".Z").value_or(emptyDoubles);

            auto iterations = reader.getOptional<vector<int>>(keyBase + ".Iterations").value_or(emptyInts);
            auto times = reader.getOptional<vector<double>>(keyBase + ".Times").value_or(emptyDoubles);
            auto values = reader.getOptional<vector<double>>(keyBase + ".Values").value_or(emptyDoubles);

            refinement.emplace_back(type, iterations, times, values, x, y, z);
        }
    }
    meshRefinement_->setRefinementConfig(std::move(refinement));
}

void SettingsReader::readSolverSettings(Reader &reader) {
  // read all solver settings
  solverSettings_ = std::make_unique<SolverSettings>();
  readInt(reader, SFK_NUM_ASSEMBLY_THREADS, [&](int v) { solverSettings_->setnThreads(v); });

  readInt(reader, SFK_Q_SOLVER, [&](int v) { solverSettings_->setQ_Solver(v); });
  readInt(reader, SFK_Q_NEWTON_PANIC_ITER, [&](int v) { solverSettings_->setQ_Newton_Panic_Iter(v); });
  readDouble(reader, SFK_Q_NEWTON_PANIC_COEFF, [&](double v) { solverSettings_->setQ_Newton_Panic_Coeff(v); });
  readInt(reader, SFK_Q_PCG_PRECONDITIONER, [&](int v) { solverSettings_->setQ_PCG_Preconditioner(v); });
  readInt(reader, SFK_Q_PCG_MAXITER, [&](int v) { solverSettings_->setQ_PCG_Maxiter(v); });
  readDouble(reader, SFK_Q_PCG_TOLER, [&](double v) { solverSettings_->setQ_PCG_Toler(v); });
  readInt(reader, SFK_Q_GMRES_PRECONDITIONER, [&](int v) { solverSettings_->setQ_GMRES_Preconditioner(v); });
  readInt(reader, SFK_Q_GMRES_MAXITER, [&](int v) { solverSettings_->setQ_GMRES_Maxiter(v); });
  readInt(reader, SFK_Q_GMRES_RESTART, [&](int v) { solverSettings_->setQ_GMRES_Restart(v); });
  readDouble(reader, SFK_Q_GMRES_TOLER, [&](double v) { solverSettings_->setQ_GMRES_Toler(v); });

  readInt(reader, SFK_V_SOLVER, [&](int v) { solverSettings_->setV_Solver(v); });
  readInt(reader, SFK_V_PCG_PRECONDITIONER, [&](int v) { solverSettings_->setV_PCG_Preconditioner(v); });
  readInt(reader, SFK_V_PCG_MAXITER, [&](int v) { solverSettings_->setV_PCG_Maxiter(v); });
  readDouble(reader, SFK_V_PCG_TOLER, [&](double v) { solverSettings_->setV_PCG_Toler(v); });

  readInt(reader, SFK_V_GMRES_PRECONDITIONER, [&](int v) { solverSettings_->setV_GMRES_Preconditioner(v); });
  readInt(reader, SFK_V_GMRES_MAXITER, [&](int v) { solverSettings_->setV_GMRES_Maxiter(v); });
  readInt(reader, SFK_V_GMRES_RESTART, [&](int v) { solverSettings_->setV_GMRES_Restart(v); });
  readDouble(reader, SFK_V_GMRES_TOLER, [&](double v) { solverSettings_->setV_GMRES_Toler(v); });
}

void SettingsReader::readAlignment(Reader &reader) {
    alignment_ = std::make_unique<Alignment>();

    if (reader.containsKeyWithPrefix("FIXLC")) {
        for (int i = 1; i <= 99; i++) {
          string keyBase = "FIXLC" + to_string(i);
          if (!reader.containsKeyWithPrefix(keyBase)) {
            continue;
          }

          auto type = reader.getValueByKey<string>(keyBase + ".Anchoring");

          auto overrideVolume = reader.getOptional<bool>(keyBase + ".overrideVolume");


          if (type == "strong") {
            std::string key = keyBase + ".Easy";
            if (reader.isValueArrayOfNumbers(key)) {
              // if all easy angles values are valid numbers, then it is an angle
              auto easyAngles = reader.getValueByKey<std::vector<double>>(key);
              double tiltAngleDegrees = easyAngles[0];
              double twistAngleDegrees = easyAngles[1];

              alignment_->addSurface( overrideVolume.has_value() ?
                                      Surface::ofStrongAnchoring(i, tiltAngleDegrees, twistAngleDegrees, overrideVolume.value()) :
                                      Surface::ofStrongAnchoring(i, tiltAngleDegrees, twistAngleDegrees));

            } else if (reader.isValueArrayOfStrings(key)) {
              // if all easy angles values are strings, then it is an expression
              auto easyAngles = reader.getValueByKey<std::vector<std::string>>(key);
              auto tiltExpression = easyAngles[0];
              auto twistExpression = easyAngles[1];

              alignment_->addSurface( overrideVolume.has_value() ?
                                      Surface::ofStrongAnchoring(i, tiltExpression, twistExpression, overrideVolume.value()) :
                                      Surface::ofStrongAnchoring(i, tiltExpression, twistExpression));
            } else {
              throw ReaderError("Failed to read easy angles for FIXLC" + std::to_string(i), fileName_.string());
            }
          } else if (type == "homeotropic") {
            // Ignore overrideVolume, as it doesnt' make sense in the strong homeotropic case
            alignment_->addSurface(Surface::ofStrongHomeotropic(i));
          } else if (type == "weak") {
            auto easyAngles = reader.getValueByKey<std::vector<double>>(keyBase + ".Easy");
            assertTrue(easyAngles.size() == 2 || easyAngles.size() == 3, keyBase + " easy angles should have 2 or 3 components, got " + std::to_string(easyAngles.size()));
            auto strength = reader.getValueByKey<double>(keyBase + ".Strength");
            auto k1 = reader.getValueByKey<double>(keyBase + ".K1");
            auto k2 = reader.getValueByKey<double>(keyBase + ".K2");

            alignment_->addSurface(overrideVolume.has_value() ?
                                  Surface::ofWeakAnchoring(i, easyAngles[0], easyAngles[1], strength, k1, k2, overrideVolume.value()) :
                                  Surface::ofWeakAnchoring(i, easyAngles[0], easyAngles[1], strength, k1, k2));

          } else if (type == "degenerate") {
            auto strength = reader.getValueByKey<double>(keyBase + ".Strength");
            if (strength >= 0) {
              alignment_->addSurface(overrideVolume.has_value() ?
                                    Surface::ofPlanarDegenerate(i, strength, overrideVolume.value()) :
                                    Surface::ofPlanarDegenerate(i, strength));
            } else {
              alignment()->addSurface(overrideVolume.has_value() ?
                                    Surface::ofWeakHomeotropic(i, -strength, overrideVolume.value()) :
                                    Surface::ofWeakHomeotropic(i, -strength));
            }
          } else if (type == "weakhomeotropic") {
            auto strength = reader.getValueByKey<double>(keyBase + ".Strength");
            alignment_->addSurface(overrideVolume.has_value() ?
                                  Surface::ofWeakHomeotropic(i, strength, overrideVolume.value()) :
                                  Surface::ofWeakHomeotropic(i, strength));
          } else if (type == "freeze") {
            alignment_->addSurface(Surface::ofFreeze(i));
          }
          else {
            throw ReaderError("Invalid anchoring type: " + type, fileName_.string() + ". This may be a typo in the settings file or it has not yet been implemented");
          }
        }
    }
}

void SettingsReader::readInitialVolumeOrientation(Reader &reader) {
  using std::string;
  const int MAX_NUM_BOXES = 100;
  // Loop over all possible boxes ad try to read
  for (int boxNum = 1; boxNum < MAX_NUM_BOXES; boxNum++) {
    string typeKey = wildcardToNum(SFK_BOX_TYPE, boxNum);
    if (reader.containsKey(typeKey)) {
      // found box with current number
      // Make keys for current box number
      string paramsKey = wildcardToNum(SFK_BOX_PARAMS, boxNum);
      string xKey = wildcardToNum(SFK_BOX_X, boxNum);
      string yKey = wildcardToNum(SFK_BOX_Y, boxNum);
      string zKey = wildcardToNum(SFK_BOX_Z, boxNum);
      string tiltKey = wildcardToNum(SFK_BOX_TILT, boxNum);
      string twistKey = wildcardToNum(SFK_BOX_TWIST, boxNum);
      // add box
      /*
      initialVolumeOrientation_->addBox(boxNum, // TODO we need a box builder to better support expressions for tilt and twist
                   reader.get<string>(typeKey, Box::DEFAULT_TYPE),
                   reader.get<vector<double>>(paramsKey, Box::DEFAULT_PARAMS),
                   reader.get<vector<double>>(xKey, Box::DEFAULT_X_Y_Z),
                   reader.get<vector<double>>(yKey, Box::DEFAULT_X_Y_Z),
                   reader.get<vector<double>>(zKey, Box::DEFAULT_X_Y_Z),
                   reader.get<vector<double>>(tiltKey, Box::DEFAULT_TILT_TWIST),
                   reader.get<vector<double>>(twistKey, Box::DEFAULT_TILT_TWIST));
      */
      auto bbuilder = BoxBuilder(boxNum);
      bbuilder.setBoxType(reader.get<string>(typeKey, Box::DEFAULT_TYPE))
              .setParams(reader.get<vector<double>>(paramsKey, Box::DEFAULT_PARAMS))
              .setX(reader.get<vector<double>>(xKey, Box::DEFAULT_X_Y_Z))
              .setY(reader.get<vector<double>>(yKey, Box::DEFAULT_X_Y_Z))
              .setZ(reader.get<vector<double>>(zKey, Box::DEFAULT_X_Y_Z));

      if (reader.isValueArray(tiltKey)) {
        bbuilder.setTilt(reader.get<vector<double>>(tiltKey, Box::DEFAULT_TILT_TWIST));
      } else { // assume it is an analytic expression
        auto tiltExpression = reader.get<string>(tiltKey, Box::DEFAULT_TILT_ANGLE_DEGREES);
        bbuilder.setTiltExpression(tiltExpression);
      }

      if (reader.isValueArray(twistKey)) {
        bbuilder.setTwist(reader.get<vector<double>>(twistKey, Box::DEFAULT_TILT_TWIST));
      } else { // assume it is an analytic expression
        auto twistExpression = reader.get<string>(twistKey, Box::DEFAULT_TWIST_ANGLE_DEGREES);
        bbuilder.setTwistExpression(twistExpression);
      }
      try {
        initialVolumeOrientation_->addBox(bbuilder.build());
      } catch (ExpressionException &e) {
        throw ReaderError("Error in expression for BOX" + std::to_string(boxNum) + ", " + e.what(), fileName_.string());
      } catch (...) {
        throw std::runtime_error("Error creating BOX" + std::to_string(boxNum));
      }
    }
  } // end for boxNum
}

// </editor-fold>