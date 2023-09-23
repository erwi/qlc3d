#include <simu.h>
#include <algorithm>
#include <cassert>
#include <omp.h>
#include "stringenum.h"
#include "util/stringutil.h"
#include <reader.h>
#include <globals.h>
#include <settings_file_keys.h>
#include <geom/vec3.h>

// Define valid enum string keys
const vector<string> Simu::VALID_END_CRITERIA = {"iterations", "time", "change"};

// NOTE: the order of these should match the order of the corresponding enum values
const vector<string> Simu::VALID_SAVE_FORMATS = {
        "lcview", "regularvtk", "regularvecmat", "dirstackz", "lcviewtxt", "csvunstructured", "vtkunstructuredasciigrid"};
const vector<string> Simu::VALID_Q_MATRIX_SOLVERS = {"Auto", "PCG", "GMRES"};
// Define values of default Simu parameters
const string Simu::DEFAULT_LOAD_Q = "";
const string Simu::DEFAULT_SAVE_DIR = "res";
const string Simu::DEFAULT_Q_MATRIX_SOLVER = "auto";
const double Simu::DEFAULT_END_VALUE = 1e-3;
const double Simu::DEFAULT_DT = 1e-9;
const double Simu::DEFAULT_TARGET_DQ = 1e-3;
const double Simu::DEFAULT_MAX_DT = 1e-4;
const double Simu::DEFAULT_MAX_ERROR = 1e-3;
const vector<double> Simu::DEFAULT_STRETCH_VECTOR = {1, 1, 1};
const vector<double> Simu::DEFAULT_DT_LIMITS = {Simu::DEFAULT_DT, Simu::DEFAULT_MAX_DT};
const vector<double> Simu::DEFAULT_DT_FUNCTION = {0.5, 0.8, 1.2, 10};
const vector<idx> Simu::DEFAULT_REGULAR_GRID_SIZE = {0, 0, 0};
const int Simu::DEFAULT_OUTPUT_ENERGY = 0;
const int Simu::DEFAULT_OUTPUT_FORMAT = SIMU_OUTPUT_FORMAT_BINARY; // TODO: Get rid of this
const int Simu::DEFAULT_SAVE_ITER = 0;
const double Simu::DEFAULT_SAVE_TIME = 0;
const int Simu::DEFAULT_NUM_ASSEMBLY_THREADS = 0;
const int Simu::DEFAULT_NUM_MATRIX_SOLVER_THREADS = 0;
const Simu::EndCriteria Simu::DEFAULT_END_CRITERION = Simu::EndCriteria::Time;

// TODO: use array of length 4 instead of pointer?
void Simu::getdtFunction(double *f) {
    f[0] = dtFunction_[0];
    f[1] = dtFunction_[1];
    f[2] = dtFunction_[2];
    f[3] = dtFunction_[3];
}

const std::vector<std::string> Simu::getSaveFormatStrings() const {
  std::vector<std::string> saveFormatStrings;
  for (auto s : getSaveFormat()) {
    saveFormatStrings.push_back(VALID_SAVE_FORMATS[s]);
  }
  return saveFormatStrings;
}

Vec3 Simu::getStretchVector() const {
    return {stretchVector_[0], stretchVector_[1], stretchVector_[2]};
}

// <editor-fold desc=SimuBuilder>
void assertTrue(bool shouldBeTrue, const std::string &msg) {
    if (!shouldBeTrue) {
        throw runtime_error(msg);
    }
}

Simu *SimuBuilder::build() const {
    return new Simu(meshFileName_, initialTimeStep_,
                    qMatrixSolver_, maxError_,
                    targetDQ_, dtLimits_,
                    dtFunction_, endCriterion_,
                    loadQ_, saveDir_,
                    endValue_, stretchVector_,
                    regularGridSize_, outputEnergy_,
                    outputFormat_, saveIter_, saveTime_,
                    saveFormat_, numAssemblyThreads_,
                    numMatrixSolverThreads_,
                    workingDir_ / saveDir_); // absolute path to result directory
}

SimuBuilder &SimuBuilder::initialTimeStep(double dt) {
    assertTrue(dt >= 0, "time step should be >= 0");
    initialTimeStep_ = dt;
    return *this;
}

SimuBuilder &SimuBuilder::qMatrixSolver(const std::string &solverName) {
    if ("auto" == solverName) {
        qMatrixSolver_ = Simu::QMatrixSolvers::Auto;
    } else if ("pcg" == solverName) {
        qMatrixSolver_ = Simu::QMatrixSolvers::PCG;
    } else if ("gmres" == solverName) {
        qMatrixSolver_ = Simu::QMatrixSolvers::GMRES;
    } else {
        throw runtime_error("valid Q solver names are [auto, pcg, gmres], got " + solverName);
    }
    return *this;
}

SimuBuilder &SimuBuilder::targetDQ(double targetDQ) {
    assertTrue(targetDQ > 0, "targetDQ should be > 0");
    targetDQ_ = targetDQ;
    return *this;
}

SimuBuilder &SimuBuilder::maxError(double maxError) {
    assertTrue(maxError > 0, "maxError should be > 0");
    maxError_ = maxError;
    return *this;
}

SimuBuilder &SimuBuilder::dtLimits(double low, double high) {
    assertTrue(low > 0, "min time step should be > 0");
    assertTrue(high > 0, "max time step should be > min time step");
    dtLimits_[0] = low;
    dtLimits_[1] = high;
    return *this;
}

SimuBuilder &SimuBuilder::dtFunction(double v1, double v2, double v3, double v4) {
    // TODO assertions
    dtFunction_[0] = v1;
    dtFunction_[1] = v2;
    dtFunction_[2] = v3;
    dtFunction_[3] = v4;
    return *this;
}

SimuBuilder &SimuBuilder::endCriterion(const string &name) {
    for (unsigned int i = 0; i < Simu::VALID_END_CRITERIA.size(); i++) {
        if (Simu::VALID_END_CRITERIA[i] == name) {
            endCriterion_ = static_cast<Simu::EndCriteria>(i);
            return *this;
        }
    }

    assertTrue(false, "invalid end criterion " + name + ", valid criteria are " + StringUtil::toString(Simu::VALID_END_CRITERIA));
    return *this; // should neve execute, but appeases CLion code analyser
}

SimuBuilder &SimuBuilder::loadQ(const std::string &loadQ) {
    loadQ_ = loadQ;
    return *this;
}

SimuBuilder &SimuBuilder::saveDir(const std::string &saveDir) {
    saveDir_ = saveDir;
    return *this;
}

SimuBuilder &SimuBuilder::endValue(double endValue) {
    assertTrue(endValue > 0, "endValue should be > 0");
    endValue_ = endValue;
    return *this;
}

SimuBuilder &SimuBuilder::stretchVector(double x, double y, double z) {
    assertTrue(x > 0, "stretch vector x should be > 0");
    assertTrue(y > 0, "stretch vector y should be > 0");
    assertTrue(z > 0, "stretch vector z should be > 0");
    stretchVector_[0] = x;
    stretchVector_[1] = y;
    stretchVector_[2] = z;
    return *this;
}

SimuBuilder &SimuBuilder::regularGridSize(size_t x, size_t y, size_t z) {
    // either all 0 or all > 0. not mixing some 0 some non-zero
    if (x > 0 || y > 0 || z > 0) {
        assertTrue(x > 0 && y > 0 && z > 0,
                   "invalid regular grid size. All values must be either 0 or non-zero");
    }
    regularGridSize_[0] = x;
    regularGridSize_[1] = y;
    regularGridSize_[2] = z;
    return *this;
}

SimuBuilder &SimuBuilder::outputEnergy(int outputEnergy) {
    outputEnergy_ = outputEnergy;
    return *this;
}

SimuBuilder &SimuBuilder::outputFormat(int outputFormat) {
    outputFormat_ = outputFormat;
    return *this;
}

SimuBuilder &SimuBuilder::saveIter(int saveIter) {
    assertTrue(saveIter >= 0, "saveIter should not be negative");
    saveIter_ = saveIter;
    return *this;
}

SimuBuilder &SimuBuilder::saveTime(double saveTime) {
    assertTrue(saveTime >= 0, "saveTime should not be negative");
    saveTime_ = saveTime;
    return *this;
}


SimuBuilder &SimuBuilder::saveFormat(const set<std::string> &saveFormats) {
    int numValidFormats = Simu::VALID_SAVE_FORMATS.size();
    for (const auto & saveFormat : saveFormats) {
        int i = 0;
        for (i = 0; i < numValidFormats; i++) {
            if (Simu::VALID_SAVE_FORMATS[i] == saveFormat) {
                saveFormat_.insert(static_cast<Simu::SaveFormats>(i));
                break;
            }
        }
        assertTrue(i < numValidFormats, "not a valid save format: " + saveFormat);
    }
    return *this;
}

SimuBuilder &SimuBuilder::numAssemblyThreads(unsigned int n) {
    numAssemblyThreads_ = n;
    return *this;
}

SimuBuilder &SimuBuilder::numMatrixSolverThreads(unsigned int n) {
    numMatrixSolverThreads_ = n;
    return *this;
}

// </editor-fold>