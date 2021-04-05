//
// Created by eero on 04/04/2021.
//
#include <settings-reader.h>
#include <settings_file_keys.h>
#include <reader.h> // reads our custom key/value based settings file format
#include <utility>
#include <cassert>
#include <set>

SettingsReader::SettingsReader(std::string fileName):
fileName_(std::move(fileName)),
simu_(nullptr) {
    read();
}

void SettingsReader::read() {
    try {
        Reader reader;
        reader.setCaseSensitivity(false);
        reader.setLowerCaseStringValues(true); // all returned string values are lower case.
        reader.readSettingsFile(fileName_);

        readSimu(reader);
        //readSimu(simu, eventList,reader);
        //readLC(lc, reader);
        //readBoxes(boxes, reader);
        //readSolverSettings(settings, reader);
        //readAlignment(alignment, reader);
        //readElectrodes( electrodes, eventList, reader);
        //readRefinement(reader, eventList);
    } catch (ReaderError &e) {
        e.printError();
        throw e;
    }
}

std::unique_ptr<Simu> SettingsReader::simu() {
    assert(simu_.get() != nullptr); // read not called
    return std::move(simu_);
}

// <editor-fold desc="Private Methods">
void SettingsReader::assertTrue(bool condition, const std::string &errorMsg) {
    if (!condition) {
        throw ReaderError(errorMsg, fileName_);
    }
}
void SettingsReader::readSimu(Reader &reader) {
    using namespace std;
    SimuBuilder builder;
    // meshName is mandatory, all other settings in Simu are optional as sensible defaults exist.
    builder.meshFileName(reader.getValueByKey<string>(SFK_MESH_NAME));

    if (auto t = reader.getOptional<double>(SFK_DT)) { builder.initialTimeStep(t.value()); }
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
// </editor-fold>