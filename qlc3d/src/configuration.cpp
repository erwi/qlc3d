#include <filesystem>
#include <settings-reader.h>
#include <cassert>
#include <memory>
#include <configuration.h>
#include <solver-settings.h>

Configuration::Configuration() :
        settingsFilePath_("./meshes/test.txt"), // default for backwards compatibility
        currentDirectory_(std::filesystem::current_path().c_str()),
        simu_(nullptr),
        solverSettings_(nullptr)
        {}

void Configuration::readSettings() {
    SettingsReader reader(settingsFile());
    simu_ = reader.simu();
    lc_ = reader.lc();
    refinement_ = reader.refinement();
    electrodes_ = reader.electrodes();
    solverSettings_ = reader.solverSettings();
}

void Configuration::solverSettings(SolverSettings *solverSettings) {
    solverSettings_ = std::shared_ptr<SolverSettings>(solverSettings);
}

std::shared_ptr<Simu> Configuration::getSimu() const {
    assert(simu_ != nullptr);
    return simu_;
}

std::shared_ptr<Electrodes> Configuration::getElectrodes() const {
    assert(electrodes_ != nullptr);
    return electrodes_;
}

std::shared_ptr<LC> Configuration::getLC() const {
    assert(lc_ != nullptr);
    return lc_;
}

std::shared_ptr<MeshRefinement> Configuration::refinement() const {
    assert(refinement_ != nullptr);
    return refinement_;
}

std::shared_ptr<SolverSettings> Configuration::getSolverSettings() const {
    assert(solverSettings_ != nullptr);
    return solverSettings_;
}