#ifndef PROJECT_QLC3D_CONFIGURATION_H
#define PROJECT_QLC3D_CONFIGURATION_H

#include <string>
#include <memory>
#include <filesystem>
class Simu;
class SolverSettings;
class InitialVolumeOrientation;

class Configuration {
    std::filesystem::path settingsFilePath_;
    // The directory of the currently running project
    std::filesystem::path currentDirectory_;

    std::shared_ptr<Simu> simu_;
    std::shared_ptr<LC> lc_;
    std::shared_ptr<MeshRefinement> refinement_;
    std::shared_ptr<Electrodes> electrodes_;
    std::shared_ptr<SolverSettings> solverSettings_;
    std::shared_ptr<Alignment> alignment_;
    std::shared_ptr<InitialVolumeOrientation> initialVolumeOrientation_;

public:
    Configuration();

    void readSettings();
    [[nodiscard]] std::shared_ptr<Simu> getSimu() const;
    [[nodiscard]] std::shared_ptr<Electrodes> getElectrodes() const;

    void simu(Simu *simu) { simu_ = std::shared_ptr<Simu>(simu); }
    void lc(LC *lc) { lc_ = std::shared_ptr<LC>(lc); }
    void solverSettings(SolverSettings *solverSettings);

    [[nodiscard]] std::shared_ptr<LC> getLC() const;
    [[nodiscard]] std::shared_ptr<MeshRefinement> refinement() const;
    [[nodiscard]] std::shared_ptr<SolverSettings> getSolverSettings() const;
    [[nodiscard]] std::shared_ptr<Alignment> getAlignment() const;
    [[nodiscard]] std::shared_ptr<InitialVolumeOrientation> getInitialVolumeOrientation() const;
    /// The file name of the current settings file
    [[nodiscard]] const std::filesystem::path &settingsFile() const { return settingsFilePath_; };
    void settingsFileName(const std::filesystem::path &path) { settingsFilePath_ = path; }

     /// The working directory of the currently running project
    [[nodiscard]] const std::filesystem::path &currentDirectory() const { return currentDirectory_; }
    void currentDirectory(const std::string &path) { currentDirectory_ = path; }
};
#endif //PROJECT_QLC3D_CONFIGURATION_H
