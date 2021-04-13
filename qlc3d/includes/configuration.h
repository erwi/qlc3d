//
// Created by eero on 02/04/2021.
//

#ifndef PROJECT_QLC3D_CONFIGURATION_H
#define PROJECT_QLC3D_CONFIGURATION_H

#include <string>
#include <memory>
class Simu;

class Configuration {
    std::string settingsFileName_;
    // The directory of the currently running project
    std::string currentDirectory_;

    std::shared_ptr<Simu> simu_;

public:
    Configuration();

    void readSettings();
    [[nodiscard]] std::shared_ptr<Simu> simu() const;
    void simu(Simu *simu) { simu_ = std::shared_ptr<Simu>(simu); }

    /// The file name of the current settings file
    [[nodiscard]] const std::string &settingsFileName() const { return settingsFileName_; };
    void settingsFileName(const std::string &name) { settingsFileName_ = name; }

     /// The working directory of the currently running project
    [[nodiscard]] const std::string &currentDirectory() const { return currentDirectory_; }
    void currentDirectory(const std::string &path) { currentDirectory_ = path; }
};
#endif //PROJECT_QLC3D_CONFIGURATION_H
