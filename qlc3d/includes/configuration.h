//
// Created by eero on 02/04/2021.
//

#ifndef PROJECT_QLC3D_CONFIGURATION_H
#define PROJECT_QLC3D_CONFIGURATION_H

#include <string>

class Configuration {
    std::string settingsFileName_;
    // The directory of the currently running project
    std::string currentDirectory_;

public:
    Configuration();

    const std::string &settingsFileName() const { return settingsFileName_; };
    void settingsFileName(const std::string &name) { settingsFileName_ = name; }

    /*!
     * The working directory of the currently running project
     */
    const std::string &currentDirectory() const { return currentDirectory_; }
    void currentDirectory(const std::string &path) { currentDirectory_ = path; }
};


#endif //PROJECT_QLC3D_CONFIGURATION_H
