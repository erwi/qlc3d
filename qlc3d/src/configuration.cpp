//
// Created by eero on 02/04/2021.
//
#include <filesystem>
#include "configuration.h"

Configuration::Configuration() :
        settingsFileName_("./meshes/test.txt"), // default for backwards compatibility
        currentDirectory_(std::filesystem::current_path().c_str()) {}