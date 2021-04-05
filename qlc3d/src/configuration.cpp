//
// Created by eero on 02/04/2021.
//
#include <filesystem>
#include <settings-reader.h>
#include "configuration.h"

Configuration::Configuration() :
        settingsFileName_("./meshes/test.txt"), // default for backwards compatibility
        currentDirectory_(std::filesystem::current_path().c_str()),
        simu_(nullptr)
        {}

void Configuration::readSettings() {
    SettingsReader reader(settingsFileName());
    simu_ = reader.simu();
}

std::shared_ptr<Simu> Configuration::simu() const {
    return simu_;
}