#include <qlc3d.h>
#include <configuration.h>
#include <simulation-container.h>
#include <util/logging.h>
#include <util/exception.h>
#include <filesystem>
#include <io/result-output.h>

namespace fs = std::filesystem;

void parseArgs(int argc, char **args, Configuration &configuration) {
    // if a working directory has been provided set to that first
    if (argc >= 3) {
        fs::path workingDir = args[2];
        if (!workingDir.is_absolute()) { RUNTIME_ERROR(fmt::format("Working directory {} must be absolute", workingDir)); }
        if (!fs::exists(workingDir) || !fs::is_directory(workingDir)) { RUNTIME_ERROR(fmt::format("Working directory {} should be an existing directory", workingDir)); }
        configuration.currentDirectory(workingDir.string());
        Log::info("Setting working directory to {}", workingDir.string());
        fs::current_path(workingDir);
    }

    // the config file name has been provided. It may be either absolute or relative. If relative convert it to absolute assuming it's in current working dir
    if (argc >= 2) {
      std::filesystem::path settingsFilePath = args[1];

      if (settingsFilePath.is_relative()) {
          settingsFilePath = fs::current_path() / settingsFilePath;
      }

      if (!std::filesystem::exists(settingsFilePath)) {
        RUNTIME_ERROR(fmt::format("The settings file {} does not exist", settingsFilePath));
      } else {
        Log::info("Using settings file {}", settingsFilePath.string());
      }
      configuration.settingsFileName(settingsFilePath);
    }
}

int runSimulation(Configuration &configuration) {
    try {
        configuration.readSettings();

        auto simu = configuration.simu();
        ResultOutput resultOutput(simu->getSaveFormat(), simu->meshName(), configuration.lc()->S0(), simu->getSaveDirAbsolutePath());

        SimulationContainer simulation(configuration, resultOutput);
        Log::clearIndent();
        Log::info("Initialising.");
        Log::incrementIndent();
        simulation.initialise();
        Log::decrementIndent();

        Log::info("Starting simulation.");
        while (simulation.hasIteration()) {
            simulation.runIteration();
            // TODO: auto state = simulation.getCurrentState() and do something with it?
        }
        simulation.postSimulationTasks();
    } catch(std::exception &e) {
        Log::error("An exception has occurred: {}", e.what());
        return 1;
    } catch(...) {
        Log::error("An error has occurred");
        return 3;
    }
    Log::info("simulation has ended without errors");
    return 0;
}

void printInfo() {
    Qlc3dInfo info;
    Log::info("qlc3d. Build date={}, build time={}, build type={}.",
              info.buildDate,
              info.buildTime,
              info.isDebug ? "DEBUG" : "RELEASE");
    // print current directory
    Log::info("Current directory: {}", std::filesystem::current_path().string());
}

int main(int argc, char **args) {
    printInfo();

    Configuration configuration;
    parseArgs(argc, args, configuration);

    return runSimulation(configuration);
}