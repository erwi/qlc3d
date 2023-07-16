#include <qlc3d.h>
#include <configuration.h>
#include <simulation-container.h>
#include <util/logging.h>
#include <util/exception.h>
#include <filesystem>
#include <io/result-output.h>


void parseArgs(int argc, char **args, Configuration &configuration) {
    // just the config file name has been provided. It is expected to be in the current directory
    if (argc >= 2) {
      std::filesystem::path settingsFilePath = std::filesystem::current_path() / args[1];
      if (!std::filesystem::exists(settingsFilePath)) {
        RUNTIME_ERROR(fmt::format("The settings file {} does not exist", settingsFilePath.string()));
      } else {
        Log::info("Using settings file {}", settingsFilePath.string());
      }
      configuration.settingsFileName(settingsFilePath);
    }

    if (argc >= 3) {
        configuration.currentDirectory(args[2]);
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