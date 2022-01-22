#include <qlc3d.h>
#include <configuration.h>
#include <simulation-container.h>
#include <util/logging.h>

void parseArgs(int argc, char **args, Configuration &configuration) {
    if (argc >= 2) {
        configuration.settingsFileName(args[1]);
    }

    if (argc >= 3) {
        configuration.currentDirectory(args[2]);
    }
}

int runSimulation(Configuration &configuration) {
    try {
        configuration.readSettings();
        SimulationContainer simulation(configuration);
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
}

int main(int argc, char **args) {
    printInfo();

    Configuration configuration;
    parseArgs(argc, args, configuration);

    return runSimulation(configuration);
}