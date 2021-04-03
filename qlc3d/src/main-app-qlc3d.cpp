#include <qlc3d.h>
#include <configuration.h>
#include <simulation-container.h>
void parseArgs(int argc, char **args, Configuration &configuration) {
    if (argc >= 2) {
        configuration.settingsFileName(args[1]);
    }

    if (argc >= 3) {
        configuration.currentDirectory(args[2]);
    }
}

int runSimulation(SimulationContainer &simulation) {
    try {
        simulation.initialise();
        while (simulation.hasIteration()) {
            simulation.runIteration();
            // TODO: auto state = simulation.getCurrentState() and do something with it?
        }
        simulation.postSimulationTasks();
    } catch(std::exception &e) {
        std::cerr << "ERROR - An error has occurred: " << e.what() << std::endl;
        return 1;
    } catch(...) {
        std::cerr << "ERROR - An error exception has occurred - bye!" << std::endl;
        return 3;
    }
    std::cout << "simulation has ended without errors" << std::endl;
    return 0;
}

int main(int argc, char **args) {
    Configuration configuration;
    parseArgs(argc, args, configuration);

    SimulationContainer simulation(configuration);
    return runSimulation(simulation);
}