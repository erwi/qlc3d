#include <qlc3d.h>
#include <configuration.h>

void parseArgs(int argc, char **args, Configuration &configuration) {
    if (argc >= 2) {
        configuration.settingsFileName(args[1]);
    }

    if (argc >= 3) {
        configuration.currentDirectory(args[2]);
    }
}

int main(int argc, char **args) {
    Configuration configuration;
    parseArgs(argc, args, configuration);

    return runQlc3d(argc, args, configuration);
}