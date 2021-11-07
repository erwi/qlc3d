
#include <catch.h>
#include <test-util.h>
#include <reader.h>

TEST_CASE("Environment variable substitution in settings file") {
    using namespace std;

    // set environment variable
    putenv(strdup("QLC3D_TEST=WOO"));

    Reader reader;
    reader.setEnvironmentVariableSubstitution(true);

    SECTION("Substitute one environment value") {
        // create temp test file
        string fileContents = "key=${QLC3D_TEST }";
        auto settingsFile = TestUtil::TemporaryFile::withContents(fileContents);

        reader.readSettingsFile(settingsFile.name());
        string value = reader.getValueByKey<string>("key");
        REQUIRE("WOO" == value);
    }

    SECTION("Invalid substitution string format") {
        string fileContents = "key=${WAA"; // fileContents with invalid formatting
        auto settingsFile = TestUtil::TemporaryFile::withContents(fileContents);

        try {
            reader.readSettingsFile(settingsFile.name());
        } catch (ReaderError &error) {
            // should throw
            REQUIRE(error.errorMessage.find("Invalid environment variable substitution format") == 0);
            return;
        }
        REQUIRE( false ); // should not reach this
    }

    SECTION("Environment variable not set") {
        string fileContents = "key=${BOO}"; // fileContents with env variable not set
        auto settingsFile = TestUtil::TemporaryFile::withContents(fileContents);

        try {
            reader.readSettingsFile(settingsFile.name());
        } catch (ReaderError &error) {
            // should throw
            REQUIRE(error.errorMessage.find("No such environment variable") == 0);
            return;
        }
        REQUIRE( false ); // should not reach this
    }
}