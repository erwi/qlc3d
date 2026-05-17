#include <catch.h>
#include <io/energy-csv-writer.h>
#include <energy/energy-result.h>
#include <test-util.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// ============================================================================
// Helper: read CSV file into rows of string values.
// ============================================================================
static std::vector<std::vector<std::string>> readCsv(const std::filesystem::path &path) {
    std::vector<std::vector<std::string>> rows;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        std::vector<std::string> cols;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            cols.push_back(cell);
        }
        rows.push_back(cols);
    }
    return rows;
}

TEST_CASE("EnergyCsvWriter: writes header and two data rows round-trip correctly", "[energy-csv-writer]") {
    // ARRANGE: two distinct EnergyResult objects
    EnergyResult r1{1.0, 2.0, 3.0, 4.0};   // elastic, thermotropic, electric, surface
    EnergyResult r2{10.0, 20.0, 30.0, 40.0};

    TestUtil::TemporaryDirectory tmpDir;
    auto csvPath = tmpDir.path() / "energy.csv";

    // ACT: write both rows
    {
        EnergyCsvWriter writer(csvPath);
        writer.write(1.5e-3, 1, r1);
        writer.write(2.5e-3, 2, r2);
    } // writer destructs and flushes

    // ASSERT: re-read and verify header + data
    auto rows = readCsv(csvPath);

    // Header row
    REQUIRE(rows.size() == 3);
    REQUIRE(rows[0][0] == "time_s");
    REQUIRE(rows[0][1] == "iteration");
    REQUIRE(rows[0][2] == "elastic_J");
    REQUIRE(rows[0][3] == "thermotropic_J");
    REQUIRE(rows[0][4] == "electric_J");
    REQUIRE(rows[0][5] == "surface_J");
    REQUIRE(rows[0][6] == "total_J");

    // First data row
    REQUIRE(std::stod(rows[1][0]) == Approx(1.5e-3).epsilon(1e-10));
    REQUIRE(std::stoi(rows[1][1]) == 1);
    REQUIRE(std::stod(rows[1][2]) == Approx(r1.elastic).epsilon(1e-10));
    REQUIRE(std::stod(rows[1][3]) == Approx(r1.thermotropic).epsilon(1e-10));
    REQUIRE(std::stod(rows[1][4]) == Approx(r1.electric).epsilon(1e-10));
    REQUIRE(std::stod(rows[1][5]) == Approx(r1.surface).epsilon(1e-10));
    REQUIRE(std::stod(rows[1][6]) == Approx(r1.total()).epsilon(1e-10));

    // Second data row
    REQUIRE(std::stod(rows[2][0]) == Approx(2.5e-3).epsilon(1e-10));
    REQUIRE(std::stoi(rows[2][1]) == 2);
    REQUIRE(std::stod(rows[2][2]) == Approx(r2.elastic).epsilon(1e-10));
    REQUIRE(std::stod(rows[2][6]) == Approx(r2.total()).epsilon(1e-10));
}


