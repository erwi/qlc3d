#include <io/energy-csv-writer.h>
#include <util/exception.h>
#include <fmt/format.h>
#include <iomanip>

EnergyCsvWriter::EnergyCsvWriter(const std::filesystem::path &path) {
    file_.open(path, std::ios::out | std::ios::trunc);
    if (!file_.is_open()) {
        RUNTIME_ERROR(fmt::format("Could not open energy CSV file: {}", path.string()));
    }
    // Write header row
    file_ << "time_s,iteration,elastic_J,thermotropic_J,electric_J,surface_J,total_J\n";
}

void EnergyCsvWriter::write(double time, int iteration, const EnergyResult &result) {
    file_ << std::setprecision(15) << std::scientific
          << time << ","
          << iteration << ","
          << result.elastic << ","
          << result.thermotropic << ","
          << result.electric << ","
          << result.surface << ","
          << result.total() << "\n";
}
