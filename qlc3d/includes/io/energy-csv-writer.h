#pragma once
#include <energy/energy-result.h>

#include <fstream>
#include <filesystem>

/**
 * @brief Writes LC free energy results to a CSV file.
 *
 * Opens a CSV file on construction and appends one row per call to @c write().
 * The header row is written automatically at construction.
 *
 * CSV columns:
 *   time_s, iteration, elastic_J, thermotropic_J, electric_J, surface_J, total_J
 *
 * Uses @c std::ofstream; the file is closed when the destructor runs.
 */
class EnergyCsvWriter {
    std::ofstream file_;

public:
    /**
     * @brief Open (or create) a CSV file and write the header row.
     * @param path  Absolute or relative path to the output CSV file.
     * @throws std::runtime_error if the file cannot be opened.
     */
    explicit EnergyCsvWriter(const std::filesystem::path &path);

    ~EnergyCsvWriter() = default;

    // Non-copyable; moveable
    EnergyCsvWriter(const EnergyCsvWriter &) = delete;
    EnergyCsvWriter &operator=(const EnergyCsvWriter &) = delete;
    EnergyCsvWriter(EnergyCsvWriter &&) = default;
    EnergyCsvWriter &operator=(EnergyCsvWriter &&) = default;

    /**
     * @brief Append one row to the CSV file.
     * @param time       Simulation time [s].
     * @param iteration  Iteration index.
     * @param result     Energy components and total [J].
     */
    void write(double time, int iteration, const EnergyResult &result);
};

