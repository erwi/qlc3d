#pragma once

/**
 * @brief Result of an LC free energy calculation.
 *
 * Pure physics result with no simulation metadata (time, iteration).
 * Time and iteration are passed separately when writing output.
 * All values are in Joules.
 */
struct EnergyResult {
    /** Total elastic distortion energy (Frank elastic terms L1..L6 combined) [J] */
    double elastic = 0.0;
    /** Landau-de Gennes bulk ordering (thermotropic) energy [J] */
    double thermotropic = 0.0;
    /** Dielectric and flexoelectric electric energy [J] */
    double electric = 0.0;
    /** Surface anchoring energy [J] */
    double surface = 0.0;

    /**
     * @brief Sum of all energy contributions.
     * @return elastic + thermotropic + electric + surface [J]
     */
    [[nodiscard]] double total() const {
        return elastic + thermotropic + electric + surface;
    }
};

