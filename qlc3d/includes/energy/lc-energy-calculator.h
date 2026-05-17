#pragma once
#include <energy/energy-result.h>

class LC;
class Geometry;
class SolutionVector;
class Alignment;

/**
 * @brief Top-level LC free energy calculator.
 *
 * Combines volume integration (elastic, thermotropic, electric) and surface
 * anchoring integration into a single @c EnergyResult.
 *
 * Usage:
 * @code
 *   LcEnergyCalculator calc;
 *   EnergyResult r = calc.calculate(lc, geom, v, q, alignment);
 * @endcode
 *
 * Writing results to a file is the responsibility of the caller (e.g. using
 * @c EnergyCsvWriter from @c io/energy-csv-writer.h).
 */
class LcEnergyCalculator {
public:
    LcEnergyCalculator() = default;

    /**
     * @brief Calculate all free energy contributions.
     *
     * @param lc        LC material parameters.
     * @param geom      Simulation geometry (mesh + coordinates).
     * @param v         Potential solution vector.
     * @param q         Q-tensor (TTensor) solution vector.
     * @param alignment Surface anchoring conditions.
     * @return          All energy components and total [J].
     */
    [[nodiscard]] EnergyResult calculate(const LC &lc, const Geometry &geom,
                                         const SolutionVector &v, const SolutionVector &q,
                                         const Alignment &alignment) const;
};
