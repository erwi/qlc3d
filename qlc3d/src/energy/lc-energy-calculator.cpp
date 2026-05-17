#include <energy/lc-energy-calculator.h>
#include <energy/lc-energy-integrator.h>
#include <lc.h>
#include <geometry.h>
#include <solutionvector.h>
#include <alignment.h>

EnergyResult LcEnergyCalculator::calculate(const LC &lc, const Geometry &geom,
                                            const SolutionVector &v, const SolutionVector &q,
                                            const Alignment &alignment) const {
    return integrateEnergy(lc, geom, v, q, alignment);
}
