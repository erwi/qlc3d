#pragma once
#include <energy/energy-result.h>

class LC;
class Geometry;
class SolutionVector;
class Alignment;

/**
 * @brief Integrate volume energy contributions over all LC tetrahedra.
 *
 * Iterates over all tetrahedral elements with material number ≤ MAT_DOMAIN7,
 * applies Gaussian quadrature (Keast8) using either linear (TET4) or quadratic
 * (TET10) shape functions depending on the mesh element type, and accumulates
 * elastic, thermotropic, and electric energy contributions.
 *
 * The returned @c EnergyResult has @c surface = 0.0; surface energy must be
 * added separately by the caller.
 *
 * @param lc    LC material parameters.
 * @param geom  Geometry containing the tetrahedral mesh and coordinates.
 * @param v     Potential solution vector (electric field source).
 * @param q     Q-tensor (TTensor) solution vector.
 * @return      Integrated volume energy components [J]. Surface is zero.
 */
EnergyResult integrateVolumeEnergy(const LC &lc, const Geometry &geom,
                                   const SolutionVector &v, const SolutionVector &q);

/**
 * @brief Integrate surface anchoring energy over all weak-anchoring triangle surfaces.
 *
 * Iterates over triangle elements, identifies which FIXLC anchoring condition applies
 * (via material number → FIXLC number), and accumulates the Rapini-Papoular energy:
 *   f_S = W·A·(q1²+…+q5²) + W·K1·(v̂₁·Q·v̂₁) + W·K2·(v̂₂·Q·v̂₂)
 * where A = (K1+K2)/(S0·6).
 *
 * Only surfaces with weak anchoring contribute; strong-anchoring surfaces are skipped.
 *
 * @param alignment  Alignment object containing anchoring conditions for each FIXLC surface.
 * @param geom       Geometry containing triangle mesh and coordinates.
 * @param q          Q-tensor (TTensor) solution vector.
 * @param S0         Equilibrium order parameter.
 * @return           Integrated surface anchoring energy [J].
 */
double integrateSurfaceEnergy(const Alignment &alignment, const Geometry &geom,
                              const SolutionVector &q, double S0);

/**
 * @brief Calculate the total LC free energy by combining volume and surface integrals.
 *
 * Calls @c integrateVolumeEnergy and @c integrateSurfaceEnergy and returns them as a
 * single @c EnergyResult with all four components populated.
 *
 * @param lc        LC material parameters.
 * @param geom      Simulation geometry (mesh + coordinates).
 * @param v         Potential solution vector.
 * @param q         Q-tensor (TTensor) solution vector.
 * @param alignment Surface anchoring conditions.
 * @return          All energy components and total [J].
 */
EnergyResult integrateEnergy(const LC &lc, const Geometry &geom,
                             const SolutionVector &v, const SolutionVector &q,
                             const Alignment &alignment);

