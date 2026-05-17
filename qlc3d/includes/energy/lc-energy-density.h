#pragma once
#include <lc.h>
#include <geom/vec3.h>

/**
 * @brief Data interpolated at a single Gauss quadrature point for energy density evaluation.
 *
 * All q1..q5 components are qlc3d TTensor DOFs. Spatial derivatives are
 * in physical units (SI metres after the MICROMETER_TO_METER conversion).
 * Ex, Ey, Ez are the electric field components E_i = -∂φ/∂x_i.
 */
struct GaussPointData {
    // TTensor components (qlc3d DOFs)
    double q1 = 0, q2 = 0, q3 = 0, q4 = 0, q5 = 0;
    // Spatial derivatives of each TTensor component
    double q1x = 0, q1y = 0, q1z = 0;
    double q2x = 0, q2y = 0, q2z = 0;
    double q3x = 0, q3y = 0, q3z = 0;
    double q4x = 0, q4y = 0, q4z = 0;
    double q5x = 0, q5y = 0, q5z = 0;
    // Electric field components E_i = -∂Φ/∂x_i
    double Ex = 0, Ey = 0, Ez = 0;
};

/**
 * @brief Material parameters required for LC energy density evaluation.
 *
 * Constructed from an @c LC object via @c lcToEnergyMaterialParams.
 */
struct EnergyMaterialParams {
    double S0;          ///< Equilibrium order parameter
    double K11;         ///< Splay elastic constant [N]
    double K22;         ///< Twist elastic constant [N]
    double K33;         ///< Bend elastic constant [N]
    double A;           ///< Thermotropic bulk coefficient [J/m³]
    double B;           ///< Thermotropic bulk coefficient [J/m³]
    double C;           ///< Thermotropic bulk coefficient [J/m³]
    double eps_par;     ///< Dielectric permittivity parallel to director (dimensionless)
    double eps_per;     ///< Dielectric permittivity perpendicular to director (dimensionless)
    double e1;          ///< Flexoelectric coefficient [C/m]
    double e3;          ///< Flexoelectric coefficient [C/m]
    double q0;          ///< Chirality wavenumber 2π/p0 [1/m]; 0 if non-chiral
};

/**
 * @brief Elastic distortion energy density at a single Gauss point.
 *
 * Computes the total Frank elastic energy density combining splay (K11), twist (K22)
 * and bend (K33) contributions using the G1, G2, G4, G6 invariants in TTensor form.
 * Note: The G3 (K24) saddle-splay invariant is present in the underlying equations
 * but is not included here because it did not appear in the original energy.cpp
 * splay/twist/bend decomposition.
 *
 * @param p  Gauss-point data (TTensor values and their spatial derivatives).
 * @param m  Material parameters.
 * @return   Elastic energy density [J/m³].
 */
double elasticEnergyDensity(const GaussPointData &p, const EnergyMaterialParams &m);

/**
 * @brief Thermotropic (Landau-de Gennes) bulk energy density at a single Gauss point.
 *
 * Returns the raw f_th without any ground-state offset. The value is negative at
 * the equilibrium state and is not zero there. Use the returned values to monitor
 * relative changes over time (e.g. that total energy decreases), not as an absolute
 * measure from the ground state.
 *
 * @param p  Gauss-point data (TTensor values; derivatives not used).
 * @param m  Material parameters (A, B, C, S0).
 * @return   Thermotropic energy density [J/m³].
 */
double thermotropicEnergyDensity(const GaussPointData &p, const EnergyMaterialParams &m);

/**
 * @brief Electric (dielectric + flexoelectric) energy density at a single Gauss point.
 *
 * Computes -f_E = -(dielectric + flexoelectric) contribution.
 * The electric field components are taken from @p p as E_i = -∂φ/∂x_i.
 *
 * @param p  Gauss-point data (TTensor values, derivatives, and electric field).
 * @param m  Material parameters (eps_par, eps_per, e1, e3, S0).
 * @return   Electric energy density [J/m³].
 */
double electricEnergyDensity(const GaussPointData &p, const EnergyMaterialParams &m);

/**
 * @brief Construct @c EnergyMaterialParams from an @c LC material object.
 *
 * @param lc  The liquid crystal material parameters.
 * @return    Corresponding @c EnergyMaterialParams struct.
 */
EnergyMaterialParams lcToEnergyMaterialParams(const LC &lc);

/**
 * @brief Rapini-Papoular surface anchoring energy density at a single Gauss point.
 *
 * Evaluates:
 *   f_S = W·A·(q1²+…+q5²)  +  W·K1·(v̂₁·Q·v̂₁)  +  W·K2·(v̂₂·Q·v̂₂)
 * where A = (K1+K2)/(S0·6) and the v̂·Q·v̂ inner product is expressed in
 * TTensor basis form.
 *
 * The returned value is the raw anchoring energy density with no ground-state offset
 * applied. It is negative at the preferred orientation and is not zero there. Use
 * the returned values to monitor relative changes over time, not as an absolute
 * measure from the ground state.
 *
 * @param q1..q5  TTensor DOF values at the Gauss point.
 * @param v1      First principal anchoring axis vector.
 * @param v2      Second principal anchoring axis vector (surface normal for homeotropic).
 * @param W       Anchoring strength [J/m²]; negative for WeakHomeotropic convention.
 * @param K1      First anchoring coefficient (dimensionless).
 * @param K2      Second anchoring coefficient (dimensionless).
 * @param S0      Equilibrium order parameter.
 * @return        Surface energy density [J/m²].
 */
double surfaceEnergyDensity(double q1, double q2, double q3, double q4, double q5,
                            const Vec3 &v1, const Vec3 &v2,
                            double W, double K1, double K2, double S0);

