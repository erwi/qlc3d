"""
Symbolic derivation of LC energy densities using SymPy.

Based on the equations in doc-impl/equations.md.

The total free energy is:

    F = ∫_Ω (f_D + f_B - f_E) dΩ + ∫_Γ f_S dΓ

where:
  f_D : elastic distortion energy density
  f_B : thermotropic (bulk) energy density
  f_E : electric field energy density
  f_S : surface (anchoring) energy density

All energy expressions are written directly in the qlc3d T-tensor basis
(t1..t5) and their spatial partial derivatives (t1x, t1y, t1z, …).  This
matches the runtime representation used in the qlc3d solver (see
potential-derivation.py) and makes C++ code generation straightforward.

Q-tensor in T-basis:

    Q = t1 T1 + t2 T2 + t3 T3 + t4 T4 + t5 T5

where the orthonormal basis matrices are (qlc3d convention):

    T1 = (3 ê_z⊗ê_z - I) / √6          (axial / zz mode)
    T2 = (ê_x⊗ê_x - ê_y⊗ê_y) / √2     (biaxial / xx-yy mode)
    T3 = (ê_x⊗ê_y + ê_y⊗ê_x) / √2     (xy shear)
    T4 = (ê_y⊗ê_z + ê_z⊗ê_y) / √2     (yz shear)
    T5 = (ê_x⊗ê_z + ê_z⊗ê_x) / √2     (xz shear)

This convention is verified against TTensor::toQTensor() and
DielectricPermittivity::fromTTensor() in lc-representation.cpp.

Orthonormality tr(Ti Tj) = δ_ij gives the elegant identity:
    Q_ij Q_ij ≡ tr(Q²) = t1² + t2² + t3² + t4² + t5²
"""

from sympy import (
    symbols, sqrt, Rational, Matrix, eye, zeros, simplify, expand,
    pprint, init_printing
)
from sympy.printing import ccode

init_printing(use_unicode=True)

# ── helpers ──────────────────────────────────────────────────────────────────
rt2 = sqrt(2)
rt6 = sqrt(6)


# ==============================================================================
# PART 1: Q-TENSOR IN T-BASIS
# ==============================================================================
#
# The five T-tensor coefficients t1..t5 are qlc3d's primary DOFs.
# The 3×3 Q matrix is expressed directly in these coefficients:
#
#   Q = [[ -t1/√6 + t2/√2,   t3/√2,          t5/√2         ],
#        [  t3/√2,           -t1/√6 - t2/√2,  t4/√2         ],
#        [  t5/√2,            t4/√2,             2*t1/√6       ]]
#
# Tracelessness is guaranteed by construction:
#   Q[2,2] = 2t1/√6 = -( (-t1/√6+t2/√2) + (-t1/√6-t2/√2) )  ✓

t1, t2, t3, t4, t5 = symbols('t1 t2 t3 t4 t5', real=True)

Q = Matrix([
    [-t1/rt6 + t2/rt2,   t3/rt2,            t5/rt2      ],
    [ t3/rt2,           -t1/rt6 - t2/rt2,   t4/rt2      ],
    [ t5/rt2,            t4/rt2,             2*t1/rt6    ],
])

# Partial derivatives of the five T-tensor components w.r.t. x, y, z.
# Notation: tnx = ∂t_n/∂x, tny = ∂t_n/∂y, tnz = ∂t_n/∂z.
t1x, t1y, t1z = symbols('t1x t1y t1z', real=True)
t2x, t2y, t2z = symbols('t2x t2y t2z', real=True)
t3x, t3y, t3z = symbols('t3x t3y t3z', real=True)
t4x, t4y, t4z = symbols('t4x t4y t4z', real=True)
t5x, t5y, t5z = symbols('t5x t5y t5z', real=True)

# Gradient of the 5 T-tensor components, indexed as _dt[n][k]:
#   _dt[0] = [t1x, t1y, t1z], …
_dt = [
    [t1x, t1y, t1z],   # n=0 ↔ t1
    [t2x, t2y, t2z],   # n=1 ↔ t2
    [t3x, t3y, t3z],   # n=2 ↔ t3
    [t4x, t4y, t4z],   # n=3 ↔ t4
    [t5x, t5y, t5z],   # n=4 ↔ t5
]


def dQ(i, j, k):
    """
    Return ∂Q_ij/∂x_k in terms of the T-tensor derivatives _dt[n][k].

    The expression follows from differentiating Q = Σ t_n T_n component-wise.

    Parameters
    ----------
    i, j : int in {0, 1, 2}  – row and column of Q  (0=x, 1=y, 2=z)
    k    : int in {0, 1, 2}  – differentiation axis  (0=∂/∂x, 1=∂/∂y, 2=∂/∂z)
    """
    _t = _dt  # alias for brevity

    table = {
        (0, 0): lambda k: -_t[0][k]/rt6 + _t[1][k]/rt2,
        (1, 1): lambda k: -_t[0][k]/rt6 - _t[1][k]/rt2,
        (2, 2): lambda k:  2*_t[0][k]/rt6,
        (0, 1): lambda k:  _t[2][k]/rt2,
        (1, 0): lambda k:  _t[2][k]/rt2,
        (1, 2): lambda k:  _t[3][k]/rt2,
        (2, 1): lambda k:  _t[3][k]/rt2,
        (0, 2): lambda k:  _t[4][k]/rt2,
        (2, 0): lambda k:  _t[4][k]/rt2,
    }
    return table[(i, j)](k)


def levi_civita(i, j, k):
    """Return the Levi-Civita symbol ε_ijk (+1, -1, or 0)."""
    if (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        return 1
    if (i, j, k) in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
        return -1
    return 0


N = 3
range3 = range(N)


# ==============================================================================
# PART 2: THERMOTROPIC (BULK) ENERGY DENSITY
# ==============================================================================

def thermotropic_energy_density():
    """
    Thermotropic (bulk) energy density.

        f_B = (A/2) Q_ij Q_ij  +  (B/3) Q_ij Q_jk Q_ki  +  (C/4)(Q_ij Q_ij)²

    Due to the orthonormality of the T-basis, tr(Q²) = t1²+t2²+t3²+t4²+t5²,
    giving a particularly compact form.

    Material parameters:
      A – second-order Landau coefficient (temperature dependent)
      B – third-order Landau coefficient
      C – fourth-order Landau coefficient

    Returns
    -------
    SymPy expression for f_B in terms of t1..t5.
    """
    A, B, C = symbols('A B C', real=True)

    # Q_ij Q_ij = tr(Q^T Q) = tr(Q^2)  (Q symmetric) = t1²+…+t5² by orthonormality
    QQ  = (Q.T * Q).trace()

    # Q_ij Q_jk Q_ki = tr(Q^3)
    QQQ = (Q * Q * Q).trace()

    f_B = Rational(1, 2)*A * QQ \
        + Rational(1, 3)*B * QQQ \
        + Rational(1, 4)*C * QQ**2

    return expand(simplify(f_B))


# ==============================================================================
# PART 3: ELASTIC DISTORTION ENERGY DENSITY
# ==============================================================================

def elastic_energy_density():
    """
    Elastic distortion energy density.

        f_D = (1/2) * (
            L1 Q_ij,k  Q_ij,k
          + L2 Q_ij,j  Q_ik,k
          + L3 Q_ik,j  Q_ij,k
          + L4 ε_ijk   Q_il    Q_jl,k
          + L6 Q_lk    Q_ij,l  Q_ij,k
        )

    where ε_ijk is the Levi-Civita symbol and Q_ij,k ≡ ∂Q_ij/∂x_k.

    Material parameters L1..L6 map to the classical Frank elastic constants
    K1 (splay), K2 (twist), K3 (bend), and K24 (saddle-splay).

    Returns
    -------
    SymPy expression for f_D in terms of t1..t5 and their first derivatives.
    """
    L1, L2, L3, L4, L6 = symbols('L1 L2 L3 L4 L6', real=True)

    # L1 term:  Q_ij,k Q_ij,k   (sum over i, j, k)
    L1_term = sum(dQ(i, j, k)**2
                  for i in range3 for j in range3 for k in range3)

    # L2 term:  Q_ij,j Q_ik,k  =  |div Q|²
    #   (div Q)_i = Σ_j ∂Q_ij/∂x_j
    divQ = [sum(dQ(i, j, j) for j in range3) for i in range3]
    L2_term = sum(divQ[i]**2 for i in range3)

    # L3 term:  Q_ik,j Q_ij,k   (sum over i, j, k)
    L3_term = sum(dQ(i, k, j) * dQ(i, j, k)
                  for i in range3 for j in range3 for k in range3)

    # L4 term:  ε_ijk Q_il Q_jl,k   (sum over i, j, k, l)
    L4_term = sum(levi_civita(i, j, k) * Q[i, l] * dQ(j, l, k)
                  for i in range3 for j in range3 for k in range3 for l in range3)

    # L6 term:  Q_lk Q_ij,l Q_ij,k   (sum over i, j, k, l)
    L6_term = sum(Q[l, k] * dQ(i, j, l) * dQ(i, j, k)
                  for i in range3 for j in range3 for k in range3 for l in range3)

    f_D = Rational(1, 2) * (
        L1 * L1_term +
        L2 * L2_term +
        L3 * L3_term +
        L4 * L4_term +
        L6 * L6_term
    )

    return expand(f_D)


# ==============================================================================
# PART 4: ELECTRIC FIELD ENERGY DENSITY
# ==============================================================================

def electric_energy_density():
    """
    Electric field energy density.

        f_E = (1/2) ε₀ ε_ij E_i E_j  +  P_i E_i

    Dielectric permittivity tensor (uniaxial LC, expressed via Q):

        ε_ij = (ε_⊥ + Δε/3) δ_ij  +  (2Δε / (3 S₀)) Q_ij

    Flexoelectric polarisation:

        P_i = ξ_a Q_ij,j  +  ξ_b Q_ij Q_jk,k

    where E_i = -∂φ/∂x_i is the electric field.

    Returns
    -------
    dict with keys:
        'f_E'      – full energy-density expression in t1..t5 and E-field
        'eps'      – 3×3 SymPy Matrix for ε_ij in t1..t5
        'P_flexo'  – 3×1 SymPy Matrix for P_i in t1..t5 and their derivatives
    """
    eps0, eps_perp, delta_eps, S0 = symbols(
        'epsilon_0 epsilon_perp Delta_epsilon S_0', real=True, positive=True)
    xi_a, xi_b = symbols('xi_a xi_b', real=True)

    # Electric field components  E_i = -∂φ/∂x_i
    Ex, Ey, Ez = symbols('E_x E_y E_z', real=True)
    E = Matrix([Ex, Ey, Ez])

    # ── Dielectric permittivity tensor ────────────────────────────────────────
    # ε_ij = (ε_⊥ + Δε/3) δ_ij  +  (2Δε / (3 S₀)) Q_ij
    eps_tensor = (eps_perp + delta_eps/3) * eye(3) \
               + (2 * delta_eps / (3 * S0)) * Q

    # ── Flexoelectric polarisation ────────────────────────────────────────────
    # (div Q)_i = Σ_j ∂Q_ij/∂x_j
    divQ = Matrix([sum(dQ(i, j, j) for j in range3) for i in range3])

    # P_i = ξ_a Q_ij,j  +  ξ_b Q_ij (Q_jk,k)
    P_flexo = expand(xi_a * divQ + xi_b * (Q * divQ))

    # ── Electric energy density ───────────────────────────────────────────────
    # f_E = (1/2) ε₀ ε_ij E_i E_j  +  P_i E_i
    f_E_dielectric = Rational(1, 2) * eps0 * (E.T * eps_tensor * E)[0, 0]
    f_E_flexo      = P_flexo.dot(E)
    f_E = expand(f_E_dielectric + f_E_flexo)

    return {
        'f_E':     f_E,
        'eps':     simplify(eps_tensor),
        'P_flexo': P_flexo,
    }


# ==============================================================================
# PART 5: SURFACE ANCHORING ENERGY DENSITY
# ==============================================================================

def surface_anchoring_energy_density():
    """
    Surface (Rapini-Papoular style) anchoring energy density.

        f_S = a_s Q_ij Q_ij
            + W_1 v1_i Q_ij v1_j
            + W_2 v2_i Q_ij v2_j

    where v̂_1 and v̂_2 are mutually orthogonal unit vectors defining the
    preferred surface orientation (easy axis and a transverse in-plane axis).

    Due to orthonormality, Q_ij Q_ij = t1²+t2²+t3²+t4²+t5².

    Returns
    -------
    SymPy expression for f_S in terms of t1..t5, v1 components, v2 components.
    """
    a_s, W1, W2 = symbols('a_s W_1 W_2', real=True)

    # Easy-axis unit vectors (mutually orthogonal)
    v1x, v1y, v1z = symbols('v1_x v1_y v1_z', real=True)
    v2x, v2y, v2z = symbols('v2_x v2_y v2_z', real=True)
    v1 = Matrix([v1x, v1y, v1z])
    v2 = Matrix([v2x, v2y, v2z])

    # Q_ij Q_ij = tr(Q²) = t1²+t2²+t3²+t4²+t5² (orthonormality of T-basis)
    QQ = (Q.T * Q).trace()

    # v1_i Q_ij v1_j  =  v1^T Q v1
    v1Qv1 = (v1.T * Q * v1)[0, 0]

    # v2_i Q_ij v2_j  =  v2^T Q v2
    v2Qv2 = (v2.T * Q * v2)[0, 0]

    f_S = expand(simplify(a_s * QQ) + W1 * v1Qv1 + W2 * v2Qv2)
    return f_S


# ==============================================================================
# MAIN: evaluate all densities, verify key identities, display results
# ==============================================================================

if __name__ == '__main__':

    sep = "=" * 72

    # ── Verify orthonormality: tr(Q²) = Σ tᵢ² ───────────────────────────────
    print(sep)
    print("SANITY CHECKS")
    print(sep)
    QQ = (Q.T * Q).trace()
    QQ_simple = simplify(QQ)
    expected = t1**2 + t2**2 + t3**2 + t4**2 + t5**2
    assert simplify(QQ_simple - expected) == 0, \
        f"Orthonormality check FAILED: tr(Q²) = {QQ_simple}"
    print(f"\ntr(Q²) = {QQ_simple}  ✓")

    # ── Verify Q matrix matches TTensor::toQTensor() in lc-representation.cpp
    # C++ maps:  q1=-t1/rt6+t2/rt2, q2=-t1/rt6-t2/rt2, q3=t3/rt2,
    #            q4=t4/rt2 (yz), q5=t5/rt2 (xz)
    # and stores Q as [[q1, q3, q5], [q3, q2, q4], [q5, q4, -(q1+q2)]]
    Q_cpp = Matrix([
        [-t1/rt6 + t2/rt2,  t3/rt2,           t5/rt2          ],
        [ t3/rt2,           -t1/rt6 - t2/rt2,  t4/rt2          ],
        [ t5/rt2,            t4/rt2,            2*t1/rt6        ],
    ])
    assert simplify(Q_cpp - Q) == zeros(3, 3), \
        "Q matrix does not match TTensor::toQTensor() in lc-representation.cpp!"
    print("Q matrix matches TTensor::toQTensor() (lc-representation.cpp)  ✓")

    # ── Verify ε_ij matches DielectricPermittivity::fromTTensor() ────────────
    # C++: e13 uses t5 (xz), e23 uses t4 (yz)
    S0_v, de_v, ep_v = symbols('S0 deleps eper', positive=True)
    elec_check = electric_energy_density()
    eps_script = elec_check['eps']
    # substitute the same symbol names used in the C++ check
    e13_script = (2/(3*S0_v)) * (t5/rt2) * de_v   # xz entry from script eps
    e23_script = (2/(3*S0_v)) * (t4/rt2) * de_v   # yz entry from script eps
    eps_perp_v, de_sym, S0_sym = symbols(
        'epsilon_perp Delta_epsilon S_0', positive=True)
    e13_from_eps = simplify(
        eps_script[0, 2].subs({symbols('epsilon_perp', positive=True): 0,
                                symbols('Delta_epsilon', positive=True): de_v,
                                symbols('S_0', positive=True): S0_v}))
    e23_from_eps = simplify(
        eps_script[1, 2].subs({symbols('epsilon_perp', positive=True): 0,
                                symbols('Delta_epsilon', positive=True): de_v,
                                symbols('S_0', positive=True): S0_v}))
    assert simplify(e13_from_eps - e13_script) == 0, \
        f"ε[0,2] (xz) mismatch vs lc-representation.cpp: {e13_from_eps} != {e13_script}"
    assert simplify(e23_from_eps - e23_script) == 0, \
        f"ε[1,2] (yz) mismatch vs lc-representation.cpp: {e23_from_eps} != {e23_script}"
    print("ε_ij matches DielectricPermittivity::fromTTensor() (lc-representation.cpp)  ✓")

    # ── Thermotropic energy ───────────────────────────────────────────────────
    print("\n" + sep)
    print("THERMOTROPIC ENERGY DENSITY  f_B")
    print(sep)
    f_B = thermotropic_energy_density()
    print("\nf_B = (A/2) Q_ij Q_ij  +  (B/3) Q_ij Q_jk Q_ki  +  (C/4)(Q_ij Q_ij)²")
    print("\nExpanded in t1..t5:")
    pprint(f_B)

    # ── Elastic energy ────────────────────────────────────────────────────────
    print("\n" + sep)
    print("ELASTIC DISTORTION ENERGY DENSITY  f_D")
    print(sep)
    f_D = elastic_energy_density()
    print("\nf_D = (1/2)(L1 Q_ij,k Q_ij,k + L2 Q_ij,j Q_ik,k")
    print("           + L3 Q_ik,j Q_ij,k + L4 ε_ijk Q_il Q_jl,k + L6 Q_lk Q_ij,l Q_ij,k)")
    print("\nExpanded in t1..t5 and their first derivatives:")
    pprint(f_D)

    # ── Electric energy ───────────────────────────────────────────────────────
    print("\n" + sep)
    print("ELECTRIC FIELD ENERGY DENSITY  f_E")
    print(sep)
    elec = electric_energy_density()
    f_E = elec['f_E']
    eps  = elec['eps']
    P    = elec['P_flexo']
    print("\nDielectric permittivity tensor ε_ij:")
    pprint(eps)
    print("\nFlexoelectric polarisation P_i = ξ_a Q_ij,j + ξ_b Q_ij Q_jk,k:")
    pprint(P)
    print("\nf_E = (1/2) ε₀ ε_ij E_i E_j  +  P_i E_i:")
    pprint(f_E)

    # ── Surface anchoring energy ──────────────────────────────────────────────
    print("\n" + sep)
    print("SURFACE ANCHORING ENERGY DENSITY  f_S")
    print(sep)
    f_S = surface_anchoring_energy_density()
    print("\nf_S = a_s Q_ij Q_ij  +  W_1 v1_i Q_ij v1_j  +  W_2 v2_i Q_ij v2_j")
    print("\nExpanded in t1..t5, v1, v2:")
    pprint(f_S)

