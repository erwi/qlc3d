"""
Symbolic derivation and C++ code generation for the FEM electric potential solver
in qlc3d.

═══════════════════════════════════════════════════════════════════════════════
PROBLEM STATEMENT
═══════════════════════════════════════════════════════════════════════════════

We want to find the electric potential V(x,y,z) inside a liquid crystal cell.
Assuming no free charges, the governing equation is Gauss' law:

    ∇·D = 0

where D = ε₀ ε_r E = -ε₀ ε_r ∇V is the electric displacement field.
Substituting, with ε ≡ ε_r (relative permittivity, a 3×3 tensor for an
anisotropic medium like LC):

    ∇·(ε·∇V) = 0   on Ω

with Dirichlet boundary conditions V = V_D on ∂Ω_D (electrodes).

This is the anisotropic Laplace equation (Poisson's equation with zero source).
The dielectric anisotropy is caused by the orientational order of the liquid
crystal molecules, captured by the Q-tensor order parameter.

Script sections:
  1. Q-tensor – definition and physical interpretation
  2. T-tensor – internal (orthonormal) re-parametrisation used in qlc3d
  3. Dielectric permittivity tensor ε from Q-tensor
  4. Weak (variational) form – symbolic derivation via product rule and
      integration by parts (divergence theorem); flexoelectric charge
      density and boundary terms are retained explicitly
  5. FEM discretisation – local element stiffness matrix K^e_ij plus the
      flexoelectric load vector and surface flux terms derived by
      substituting the Galerkin ansatz into the weak form; Neumann
      boundary stiffness at the LC–vacuum interface derived from the IBP
      surface term and the unmeshed vacuum exterior contribution
  6. C++ code generation via sympy.printing.ccode

Notes on conventions used in qlc3d:
  * Surface normals at Neumann boundaries point INWARD (into the LC
    element), i.e. n_qlc3d = -n_math where n_math is the outward normal
    used in classical mathematics / integration-by-parts formulas.  This
    convention is shared with the surface-anchoring calculation.
  * The Neumann boundaries represent the external surfaces of the LC mesh
    that face the surrounding unmodelled medium (vacuum, ε_vac = I).
    Dirichlet boundaries (electrodes) are handled separately.
"""

from sympy import (
    symbols, sqrt, Rational, Matrix, eye, simplify, expand,
    pprint, init_printing, Function, Eq, Integral
)
from sympy.printing import ccode

init_printing(use_unicode=True)

# ── helpers ─────────────────────────────────────────────────────────────────
rt2 = sqrt(2)
rt6 = sqrt(6)

# ==============================================================================
# PART 1: THE Q-TENSOR
# ==============================================================================

# Scalar order parameter S and director components (unit vector n).
nx, ny, nz, S = symbols('n_x n_y n_z S', real=True)

# The Q-tensor is the traceless symmetric second-rank tensor
#   Q_ij = (S/2) * (3 n_i n_j  −  δ_ij)
# It vanishes in the isotropic phase (S=0) and captures orientational order.
# For a 3D system it has 5 independent components (symmetric, traceless).

# qlc3d stores Q via 5 independent components (q1..q5):
#   Q = [[q1,  q3,  q5 ],
#        [q3,  q2,  q4 ],
#        [q5,  q4,  -(q1+q2)]]   ← tracelessness enforces Q_33 = -(Q_11+Q_22)

q1, q2, q3, q4, q5 = symbols('q1 q2 q3 q4 q5', real=True)

Q = Matrix([
    [q1,        q3,        q5       ],
    [q3,        q2,        q4       ],
    [q5,        q4,        -(q1+q2) ],
])

# Components from the director:
#   q1 = (S/2)*(3*nx^2-1),  q2 = (S/2)*(3*ny^2-1),  q3 = (S/2)*(3*nx*ny)
#   q4 = (S/2)*(3*ny*nz),   q5 = (S/2)*(3*nx*nz)

# ==============================================================================
# PART 2: T-TENSOR (qlc3d internal representation)
# ==============================================================================

# The T-tensor is an orthonormal re-parametrisation of the five Q components,
# choosing basis matrices that diagonalise the most common LC configurations.
#
# Forward Q → T:
#   t1 = −(q1+q2)*√6/2       (zz-vs-xx+yy mode)
#   t2 =  (q1−q2)/2  * √2    (xx-vs-yy mode)
#   t3 =  q3 * √2            (xy shear)
#   t4 =  q4 * √2            (yz shear)
#   t5 =  q5 * √2            (xz shear)
#
# Inverse T → Q:
#   q1 = −t1/√6 + t2/√2,  q2 = −t1/√6 − t2/√2
#   q3 = t3/√2,            q4 = t4/√2,   q5 = t5/√2

t1, t2, t3, t4, t5 = symbols('t1 t2 t3 t4 t5', real=True)

q1_from_T = -t1/rt6 + t2/rt2
q2_from_T = -t1/rt6 - t2/rt2
q3_from_T =  t3/rt2
q4_from_T =  t4/rt2
q5_from_T =  t5/rt2

# ==============================================================================
# PART 3: DIELECTRIC PERMITTIVITY TENSOR FROM Q-TENSOR
# ==============================================================================

# For a uniaxial LC the permittivity tensor is:
#   ε_ij = ε_perp δ_ij + Δε n_i n_j
# where Δε = ε_∥ − ε_⊥ is the dielectric anisotropy.
#
# The outer product n⊗n is related to the Q-tensor by:
#   Q_ij = (S/2)*(3 n_i n_j − δ_ij)  →  n_i n_j = 2Q_ij/(3S₀) + δ_ij/3
#
# Substituting:
#   ε_ij = ε_perp δ_ij + Δε * (2Q_ij/(3S₀) + δ_ij/3)
#         = (ε_perp + Δε/3) δ_ij + (2Δε)/(3S₀) Q_ij

S0, deleps, eps_perp = symbols('S_0 Delta_epsilon epsilon_perp',
                                real=True, positive=True)
I3 = eye(3)

# Build the permittivity tensor in Q-component form
eps_tensor_Q = eps_perp * I3 + deleps * (Rational(2, 3) / S0 * Q + I3 / 3)
eps_tensor_Q = simplify(eps_tensor_Q)

# Extract the 6 independent symmetric components (exx, eyy, ezz, exy, exz, eyz)
exx_q = eps_tensor_Q[0, 0]
eyy_q = eps_tensor_Q[1, 1]
ezz_q = eps_tensor_Q[2, 2]
exy_q = eps_tensor_Q[0, 1]
exz_q = eps_tensor_Q[0, 2]
eyz_q = eps_tensor_Q[1, 2]

# Express in T-tensor components (the representation used in qlc3d at runtime)
sub_q_to_t = {q1: q1_from_T, q2: q2_from_T, q3: q3_from_T,
              q4: q4_from_T, q5: q5_from_T}

exx_T = simplify(exx_q.subs(sub_q_to_t))
eyy_T = simplify(eyy_q.subs(sub_q_to_t))
ezz_T = simplify(ezz_q.subs(sub_q_to_t))
exy_T = simplify(exy_q.subs(sub_q_to_t))
exz_T = simplify(exz_q.subs(sub_q_to_t))
eyz_T = simplify(eyz_q.subs(sub_q_to_t))

# ==============================================================================
# PART 4: WEAK (VARIATIONAL) FORM – SYMBOLIC DERIVATION
# ==============================================================================

# Set up 3D spatial coordinates and treat V, w as generic smooth functions.
xc, yc, zc = symbols('x y z', real=True)
V = Function('V')(xc, yc, zc)   # electric potential (unknown)
w = Function('w')(xc, yc, zc)   # test function (vanishes on Dirichlet boundary)

# Generic flexoelectric polarisation field and its divergence.
Px_fun = Function('P_x')(xc, yc, zc)
Py_fun = Function('P_y')(xc, yc, zc)
Pz_fun = Function('P_z')(xc, yc, zc)
P_fun = Matrix([Px_fun, Py_fun, Pz_fun])

# Gradient column vectors  ∇V and ∇w
grad_V = Matrix([V.diff(xc), V.diff(yc), V.diff(zc)])
grad_w = Matrix([w.diff(xc), w.diff(yc), w.diff(zc)])

# Use generic (piecewise-constant within an element) permittivity ε for this
# derivation; components will be substituted from the Q/T-tensor in Part 6.
exx_s, eyy_s, ezz_s, exy_s, exz_s, eyz_s = symbols(
    'eps_xx eps_yy eps_zz eps_xy eps_xz eps_yz', real=True)
eps_deriv = Matrix([
    [exx_s, exy_s, exz_s],
    [exy_s, eyy_s, eyz_s],
    [exz_s, eyz_s, ezz_s],
])

# ── Step 0: the governing PDE ────────────────────────────────────────────────
# Flux field  D = ε · ∇V - P
eps_grad_V = eps_deriv * grad_V          # 3×1 column vector
div_eps_grad_V = eps_grad_V[0].diff(xc) + eps_grad_V[1].diff(yc) + eps_grad_V[2].diff(zc)
div_P = Px_fun.diff(xc) + Py_fun.diff(yc) + Pz_fun.diff(zc)
rho_flexo = div_P
D = eps_grad_V - P_fun

# Divergence  ∇·D  (with ε treated as constant inside each element)
div_D = D[0].diff(xc) + D[1].diff(yc) + D[2].diff(zc)

# Boundary normal used to keep the surface term explicit.
nx_b, ny_b, nz_b = symbols('n_x n_y n_z', real=True)
n_boundary = Matrix([nx_b, ny_b, nz_b])
surface_flux = expand(D.dot(n_boundary))

print("=" * 72)
print("PART 4: WEAK (VARIATIONAL) FORM")
print("=" * 72)
print("\nGoverning PDE:  ∇·(ε·∇V) = ρ_flexo = ∇·P")
print("Equivalent conservation form:  ∇·(ε·∇V - P) = 0")
print("Expanded PDE:")
pprint(Eq(expand(div_eps_grad_V), rho_flexo))

# ── Step 1: weighted residual ────────────────────────────────────────────────
# Multiply PDE by test function w and integrate over domain Ω:
#   ∫_Ω  w · ∇·(ε·∇V) dΩ = ∫_Ω w · ρ_flexo dΩ
weighted_residual = w * (div_eps_grad_V - rho_flexo)

print("\nStep 1 – weighted residual  w · [∇·(ε·∇V) - ρ_flexo] :")
pprint(expand(weighted_residual))

# ── Step 2: product / Leibniz rule ───────────────────────────────────────────
# Identity:  w · ∇·F  =  ∇·(w F)  −  ∇w · F
# Compute both sides symbolically and verify the identity for the
# electric flux term and the flexoelectric charge-density term.

wepsV = w * eps_grad_V
div_wepsV = wepsV[0].diff(xc) + wepsV[1].diff(yc) + wepsV[2].diff(zc)
gradw_dot_epsV = grad_w.dot(eps_grad_V)

wP = w * P_fun
div_wP = wP[0].diff(xc) + wP[1].diff(yc) + wP[2].diff(zc)
gradw_dot_P = grad_w.dot(P_fun)

# LHS =  w · ∇·(ε∇V) ;  RHS = ∇·(w ε∇V) − ∇w·(ε∇V)
identity_residual_eps = simplify(div_wepsV - gradw_dot_epsV - w * div_eps_grad_V)
assert identity_residual_eps == 0, \
    f"Product-rule identity w·∇·(ε∇V) = ∇·(w ε∇V) − ∇w·(ε∇V) FAILED: {identity_residual_eps}"

# LHS =  w · ∇·P  ;  RHS = ∇·(wP) − ∇w·P
identity_residual_P = simplify(div_wP - gradw_dot_P - w * div_P)
assert identity_residual_P == 0, \
    f"Product-rule identity w·∇·P = ∇·(wP) − ∇w·P FAILED: {identity_residual_P}"

print("\nStep 2 – product rule identities verified ✓")
print("  w·∇·(ε∇V) = ∇·(w ε∇V) − ∇w·(ε∇V)")
print("  w·∇·P     = ∇·(w P)     − ∇w·P")

# ── Step 3: split the volume integral ────────────────────────────────────────
# ∫_Ω w·∇·(ε∇V) dΩ  =  ∫_Ω w·ρ_flexo dΩ

print("\nStep 3 – split volume integral:")
print("  ∫_Ω w·∇·(ε∇V) dΩ  =  ∫_Ω w·ρ_flexo dΩ")
print("  with ρ_flexo = ∇·P")

# ── Step 4: divergence theorem on the first term ─────────────────────────────
# ∫_Ω ∇·(w ε∇V) dΩ  =  ∮_∂Ω w(ε∇V)·n dS
# ∫_Ω ∇·(w P) dΩ      =  ∮_∂Ω w P·n dS
print("\nStep 4 – divergence theorem:")
print("  ∫_Ω ∇·(w ε∇V) dΩ  =  ∮_∂Ω w(ε∇V)·n dS")
print("  ∫_Ω ∇·(w P) dΩ      =  ∮_∂Ω w P·n dS")

# ── Step 5: apply Dirichlet boundary condition w = 0 on ∂Ω_D ────────────────
print("\nStep 5 – Dirichlet BC w = 0 on ∂Ω_D:")
print("  the surface term vanishes only on ∂Ω_D; natural boundaries keep it")
print("  ∮_∂Ω w((ε·∇V - P)·n) dS remains on the non-Dirichlet part")
print("\nWeak form:")
print("  ∫_Ω ∇w · (ε·∇V - P) dΩ  =  ∮_∂Ω w((ε·∇V - P)·n) dS")
print("  equivalently,  ∫_Ω ∇w·(ε·∇V) dΩ = ∫_Ω ∇w·P dΩ + boundary terms")

# ── The integrand of the weak bilinear form ───────────────────────────────────
weak_integrand_sym = expand(gradw_dot_epsV)
print("\nWeak-form integrand  ∇w · (ε · ∇V)  (fully expanded) =")
pprint(weak_integrand_sym)

polarization_integrand_sym = expand(gradw_dot_P)
print("\nFlexoelectric source integrand  ∇w · P  (fully expanded) =")
pprint(polarization_integrand_sym)

boundary_flux_sym = expand(surface_flux)
print("\nBoundary flux integrand  (ε·∇V - P)·n  =")
pprint(boundary_flux_sym)

# ==============================================================================
# PART 5: FEM DISCRETISATION – LOCAL STIFFNESS MATRIX DERIVED FROM WEAK FORM
# ==============================================================================

print("\n" + "=" * 72)
print("PART 5: FEM DISCRETISATION")
print("=" * 72)

# Galerkin ansatz:
#   V(x) ≈ Σ_j  Vj · Nj(x,y,z)     (sum over element nodes j)
#   w    =  Ni(x,y,z)               (test function = i-th shape function)
#
# For the bilinear form we consider one j-contribution at a time:
#   ∂V/∂x → Njx · Vj,   ∂V/∂y → Njy · Vj,   ∂V/∂z → Njz · Vj
#   ∂w/∂x → Nix,        ∂w/∂y → Niy,         ∂w/∂z → Niz
# where Nix = ∂Ni/∂x etc. are evaluated at a quadrature point.

# Shape-function gradient symbols at a quadrature point
Nix, Niy, Niz = symbols('Nix Niy Niz', real=True)   # ∂Ni/∂x, ∂Ni/∂y, ∂Ni/∂z
Njx, Njy, Njz = symbols('Njx Njy Njz', real=True)   # ∂Nj/∂x, ∂Nj/∂y, ∂Nj/∂z
Vj = symbols('Vj', real=True)                        # nodal potential at j

# Substitute the Galerkin ansatz into the weak-form integrand.
# SymPy treats V.diff(xc) as Derivative(V(x,y,z), x), which subs() handles.
K_fem_with_Vj = weak_integrand_sym.subs({
    V.diff(xc): Njx * Vj,
    V.diff(yc): Njy * Vj,
    V.diff(zc): Njz * Vj,
    w.diff(xc): Nix,
    w.diff(yc): Niy,
    w.diff(zc): Niz,
})

# The expression is linear in Vj; factor it out to obtain the bilinear form.
K_integrand_derived = expand(K_fem_with_Vj / Vj)
print("\nGalerkin ansatz substituted into weak-form integrand.")
print("Local stiffness integrand  K^e_ij = ∇Ni · ε · ∇Nj  =")
pprint(K_integrand_derived)

# Verify K_integrand equals the manually written result from the original
# formulation  ∇Ni · ε_mat · ∇Nj  (as a sanity check).
grad_Ni = Matrix([Nix, Niy, Niz])
grad_Nj = Matrix([Njx, Njy, Njz])
eps_check = Matrix([
    [exx_s, exy_s, exz_s],
    [exy_s, eyy_s, eyz_s],
    [exz_s, eyz_s, ezz_s],
])
K_direct = expand((grad_Ni.T * eps_check * grad_Nj)[0, 0])
assert simplify(K_integrand_derived - K_direct) == 0, \
    "Derived K_integrand does not match direct ∇Ni·ε·∇Nj computation!"
print("\nCross-check ∇Ni · ε · ∇Nj (direct) vs. derived from weak form : match ✓")

# Flexoelectric polarisation is now written explicitly in the Q-tensor basis.
xi_a, xi_b = symbols('xi_a xi_b', real=True)

q1x, q2x, q3x, q4x, q5x = symbols('q1x q2x q3x q4x q5x', real=True)
q1y, q2y, q3y, q4y, q5y = symbols('q1y q2y q3y q4y q5y', real=True)
q1z, q2z, q3z, q4z, q5z = symbols('q1z q2z q3z q4z q5z', real=True)

divQ = Matrix([
    q1x + q3y + q5z,
    q3x + q2y + q4z,
    q5x + q4y - q1z - q2z,
])

P_flexo = expand(xi_a * divQ + xi_b * (Q * divQ))

P_load_integrand = expand(grad_Ni.dot(P_flexo))

nx_s, ny_s, nz_s = symbols('n_x n_y n_z', real=True)
n_surface = Matrix([nx_s, ny_s, nz_s])
P_surface_flux = expand(P_flexo.dot(n_surface))

print("\nFlexoelectric polarisation from the Q-tensor:")
print("  P_i = ξ_a Q_ij,j + ξ_b Q_ij Q_jk,k")
pprint(P_flexo)

print("\nFlexoelectric volume source integrand  ∇Ni · P =")
pprint(P_load_integrand)

print("\nFlexoelectric boundary flux integrand  P·n =")
pprint(P_surface_flux)

# ==============================================================================
# NEUMANN BOUNDARY STIFFNESS – LC / vacuum interface
# ==============================================================================
#
# The Neumann boundaries are the external surfaces of the LC mesh that face
# the surrounding unmodelled medium (vacuum, ε_vac = I).
#
# In qlc3d the surface normals at these boundaries point INWARD (into the LC
# element); call this n_in.  The outward-pointing normal used in standard
# integration-by-parts is n_out = -n_in.
#
# Two contributions combine at each Neumann surface Γ_N:
#
#  (a) IBP surface term of the LC ε volume integral at Γ_N
#      Rewriting ∫_{Ω_LC} ∇·(w ε ∇V) dΩ via the divergence theorem:
#        ∮_{∂Ω_LC} w (ε ∇V)·n_out dS
#      At Γ_N with n_out = -n_in this contributes:
#        -∮_{Γ_N} w (ε ∇V)·n_in dS
#
#  (b) Unmeshed vacuum exterior contribution
#      The exterior satisfies ∇²V = 0 with zero normal flux at the far
#      boundary.  Its weak-form volume integral ∫_{Ω_ext} ∇w·∇V dΩ
#      reduces (by IBP and the two boundary conditions) to a surface
#      integral at Γ_N.  The outward normal of the exterior at Γ_N equals
#      n_in (it points away from the exterior, i.e. into the LC), so:
#        ∫_{Ω_ext} ∇w · ∇V dΩ = +∮_{Γ_N} w (∇V)·n_in dS
#
# Adding (a) + (b) and moving to the bilinear-form LHS gives the Neumann
# stiffness integrand:
#
#   K_N(i,j) = Ni · ((ε − I) · ∇Nj) · n_in
#
# Physical check: substituting back via IBP of the full volume integral
# shows the combined system enforces  ∂V/∂n_in = 0  at Γ_N, which is the
# zero-outward-flux (homogeneous Neumann) condition:
#   D·n_out = (ε ∇V)·n_out = 0.
#
# The minus sign in the code (lK -= ...) is the standard sign convention
# where lK stores the NEGATIVE of the stiffness contribution.

print("\n" + "=" * 72)
print("NEUMANN BOUNDARY STIFFNESS – LC/vacuum interface")
print("=" * 72)
print("\nn_in points INTO the LC element (qlc3d sign convention).")
print("n_out = -n_in  (standard outward normal used in IBP).")
print("\n(a) IBP surface term from ε volume at Γ_N : -∮ w (ε·∇V)·n_in dS")
print("(b) Unmeshed vacuum exterior (∇²V=0, zero far-BC): +∮ w (∇V)·n_in dS")
print("\nCombined Neumann stiffness integrand:")
print("  K_N(i,j) = Ni · ((ε − I) · ∇Nj) · n_in")
print("Code convention (lK -= ...): lK(i,j) -= mul * Ni · ((ε−I)·∇Nj)·n_in")

# Ni_s  = shape-function VALUE at the boundary quadrature point (not a gradient)
Ni_s = symbols('Ni', real=True)
eps_minus_I = eps_check - eye(3)              # uses exx_s..eyz_s symbols
grad_Nj_vec = Matrix([Njx, Njy, Njz])
n_in_vec    = Matrix([nx_b, ny_b, nz_b])      # inward normal (qlc3d convention)

K_neumann_integrand = expand(Ni_s * (eps_minus_I * grad_Nj_vec).dot(n_in_vec))
print("\nNeumann stiffness integrand  Ni · ((ε−I) · ∇Nj) · n_in =")
pprint(K_neumann_integrand)

# Verify the symbolic expression matches the explicit form found in the C++
# implementation (localKLNeumann):
#   Ni * (((exx-1)*Njx + exy*Njy + exz*Njz)*nx
#        + (exy*Njx + (eyy-1)*Njy + eyz*Njz)*ny
#        + (exz*Njx + eyz*Njy + (ezz-1)*Njz)*nz)
K_neumann_code_check = Ni_s * (
    ((exx_s - 1)*Njx + exy_s*Njy + exz_s*Njz) * nx_b
  + (exy_s*Njx + (eyy_s - 1)*Njy + eyz_s*Njz) * ny_b
  + (exz_s*Njx + eyz_s*Njy + (ezz_s - 1)*Njz) * nz_b
)
assert simplify(K_neumann_integrand - expand(K_neumann_code_check)) == 0, \
    "Neumann stiffness integrand does not match C++ code formula!"
print("\nCross-check vs. C++ code formula (localKLNeumann): match ✓")

# Rename generic ε symbols to the short names used in the C++ generation below.
exx, eyy, ezz, exy, exz, eyz = symbols('exx eyy ezz exy exz eyz', real=True)
K_integrand = K_integrand_derived.subs({
    exx_s: exx, eyy_s: eyy, ezz_s: ezz,
    exy_s: exy, exz_s: exz, eyz_s: eyz,
})

# ==============================================================================
# PART 6: C++ CODE GENERATION
# ==============================================================================

# Rename symbols to C++ variable names matching the qlc3d shape-function API
# (s.Nx(i) for ∂N_i/∂x, etc.)
Nix_c, Niy_c, Niz_c = symbols('Nix Niy Niz')
Njx_c, Njy_c, Njz_c = symbols('Njx Njy Njz')

K_integrand_cpp = K_integrand.subs({
    Nix: Nix_c, Niy: Niy_c, Niz: Niz_c,
    Njx: Njx_c, Njy: Njy_c, Njz: Njz_c,
})

P_load_integrand_cpp = P_load_integrand.subs({
    Nix: Nix_c, Niy: Niy_c, Niz: Niz_c,
})

P_surface_flux_cpp = P_surface_flux.subs({
    nx_s: nx_b, ny_s: ny_b, nz_s: nz_b,
})

# Neumann stiffness integrand with C++ variable names.
# n_x, n_y, n_z correspond to n.x(), n.y(), n.z() in localKLNeumann.
# Ni corresponds to shapes.N(i); Njx/Njy/Njz to shapes.Nx/Ny/Nz(j).
K_neumann_integrand_cpp = K_neumann_integrand.subs({
    exx_s: exx, eyy_s: eyy, ezz_s: ezz,
    exy_s: exy, exz_s: exz, eyz_s: eyz,
    Njx: Njx_c, Njy: Njy_c, Njz: Njz_c,
    # nx_b/ny_b/nz_b keep their 'n_x'/'n_y'/'n_z' symbol names (= n.x() etc.)
    # Ni_s keeps its 'Ni' symbol name
})

# Dielectric permittivity components in T-tensor form, renamed to C++ variables
t1_c, t2_c, t3_c, t4_c, t5_c = symbols('t.t1() t.t2() t.t3() t.t4() t.t5()')
S0_c, deleps_c, eps_perp_c = symbols('S0 deleps eper_lc')
xi_a_c, xi_b_c = symbols('efe efe2')

q1x_c, q2x_c, q3x_c, q4x_c, q5x_c = symbols('q1x q2x q3x q4x q5x')
q1y_c, q2y_c, q3y_c, q4y_c, q5y_c = symbols('q1y q2y q3y q4y q5y')
q1z_c, q2z_c, q3z_c, q4z_c, q5z_c = symbols('q1z q2z q3z q4z q5z')

eps_components_cpp = {
    'exx': exx_T.subs({t1: t1_c, t2: t2_c, S0: S0_c, deleps: deleps_c, eps_perp: eps_perp_c}),
    'eyy': eyy_T.subs({t1: t1_c, t2: t2_c, S0: S0_c, deleps: deleps_c, eps_perp: eps_perp_c}),
    'ezz': ezz_T.subs({t1: t1_c,            S0: S0_c, deleps: deleps_c, eps_perp: eps_perp_c}),
    'exy': exy_T.subs({t3: t3_c,            S0: S0_c, deleps: deleps_c}),
    'exz': exz_T.subs({t5: t5_c,            S0: S0_c, deleps: deleps_c}),
    'eyz': eyz_T.subs({t4: t4_c,            S0: S0_c, deleps: deleps_c}),
}

P_components_cpp = {
    'Px': P_flexo[0].subs({
        xi_a: xi_a_c, xi_b: xi_b_c,
        q1x: q1x_c, q2x: q2x_c, q3x: q3x_c, q4x: q4x_c, q5x: q5x_c,
        q1y: q1y_c, q2y: q2y_c, q3y: q3y_c, q4y: q4y_c, q5y: q5y_c,
        q1z: q1z_c, q2z: q2z_c, q3z: q3z_c, q4z: q4z_c, q5z: q5z_c,
    }),
    'Py': P_flexo[1].subs({
        xi_a: xi_a_c, xi_b: xi_b_c,
        q1x: q1x_c, q2x: q2x_c, q3x: q3x_c, q4x: q4x_c, q5x: q5x_c,
        q1y: q1y_c, q2y: q2y_c, q3y: q3y_c, q4y: q4y_c, q5y: q5y_c,
        q1z: q1z_c, q2z: q2z_c, q3z: q3z_c, q4z: q4z_c, q5z: q5z_c,
    }),
    'Pz': P_flexo[2].subs({
        xi_a: xi_a_c, xi_b: xi_b_c,
        q1x: q1x_c, q2x: q2x_c, q3x: q3x_c, q4x: q4x_c, q5x: q5x_c,
        q1y: q1y_c, q2y: q2y_c, q3y: q3y_c, q4y: q4y_c, q5y: q5y_c,
        q1z: q1z_c, q2z: q2z_c, q3z: q3z_c, q4z: q4z_c, q5z: q5z_c,
    }),
}

# ── Print the generated C++ ──────────────────────────────────────────────────
print("// -----------------------------------------------------------------------")
print("// Dielectric permittivity tensor components from T-tensor (LC elements)")
print("// Derived from: ε_ij = (ε_⊥ + Δε/3)δ_ij + (2Δε/3S₀) Q_ij")
print("// -----------------------------------------------------------------------")
for var, expr in eps_components_cpp.items():
    print(f"double {var} = {ccode(expr)};")

print()
print("// -----------------------------------------------------------------------")
print("// Flexoelectric polarisation and charge-density terms")
print("// P_i = ξ_a Q_ij,j + ξ_b Q_ij Q_jk,k")
print("// -----------------------------------------------------------------------")
for var, expr in P_components_cpp.items():
    print(f"double {var} = {ccode(expr)};")

print()
print("// -----------------------------------------------------------------------")
print("// Local element stiffness matrix integrand: K^e_ij += mul * (∇Nᵢ · ε · ∇Nⱼ)")
print("// mul = quadrature_weight * jacobian_determinant")
print("// lK(i,j) uses the sign convention lK = −K^e (residual assembly).")
print("// -----------------------------------------------------------------------")
print("for (unsigned int i = 0; i < nodesPerTet; i++) {")
print("  double Nix = s.Nx(i), Niy = s.Ny(i), Niz = s.Nz(i);")
print("  for (unsigned int j = 0; j < nodesPerTet; j++) {")
print("    double Njx = s.Nx(j), Njy = s.Ny(j), Njz = s.Nz(j);")
print(f"    lK(i, j) -= mul * ({ccode(K_integrand_cpp)});")
print("  }")
print("}")

print()
print("// -----------------------------------------------------------------------")
print("// Flexoelectric volume source and boundary flux terms")
print("// Volume term:  ∫_Ω ∇Nᵢ · P dΩ")
print("// Boundary term: ∮_∂Ω Nᵢ (P·n) dS")
print("// The sign on the assembled RHS depends on the residual convention.")
print("// -----------------------------------------------------------------------")
print("for (unsigned int i = 0; i < nodesPerTet; i++) {")
print("  double Nix = s.Nx(i), Niy = s.Ny(i), Niz = s.Nz(i);")
print(f"  lL(i) += mul * ({ccode(P_load_integrand_cpp)});")
print("  // For a boundary quadrature rule, multiply N(i) by the surface flux P·n.")
print(f"  // lL(i) += mul * Ni * ({ccode(P_surface_flux_cpp)});")
print("}")

print()
print("// -----------------------------------------------------------------------")
print("// Neumann boundary stiffness: LC/vacuum interface  (localKLNeumann)")
print("// Integrand: Nᵢ · ((ε − I) · ∇Nⱼ) · n_in")
print("// n_in  = surface normal pointing INTO the LC element (qlc3d convention).")
print("// mul   = quadrature_weight * triangle_jacobian_determinant")
print("//")
print("// Derivation: two contributions combine at Γ_N (LC–vacuum surface):")
print("//  (a) IBP surface term of the ε volume integral at Γ_N:")
print("//        -∮ Nᵢ (ε·∇Nⱼ)·n_in dS   (n_out = -n_in)")
print("//  (b) Unmeshed vacuum exterior (∇²V=0, zero far-BC):")
print("//        +∮ Nᵢ (∇Nⱼ·n_in) dS")
print("//  Combined → Nᵢ · ((ε−I) · ∇Nⱼ) · n_in")
print("//")
print("// Together with the volume term this enforces ∂V/∂n_out = 0 at Γ_N")
print("// (zero outward electric-displacement flux, homogeneous Neumann BC).")
print("// -----------------------------------------------------------------------")
print("for (unsigned int i = 0; i < nodesPerTet; i++) {")
print("  double Ni = shapes.N(i);")
print("  for (unsigned int j = 0; j < nodesPerTet; j++) {")
print("    double Njx = shapes.Nx(j), Njy = shapes.Ny(j), Njz = shapes.Nz(j);")
print(f"    lK(i, j) -= mul * ({ccode(K_neumann_integrand_cpp)});")
print("  }")
print("}")

