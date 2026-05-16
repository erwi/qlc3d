# Equations

## Notation

This document uses Einstein summation convention with tensor notation. Indices range from 1 to 3, representing the x, y, and z coordinates. A comma in subscripts denotes partial differentiation with respect to the following index.


## Electric Potential

The Poisson equation:

$$
(\varepsilon_{ij} \phi_{,j})_{,i} = \rho
$$

where $\phi$ is the electric potential and $\rho$ is a charge density. This follows from the condition $D_{i,i} = 0$ (no free charges), where the displacement field is $D_i = \varepsilon_{ij} E_j + P_i$, with $E_i = -\phi_{,i}$. Expanding gives the equivalent divergence-free form:

$$
(\varepsilon_{ij} \phi_{,j} - P_i)_{,i} = 0
$$

so the effective bound charge density due to the flexoelectric polarisation is:

$$
\rho = P_{i,i}
$$

## LC

The total energy $F$ is given by 
$$
F = \int _\Omega f_D + f_B - f_E \, d\Omega + 
\int_ \Gamma f_S \, d \Gamma
$$

where $f_D$ is the elastic distortion energy density, $f_B$ is the bulk or thermotropic energy density, $f_E$ is the external (electric) field-induced energy density, and $f_S$ is the surface or anchoring energy density.

### Q-tensor
Written in limit of uniaxial order with with director $\mathbf{n}$
$$
Q_{ij}=\frac{S}{2}(3n_i n_j - \delta_{ij})
$$

Symmetry and tracelessness must be maintained by writing the tensor in a five dimensional subspace as:
$$
\mathbf{Q} = \sum_{i=1}^5q_i\mathbf{T}_i
$$ 

where 
$$
\mathbf{T}_1 = (3\hat{\mathbf{e}}_z \otimes \hat{\mathbf{e}}_z - \mathbf{I}) / \sqrt{6}
$$

$$
\mathbf{T}_2 = (\hat{\mathbf{e}}_x \otimes \hat{\mathbf{e}}_x
    - \hat{\mathbf{e}}_y \otimes \hat{\mathbf{e}}_y
) / \sqrt{2}
$$

$$
\mathbf{T}_3 = (\hat{\mathbf{e}}_x \otimes \hat{\mathbf{e}}_y
    + \hat{\mathbf{e}}_y \otimes \hat{\mathbf{e}}_x
) / \sqrt{2}
$$

$$
\mathbf{T}_4 = (\hat{\mathbf{e}}_y \otimes \hat{\mathbf{e}}_z
    + \hat{\mathbf{e}}_z \otimes \hat{\mathbf{e}}_y
) / \sqrt{2}
$$

$$
\mathbf{T}_5 = (\hat{\mathbf{e}}_x \otimes \hat{\mathbf{e}}_z
    + \hat{\mathbf{e}}_z \otimes \hat{\mathbf{e}}_x
) / \sqrt{2}
$$

where the $\hat{\mathbf{e}}$ are unit vectors in $x$, $y$, $z$ directions.


### Elastic Energy

$$
f_D = \frac{1}{2}(
L_1 Q_{ij,k} Q_{ij,k} + 
L_2 Q_{ij,j}Q_{ik,k} + 
L_3 Q_{ik,j}Q_{ij,k} + 
L_4 \sigma_{ijk} Q_{il} Q_{jl,k} + 
L_6 Q_{lk} Q_{ij,l} Q_{ij,k}
)
$$


### Thermotropic Energy

$$
f_B = \frac{A}{2} Q_{ij} Q_{ij} + \frac{B}{3} Q_{ij} Q_{jk} Q_{ki} + \frac{C}{4} (Q_{ij} Q_{ij})^2
$$

### Electric Field Energy
$$
f_E = \frac{1}{2} \varepsilon_0 \varepsilon_{ij} E_i E_j + P_i E_i
$$

where $E_i = -\phi_{,i}$ is the electric field and $P_i$ is the flexoelectric polarisation.

#### Dielectric anisotropy
$$
\varepsilon_{ij} = \varepsilon_{\perp} \delta_{ij} + \Delta\varepsilon \left(\frac{2}{3S_0} Q_{ij} + \frac{1}{3}\delta_{ij}\right)
$$

#### Flexoelectric Polarisation 

$$
P_i = \xi_a Q_{ij,j} + \xi_b Q_{ij}Q_{jk,k}
$$


### Surface Anchoring Energy
A Rapini-Papoular style anchoring energy density


$$
f_S = a_s Q_{ij} Q_{ij} +
W_1 v_{1i} Q_{ij} v_{1j} +
W_2 v_{2i} Q_{ij} v_{2j}
$$

where $\hat{v}_1$ and $\hat{v}_2$ are mutually orthogonal unit vectors