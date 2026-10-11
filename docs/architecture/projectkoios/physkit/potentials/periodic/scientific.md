# Periodic-potential scientific contract

## Fractional-coordinate convention

For dimension $d$, fractional coordinates are
$u=(u_1,\ldots,u_d)$. A lattice supplied by another package may map them to a
physical position through $r=A u$. Periodicity is expressed without selecting a
particular lattice geometry:

$$
V(u+n)=V(u),\qquad n\in\mathbb Z^d.
$$

The Fourier convention is

$$
V(u)=\sum_{m\in\mathbb Z^d}V_m e^{i2\pi m\cdot u}.
$$

For a real scalar potential, the coefficients obey
$V_{-m}=V_m^*$.

## Cosine-product family

The implemented family is

$$
V(u)=A\prod_{j=1}^{d}
\cos\!\left(2\pi h_j u_j+\phi_j\right),
$$

where every $h_j$ is a positive integer. Its exact finite support contains the
$2^d$ modes

$$
m=(s_1h_1,\ldots,s_dh_d),\qquad s_j\in\{-1,+1\},
$$

with coefficients

$$
V_m=\frac{A}{2^d}
\exp\!\left(i\sum_j s_j\phi_j\right).
$$

`CosinePotential1D` is the $d=1$ case. `CheckerboardPotential2D` is the $d=2$
case and reproduces the four $A/4$ couplings used by the checkerboard Bloch
model when phases vanish and both harmonics are one.
`CosineProductPotential3D` supplies the corresponding eight-coupling field; it
does not by itself define a crystal, particle, Hamiltonian, or material.

## Scope and interpretation

These models are idealized mathematical potentials. An energy-valued amplitude
does not establish empirical fidelity to a material. Unitless amplitudes are
valid only inside an explicitly stated nondimensional convention owned by the
consumer.

The potential library makes no assertion about electron filling, spin,
interactions, density of states, temperature, chemical potential, or observable
response. Those require separate physical models and validation.

## Validation status

The analytical Fourier identities and selected sampled values are covered by
software and numerical tests. No claim is made here of experimental validation,
uncertainty quantification, or suitability for a named material.

## Navigation

- [Package index](index.md)
- [Schematics](schematics.md)
- [Numeric contract](numeric.md)
