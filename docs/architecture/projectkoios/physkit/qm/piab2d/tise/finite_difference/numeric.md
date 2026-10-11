# Oblique PIAB2D finite-difference numerical contract

## Grid and indexing

For full-grid point counts `(N0, N1)`, fractional spacings are

$$
h_i=\frac{1}{N_i-1}.
$$

Boundary points at zero and one carry the declared zero Dirichlet values and are
not unknowns. The represented interior shape is `(N0 - 2, N1 - 2)`. The second
fractional index varies fastest under C-order flattening.

## Sparse operator

For $g^{-1}=(A^{\mathsf T}A)^{-1}$, the geometric operator is

$$
\Delta_h=g^{11}D_{11}+g^{22}D_{22}+2g^{12}D_{12}.
$$

Pure terms use centered second differences. The mixed term is the Kronecker
product of centered first-difference matrices. Missing edge targets multiply
zero boundary values. The Hamiltonian

$$
H_h=-\frac{\hbar^2}{2m}\Delta_h
$$

is retained as canonical real CSR storage. Physical models convert the compound
kinetic unit into the selected unit-system energy; nondimensional models use
their declared numerical $\hbar$ and mass directly.

## Eigensolving and normalization

A sparse smallest-algebraic eigensolve is used when the requested count is less
than the represented dimension minus one. Near-complete spectra use an explicit
dense Hermitian solve. Returned energies are nondecreasing.

The constant cell Jacobian gives each interior unknown the quadrature weight

$$
w_A=|\det A|h_0h_1.
$$

Eigenfunction columns are rescaled so

$$
w_A\sum_p|\psi_{pj}|^2=1.
$$

## Verification evidence

Tests cover argument failures, canonical sparse Hermiticity, interior coordinate
bounds, area-weighted normalization, the exact discrete tensor-product sine
spectrum for an orthogonal rectangle, and second-order ground-energy refinement
against the rectangular continuum result. The maintained laboratory separately
demonstrates the mixed stencil and physical-parallelogram visualization.

These bounded tests do not prove accuracy for every skew angle, conditioning,
resolution, eigensolver configuration, or high-energy state.
