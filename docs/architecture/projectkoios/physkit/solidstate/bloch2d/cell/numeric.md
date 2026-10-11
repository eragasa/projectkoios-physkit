# `BlochCell2D` numeric contract

## Direct basis and numerical rank

Let the ordered direct vectors be the columns of

$$
A=[\mathbf a_1\;\mathbf a_2]\in\mathbb R^{2\times2}.
$$

The cell uses the dimensionless normalized column area

$$
\nu(A)=\left|\det\left[
\frac{\mathbf a_1}{\lVert\mathbf a_1\rVert_2}\;
\frac{\mathbf a_2}{\lVert\mathbf a_2\rVert_2}
\right]\right|.
$$

Each column is pre-scaled by its largest absolute component before norm
evaluation. This avoids unnecessary overflow or underflow. The maintained rank
criterion is

$$
\nu(A)>10^{-8}.
$$

Consequently, multiplying either column by a nonzero scalar does not change the
rank decision. This is a binary64 numerical-independence criterion, not an exact
symbolic determinant test.

## Derived geometry

The positive displayed-cell area is

$$
\Omega=|\det A|.
$$

The reciprocal basis is

$$
B=2\pi A^{-\mathsf T},
$$

and satisfies

$$
A^{\mathsf T}B=2\pi I.
$$

The covariant direct metric is

$$
G=A^{\mathsf T}A.
$$

For `UnitSystem.METAL`, `A`, `Omega`, `B`, and `G` use ångström,
ångström-squared, inverse-ångström, and ångström-squared units respectively.
For `UnitSystem.NONDIMENSIONAL`, all four quantities use `Unitless`.

## Scale-invariant square classification

The checkerboard model imposes a stronger shape constraint. To keep that
classification independent of a uniform basis scale, it evaluates

$$
\widehat A=\frac{A}{\max_{ij}|A_{ij}|},\qquad
\widehat G=\widehat A^{\mathsf T}\widehat A,
$$

then normalizes once more by the largest diagonal entry,

$$
C=\frac{\widehat G}{\max_i \widehat G_{ii}}.
$$

The cell is accepted as square when `C` agrees with the identity with zero
relative tolerance and absolute tolerance `2e-13`. This simultaneously checks
orthogonality and equal vector lengths without applying an absolute tolerance
to dimensional squared lengths.

## Hamiltonian conventions

For the metal convention, the plane-wave kinetic term is evaluated in the
selected numerical scale and converted to electron-volts:

$$
E_{\mathrm{kin}}(\mathbf q)=
\frac{\hbar^2}{2m_e}|\mathbf q|^2.
$$

For the nondimensional convention,

$$
\hbar=1,\qquad m_e=1,
$$

so

$$
E_{\mathrm{kin}}(\mathbf q)=\frac12|\mathbf q|^2.
$$

The direct basis, reciprocal basis, wave vectors, potential amplitude, matrix,
and returned energies must all use the same selected convention.

## Sparse Hamiltonian representation and solver dispatch

For reciprocal cutoff $-M\leq m_1,m_2\leq M$, the basis dimension is

$$
N=(2M+1)^2.
$$

The kinetic term is diagonal. Each reciprocal mode couples to at most four
other modes through transfers $(\pm1,\pm1)$. The constructor therefore assembles
only the diagonal and present four-neighbor transfers and canonicalizes them as
an immutable complex CSR quantity. Stored entries scale as $O(N)$ rather than
the $O(N^2)$ inventory of a dense matrix.

The band solver inspects the represented operator. For a sparse matrix and a
requested count $k<N-1$, it uses `scipy.sparse.linalg.eigsh` with the
smallest-algebraic selection and a deterministic initial vector, then sorts the
returned eigenvalues. When $k\geq N-1$, it materializes an explicit dense copy
and uses the complete Hermitian solver because ARPACK does not support that
selection range. Sparse storage remains the default representation in either
case.

## Numerical evidence

Mapped tests verify scale-independent rank rejection, direct/reciprocal duality,
area and metric identities, rotated-square acceptance, nonsquare rejection over
multiple scales, canonical sparse storage, automatic sparse partial-spectrum
dispatch, metal Fourier couplings, and the unitless free-particle kinetic
diagonal. These checks do not bound conditioning for every accepted basis or
establish convergence with reciprocal cutoff.

## Navigation

- [Module index](index.md)
- [Schematics](schematics.md)
- [Scientific page](scientific.md)
