# Periodic-potential numeric contract

## Coordinates and mode identifiers

Real-space evaluation accepts a nonempty `MatrixQuantity` with shape
`(sample_count, dimension)` and a `Unitless` unit. Coordinates are fractional
unit-cell coordinates. Values outside `[0, 1)` are accepted because the
represented functions are periodic.

Fourier evaluation accepts a nonempty immutable tuple of integer mode tuples.
Every mode has exactly `dimension` built-in integer components. Integer indices
avoid floating-point tolerance decisions at the potential boundary.

## Quantities and storage

Potential amplitudes are finite `ScalarQuantity` values. A physical amplitude
must have energy dimensionality; an explicitly nondimensional model uses
`Unitless`. Sampled values are returned as `VectorQuantity` and coefficients as
`ComplexVectorQuantity` in the unchanged amplitude unit.

The sampled converter returns `SparseMatrixQuantity`; the Fourier converter
returns `ComplexSparseMatrixQuantity`. Both use canonical immutable CSR storage.
Zero entries are eliminated. Conversion never constructs a Hamiltonian and
never densifies an existing sparse operator.

## Fourier operator construction

For ordered reciprocal modes $m_r$ and $m_c$, the represented multiplication
operator is

$$
[V]_{rc}=V_{m_r-m_c}.
$$

A potential with declared finite Fourier support is assembled in
$O(Ns)$ candidate operations for $N$ represented modes and support size $s$.
The cosine, checkerboard, and three-dimensional cosine-product models have
support sizes 2, 4, and 8 respectively. A future potential with unknown or
infinite support may use the explicit dense-coefficient fallback before the
result is canonicalized as CSR.

## Floating-point behavior

Real-space values use NumPy binary64 trigonometric evaluation. Fourier support
selection is exact in integer-index space; it does not use `isclose` against
physical reciprocal vectors. Arbitrary phases may produce complex coefficients,
but a real potential produces conjugate symmetry $V_{-m}=V_m^*$ up to ordinary
binary64 roundoff.

## Verification boundary

Tests verify known sampled values, periodic repetition, exact support, Fourier
coefficients, Hermitian operator assembly, units, and CSR sparsity. Those tests
establish represented numerical behavior, not material-specific validation or
convergence of any downstream Hamiltonian discretization.

## Navigation

- [Package index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientific.md)
