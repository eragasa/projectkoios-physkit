# `projectkoios.physkit.solidstate.bloch2d.checkerboard`

This package owns a bare electron in the prescribed square-cell potential

$$
V(\mathbf r)=V_0\cos(\mathbf b_1\cdot\mathbf r)
                  \cos(\mathbf b_2\cdot\mathbf r).
$$

`CheckerboardPotentialBloch2D` requires a [`BlochCell2D`](../cell/index.md)
with an orthogonal, equal-length direct basis. The supplied square is a
displayed conventional cell: this product potential has an additional
half-diagonal translation symmetry. The model supports `UnitSystem.METAL` and
the explicitly unitless `UnitSystem.NONDIMENSIONAL` convention; it rejects
`UnitSystem.SI`.

The model constructs a reciprocal plane-wave basis with the square cutoff
$-M\leq m_1,m_2\leq M$. The represented Hamiltonian contains free-electron
kinetic energies and the four checkerboard Fourier coefficients
$V_{(\pm1,\pm1)}=V_0/4$. In the nondimensional convention, $\hbar=1$ and the
bare-electron mass is $m_e=1$.

`CheckerboardPotentialBloch2DHamiltonianConstructor` returns the immutable
canonical complex CSR Hamiltonian and reciprocal-mode ordering for one Bloch
vector. Sparse storage is the default representation. For a partial spectrum,
`CheckerboardPotentialBloch2DBandSolver` detects the sparse operator and uses
the sparse Hermitian eigensolver. A request for the complete or all-but-one
spectrum crosses an explicit dense boundary because ARPACK requires a smaller
partial selection.

The square-cell test is invariant under uniform scaling of the represented
basis. It classifies shape from a normalized Gram matrix rather than applying a
dimensional absolute tolerance to squared lengths.

The package owns the idealized represented Hamiltonian and its numerical band
solve. It does not assert a real-material potential, select an electronic
filling, define density-of-states broadening, or own finite-temperature carrier
statistics.

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/bloch2d/checkerboard/model.py`
- Tests: `tests/projectkoios/physkit/solidstate/bloch2d/checkerboard/`
