# `projectkoios.physkit.solidstate.bloch2d.checkerboard`

This package owns a bare electron in the prescribed square-cell potential

$$
V(\mathbf r)=V_0\cos(\mathbf b_1\cdot\mathbf r)
                  \cos(\mathbf b_2\cdot\mathbf r).
$$

`CheckerboardPotentialBloch2D` requires an orthogonal, equal-length primitive
basis. It constructs a reciprocal plane-wave basis with the square cutoff
$-M\leq m_1,m_2\leq M$. The represented Hamiltonian contains exact
free-electron kinetic energies in the selected PhysKit unit system and the four
checkerboard Fourier coefficients $V_{(\pm1,\pm1)}=V_0/4$.

`CheckerboardPotentialBloch2DHamiltonianConstructor` returns the immutable
complex Hermitian matrix and canonical reciprocal-mode ordering for one Bloch
vector. `CheckerboardPotentialBloch2DBandSolver` diagonalizes that
representation over an identified Bloch-vector array and returns the requested
ordered energies.

The package owns the idealized represented Hamiltonian and its numerical band
solve. It does not assert a real-material potential, select an electronic
filling, define density-of-states broadening, or own finite-temperature carrier
statistics.
