# Package `projectkoios.physkit.potentials.periodic2d`

This package owns `CheckerboardPotential2D`, the fractional-coordinate field

$$
V(u_1,u_2)=A\prod_{j=1}^{2}\cos(2\pi h_j u_j+\phi_j).
$$

The fundamental zero-phase model has four equal coefficients $A/4$ at
$(\pm1,\pm1)$. The maintained checkerboard Bloch model obtains its potential
matrix through the reciprocal-mode operator converter rather than embedding
those couplings in its Hamiltonian constructor.

- Code:
  `src/python/projectkoios/physkit/potentials/periodic2d/checkerboard.py`
- Tests:
  `tests/projectkoios/physkit/potentials/periodic2d/checkerboard/`
- [Shared periodic-potential contracts](../periodic/index.md)
