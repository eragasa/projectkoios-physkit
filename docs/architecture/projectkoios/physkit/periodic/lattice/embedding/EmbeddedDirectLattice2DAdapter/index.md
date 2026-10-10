# Class `EmbeddedDirectLattice2DAdapter`

A stateless `DataObjectActionizer` constructing an ordered ambient orthonormal
frame, tangent projector, induced metric and inverse, ambient dual reciprocal
basis, and canonical intrinsic `2×2` basis. The intrinsic basis is passed to the
unchanged nominal `DirectLattice2D` and `DirectLatticeMetric` contracts.

The action performs geometry adaptation only. It owns no finite-difference grid,
Bloch twist, Hamiltonian coefficient, atomistic slab, or material interpretation.

- [Module index](../index.md)
- [Numerical equations](../numerical/index.md)
- Code: `src/python/projectkoios/physkit/periodic/lattice/embedding/__init__.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/embedding/test__EmbeddedDirectLattice2DAdapter__action.py`
