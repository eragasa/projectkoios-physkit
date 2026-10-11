# Module `projectkoios.physkit.periodic.lattice.lattice2d`

## Direct-lattice contract

`DirectLattice2D` owns two finite intrinsic two-dimensional primitive vectors as
immutable matrix columns. Linear independence is evaluated from the
dimensionless determinant of column-normalized vectors with absolute tolerance
`1e-8`. The decision is therefore invariant under a common or independent
nonzero rescaling of the primitive vectors and admits well-conditioned bases
represented in small physical scales such as meters.

The class represents a square intrinsic `2×2` basis. An embedding of a 2D cell
in three Cartesian coordinates is not part of this contract.

## Navigation

- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/lattice2d.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/lattice2d/test__DirectLattice2D__*.py`
