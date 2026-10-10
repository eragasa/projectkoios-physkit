# Class `EmbeddedDirectLattice2D`

A frozen, slotted `DataObject` owning one immutable `MatrixQuantity` with shape
`(3,2)`. Columns are ordered ambient primitive vectors. The unit is physical
length or `Unitless`; column-rank validation uses normalized directions and is
therefore scale invariant.

The record preserves ambient orientation and does not construct a third slab
vector.

- [Module index](../index.md)
- Code: `src/python/projectkoios/physkit/periodic/lattice/embedding/__init__.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/embedding/test__EmbeddedDirectLattice2D__init.py`
