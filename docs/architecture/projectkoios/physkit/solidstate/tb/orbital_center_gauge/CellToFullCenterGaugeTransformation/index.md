# Class `CellToFullCenterGaugeTransformation`

## Current contract

Correlates one exact gauge definition with its immutable nondimensional complex diagonal matrix.

The public boundary accepts a validated transformation definition and one square `ComplexMatrixQuantity`. The matrix shape must equal `(basis_size, basis_size)` and its unit must be the nondimensional unit-system view.

The central represented relation is

$$
D(q)=\operatorname{diag}\!\left(e^{+2\pi i q\cdot\tau_0},\ldots,e^{+2\pi i q\cdot\tau_{N_b-1}}\right).
$$

## Ownership boundary

The result retains a basis-change representation. It does not claim that a matrix subsequently transformed with it is scientifically admissible. It does not apply the transformation, diagonalize operators, or verify spectral equivalence.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/orbital_center_gauge.py::CellToFullCenterGaugeTransformation`
- Tests: `tests/projectkoios/physkit/solidstate/tb/orbital_center_gauge/`
- [Owning module](../index.md)
