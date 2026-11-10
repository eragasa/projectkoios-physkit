# Class `CellToFullCenterGaugeTransformationDefinition`

## Current contract

Retains the exact ordered fractional orbital centers and one fractional reciprocal coordinate that define a cell-to-full-center gauge transformation.

The public boundary accepts an immutable nondimensional center matrix of shape `(basis_size, 3)` and an immutable nondimensional wave vector of shape `(3,)`. Validation is nominal, shape-exact, and unit-exact. It performs no center sorting, basis inference, dimensional projection, or reciprocal-axis permutation.

The central represented relation is

$$
\theta_a=q\cdot\tau_a,\qquad a=0,\ldots,N_b-1.
$$

## Ownership boundary

The record defines representation data only. Retaining centers and a wave vector together does not establish that they belong to the same independently supplied basis. It does not construct phases, transform a matrix, choose centers, or define reciprocal sewing.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/orbital_center_gauge.py::CellToFullCenterGaugeTransformationDefinition`
- Tests: `tests/projectkoios/physkit/solidstate/tb/orbital_center_gauge/`
- [Owning module](../index.md)
