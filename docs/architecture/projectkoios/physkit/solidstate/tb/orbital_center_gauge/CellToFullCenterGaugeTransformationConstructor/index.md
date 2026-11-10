# Class `CellToFullCenterGaugeTransformationConstructor`

## Current contract

Constructs the positive-phase diagonal orbital-center gauge matrix from one validated definition.

The public boundary accepts one keyword-only `CellToFullCenterGaugeTransformationDefinition`. The action performs a rowwise center-wave-vector contraction, vectorized complex exponential, and diagonal embedding while preserving input order.

The central represented relation is

$$
D_{aa}(q)=e^{+2\pi i q\cdot\tau_a}.
$$

## Ownership boundary

The positive exponent is part of the representation convention. The constructor changes no physical model and performs no scientific comparison. It does not construct a Bloch Hamiltonian, apply left/right conjugation, or infer sewing matrices.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/orbital_center_gauge.py::CellToFullCenterGaugeTransformationConstructor`
- Tests: `tests/projectkoios/physkit/solidstate/tb/orbital_center_gauge/test__CellToFullCenterGaugeTransformationConstructor__execute.py`
- [Owning module](../index.md)
