# Module `projectkoios.physkit.solidstate.tb.orbital_center_gauge`

This module owns the generic diagonal transformation between ordinary-cell Bloch
coefficients and a representation carrying explicit fractional orbital centers.

For ordered centers $\tau_a$ and fractional reciprocal coordinate $q$,

$$
D_{aa}(q)=\exp(+2\pi i\,q\cdot\tau_a),\qquad
H_{\mathrm{full}}(q)=D(q)^\dagger H_{\mathrm{cell}}(q)D(q).
$$

The result retains the exact ordered center matrix, wave-vector object, and immutable
nondimensional complex transformation matrix. Public records and constructor inputs
require named arguments.

## Exclusions

The module does not choose orbital centers, infer a basis definition, assemble hopping
coefficients, define reciprocal sewing policy, or judge scientific admissibility.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/orbital_center_gauge.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/orbital_center_gauge/`
