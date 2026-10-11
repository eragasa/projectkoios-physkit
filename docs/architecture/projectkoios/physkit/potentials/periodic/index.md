# Package `projectkoios.physkit.potentials.periodic`

This package owns representation-independent periodic scalar potentials on
fractional unit-cell coordinates and explicit conversions to represented
multiplication operators.

`BasePeriodicPotential` defines three public facts: coordinate dimension,
energy unit, and exact finite Fourier support when known. Its evaluation façade
returns typed results for either sampled values or reciprocal-mode Fourier
coefficients. Fractional coordinates are always unitless; potential values and
coefficients retain the amplitude unit, which may be a physical energy unit or
`Unitless` for an explicitly nondimensional model.

Two converters keep field models separate from numerical representations:

- `PeriodicPotentialToSampledDiagonalOperatorConverter` constructs a canonical
  real CSR multiplication operator on supplied fractional sample coordinates.
- `PeriodicPotentialToFourierOperatorConverter` constructs the canonical complex
  CSR convolution operator in an ordered reciprocal-mode basis.

The initial concrete families are:

- [`CosinePotential1D`](../periodic1d/index.md),
- [`CheckerboardPotential2D`](../periodic2d/index.md), and
- [`CosineProductPotential3D`](../periodic3d/index.md).

The shared `periodic` package owns contracts and converters. Concrete dimensions
use the sibling package names `periodic1d`, `periodic2d`, and `periodic3d`.
These are implemented model packages, not empty dimensional placeholders.

The package does not own lattice vectors, the mapping from fractional to
physical coordinates, kinetic energy, boundary-condition discretization,
Hamiltonian composition, or spectral solving.

## Local mapping

- Base contracts:
  `src/python/projectkoios/physkit/potentials/periodic/base.py`
- Operator conversions:
  `src/python/projectkoios/physkit/potentials/periodic/operators.py`
- Tests: `tests/projectkoios/physkit/potentials/periodic/`

## Navigation

- [Schematics](schematics.md)
- [Numeric contract](numeric.md)
- [Scientific contract](scientific.md)
