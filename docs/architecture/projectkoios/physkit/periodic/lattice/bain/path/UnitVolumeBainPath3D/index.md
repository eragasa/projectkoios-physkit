# Class `projectkoios.physkit.periodic.lattice.bain.path.UnitVolumeBainPath3D`

## Current contract

The stateless action evaluates one explicit unit-volume FCC-to-BCC Bain path.

- `evaluate(...)` accepts exactly one finite built-in float in `[0, 1]`.
- Parameter zero has axial ratio `sqrt(2)` and parameter one has axial ratio one.
- Axial-ratio interpolation is linear in the path parameter.
- Scale restoration preserves the represented primitive volume.
- Every call returns a new correlated `UnitVolumeBainPathPoint3D`.

## Mathematical definition

The path uses one body-centered-tetragonal scale convention.

- For axial ratio `r`, the transverse scale is `(2 / r)^(1 / 3)`.
- The longitudinal scale is `r` times the transverse scale.
- The primitive determinant is `transverse^2 * longitudinal / 2 = 1`.
- Basis vectors and their order are part of the public contract.
- No extrapolation outside the endpoint interval is defined.

## Evidence limits

Evaluation establishes represented geometry only.

- Dense sampled tests verify the binary64 volume invariant.
- Endpoint labels identify crystallographic lattices under the declared basis
  convention.
- The action does not model kinetics, barriers, stress, temperature, or a material
  transformation mechanism.
- Path sampling and scientific acceptance belong to the consumer.

## Local mapping

Source and tests remain adjacent to the module owner.

- Code: `src/python/projectkoios/physkit/periodic/lattice/bain/path.py::UnitVolumeBainPath3D`
- Tests: `tests/projectkoios/physkit/periodic/lattice/bain/path/test__UnitVolumeBainPath3D__evaluate.py`
- [Parent module](../index.md)
