# Package `projectkoios.physkit.solidstate.tb.slater_koster`

This package owns reusable spinless orthogonal homopolar Slater–Koster spatial-basis,
bond, channel, and two-center matrix mechanics.

## Current responsibility

- Represent ordered real-cubic spatial orbitals with explicit center and radial-channel
  labels.
- Admit the exact homopolar `(s, p_x, p_y, p_z, s*)` specialization, with `s*`
  represented as radial channel one of the `s` component.
- Keep ordinary-cell displacement separate from unit-bearing Cartesian bond geometry.
- Represent seven explicit two-center energy channels in a stable order.
- Construct the real five-by-five hopping matrix as a linear combination of
  nondimensional angular channel matrices.
- Preserve `T(-d) = T(d)^T` for real homopolar integrals.
- Require named arguments for public records and action inputs.

## Exclusions

This tranche does not own onsite terms, neighbor-shell selection, material selection,
fitting, residuals, sample weights, evidence, acceptance, spin, spin–orbit coupling,
overlap matrices, heteropolar bonds, or `f` orbitals.

## Public capability groups

- Basis: `basis.py`
- Bond geometry: `bond.py`
- Angular channels: `channels.py`
- Explicit parameter groups: `parameters.py`
- Two-center matrix construction: `sp3s_star.py`

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/`
- Provenance: [`docs/provenance/slater-koster-source-mapping.md`](../../../../../../provenance/slater-koster-source-mapping.md)
