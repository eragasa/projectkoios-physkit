# Package `projectkoios.physkit.solidstate.tb.slater_koster`

This package owns reusable spinless orthogonal homopolar Slater–Koster spatial-basis,
bond, channel, two-center matrix, ordinary-cell assembly, linear parameterization, and
reciprocal-sample evaluation mechanics.

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
- Embed directed five-orbital representatives into an exact ordered primitive-cell
  basis and generate their Hermitian reverses.
- Accumulate and sort dense coefficients by exact ordinary-cell displacement.
- Construct one nondimensional hopping model for every explicit parameter-group and
  channel pair.
- Evaluate the compatible linear model family while retaining group, channel, and
  reciprocal-sample order.
- Require named arguments for public records and action inputs.

## Exclusions

This tranche does not own onsite terms, neighbor-shell selection, material selection,
parameter fitting, design-matrix flattening, residuals, sample selection or weights,
evidence, acceptance, spin, spin–orbit coupling, overlap matrices, heteropolar bonds, or
`f` orbitals.

## Module architecture

Every maintained source module has its own four-page current-state contract.

- [`basis.py`](basis/index.md): exact ordered centers and spatial coordinates.
- [`bond.py`](bond/index.md): directed ordinary-cell bonds and physical geometry.
- [`channels.py`](channels/index.md): seven nondimensional angular matrices.
- [`parameters.py`](parameters/index.md): explicit integral-sharing groups.
- [`sp3s_star.py`](sp3s_star/index.md): energy-valued two-center matrices.
- [`assembly.py`](assembly/index.md): dense ordinary-cell coefficient assembly.
- [`linear_parameterization.py`](linear_parameterization/index.md): nondimensional
  group/channel hopping models.
- [`linear_bloch.py`](linear_bloch/index.md): ordered reciprocal evaluation of the
  linear model family.

The package façade in `__init__.py` deliberately aggregates this cohesive specialization
without re-exporting it through broader parent packages.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/`
- Provenance: [`docs/provenance/slater-koster-source-mapping.md`](../../../../../../provenance/slater-koster-source-mapping.md)
