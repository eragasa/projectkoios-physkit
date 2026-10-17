# Package `projectkoios.physkit`

## Current responsibility

The package owns reusable computational-physics representations and models.

- Owned surfaces include mathematical representations, physical units, numerical
  methods, crystallography, and pedagogical helpers.
- Optimization engines and generic backend arithmetic remain with their numerical
  owners.
- Campaign orchestration, application workflow state, and material-specific acceptance
  remain outside PhysKit.

## Public contract

The distribution contributes one portion of the implicit `projectkoios` namespace.

- `projectkoios-physkit` must not install a top-level `projectkoios/__init__.py`.
- Defining modules are supported import locations; package façades aggregate only
  deliberate public identities.
- `src/python/projectkoios/physkit/` is the canonical packaged source root.
- Any optional `src/physkit/` Python source must mirror the same relative path under the
  canonical tree.
- Caches and workspace artifacts do not participate in the source-mirror contract.

## Navigation

Architecture pages follow package-level ownership.

- [Package `projectkoios.physkit.core`](core/index.md)
- [Package `projectkoios.physkit.deposition`](deposition/index.md)
- [Package `projectkoios.physkit.dos`](dos/index.md)
- [Package `projectkoios.physkit.dos1d`](dos1d/index.md)
- [Package `projectkoios.physkit.dos2d`](dos2d/index.md)
- [Package `projectkoios.physkit.dos3d`](dos3d/index.md)
- [Package `projectkoios.physkit.math`](math/index.md)
- [Package `projectkoios.physkit.materials`](materials/index.md)
- [Package `projectkoios.physkit.mechanics`](mechanics/index.md)
- [Package `projectkoios.physkit.numerics`](numerics/index.md)
- [Package `projectkoios.physkit.periodic`](periodic/index.md)
- [Package `projectkoios.physkit.potentials`](potentials/index.md)
- [Package `projectkoios.physkit.plasmas`](plasmas/index.md)
- [Package `projectkoios.physkit.qm`](qm/index.md)
- [Package `projectkoios.physkit.solidstate`](solidstate/index.md)
- [Package `projectkoios.physkit.thermal`](thermal/index.md)
- [Package `projectkoios.physkit.units`](units/index.md)
- [Package `projectkoios.physkit.waves`](waves/index.md)
- [Finite-periodic extraction provenance](../../../provenance/finite-periodic-lattice-source-mapping.md)

## Local mapping

Maintained source, packaging, and verification remain explicit.

- Code: `src/python/projectkoios/physkit/__init__.py`
- Packaging: `pyproject.toml`
- Tests: `tests/projectkoios/physkit/`

## Evidence

Mapped tests and wheel checks establish implementation conformance.

- Complete-suite and isolated-install checks exercise the packaged namespace.
- Passing checks do not establish numerical or scientific validation for every model.
