# Package `projectkoios.physkit`

## Current responsibility

`projectkoios.physkit` owns reusable computational-physics models,
mathematical representations, numerical methods, unit-aware quantities, and
pedagogical helpers. It does not own campaign orchestration, application
workflow state, optimization engines, or material-specific acceptance policy.

## Public contract

The distribution `projectkoios-physkit` contributes the
`projectkoios.physkit` portion of the PEP 420 `projectkoios` namespace. The
repository does not install a top-level `projectkoios/__init__.py`. Defining
modules are the default supported import locations; parent packages do not
indiscriminately aggregate implementation classes.

`src/python/projectkoios/physkit/` is the canonical and only packaged source
root. An optional `src/physkit/` tree may exist, but every Python source path in
it must be mirrored at the same relative path in the canonical tree. The legacy
tree does not define a supported top-level import or a second implementation
owner; caches and non-source workspace artifacts do not participate in the
mirror contract.

## Navigation

- [Package `projectkoios.physkit.core`](core/index.md)
- [Package `projectkoios.physkit.deposition`](deposition/index.md)
- [Package `projectkoios.physkit.math`](math/index.md)
- [Package `projectkoios.physkit.materials`](materials/index.md)
- [Package `projectkoios.physkit.mechanics`](mechanics/index.md)
- [Package `projectkoios.physkit.numerics`](numerics/index.md)
- [Package `projectkoios.physkit.periodic`](periodic/index.md)
- [Package `projectkoios.physkit.plasmas`](plasmas/index.md)
- [Package `projectkoios.physkit.qm`](qm/index.md)
- [Package `projectkoios.physkit.solidstate`](solidstate/index.md)
- [Package `projectkoios.physkit.thermal`](thermal/index.md)
- [Package `projectkoios.physkit.waves`](waves/index.md)
- [Finite-periodic extraction provenance](../../../provenance/finite-periodic-lattice-source-mapping.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/__init__.py`
- Packaging: `pyproject.toml`
- Tests: `tests/projectkoios/physkit/`

## Evidence

Implementation conformance is supported by the complete test suite and wheel
import checks. This result does not establish numerical or scientific
validation for every model in the package.
