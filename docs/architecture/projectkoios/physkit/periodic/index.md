# Package `projectkoios.physkit.periodic`

## Current responsibility

The package owns reusable periodic-system representations and mechanics.

- Physical geometry includes direct and reciprocal lattices, reciprocal bases and
  modes, and Brillouin-zone paths.
- Represented mechanics include finite periodic integer domains, boundary twists, and
  nonorthogonal-cell periodic finite differences.
- Immutable unit cells retain deterministic serialization contracts.

## Public contract

Physical geometry and represented integer domains retain distinct identities.

- Physical direct-lattice vectors are not finite periodic integer indices.
- `FinitePeriodicDomain` carries no physical vectors, units, or unit-cell ownership.
- Reciprocal-mode inventories preserve explicit caller-supplied ordering.

## Navigation

Architecture pages map cohesive periodic ownership.

- [Module `projectkoios.physkit.periodic.basis`](basis/index.md)
- [Package `projectkoios.physkit.periodic.lattice`](lattice/index.md)
- [Package `projectkoios.physkit.periodic.finite_difference`](finite_difference/index.md)
- [Package `projectkoios.physkit.periodic.unit_cell`](unit_cell/index.md)
- [Package `projectkoios.physkit.periodic.pbc1d`](pbc1d/index.md)
- [Parent package](../index.md)

## Local mapping

Maintained source and tests mirror periodic ownership.

- Code: `src/python/projectkoios/physkit/periodic/__init__.py`
- Tests: `tests/projectkoios/physkit/periodic/`

## Evidence

Mapped package tests establish implementation conformance.

- Lattice tests bind the supported direct and reciprocal identities.
- Basis and mode tests bind canonical reciprocal-vector composition.
- Exact represented agreement does not establish physical alignment or scientific
  validation.
