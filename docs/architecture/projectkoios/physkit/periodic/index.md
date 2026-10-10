# Package `projectkoios.physkit.periodic`

## Current responsibility

The package owns reusable periodic-system representations, including physical
direct and reciprocal lattices, reciprocal bases and modes, Brillouin-zone
paths, finite periodic integer domains, boundary-twist mechanics, and immutable
unit cells with deterministic serialization.

## Public contract

Physical direct-lattice vectors remain distinct from finite periodic integer
index domains. `FinitePeriodicDomain` represents only a one-, two-, or
three-dimensional integer quotient domain and carries no physical vectors,
units, or unit-cell ownership.

## Navigation

- [Package `projectkoios.physkit.periodic.lattice`](lattice/index.md)
- [Package `projectkoios.physkit.periodic.unit_cell`](unit_cell/index.md)
- [Package `projectkoios.physkit.periodic.pbc1d`](pbc1d/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/__init__.py`
- Tests: `tests/projectkoios/physkit/periodic/`

## Evidence

Implementation conformance is supported by the mapped periodic-package tests.
Exact represented agreement does not establish physical alignment or scientific
validation.
