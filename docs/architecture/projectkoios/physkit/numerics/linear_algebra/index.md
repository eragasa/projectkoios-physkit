# Package `projectkoios.physkit.numerics.linear_algebra`

## Current responsibility

The package owns reusable unit-free linear-algebra calculations, sparse
separable-operator assembly, and numerical predicates that are independent of
lattice or model policy.

## Public contract

Callers supply plain numerical matrices and explicit tolerances where a
numerical predicate requires policy. Defining modules are the supported import
locations.

## Navigation

- [Module `matrix_properties`](matrix_properties/index.md)
- [Module `kronecker`](kronecker/index.md)
- [Module `volume`](volume/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/numerics/linear_algebra/`
- Tests: `tests/projectkoios/physkit/numerics/linear_algebra/`

## Evidence

Mapped tests verify sparse Kronecker-sum shape, ordering, dtype, and spectrum;
scale-invariant normalized-volume behavior; and dense/sparse Hermiticity
predicates with explicit tolerances. They establish numerical software behavior
only.
