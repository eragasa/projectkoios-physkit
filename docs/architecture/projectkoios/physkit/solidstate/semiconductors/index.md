# `projectkoios.physkit.solidstate.semiconductors`

This package owns semiconductor-specific reduced physical models.

- [`mobility`](mobility/index.md) owns empirical lattice and ionized-impurity
  mobility contributions and a relaxation-time trajectory sampler.
- [`junction`](junction/index.md) owns a nondegenerate abrupt p–n-junction
  depletion model and its correlated equilibrium results.

Dimension-specific effective-mass ownership is explicit:

- `effmass_1d` is reserved for a future one-dimensional contract and currently
  provides no maintained model.
- [`effmass_2d`](effmass_2d/index.md) owns the maintained two-dimensional
  periodic primitive-cell methodological model.
- [`effmass_3d`](effmass_3d/index.md) owns the maintained three-dimensional
  periodic primitive-cell reference model.

The reserved `effmass_1d` package establishes an authorized naming boundary
only. It does not claim an implemented or verified capability.
