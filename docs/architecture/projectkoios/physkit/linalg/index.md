# Package `projectkoios.physkit.linalg`

## Responsibility

This structural package organizes generic finite-matrix verification, mathematical
matrix operations, and eigensystem contracts. It defines no numerical behavior and
exports no module-owned or subpackage-owned identities from `__init__.py`.

## Canonical ownership paths

- [`matrices`](matrices/index.md) owns raw-or-quantity matrix verification and explicit
  canonical dense conversion.
- [`matrix`](matrix/index.md) organizes mathematical matrix operation families.
- [`matrix.norm`](matrix/norm/index.md) owns quantity-aware matrix norms.
- [`eigen`](eigen/index.md) owns eigensystem and spectral-decomposition contracts.

Callers import public identities from the narrowest owning module or integration
package. Parent-package convenience re-exports are not part of the contract.

## Representation boundaries

`matrices` accepts raw NumPy dense arrays, raw SciPy CSR values, and canonical PhysKit
matrix quantities because verification owns raw-to-quantity canonicalization.

`matrix.norm` accepts only canonical real/complex dense/sparse matrix quantities. It
delegates unit-free arithmetic to exact Optimization owners and restores the exact input
unit.

`eigen` accepts its documented unit-aware represented matrices and retains correlated
eigenvalue, eigenvector, and residual evidence. Matrix verification, matrix norms, and
spectral decomposition remain separate ownership contracts.

## Package initializer

`projectkoios.physkit.linalg.__all__` is empty. The initializer performs no eager imports,
wrapper construction, aliasing, dispatch, or compatibility behavior.

## Local mapping

- Package initializer: `src/python/projectkoios/physkit/linalg/__init__.py`
- Tests: `tests/projectkoios/physkit/linalg/test__linalg_public_api.py`

## Navigation

- [Module `matrices`](matrices/index.md)
- [Package `matrix`](matrix/index.md)
- [Package `eigen`](eigen/index.md)
- [Parent package](../index.md)
