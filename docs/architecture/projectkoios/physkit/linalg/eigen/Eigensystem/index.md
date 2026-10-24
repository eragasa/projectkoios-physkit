# `Eigensystem`

## Responsibility

`Eigensystem` is the public immutable evidence record returned by represented
self-adjoint spectral decompositions. It retains the exact matrix, ordered real
eigenvalues, unitless column eigenvectors, per-column residual norms, maximum
residual, and derived maximum/selected rank properties.

`max_rank` is exactly the square matrix dimension N and is never a
tolerance-dependent numerical matrix rank. `n_rank` is the retained eigenpair
count. `is_complete` reports whether `n_rank == max_rank`.

## Structurally admitted representations

The record accepts real or complex dense/CSR matrix quantities and real or
complex dense eigenvector quantities. Every matrix/eigenvector family
combination is structurally admitted: a real self-adjoint matrix may use
complex-phased eigenvectors, while some complex-Hermitian matrices may have
real eigenvectors. Eigenvalues and residual norms remain real quantities, and
eigenvectors must be Unitless.

Construction validates nominal types, a nonempty square exactly Hermitian
matrix, dimensions, units, nondecreasing eigenvalue order, conjugate
orthonormality, nonnegative residual norms, and exact maximum-residual
correlation. It deliberately does **not** verify
`matrix @ eigenvectors == eigenvectors * eigenvalues`. Residual quantities are
retained supplied evidence. `SpectralDecomposition` owns their computation;
direct record construction does not establish that the supplied vectors are
eigenvectors or that the supplied residual evidence was independently derived.

The generic record does not claim that a complex decomposition algorithm is
implemented. The public façade continues to report exact complex-Hermitian
matrices as unsupported.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/eigen/eigensystem.py::Eigensystem`
- Constructor tests: `tests/projectkoios/physkit/linalg/eigen/test__Eigensystem__init__.py`
- Nominal results ownership: `tests/projectkoios/physkit/core/results/test__ResultsObject__inheritance.py`

## Navigation

- [Package](../index.md)
