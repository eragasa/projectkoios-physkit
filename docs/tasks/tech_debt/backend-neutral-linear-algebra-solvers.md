# Backend-neutral linear-algebra solver boundary

**Status:** Open technical debt

## Problem

PhysKit currently owns both dense and sparse represented operators, but several
model layers call backend APIs directly. Typical examples include
`numpy.linalg.eigh` or `numpy.linalg.eigvalsh` for dense arrays and
`scipy.sparse.linalg.eigsh` for sparse partial spectra. This makes solver
selection, ordering, fallback, convergence reporting, and densification policy
model-specific.

The checkerboard Bloch model now correctly defaults to a canonical sparse
Hamiltonian and detects that representation when selecting a partial sparse
solve. Its local dispatch is intentionally bounded, but it should eventually
consume a reusable linear-algebra-neutral contract rather than own NumPy/SciPy
branching itself.

## Desired ownership

Add a reusable numerical boundary under
`projectkoios.physkit.numerics.eigenproblems` with nominal request and result
owners for real-symmetric and complex-Hermitian operators. The public numerical
contract should describe:

- complete versus lowest-`k` spectrum requests;
- eigenvalues-only versus correlated eigenpair results;
- deterministic ascending ordering;
- real-symmetric and complex-Hermitian validation;
- dense and sparse represented inputs without structural protocols;
- explicit prohibition of accidental densification;
- supported complete-spectrum fallbacks;
- deterministic initialization where an iterative backend needs it; and
- convergence or backend failure surfaced without model-specific translation.

Backend adapters should implement that contract:

1. a NumPy dense implementation for complete or selected dense spectra;
2. a SciPy sparse implementation for proper partial spectra; and
3. an explicitly named dense fallback policy for sparse requests that cannot be
   completed by the maintained sparse algorithm.

Model code should choose the mathematical request. Representation-aware solver
dispatch should occur at the numerical boundary.

## Representation constraints

- `MatrixQuantity` and `ComplexMatrixQuantity` remain explicit dense owners.
- `SparseMatrixQuantity` and `ComplexSparseMatrixQuantity` remain explicit
  canonical CSR owners.
- Backend-neutral does not mean representation-neutral storage. Sparse input
  must remain sparse unless the selected, documented policy authorizes an
  explicit dense boundary.
- Unit metadata stays outside the numerical backend. Model code unwraps
  magnitudes and reattaches units to correlated results.
- Sparse and dense implementations must represent the same mathematical
  operator and eigenvalue ordering for their common supported cases.

## Migration inventory

Before implementation, inventory every direct use of:

- `numpy.linalg.eigh` and `numpy.linalg.eigvalsh`;
- `scipy.linalg` complete-spectrum routines;
- `scipy.sparse.linalg.eigsh` and related iterative solvers; and
- implicit or explicit sparse-to-dense conversion before an eigenproblem.

Classify each caller by Hermitian class, scalar type, full/partial selection,
required eigenvectors, degeneracy behavior, matrix size, and current failure
contract. Migrate one connected model family at a time rather than replacing
all calls mechanically.

The checkerboard Bloch model is the first candidate consumer after the neutral
complex-Hermitian contract exists. The existing real-symmetric eigenproblem
implementation should be reconciled with the shared request/result vocabulary
without losing its tridiagonal complete-sparse specialization.

## Required verification

Completion requires tests that establish:

- dense NumPy and sparse SciPy agreement on the same real-symmetric and
  complex-Hermitian operators;
- deterministic ascending eigenvalue ordering;
- eigenvector phase/sign canonicalization where eigenvectors are returned;
- behavior at `k=1`, an ordinary partial selection, `k=N-1`, and `k=N`;
- no dense conversion on maintained sparse partial paths;
- explicit dense fallback only for declared cases;
- stable handling or explicit rejection of degenerate spectra;
- propagation of iterative nonconvergence; and
- preservation of unit attachment at model boundaries.

Run focused numerical tests, every migrated model family, the full suite,
changed-module compilation, Ruff, wheel build, and `git diff --check`.

## Documentation requirements

Update the owning numerical architecture index plus schematic, numeric, and
scientific views. Each migrated model page must identify its mathematical
request, representation, selected backend policy, dense-boundary behavior, and
evidence limitations.

## Non-goals

This task does not select GPU, distributed, or out-of-core backends; promise
bitwise equality between dense and iterative algorithms; define
material-specific convergence; or treat backend agreement as scientific
validation.

## Completion condition

The debt is resolved when model layers no longer branch directly on NumPy versus
SciPy eigenvalue APIs, every migrated sparse path has an auditable densification
policy, and backend-neutral results preserve the existing mathematical, unit,
and evidence contracts.
