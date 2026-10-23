# KSD real-symmetric linear-algebra extraction lineage

## Donor snapshot

- Read-only donor checkout:
  `/Users/eugene/worktrees/ksdft2effmass-source-package-ownership-cleanup`
- Donor revision: `cf36da4387062c4daad87e6a0d77bf336551be67`
- `python/src/ksdft2effmass/linalg/matrices.py`
  - Git blob: `318256de3ccd46d8c973febe4f81b0e7597a7d67`
  - SHA-256: `c704dd74ca167e5ba48bed174ad2862a789d3100b130edcff25ccbf4b29f49ee`
- `python/src/ksdft2effmass/operators/eigenpairs.py`
  - Git blob: `2bc255b247b24d51fd2993d222827ec0cd3dfa78`
  - SHA-256: `0d541b6d25cbdd46857825ff71e081ed78fe704877f6a3df9cfb81daa130d1de`

Both donor paths were clean relative to the stated revision during extraction.
The donor checkout remained read-only.

## Exact ownership mapping

| Donor responsibility | PhysKit owner |
| --- | --- |
| Dense/CSR input canonicalization and verified quantity union | `projectkoios.physkit.linalg.matrices` |
| Sparse-safe Hermiticity predicate | `projectkoios.physkit.numerics.linear_algebra.matrix_properties` |
| Unit-aware real-symmetric result | `projectkoios.physkit.linalg.eigenproblems.real_symmetric.RealSymmetricEigendecomposition` |
| Unit-aware complete/lowest decomposition action | `projectkoios.physkit.linalg.eigenproblems.real_symmetric.RealSymmetricEigendecomposer` |
| Unit-free eigensolver algorithm and arrays | Existing `projectkoios.physkit.numerics.eigenproblems.real_symmetric` |
| Dense and CSR quantity representations | Existing `projectkoios.physkit.units.quantities` |

## Adaptation notes

- Generic linear algebra moved to a top-level PhysKit `linalg` boundary that can
  consume raw NumPy/SciPy representations or canonical PhysKit quantities.
- The Hermiticity calculation was separated into a unit-free numerical kernel.
  Sparse comparison remains sparse and does not materialize a dense matrix.
- The represented eigendecomposition is distinct from the unit-free backend
  result. It retains the exact operator quantity, restores the operator unit on
  eigenvalues, and owns unitless eigenvectors.
- PhysKit canonical names replace KSD-specific quantity names. No wrapper
  subclasses or KSD-to-PhysKit conversion records were introduced.
- PIAB1D composes this represented decomposition rather than storing a separate
  Hamiltonian and bare numerical eigensystem.

## Deliberately excluded policy

PhysKit does not import KSD's `HermitianHamiltonian`. Spectral subspace
selection, compression, residual analysis, acceptance decisions, verification,
serialization, provenance records, and workflow state remain outside this
migration.

Generic matrix Hermiticity supports complex input, but no complex-Hermitian
eigendecomposition was extracted or claimed.

## Verification boundary

Mapped tests establish immutable representation, identity preservation,
unit restoration, exact public input validation, sparse-safe Hermiticity, and
real-symmetric numerical behavior. They do not establish scientific acceptance
of a downstream reduced model or complex-Hermitian solver.
