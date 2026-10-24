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
- `python/src/ksdft2effmass/analysis/operators/eigenpair_diagnostics.py`
  - Git blob: `2aa7e6399ff599bc536699faa81c094eb60786a7`
  - SHA-256: `9f56013d9a0a8da0836d7daa64f172eb429a8150d8030a5957c249ed10cf2eef`

All donor paths were clean relative to the stated revision during extraction.
The donor checkout remained read-only.

## Exact ownership mapping

| Donor responsibility | PhysKit owner |
| --- | --- |
| Dense/CSR input canonicalization and verified quantity union | `projectkoios.physkit.linalg.matrices` |
| Sparse-safe Hermiticity predicate | `projectkoios.physkit.numerics.linear_algebra.matrix_properties` |
| Eager represented decomposition façade, nominal contract, and immutable eigensystem | `projectkoios.physkit.linalg.eigen` |
| Bounded real/complex matrix-family dispatch | `projectkoios.physkit.linalg.eigen.factory` |
| Eager complete real-symmetric operation and retained-prefix selection | Internal `projectkoios.physkit.linalg.eigen._real_symmetric` |
| Unit-free column residual norms | `projectkoios.physkit.numerics.eigenproblems.residuals.eigenpair_residual_norms` |
| Unit-free eigensolver algorithm and arrays | Existing `projectkoios.physkit.numerics.eigenproblems.real_symmetric` |
| Dense and CSR quantity representations | Existing `projectkoios.physkit.units.quantities` |

## Adaptation notes

- Generic linear algebra moved to a top-level PhysKit `linalg` boundary that can
  consume raw NumPy/SciPy representations or canonical PhysKit quantities.
- The Hermiticity calculation was separated into a unit-free numerical kernel.
  Sparse comparison remains sparse and does not materialize a dense matrix.
- `SpectralDecomposition(matrix)` eagerly computes and retains the complete
  N-state public `Eigensystem`. `eigensystem(n_rank=None)` selects only from
  that retained evidence: `None`, zero, and N return the full object by
  identity, while positive smaller ranks return immutable lowest prefixes.
  `max_rank` is matrix dimension N and never a tolerance-dependent numerical
  rank. All paths retain the exact matrix quantity, restore the matrix unit on
  eigenvalues and residual evidence, and own unitless eigenvectors.
- Generic residual algebra was separated from KSD's diagnostic policy. The
  numerical kernel computes unscaled column 2-norms without densifying sparse
  matrices; the represented spectral operation restores the matrix unit. The
  generic public `Eigensystem` retains supplied residual evidence and does not
  directly verify the eigen-equation or derive residuals. Tolerance and
  pass/fail selection remain outside PhysKit decomposition.
- PhysKit canonical names replace KSD-specific quantity names. No wrapper
  subclasses or KSD-to-PhysKit conversion records were introduced.
- PIAB1D uses only the public `SpectralDecomposition` façade and exact public
  `Eigensystem` rather than importing the internal real implementation or
  storing a separate Hamiltonian and bare numerical eigensystem.

## Deliberately excluded policy

PhysKit does not import KSD's `HermitianHamiltonian`. Spectral subspace
selection, compression, scaled/reference-vector residual diagnostics,
domain-specific acceptance decisions, verification, serialization, provenance
records, and workflow state remain outside this migration.

Generic matrix Hermiticity supports complex input, but no complex-Hermitian
eigendecomposition was extracted or claimed.

## Verification boundary

Mapped tests establish immutable representation, identity preservation,
unit restoration, exact public input validation, sparse-safe Hermiticity,
sparse-safe eigenpair residuals, and real-symmetric numerical behavior. They do
not establish scientific acceptance
of a downstream reduced model or complex-Hermitian solver.
