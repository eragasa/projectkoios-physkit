# PIAB2D numerical laboratories and prototypes

- [`oblique-periodic-free-particle.ipynb`](oblique-periodic-free-particle.ipynb)
  is the maintained package-backed laboratory for a nondimensional free
  particle on an oblique periodic primitive cell. It demonstrates the complete
  inverse-metric stencil, Bloch quotient seams, sparse eigensolving, exact
  discrete plane-wave action, second-order refinement, Jacobian-aware
  normalization, and physical-parallelogram density rendering.
- [`oblique-dirichlet-free-particle.ipynb`](oblique-dirichlet-free-particle.ipynb)
  is the maintained package-backed laboratory for a free particle bounded by
  four hard-wall parallelogram edges. It demonstrates interior-only Dirichlet
  unknowns, sparse inverse-metric finite differences, area normalization,
  physical-space states, and rectangular-limit refinement.
- [`oblique-periodic-boundary-prototype.ipynb`](oblique-periodic-boundary-prototype.ipynb)
  preserves the former `notebooks/solidstate/piab2d_pbc.ipynb` byte-for-byte.
The historical oblique periodic prototype remains present. Its maintained
successor is `oblique-periodic-free-particle.ipynb`. The stored prototype stops
at a syntactically invalid constructor annotation, so its saved outputs cannot
be reproduced from the stored source.

The prototype's active code constructs centered Dirichlet-style operators
without seam couplings and maps an endpoint-included parameter grid through an
oblique direct basis. Its unused periodic helper does not establish periodic
boundary behavior. The maintained successor instead consumes the package-owned
endpoint-excluded Bloch-periodic metric Laplacian.

The preserved goals include:

- metric-aware differentiation in oblique coordinates;
- a two-dimensional free-particle eigenproblem;
- physical-coordinate mapping of represented states; and
- probability-density visualization on a parallelogram.

The generalized prototype was also reviewed individually. Clean sequential
execution stops before its first reciprocal-space calculation because `B` is
undefined, and a later prose cell stores a `SyntaxError`. Its first route is a
bounded homogeneous-Dirichlet curvilinear problem, but it uses unit spacing for
50 points on a declared unit interval and omits cell-Jacobian state
normalization. Its second route explores periodic reciprocal plane waves, with
later cells mixing incompatible array lengths and direct bases.

The generalized prototype's periodic reciprocal-space content is superseded by
the maintained `PeriodicFreeParticle2D` methodological review and the
nondimensional oblique periodic laboratory. Its scientifically distinct bounded
route is represented by `ObliquePiab2D`, its sparse finite-difference solver,
and the maintained oblique Dirichlet laboratory. The prototype was removed
after final disposition; exact identity, recovery data, saved-figure review, and
destination mapping are recorded in
`docs/provenance/generalized-oblique-coordinate-prototype-review.md`.

The remaining oblique periodic prototype is historical. Use
`oblique-periodic-free-particle.ipynb` for the maintained calculation.
