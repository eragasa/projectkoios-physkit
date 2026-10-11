# PIAB2D numerical laboratories and prototypes

- [`oblique-periodic-free-particle.ipynb`](oblique-periodic-free-particle.ipynb)
  is the maintained package-backed laboratory for a nondimensional free
  particle on an oblique periodic primitive cell. It demonstrates the complete
  inverse-metric stencil, Bloch quotient seams, sparse eigensolving, exact
  discrete plane-wave action, second-order refinement, Jacobian-aware
  normalization, and physical-parallelogram density rendering.
- [`oblique-periodic-boundary-prototype.ipynb`](oblique-periodic-boundary-prototype.ipynb)
  preserves the former `notebooks/solidstate/piab2d_pbc.ipynb` byte-for-byte.
- [`generalized-oblique-coordinate-prototype.ipynb`](generalized-oblique-coordinate-prototype.ipynb)
  preserves the former `notebooks/scratch/generalized_piab2d.ipynb`
  byte-for-byte.

The oblique periodic prototype remains present pending final disposition, but
its scientific and pedagogical intent now has a maintained successor. The
stored prototype does not execute cleanly because one constructor annotation is
syntactically invalid. Its saved outputs are therefore not reproducible from
the stored source.

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

The generalized prototype adds derivation of the metric tensor, a mixed-
derivative curvilinear Laplacian, reciprocal-basis construction, plane-wave
energies, real-space reconstruction, and both physical-space and reciprocal-
space plots. It also contains repeated experimental routes and a saved error.
Those layers are preserved because they document distinct numerical ideas and
presentation attempts, not because every route is correct.

The two prototypes are related but not interchangeable. One expresses an
intended periodic-boundary PIAB2D path; the other explores both bounded
curvilinear finite differences and periodic reciprocal-space plane waves. A
future consolidation must first identify each route's state space, boundary
condition, endpoint convention, basis normalization, and error source.

These are preserved exploratory prototypes, not maintained or validated
laboratories. They must not be cited as evidence of a periodic or generalized
PIAB2D solver until the boundary conditions, endpoint treatment, spacing, mixed
derivative, state normalization, and independent numerical checks are corrected
explicitly.
