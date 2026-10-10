# Package `projectkoios.physkit.periodic.finite_difference`

## Current responsibility

This package owns endpoint-excluded finite-difference grids and geometric
Bloch-periodic Laplacians on intrinsic two-dimensional and full
three-dimensional primitive cells. It composes direct-lattice metric geometry,
finite periodic integer domains, boundary-twist reductions, and immutable
complex sparse matrix quantities.

## Public contract

The implementation is bounded to two and three dimensions. Direct primitive
vectors are stored as the columns of a square `2×2` or `3×3` basis. A 2D basis
embedded as a `3×2` matrix requires a separate explicit coordinate map and is
not accepted here. Grid arrays use C-order flattening with the last fractional
axis varying fastest.

The represented operator is the geometric Laplacian only. This package does not
apply `-hbar²/(2m)`, an effective-mass tensor, a potential, or any other
Hamiltonian scaling.

## Architecture views

- [Schematic](schematic/index.md): Mermaid component, construction-sequence,
  indexing, phase, and ownership diagrams.
- [Numerical](numerical/index.md): coordinate transformation, full metric
  stencil, Bloch seam, discrete symbol, convergence, and numerical evidence.
- [Scientific](scientific/index.md): assumptions, intended interpretation,
  explicit model separations, limitations, and evidence boundaries.

## Module navigation

- [Module `grid`](grid/index.md)
- [Module `laplacian`](laplacian/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/finite_difference/`
- Tests: `tests/projectkoios/physkit/periodic/finite_difference/`

## Evidence boundary

Mapped tests establish the declared data invariants, sparse representation,
finite plane-wave identities, and second-order numerical convergence. They do
not establish material-specific scientific validity, uncertainty
quantification, or human acceptance.
