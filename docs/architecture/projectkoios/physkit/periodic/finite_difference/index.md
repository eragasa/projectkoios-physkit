# Package `projectkoios.physkit.periodic.finite_difference`

## Current responsibility

This package owns endpoint-excluded finite-difference grids and geometric
Bloch-periodic Laplacians on intrinsic two-dimensional and full
three-dimensional primitive cells. It composes direct-lattice metric geometry,
finite periodic integer domains, boundary-twist reductions, and immutable
complex sparse matrix quantities.

## Public contract

The implementation is bounded to two and three intrinsic dimensions. Direct
primitive vectors consumed by the grid remain the columns of a square `2×2` or
`3×3` basis. A flat 2D lattice supplied as an ambient `3×2` basis must first use
the explicit [`periodic.lattice.embedding`](../lattice/embedding/index.md)
adapter; the grid does not silently discard orientation or reinterpret the plane
as a full slab. Grid arrays use C-order flattening with the last fractional axis
varying fastest.

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
- Laboratory: `notebooks/qm/piab2d/tise/numerical/oblique-periodic-free-particle.ipynb`

## Evidence boundary

Mapped tests establish the declared data invariants, sparse representation,
finite plane-wave identities, and second-order numerical convergence. They do
not establish material-specific scientific validity, uncertainty
quantification, or human acceptance.
