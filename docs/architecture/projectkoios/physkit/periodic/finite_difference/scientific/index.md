# Periodic finite-difference scientific scope

## Represented scientific statement

This capability represents the scalar geometric Laplacian on a primitive
periodic cell with a constant, possibly nonorthogonal basis. Its coordinate
identity is

$$
\nabla_{\mathbf r}^{2}
=\sum_{i,j}(A^{\mathsf T}A)^{-1}_{ij}
\frac{\partial^2}{\partial s_i\partial s_j}.
$$

The numerical representation is a centered second-order finite-difference
approximation on a uniform endpoint-excluded fractional grid with a specified
Bloch twist.

## Assumptions

- The primitive basis is constant over the represented cell.
- The 2D case uses intrinsic in-plane coordinates and a square `2×2` basis.
- The 3D case uses a full-rank square `3×3` basis.
- The scalar field is sufficiently smooth for centered second-order differences
  to have their expected asymptotic behavior.
- Every fractional-axis extent is at least three.
- The caller supplies basis magnitudes in the declared `PhysicalUnit` or
  explicitly declares a nondimensional basis with `Unitless`.
- Bloch twists are expressed in turns along direct primitive axes, not as
  Cartesian inverse-length vectors.

## Explicit separations

### Geometry versus kinetic energy

The result represents `Delta`, not a Hamiltonian. For an isotropic scalar mass,
a consuming quantum model may form

$$
T=-\frac{\hbar^2}{2m}\Delta,
$$

but the coefficient, its units, and the energy scale remain outside this
package.

### Geometry versus an effective-mass tensor

A Cartesian anisotropic effective-mass model represents a constitutive operator
such as

$$
-\frac{\hbar^2}{2}\nabla\cdot M^{-1}\nabla.
$$

That operator is not obtained by multiplying this scalar Laplacian by a matrix.
It requires an explicitly transformed constitutive tensor and a separately
owned contract.

### Intrinsic 2D versus embedded 2D

A `3×2` basis embedding a flat sheet in three-dimensional Cartesian space is
not silently projected by the grid. Such input must use the explicit
[`periodic.lattice.embedding`](../../lattice/embedding/index.md) adapter, which
retains ambient orientation, frame, projector, and dual-basis data before
supplying the intrinsic metric. A full `3×3` slab supercell remains distinct.

## Units and physical interpretation

`DirectLatticeMetric` owns finite numerical basis and metric arrays without unit
metadata. `PeriodicFiniteDifferenceGrid.direct_basis_unit` declares how those
basis magnitudes are interpreted. A physical basis produces an inverse-length-
squared matrix unit; a nondimensional basis produces `Unitless`.

The unit declaration does not independently prove that caller-supplied
magnitudes use that unit. It is an explicit package-boundary assertion checked
for dimensional compatibility, not provenance validation.

## Evidence classification

### Implementation conformance

Supported by typed-record failure tests, exact index/phase cases, sparse storage
checks, and architecture-to-source mapping.

### Numerical verification

Supported for the tested orthogonal, skew-2D, and fcc-primitive 3D cases by
analytic discrete symbols, Hermiticity, null-space behavior, and observed
second-order convergence.

### Scientific validation

Not established. The maintained nondimensional oblique-cell laboratory is a
software and numerical demonstration, not comparison with material-specific
experiments, higher-fidelity electronic-structure calculations, or an
independently validated effective-mass model.

### Uncertainty quantification

Not established. The capability does not estimate discretization uncertainty,
cell-parameter uncertainty, floating-point conditioning error, or model-form
uncertainty for a consuming calculation.

### Human acceptance

Not established by tests, documentation, or the existence of a pinnable
revision. A consuming project retains responsibility for its declared-use
review and acceptance criteria.

## Intended use and non-use

Appropriate uses include reusable geometric operator construction, numerical
method verification, and composition into separately owned isotropic kinetic
models. Unsupported claims include material band-structure validation,
anisotropic effective-mass equivalence, equivalence to a full slab-supercell
model, or production convergence at an undeclared grid resolution.

## Navigation

- [Package index](../index.md)
- [Schematic](../schematic/index.md)
- [Numerical page](../numerical/index.md)
