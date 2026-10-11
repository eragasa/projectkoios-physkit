# Oblique PIAB2D finite-difference scientific scope

## Represented model

The model is one noninteracting scalar particle confined to the parallelogram

$$
\mathbf r=A\mathbf s,\qquad \mathbf s\in[0,1]^2,
$$

with homogeneous Dirichlet data on every edge. The direct basis is constant,
the mass is positive and isotropic, and the interior potential is zero.

This boundary-value problem is not a periodic primitive cell. Opposite edges
are not identified, no Bloch phase is applied, and the constant function is not
an admissible state.

## Assumptions and non-use

- The direct basis has two finite linearly independent Cartesian columns.
- The coordinate map is affine, so the metric and Jacobian are constant.
- Centered second-order differences are appropriate for the represented states.
- `UnitSystem` values declare one coherent numerical scale.
- The contract covers a constant affine metric, zero boundary data, a scalar
  free particle, and isotropic mass.

## Evidence classification

### Implementation conformance

Nominal type, unit, geometry, sparse-storage, shape, ordering, and normalization
checks are covered by maintained tests.

### Numerical verification

The zero-skew limit is checked against the exact discrete sine spectrum, and its
ground energy exhibits second-order refinement toward the rectangular
continuum value. Sparse Hermiticity and area-weighted normalization are checked
for a skew cell.

### Scientific validation

Scope: mathematical free-particle boundary-value problem. Experimental and
material-specific comparison is outside this package.

### Uncertainty quantification

Output: represented eigenpairs without error bars. Grid refinement supplies the
maintained discretization-error evidence.

### Acceptance

Owner: the consuming project for its declared use.
