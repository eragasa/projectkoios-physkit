# Embedded 2D lattice scientific scope

## Represented system

This capability represents a flat two-dimensional translation lattice whose
primitive vectors are embedded in three-dimensional Cartesian space. It retains
ambient orientation and nonorthogonality while adapting the induced metric to an
intrinsic two-coordinate representation.

Examples include an idealized planar crystal sheet described in a laboratory
Cartesian frame. Atomic thickness, vacuum spacing, surface relaxation, and
out-of-plane electronic structure are not inferred by this geometric record.

## Assumptions

- The two ambient primitive vectors are finite and numerically full column rank.
- The periodic domain is intrinsically two-dimensional and flat.
- Both primitive vectors share one declared physical length unit or are
  explicitly nondimensional.
- Ambient transformations used for invariance statements are orthogonal
  rotations or reflections.
- Downstream finite differences act on fractional in-plane coordinates only.
- Any ambient vector passed to a method named `project_...` may contain a normal
  component that is deliberately projected away.

## Embedded plane versus full slab supercell

An embedded 2D plane uses

$$
A_{\mathrm{plane}}\in\mathbb R^{3\times2}
$$

and has two independent periodic translations. A full slab or supercell uses

$$
A_{\mathrm{slab}}\in\mathbb R^{3\times3}
$$

and has three lattice coordinates, even when one vector primarily represents
vacuum separation.

This adapter never invents a third vector, vacuum length, boundary condition, or
periodic coordinate. Conversion from an embedded plane to a slab is a separate
modeling operation with separate scientific assumptions.

## Preserved orientation and intrinsic invariance

The original ambient basis is retained exactly as an immutable matrix quantity.
The frame `Q`, projector `P`, and ambient dual basis `B` make orientation and
coordinate mappings inspectable. The intrinsic metric and scalar geometric
Laplacian are invariant under ambient rotation and reflection because they
depend on `A^T A`.

This invariance does not imply that every ambient physical observable is
rotation/reflection invariant. Vector, pseudovector, spin, polarization,
out-of-plane field, and handedness-sensitive models may require the retained
ambient frame and their own transformation laws.

## Separation from physical models

### Geometric Laplacian

The adapter enables reuse of the scalar geometric operator

$$
\Delta=\sum_{i,j}(G^{-1})_{ij}\partial_i\partial_j.
$$

It does not apply electron mass, an effective-mass tensor, `hbar`, a potential,
or an energy scale.

### Atomic structure

The direct basis does not define atomic positions, species, layer thickness,
buckling, defects, or surface termination. Those belong to explicit unit-cell
or material models.

### Electromagnetic and constitutive tensors

Ambient fields and anisotropic constitutive tensors require explicit projection
or pullback rules. They must not be inferred from the scalar metric adapter.

### Serialization versus scientific provenance

The domain codec preserves the defining ambient basis and unit under a strict,
deterministic schema. A successful byte or value round trip establishes wire
conformance only. It does not establish the experimental provenance, relaxation
state, material identity, or scientific suitability of the encoded basis.
Derived geometry is recomputed rather than trusted from serialized caches.

## Evidence classification

### Implementation conformance

Supported by nominal type, unit, shape, rank, immutable-storage, mapping, and
correlated-result tests.

### Numerical verification

Supported for the tested bases by duality, projector, metric, coordinate-map,
orthogonal-transformation, and represented-Laplacian equivalence checks.

### Scientific validation

Not established. The tests do not compare with a material-specific slab model,
experimental data, or an independently validated electronic-structure result.

### Uncertainty quantification

Not established. The capability does not propagate lattice-parameter,
orientation, conditioning, discretization, or model-form uncertainty.

### Human acceptance

Not established by passing software or numerical checks. Consuming projects
retain responsibility for declared-use review and scientific acceptance.

## Navigation

- [Module index](../index.md)
- [Schematic](../schematic/index.md)
- [Numerical page](../numerical/index.md)
