# Embedded 2D lattice numerical contract

## Ambient direct basis and induced metric

Let the ordered primitive vectors be the columns of a full-column-rank matrix

$$
A=\begin{bmatrix}\mathbf a_1&\mathbf a_2\end{bmatrix}
\in\mathbb R^{3\times2}.
$$

Fractional coordinates remain two-dimensional:

$$
\mathbf r=A\mathbf s,
\qquad \mathbf s\in\mathbb R^2,
\qquad \mathbf r\in\mathbb R^3.
$$

The induced metric and inverse metric are

$$
G=A^{\mathsf T}A,
\qquad G^{-1}=(A^{\mathsf T}A)^{-1}.
$$

`EmbeddedDirectLattice2D` requires shape `(3,2)`, finite entries, one common
physical length unit or `Unitless`, and a scale-invariant normalized area above
`1e-8`.

## Ordered orthonormal frame

The adapter constructs

$$
\mathbf q_1=\frac{\mathbf a_1}{\lVert\mathbf a_1\rVert_2},
$$

$$
\mathbf q_2=
\frac{(I-\mathbf q_1\mathbf q_1^{\mathsf T})\mathbf a_2}
{\lVert(I-\mathbf q_1\mathbf q_1^{\mathsf T})\mathbf a_2\rVert_2},
$$

and

$$
Q=\begin{bmatrix}\mathbf q_1&\mathbf q_2\end{bmatrix},
\qquad Q^{\mathsf T}Q=I_2.
$$

The tangent-plane projector is

$$
P=QQ^{\mathsf T},
\qquad P^2=P,
\qquad P^{\mathsf T}=P.
$$

The ordered Gram--Schmidt convention is deterministic for an ordered basis and
is equivariant under any ambient orthogonal operation `O`, including
reflections:

$$
A'=OA \Longrightarrow Q'=OQ.
$$

Consequently `P'=OPO^T` while the induced metric remains unchanged.

## Intrinsic square adapter

The intrinsic basis is

$$
C=Q^{\mathsf T}A\in\mathbb R^{2\times2}.
$$

It has the same metric:

$$
C^{\mathsf T}C
=A^{\mathsf T}QQ^{\mathsf T}A
=A^{\mathsf T}A
=G,
$$

because both columns of `A` lie in the range of `Q`. The adapter constructs a
nominal `DirectLattice2D` from `C`, followed by the existing
`DirectLatticeMetric`. No rectangular matrix is passed into or accepted by the
existing square metric contract.

The nonorthogonal finite-difference Laplacian therefore uses the same

$$
\Delta=\sum_{i,j=1}^{2}(G^{-1})_{ij}
\frac{\partial^2}{\partial s_i\partial s_j}
$$

whether geometry arrives through the ambient adapter or an explicit intrinsic
basis with metric `G`.

## Ambient dual reciprocal basis

The in-plane ambient dual basis is

$$
B=2\pi A G^{-1}\in\mathbb R^{3\times2}.
$$

It satisfies

$$
A^{\mathsf T}B=2\pi I_2.
$$

For an ambient row vector `r`, the fractional coefficients of its orthogonal
in-plane projection are

$$
\mathbf s^{\mathsf T}
=\frac{1}{2\pi}\mathbf r^{\mathsf T}B.
$$

An out-of-plane component is projected away explicitly by the method named
`project_ambient_to_fractional`; it is not silently interpreted as a third
periodic coordinate.

## Unit transformations

If `A` has unit `L`, then

| Value | Unit |
| --- | --- |
| `C` | `L` |
| `G` | `L²` |
| `G^{-1}` | `L⁻²` |
| `B` | `L⁻¹` |
| `Q`, `P` | unitless |
| adapted geometric Laplacian | `L⁻²` |

A `Unitless` basis gives unitless values throughout. Changing only the physical
length unit rescales numerical magnitudes but yields the same converted geometry
and represented operator.

## Strict wire reconstruction

Schema version 1 serializes only `A`, its column convention, and its unit. The
six basis numbers are finite standard JSON numbers. Decoding reconstructs the
exact nominal `EmbeddedDirectLattice2D`; subsequent adaptation recomputes

$$
G=A^{\mathsf T}A,\qquad
G^{-1},\qquad
B=2\pi A G^{-1},\qquad
Q,\qquad
P=QQ^{\mathsf T}.
$$

No derived numeric array is admitted from the wire. This prevents inconsistent
or stale cached geometry from bypassing the correlated-result invariants.

## Numerical evidence

Mapped tests verify:

- scale-invariant full-column-rank decisions, including meter-scale vectors;
- `G`, `G^{-1}`, `Q`, `P`, and `A^T B=2πI` identities;
- fractional-to-ambient round trips and explicit out-of-plane projection;
- intrinsic-to-ambient frame round trips;
- rotation and reflection invariance of metrics and geometric Laplacians;
- equality between embedded and explicit intrinsic finite-difference routes at
  nonzero Bloch twist;
- physical length-unit conversion and nondimensional support;
- immutable quantity storage;
- golden deterministic schema-v1 JSON and exact nominal-type round trips;
- derived-geometry reconstruction from the decoded defining basis; and
- fail-closed shape, rank, unit, coordinate-shape, type, finiteness, schema,
  duplicate-field, unknown-field, missing-field, and resource-bound checks.

These checks use binary64 tolerances scaled to represented magnitudes. They do
not prove adequacy for arbitrarily ill-conditioned cells or establish a grid
resolution for a downstream scientific calculation.

## Navigation

- [Module index](../index.md)
- [Schematic](../schematic/index.md)
- [Scientific page](../scientific/index.md)
