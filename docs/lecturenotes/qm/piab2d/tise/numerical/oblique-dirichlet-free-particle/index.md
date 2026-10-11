# Free particle in an oblique hard-wall box

This model confines a free particle to a parallelogram. All four physical edges
are hard walls with homogeneous Dirichlet data. Opposite edges are not
identified, so this problem is distinct from a free particle on a periodic
oblique primitive cell.

## Affine coordinates

Let the columns of $A$ be the parallelogram edges and map the unit fractional
square to physical space by

$$
\mathbf r=A\mathbf s,
\qquad \mathbf s\in[0,1]^2.
$$

For constant metric $g=A^{\mathsf T}A$,

$$
\nabla^2
=g^{11}\partial_1^2+g^{22}\partial_2^2
+2g^{12}\partial_1\partial_2.
$$

A skew cell has $g^{12}\ne0$, so independent one-dimensional sine products are
not generally eigenfunctions. They become the familiar rectangular-box basis
only when the metric is diagonal.

## Homogeneous-Dirichlet grid

A full fractional grid contains both boundaries with spacing

$$
h_i=\frac{1}{N_i-1}.
$$

Only points with indices $1,\ldots,N_i-2$ are unknowns. A centered stencil that
would touch index zero or $N_i-1$ instead uses the declared zero boundary value.
This is not the endpoint-excluded quotient grid used for periodic boundaries.

With $\hbar=m=1$, the represented free-particle Hamiltonian is

$$
H_h=-\frac12\Delta_h.
$$

The pure and mixed centered differences are second-order approximations for
sufficiently smooth states.

## Physical normalization

The affine Jacobian is $|\det A|$. Every represented interior point has equal
physical area weight

$$
w_A=|\det A|h_1h_2.
$$

Returned eigenfunction columns satisfy

$$
w_A\sum_p|\psi_p|^2=1.
$$

Boundary values can be restored as zero rows and columns for physical-space
visualization. The full grid must then be mapped through $A$ so the probability
density appears on the physical parallelogram rather than an unrelated
Cartesian rectangle.

## Numerical checks and limitations

Sparse symmetry checks establish Hermiticity of the represented Hamiltonian.
The orthogonal limit can be compared with both the exact discrete sine spectrum
and the continuum rectangular-box energies. Grid refinement measures numerical
convergence but does not supply an error bar for an arbitrary skew cell.

This model is not a material-specific quantum dot, a periodic band calculation,
an interacting system, or scientific validation against experiment.

See the [computational laboratory](../../../../../../notebooks/qm/piab2d/tise/numerical/oblique-dirichlet-free-particle.ipynb).
