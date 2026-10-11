# Free particle on an oblique periodic cell

This model represents a nondimensional free particle on one two-dimensional
primitive cell with Bloch-periodic boundary conditions. It is not a hard-wall
particle in a box: opposite cell edges are identified up to declared Bloch
phases.

## Fractional coordinates and metric

Let the direct primitive vectors be the columns of $A$ and write

$$
\mathbf r=A\mathbf s,
\qquad
\mathbf s=(s_1,s_2).
$$

For the constant direct metric

$$
g=A^{\mathsf T}A,
\qquad
g^{-1}=(A^{\mathsf T}A)^{-1},
$$

the Cartesian Laplacian becomes

$$
\nabla^2
=g^{11}\partial_1^2+g^{22}\partial_2^2
+2g^{12}\partial_1\partial_2.
$$

The mixed derivative is required whenever the inverse metric has a nonzero
off-diagonal entry. Mapping a rectangular fractional grid through $A$ produces
the physical parallelogram; it does not change the topology into a hard-wall
box.

## Periodic grid and Bloch twists

A regular represented grid contains unique fractional coordinates

$$
s_{i,n}=\frac{n}{N_i},
\qquad n=0,\ldots,N_i-1.
$$

The duplicated endpoint $s_i=1$ is excluded. A twist lift
$\boldsymbol\kappa$ in turns imposes

$$
\psi(\mathbf s+\mathbf q)
=\exp(2\pi i\,\boldsymbol\kappa\cdot\mathbf q)\psi(\mathbf s),
\qquad \mathbf q\in\mathbb Z^2.
$$

Pure and mixed centered stencils crossing one or both seams receive the complete
product phase. Merely defining an unused periodic helper does not make a
Dirichlet stencil periodic.

## Nondimensional Hamiltonian

With $\hbar=m=1$, the represented Hamiltonian is

$$
H_h=-\frac12\Delta_h.
$$

At zero twist, the constant vector is an exact null mode of the geometric
Laplacian. At general twist, sampled plane waves are

$$
\psi_{\mathbf m\boldsymbol\kappa}(\mathbf n)
=\exp\left[
2\pi i\sum_j\frac{(m_j+\kappa_j)n_j}{N_j}
\right].
$$

Their continuum energies are

$$
E_{\mathbf m}(\boldsymbol\kappa)
=\frac12[2\pi(\mathbf m+\boldsymbol\kappa)]^{\mathsf T}
 g^{-1}[2\pi(\mathbf m+\boldsymbol\kappa)].
$$

The centered finite-difference symbol differs from this expression by
second-order grid error for fixed mode and cell geometry.

## State normalization and visualization

For primitive-cell area $\Omega=|\det A|$, each uniform fractional sample has
physical area weight

$$
\Delta A=\frac{\Omega}{N_1N_2}.
$$

A represented state is normalized using

$$
\sum_{n_1,n_2}|\psi_{n_1n_2}|^2\Delta A=1.
$$

Probability density must be drawn at Cartesian positions $A\mathbf s$ on the
physical parallelogram. A single free-particle plane wave has uniform density;
nonuniform density from an arbitrary combination inside a degenerate eigenspace
is basis-dependent.

## Evidence boundary

Hermiticity, the zero-twist constant mode, exact discrete plane-wave action, and
observed second-order refinement verify the represented operator. The model
scope is a noninteracting scalar free particle on a constant-metric periodic
cell.

See the [computational laboratory](../../../../../../notebooks/qm/piab2d/tise/numerical/oblique-periodic-free-particle.ipynb).
