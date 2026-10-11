# Periodic finite-difference numerical contract

## Coordinate transformation

Let direct primitive vectors be the columns of `A`, with Cartesian and
fractional coordinates related by

$$
\mathbf r=A\mathbf s.
$$

The direct metric and inverse metric are

$$
g=A^{\mathsf T}A,
\qquad
g^{-1}=(A^{\mathsf T}A)^{-1}.
$$

For a constant primitive basis, the Cartesian Laplacian is

$$
\nabla_{\mathbf r}^{2}
=\sum_{i,j}g^{ij}\frac{\partial^2}{\partial s_i\partial s_j}
=\sum_i g^{ii}\partial_i^2
+2\sum_{i<j}g^{ij}\partial_i\partial_j.
$$

The factor two in the second expression is required because the inverse metric
is symmetric and the original double sum contains both `(i,j)` and `(j,i)`.

## Endpoint-excluded grid

For axis extent `N_i`, the represented fractional points and spacing are

$$
s_{i,n}=\frac{n}{N_i},
\qquad n=0,\ldots,N_i-1,
\qquad h_i=\frac{1}{N_i}.
$$

Every extent is at least three so positive and negative centered neighbors are
distinct periodic sites.

## Centered stencil

The pure derivative is

$$
(D_{ii}f)_{\mathbf n}
=\frac{f_{\mathbf n+\mathbf e_i}-2f_{\mathbf n}
+f_{\mathbf n-\mathbf e_i}}{h_i^2}.
$$

For `i != j`, the mixed derivative is

$$
(D_{ij}f)_{\mathbf n}
=\frac{
 f_{\mathbf n+\mathbf e_i+\mathbf e_j}
-f_{\mathbf n+\mathbf e_i-\mathbf e_j}
-f_{\mathbf n-\mathbf e_i+\mathbf e_j}
+f_{\mathbf n-\mathbf e_i-\mathbf e_j}}
 {4h_i h_j}.
$$

The complete represented operator is

$$
\Delta_h=\sum_i g^{ii}D_{ii}+2\sum_{i<j}g^{ij}D_{ij}.
$$

Consequently, an axial neighbor has coefficient `g^{ii}/h_i^2`, while a mixed
neighbor with signs `sigma_i,sigma_j` has coefficient

$$
\frac{g^{ij}\sigma_i\sigma_j}{2h_i h_j}.
$$

In 2D the maximum inventory is nine locations per row. In 3D it is nineteen:
one center, six axial neighbors, and twelve pairwise-mixed neighbors. Exact-zero
mixed coefficients are omitted.

## Bloch quotient seam

The maintained convention is

$$
\psi(\mathbf s+\mathbf q)
=\exp(2\pi i\,\boldsymbol\kappa\cdot\mathbf q)\psi(\mathbf s),
\qquad \mathbf q\in\mathbb Z^d.
$$

If a raw target resolves as

$$
\mathbf n+\boldsymbol\delta=\bar{\mathbf n}+N\mathbf q,
$$

then its matrix coefficient is multiplied by

$$
\exp(2\pi i\,\boldsymbol\kappa\cdot\mathbf q)
=\prod_a\exp(2\pi i\kappa_a q_a).
$$

The constructor uses the representative of `kappa` in `[0,1)` but retains the
unreduced lift and quotient in the result's `TwistFiber`.

## Discrete plane-wave symbol

For integer mode `m` and

$$
\theta_i=\frac{2\pi(m_i+\kappa_i)}{N_i},
$$

a compatible sampled plane wave is an eigenvector with eigenvalue

$$
\lambda_h=
-\sum_i\frac{4g^{ii}}{h_i^2}\sin^2\!\left(\frac{\theta_i}{2}\right)
-2\sum_{i<j}\frac{g^{ij}}{h_i h_j}
\sin(\theta_i)\sin(\theta_j).
$$

As every `h_i` tends to zero for a fixed physical mode,

$$
\lambda_h
=-(2\pi(\mathbf m+\boldsymbol\kappa))^{\mathsf T}
 g^{-1}(2\pi(\mathbf m+\boldsymbol\kappa))
+O(\max_i h_i^2).
$$

## Numerical evidence

Mapped tests verify:

- C-order row/column identities and signed multi-axis phase products;
- Hermiticity for nonzero twists and skew cells;
- the constant null mode at zero twist;
- exact reduction to an orthogonal Kronecker sum when off-diagonal inverse-metric
  entries vanish;
- the discrete symbol in skew 2D and an fcc primitive 3D cell;
- observed second-order plane-wave convergence;
- physical-unit conversion equivalence; and
- canonical, immutable CSR storage.

The maintained oblique periodic free-particle laboratory independently composes
$H=-\Delta_h/2$, verifies the zero-twist constant mode, evaluates a nonzero-
twist discrete plane-wave identity, and displays second-order refinement and
Jacobian-aware state normalization.

These are bounded binary64 checks. They do not prove stability or error bounds
for every cell conditioning, resolution, function class, or downstream solver.

## Navigation

- [Package index](../index.md)
- [Schematic](../schematic/index.md)
- [Scientific page](../scientific/index.md)
