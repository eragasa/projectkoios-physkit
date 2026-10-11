# One-dimensional free-electron Bloch representation

The physical free electron on the infinite line has continuous wave number $q$
and dispersion

$$
E(q)=\frac{\hbar^2q^2}{2m_e}.
$$

A reference-cell Bloch representation writes

$$
q_n=k+G_n,
\qquad
G_n=\frac{2\pi n}{a},
\qquad
k\in[-\pi/a,\pi/a].
$$

The resulting folded branches are

$$
E_n(k)=\frac{\hbar^2}{2m_e}
\left(k+\frac{2\pi n}{a}\right)^2.
$$

These branches are a representation of the free-electron parabola. They are not
caused by bonding, a periodic potential, or lattice scattering.

## Fixed Bloch fibers and the ring interpretation

For fixed $k$, states on the reference cell satisfy

$$
\psi(x+a)=e^{ika}\psi(x).
$$

At $k=0$, this fixed fiber is mathematically equivalent to a particle on a
periodic ring. At nonzero $k$, it is equivalent to a twisted ring. The complete
family of $k$ fibers across the Brillouin zone represents the free electron on
the infinite line; the physical model is not restricted to one literal ring.

## `pbc1d` numerical representation

`projectkoios.physkit.periodic.pbc1d` independently owns the endpoint-excluded
grid and centered finite-difference Laplacian. The positive seam has phase
$e^{ika}$ and the reverse seam has its complex conjugate. This makes the
represented Laplacian Hermitian.

For spacing $\Delta x$, its exact represented symbol is

$$
\lambda_{\Delta x}(q)
=-\frac{4}{\Delta x^2}
\sin^2\left(\frac{q\Delta x}{2}\right).
$$

The corresponding free-electron Hamiltonian is

$$
H_{\Delta x}(k)=-\frac{\hbar^2}{2m_e}L_{\Delta x}(k).
$$

Its difference from the continuum dispersion is discretization error, not
free-electron model error.

## Sampled states and reciprocal content

A grid-normalized physical state satisfies

$$
\sum_j |\psi_k(x_j)|^2\,\Delta x=1,
$$

so its amplitude and probability density carry inverse-square-root-length and
inverse-length units respectively. They should not share one unlabeled ordinate.
The phase of the twisted state contains the explicit Bloch carrier.

To inspect integer reciprocal-cell components, first form

$$
u_k(x_j)=e^{-ikx_j}\psi_k(x_j)
$$

and apply the discrete Fourier transform to $u_k$, not directly to the twisted
full state. Otherwise the noninteger Bloch carrier produces spectral leakage
that can be mistaken for reciprocal-mode mixing.

At an exact or nearly degenerate energy, an eigensolver may return any
orthonormal basis for the represented eigenspace. Real-space or Fourier plots of
one eigenvector are therefore basis-dependent; subspace-level conclusions are
the stable objects.

## Evidence boundary

The exact and finite-difference routes provide software and numerical
verification for the stated idealized model. They do not include a crystal
potential, an effective semiconductor mass, electron interactions, or DFT.

See the [computational laboratory](../../../../../notebooks/solidstate/bloch1d/freeelectron/free-electron-bloch-bands.ipynb).
