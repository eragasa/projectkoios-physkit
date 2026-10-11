# One-dimensional cosine-potential Bloch bands

The model adds one prescribed periodic Fourier component to the bare
free-electron Hamiltonian:

$$
H=-\frac{\hbar^2}{2m_e}\frac{d^2}{dx^2}
+V_0\cos\left(\frac{2\pi x}{a}\right).
$$

The reference-cell state satisfies

$$
\psi(x+a)=e^{ika}\psi(x).
$$

The generic endpoint-excluded grid and seam phases are owned by `pbc1d`; the
cosine model owns the potential and band solve.

## Free-electron limit

At $V_0=0$, the represented Hamiltonian is exactly the maintained
free-electron finite-difference Hamiltonian. The continuum reference branches
are

$$
E_n(k)=\frac{\hbar^2}{2m_e}
\left(k+\frac{2\pi n}{a}\right)^2.
$$

This exact reduction is a software-verification requirement.

## First gap

The cosine potential has reciprocal components at $G=\pm2\pi/a$. At the first
zone boundary it couples the degenerate free-electron modes. In the weak
potential limit, the first gap approaches

$$
\Delta E\approx |V_0|.
$$

A represented finite-difference calculation additionally carries grid error.
The potential effect and discretization error must remain distinct.

Changing $V_0$ to $-V_0$ translates the cosine by half a cell, so the complete
spectrum is unchanged. The spectra at $k=-\pi/a$ and $k=\pi/a$ also agree.

## One-dimensional density of states

For an explicit state degeneracy $g$, the DOS per length is

$$
g_{1D}(E)=\frac{g}{2\pi}\sum_n\int_{\mathrm{BZ}}
\delta(E-E_n(k))\,dk.
$$

Where isolated roots satisfy $E_n(k_i)=E$, the delta-function identity gives
contributions proportional to $|dE_n/dk|_{k_i}^{-1}$. Flat band regions
therefore produce large one-dimensional DOS, while an exact band gap contains
no unbroadened states.

A regular numerical mesh must contain unique representatives of the periodic
Brillouin zone. Sampling both $-\pi/a$ and $+\pi/a$ and assigning both a full
uniform weight double-counts one boundary representative. An endpoint-excluded
mesh gives every sampled eigenvalue weight

$$
\frac{g}{aN_k},
$$

so each complete represented band integrates to $g/a$ states per length.
Gaussian broadening or energy bins introduce numerical resolution scales; they
do not redefine exact band edges.

## Evidence boundary

The maintained tests verify the represented free-electron limit, Brillouin-zone
endpoint agreement, weak-potential first-gap behavior, amplitude-sign
invariance, and DOS state-count normalization. These checks establish the
idealized software and numerical model, not a material-specific semiconductor
gap. Finite-temperature analyses retained in exploratory notebooks still
require separate occupation and chemical-potential contracts.

See the [computational laboratory](../../../../../notebooks/solidstate/bloch1d/cosinepotential/cosine-potential-bands.ipynb).
