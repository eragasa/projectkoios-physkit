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

## Fixed-density chemical potential

For unique sampled energies $E_i$ with explicit state-density weights $w_i$,
the finite-temperature chemical potential satisfies

$$
n_{\mathrm{target}}=\sum_i w_i
\left[\exp\left(\frac{E_i-\mu}{k_BT}\right)+1\right]^{-1}.
$$

The inversion should use the original sampled energies and weights. Solving on
histogram centers or a broadened plotting grid adds a separate representation
error. Spin, valley, and other degeneracies belong in the declared weights; they
are not inferred from the band array.

Once $\mu$ is known, a rendered DOS $g(E)$ can be partitioned into occupied and
unoccupied parts, $g(E)f(E)$ and $g(E)[1-f(E)]$. The thermal window is

$$
-\frac{\partial f}{\partial E}
=\frac{f(E)[1-f(E)]}{k_BT},
$$

and $g(E)(-\partial f/\partial E)$ identifies states active in linearized
finite-temperature responses. These plotted transforms inherit the DOS grid and
broadening; their quadrature is not the exact sampled-state inversion.

At strictly zero temperature, a completely filled band followed by a gap does
not determine a unique chemical potential within that gap. Any reported
``Fermi energy'' convention, such as the midgap value, must therefore be stated
separately. The maintained weighted solver intentionally requires positive
thermal energy.

## Evidence boundary

The maintained tests verify the represented free-electron limit, Brillouin-zone
endpoint agreement, weak-potential first-gap behavior, amplitude-sign
invariance, DOS state-count normalization, and a weighted finite-temperature
occupation target. These checks establish the idealized software and numerical
model, not a material-specific semiconductor gap or a zero-temperature
Fermi-energy convention.

See the [computational laboratory](../../../../../notebooks/solidstate/bloch1d/cosinepotential/cosine-potential-bands.ipynb).
