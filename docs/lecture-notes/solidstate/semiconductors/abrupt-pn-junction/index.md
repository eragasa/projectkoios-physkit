# Nondegenerate abrupt p–n junction

This note describes a one-dimensional equilibrium p–n junction under the
abrupt depletion approximation. Material parameters and device coordinates use
SI units, while band energies are represented in electron-volts.

## Bulk statistics

For effective conduction- and valence-band densities of states $N_C$ and $N_V$,
the nondegenerate intrinsic carrier density is

$$
n_i=\sqrt{N_C N_V}\exp\!\left(-\frac{E_g}{2k_B T}\right),
$$

where $k_B T$ is expressed in electron-volts in the exponential. The intrinsic
level relative to the valence band is

$$
E_i-E_V=k_B T\ln\!\left(\frac{N_V}{n_i}\right).
$$

## Built-in potential and depletion widths

For uniform acceptor and donor densities $N_A$ and $N_D$,

$$
V_{\mathrm{bi}}=\frac{k_B T}{q}
\ln\!\left(\frac{N_A N_D}{n_i^2}\right).
$$

At applied potential $V$, the total depletion width is

$$
W=\sqrt{\frac{2\epsilon_s}{q}
\frac{N_A+N_D}{N_A N_D}(V_{\mathrm{bi}}-V)}.
$$

Charge neutrality partitions the width as

$$
W_n=W\frac{N_A}{N_A+N_D},\qquad
W_p=W\frac{N_D}{N_A+N_D}.
$$

The reduced model assigns fixed ionized-dopant charge inside these regions,
solves the resulting finite-difference Poisson equation with Dirichlet boundary
values, constructs nondegenerate carrier profiles, and reports separate drift
and diffusion current-density terms.

## Interpretation limits

The calculation is an equilibrium pedagogical model. It does not solve Poisson
and carrier continuity equations self-consistently, include recombination,
model graded interfaces, represent degenerate statistics, or establish
material-specific predictive validity. Numerical execution verifies the stated
software calculation only.
