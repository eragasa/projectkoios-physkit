# Density-of-states numeric contract

## Equal-weight sampling

For `N` unique full-zone wave-vector samples, every sampled band energy carries
state-density weight

$$
w=\frac{g\,\Omega_{\mathrm{BZ}}}{(2\pi)^d N},
$$

where `g` is the positive finite state degeneracy, `d` is dimension, and
$\Omega_{\mathrm{BZ}}$ is reciprocal length, area, or volume. The energy matrix
has shape `(N, band_count)`. The represented state density is therefore

$$
\frac{g\,\Omega_{\mathrm{BZ}}}{(2\pi)^d}\,\text{band_count}.
$$

Equal weights require unique samples. In particular, a regular mesh spanning a
periodic interval must use `endpoint=False`; including both opposite boundary
representatives biases the state count.

## Gaussian estimator

The broadened estimate is

$$
g_\sigma(E)=w\sum_{i,n}
\frac{\exp[-(E-E_{ni})^2/(2\sigma^2)]}{\sqrt{2\pi}\sigma}.
$$

The energy grid is strictly increasing and broadening is positive in the exact
band-energy unit. Kernel accumulation uses bounded chunks of sampled energies
to avoid allocating the complete energy-grid by state-count matrix. The result
reports both analytical represented state density and the trapezoidal integral
captured by the finite requested grid.

## Histogram estimator

Explicit strictly increasing bin edges may be nonuniform. Weighted counts are
divided by each bin width. The result distinguishes all represented states from
states included inside the requested edge range; energies outside that range
are intentionally excluded.

## Units

Physical energy matrices require a physical Brillouin-zone measure with
inverse-length power matching dimension. Outputs carry measure divided by
energy, such as `1 / (electron_volt * angstrom ** 2)` in two dimensions.
Explicitly nondimensional energies and zone measures produce `Unitless` output;
physical and unitless conventions cannot be mixed.

## Verification boundary

Tests verify dimensional unit rejection, state-count normalization, Gaussian
integration, histogram bin widths, and visualization data. They do not establish
mesh convergence, broadening convergence, or validity for a named material.

## Navigation

- [Package index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientific.md)
