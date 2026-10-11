# Density-of-states scientific contract

For bands $E_n(k)$ in $d$ reciprocal dimensions, the DOS per physical
$d$-measure is

$$
g_d(E)=\frac{g}{(2\pi)^d}\sum_n
\int_{\mathrm{BZ}}\delta(E-E_n(k))\,d^d k.
$$

Consequently, a complete band integrates to

$$
\int g_d(E)\,dE
=\frac{g\,\Omega_{\mathrm{BZ}}}{(2\pi)^d}.
$$

For a primitive direct cell of measure $\Omega_{\mathrm{cell}}$,
$\Omega_{\mathrm{BZ}}=(2\pi)^d/\Omega_{\mathrm{cell}}$, so one complete band
contains $g/\Omega_{\mathrm{cell}}$ states per direct-space measure.

The dimensional packages differ in normalization and physical output units, not
in the abstract delta-function definition:

- `dos1d` returns states per energy per length;
- `dos2d` returns states per energy per area;
- `dos3d` returns states per energy per volume.

Gaussian broadening is a declared visualization and approximation parameter. It
can fill a true spectral gap with exponentially small tails and must not be used
alone to redefine band edges. A histogram depends on bin boundaries and does
not resolve singularities below its energy-bin scale.

State degeneracy is explicit. The estimator does not infer spin, valley,
orbital, or other multiplicities. It also does not infer occupations or carrier
densities; those require a separately stated statistical model.

A normalized numerical DOS is not experimental validation. Mesh density,
energy resolution, basis convergence, model fidelity, and uncertainty remain
separate scientific responsibilities.

## Navigation

- [Package index](index.md)
- [Schematics](schematics.md)
- [Numeric contract](numeric.md)
