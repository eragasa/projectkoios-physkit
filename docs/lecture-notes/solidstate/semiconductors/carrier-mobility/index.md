# Empirical carrier mobility and relaxation-time trajectories

This note introduces a reduced carrier-mobility model for a nondegenerate
semiconductor. Mobility is represented in the conventional unit
$\mathrm{cm^2/(V\,s)}$, effective mass in kilograms, impurity concentration in
$\mathrm{cm^{-3}}$, temperature in kelvin, and electric field in
$\mathrm{V/m}$.

## Scattering contributions

The lattice-scattering contribution is represented by the power law

$$
\mu_{\mathrm{lat}}(T)
= \mu_{\mathrm{lat,ref}}
  \left(\frac{T}{T_{\mathrm{lat,ref}}}\right)^{-\alpha_{\mathrm{lat}}}.
$$

The screened ionized-impurity contribution is represented by

$$
\mu_{\mathrm{imp}}(T,N_I)
= \mu_{\mathrm{imp,ref}}
  \left(\frac{T}{T_{\mathrm{imp,ref}}}\right)^{\alpha_{\mathrm{imp}}}
  \frac{N_{I,\mathrm{ref}}}{N_I}.
$$

Matthiessen's rule combines the two channels:

$$
\frac{1}{\mu}=\frac{1}{\mu_{\mathrm{lat}}}
              +\frac{1}{\mu_{\mathrm{imp}}}.
$$

These equations are empirical teaching approximations. They omit detailed band
structure, screening corrections, carrier degeneracy, and field-dependent
velocity saturation.

## Thermal velocity and relaxation time

The three-dimensional root-mean-square thermal velocity is

$$
v_{\mathrm{th}}=\sqrt{\frac{3k_B T}{m^*}},
$$

and the relaxation-time relation is

$$
\mu=\frac{q\tau}{m^*}.
$$

The associated trajectory laboratory uses a first-order collision probability
$\Delta t/\tau$ and redraws one velocity component from a Gaussian thermal
distribution after a collision. It is a pedagogical random walk, not a
Boltzmann-transport solver or a material acceptance model.
