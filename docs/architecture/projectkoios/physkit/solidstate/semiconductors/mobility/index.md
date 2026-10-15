# `projectkoios.physkit.solidstate.semiconductors.mobility`

This module owns reusable empirical carrier-mobility calculations extracted
from the retired Sze mobility notebook.

- `CarrierMobilityModel` owns one carrier's effective mass, lattice and
  ionized-impurity reference mobilities, reference conditions, power-law
  exponents, mobility contributions, Matthiessen combination, thermal velocity,
  relaxation-time conversion, and drift velocity.
- `SemiconductorMobilityModel` groups nominal electron and hole mobility models.
- `CarrierTrajectory` stores one one-dimensional stochastic realization.
- `simulate_electron_trajectory` reproduces the source notebook's first-order
  relaxation-time random walk with an explicit random seed.

Array annotations come from
`projectkoios.numerics.numpy.typing`, the pinned Project Koios Optimization
numerical contract. Effective mass is in kilograms, mobility in
`cm^2 / (V s)`, impurity concentration in `cm^-3`, temperature in kelvin,
electric field in `V / m`, and velocity in `m / s`.

The module does not own material-specific calibration, high-field saturation,
degenerate transport, or scientific acceptance policy.
