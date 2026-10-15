# Sze semiconductor notebook migration

## Retired in-repository sources

- `notebooks/sze_ch2_01_01_mobility.ipynb`
- `notebooks/sze_ch3_01.ipynb`

These exploratory notebooks mixed reusable implementation, repeated prototypes,
plots, and stale cross-notebook excerpts. They were retired after extraction;
they are not runtime dependencies.

## Maintained package adaptations

| Retired source behavior | Maintained owner |
| --- | --- |
| Lattice mobility power law | `projectkoios.physkit.solidstate.semiconductors.mobility.CarrierMobilityModel.lattice_mobility` |
| Ionized-impurity mobility power law | `CarrierMobilityModel.impurity_mobility` |
| Matthiessen combination, thermal velocity, relaxation time, drift velocity | `CarrierMobilityModel` |
| Relaxation-time random walk | `simulate_electron_trajectory` |
| Bulk nondegenerate semiconductor properties | `projectkoios.physkit.solidstate.semiconductors.junction.SemiconductorMaterial` |
| Built-in potential and abrupt depletion widths | `AbruptPNJunction` |
| Depletion charge and finite-difference Poisson solution | `AbruptPNJunction.solve_equilibrium` |
| Band, carrier, and drift-diffusion records | `PNJunctionEquilibrium` and its owned result records |

The maintained implementations preserve the represented equations and unit
conventions while adding explicit argument validation, nominal result types,
shared array typing from the pinned Project Koios Optimization dependency, and
ordinary software tests. Repeated or incompatible prototypes were not retained.

## Maintained teaching artifacts

- `notebooks/solidstate/semiconductors/carrier-mobility.ipynb`
- `notebooks/solidstate/semiconductors/abrupt-pn-junction.ipynb`
- `docs/lecture-notes/solidstate/semiconductors/carrier-mobility/index.md`
- `docs/lecture-notes/solidstate/semiconductors/abrupt-pn-junction/index.md`

The notebooks import package behavior rather than redefining it. Successful
execution demonstrates only the recorded software calculation; it does not
establish pedagogical validation or predictive scientific validity.
