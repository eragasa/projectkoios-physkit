# `projectkoios.physkit.solidstate.semiconductors.junction`

This module owns the reusable nondegenerate abrupt-junction model extracted from
the retired Sze p–n-junction notebook.

- `SemiconductorMaterial` owns SI bulk material parameters and derived
  permittivity, thermal energy, intrinsic density, and intrinsic-level offset.
- `AbruptDopingProfile` owns uniform acceptor and donor densities.
- `AbruptPNJunction` owns the built-in-potential and depletion-width formulas,
  depletion charge, finite-difference Poisson evaluation, nondegenerate carrier
  profiles, and drift-diffusion terms.
- `DepletionWidths`, `EnergyBands`, `CarrierProfiles`, `CurrentProfiles`, and
  `PNJunctionEquilibrium` are immutable correlated result records.

Array annotations come from
`projectkoios.numerics.numpy.typing`, the pinned Project Koios Optimization
numerical contract. The model uses SI coordinates and densities and represents
band energies in electron-volts.

This is a reduced equilibrium teaching model. It is not a self-consistent
production device simulator and does not establish predictive scientific
validity for a material or junction.
