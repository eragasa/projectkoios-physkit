# One-dimensional cosine-potential Bloch laboratories

- [`cosine-potential-bands.ipynb`](cosine-potential-bands.ipynb) is the
  maintained package-backed band, first-gap, and normalized full-zone 1D DOS
  laboratory.
- [`intrinsic-semiconductor-statistics.ipynb`](intrinsic-semiconductor-statistics.ipynb)
  uses the maintained PhysKit cosine-potential solver to calculate a
  dimensionally explicit 1D intrinsic carrier state, then evaluates separate
  3D bulk nondegenerate intrinsic-semiconductor relations.
- [`exploratory-bands-dos-thermal-prototype.ipynb`](exploratory-bands-dos-thermal-prototype.ipynb)
  preserves the former
  `notebooks/scratch/solidstate/bloch_1d_cosine_potential_bands_dos.ipynb`
  byte-for-byte.

The former band-gap prototype was removed after individual review because all
of its reusable and pedagogical layers have maintained destinations. Its exact
identity and disposition are recorded under `docs/provenance/`. The extended
exploratory artifact retains finite-temperature chemical-potential,
occupied/unoccupied DOS, thermal-window, and interactive presentation layers.
Those layers remain pending separate maintained contracts; its ordinary DOS
layer is superseded by the maintained `dos1d` calculation.
