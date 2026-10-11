# One-dimensional cosine-potential Bloch laboratories

- [`cosine-potential-bands.ipynb`](cosine-potential-bands.ipynb) is the
  maintained package-backed band, first-gap, normalized full-zone 1D DOS,
  fixed-density chemical-potential, and thermal-DOS laboratory.
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
exploratory artifact remains preserved pending final disposition, but its
finite-temperature chemical-potential, occupied/unoccupied DOS, thermal-window,
thermally active DOS, and ordinary DOS layers now have maintained contracts and
package-backed examples.
