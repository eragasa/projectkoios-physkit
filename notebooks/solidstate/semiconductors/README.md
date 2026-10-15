# Semiconductor laboratories

Core transport laboratories:

- [`carrier-mobility.ipynb`](carrier-mobility.ipynb) evaluates empirical
  lattice and ionized-impurity mobility and a deterministic relaxation-time
  trajectory.
- [`abrupt-pn-junction.ipynb`](abrupt-pn-junction.ipynb) evaluates the
  nondegenerate abrupt-junction depletion model, band diagram, and carrier
  profiles.

Dimension-specific effective-mass ownership follows the package naming scheme:

- `effmass_1d` has no maintained laboratory yet;
- [`effmass_2d`](effmass_2d/README.md) contains the maintained periodic
  primitive-cell methodological review; and
- [`effmass_3d`](effmass_3d/README.md) contains the maintained
  research-facing periodic primitive-cell reference.

The absent 1D effective-mass laboratory is not an implemented capability.
The related
[cosine-potential intrinsic-statistics laboratory](../bloch1d/cosinepotential/intrinsic-semiconductor-statistics.ipynb)
uses the maintained Bloch-band solver and keeps its exact 1D carrier model
separate from the notebook's 3D bulk nondegenerate comparison. The related
[checkerboard intrinsic-statistics laboratory](../bloch2d/checkerboard/intrinsic-semiconductor-statistics.ipynb)
uses a genuine 2D Bloch Hamiltonian and reports dimensionally distinct sheet
carrier densities.
