# Full-zone density-of-states library provenance

## Extraction baseline

The DOS library was adapted from repository-owned notebook behavior at commit
`a45184c67971839738a92a208c0ae6a0fbaff945`.

| Source | Git blob | Retained contribution |
|---|---|---|
| `notebooks/solidstate/bloch1d/cosinepotential/exploratory-band-gap-prototype.ipynb` | `f4b38c1cba94e1ad8765eab60e44f54df98b12a0` | One-dimensional band-energy histogram and $dk/(2\pi)$ normalization intent |
| `notebooks/solidstate/bloch2d/checkerboard/intrinsic-semiconductor-statistics.ipynb` | `138975b03230b1d35c363c31babd3719d02c2a69` | Endpoint-excluded full-zone mesh, Gaussian broadening, spin degeneracy, and per-area state-count check |
| `notebooks/solidstate/semiconductors/effmass_3d/primitive-cell.ipynb` | `33ce85c5a1ae794baf9404b1ac85df1b0c951c7d` | Unit-aware nonorthogonal reciprocal-cell geometry and three-dimensional folded free-particle spectrum |

## Corrections and generalization

The exploratory one-dimensional histogram sampled both equivalent Brillouin-zone
endpoints and multiplied every sample by the same `dk`, overcounting the
represented periodic interval. The maintained contract instead accepts the
complete Brillouin-zone measure and divides it among a declared count of unique
equal-weight samples. Maintained regular meshes use `endpoint=False`.

Output is no longer labeled arbitrary units. Dimensional models validate
reciprocal length, area, or volume and return DOS per energy per corresponding
direct-space measure. State degeneracy is explicit. Histogram results distinguish
represented states from states inside the selected bin range, while Gaussian
results distinguish analytical represented state density from the integral
captured by the finite plotting grid.

The three-dimensional implementation is the direct dimensional form of the same
full-zone definition; no empirical material DOS was imported or inferred.

## Evidence boundary

Tests verify represented normalization, units, bin widths, Gaussian integration,
and visualization adapters. Notebook execution checks demonstrate use with the
maintained band models. This evidence does not establish mesh or broadening
convergence for every spectrum, experimental validity, or a material-specific
carrier model.
