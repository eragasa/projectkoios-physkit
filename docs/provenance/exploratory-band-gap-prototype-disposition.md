# Exploratory cosine band-gap prototype disposition

## Source identity

The reviewed artifact was
`notebooks/solidstate/bloch1d/cosinepotential/exploratory-band-gap-prototype.ipynb`,
relocated byte-for-byte from `notebooks/semphy/cosine.bloch1d.ipynb`.

- Git blob: `f4b38c1cba94e1ad8765eab60e44f54df98b12a0`
- SHA-256: `e10afcc2b0fc2c73d77ba2934e570c4477c548f59ab9f4ecb8bf599679eba8cd`
- Preservation commit: `fbf7e9968671cffc18760e9a4aadaf165843d4cc`

The preserved bytes remain recoverable from that commit and blob after removal
from the maintained notebook tree.

## Review evidence

All 12 nonempty code cells were executed in memory without error. Saved figures
were inspected without rewriting notebook outputs. Successful execution
established only that the prototype ran in the reviewed environment; it did not
validate its physical interpretation or normalization claims.

## Maintained destinations

| Prototype layer | Maintained destination |
| --- | --- |
| Cosine potential | `projectkoios.physkit.potentials.periodic1d` |
| Finite-difference Bloch Hamiltonian and bands | `projectkoios.physkit.solidstate.bloch1d.cosinepotential` |
| Potential rendering | `PeriodicPotential1DVisualizer` |
| Scalar-path band rendering | `BandStructureVisualizer` |
| Full-zone one-dimensional DOS | `projectkoios.physkit.dos1d` and `DensityOfStatesVisualizer` |
| Gap and DOS explanation | Cosine-potential lecture note and maintained computational laboratory |

The maintained successor is
`notebooks/solidstate/bloch1d/cosinepotential/cosine-potential-bands.ipynb`.
It checks the free-electron limit, endpoint spectral agreement, weak-potential
gap scaling, amplitude-sign invariance, grid refinement, and the represented DOS
state count.

## Adaptations and rejected behavior

The exploratory Brillouin-zone sampling included both periodic endpoints with
full equal weight. Its DOS histogram therefore double-counted a boundary
representative, and its density axis was presented in arbitrary units. The
maintained laboratory instead uses unique endpoint-excluded samples and reports
states per energy per length with an explicit degeneracy and state-count
integral.

Notebook-local Hamiltonian assembly and plotting were not retained as reusable
implementations. Illustrative model results were not promoted as evidence for a
real semiconductor or material-specific band gap.

## Disposition

**Removed after replacement.** Every reviewed reusable, numerical, and
pedagogical layer has a named maintained destination. Keeping a second live
notebook would duplicate behavior while preserving known normalization defects.
Removal does not authorize deletion of the separate extended thermal prototype,
whose finite-temperature and interactive layers require their own review.
