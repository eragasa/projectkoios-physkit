# Exploratory free-electron Bloch prototype review

## Source identity

The reviewed artifact is
`notebooks/solidstate/bloch1d/freeelectron/exploratory-finite-difference-prototype.ipynb`,
relocated byte-for-byte from `notebooks/semphy/freeelectron.bloch1d.ipynb`.

- Git blob: `9a5a1a14b35934061108e93b76ace6d80ad5737b`
- SHA-256: `1dd413a740d143866141fd94af04b27d5c85dd85fdff423926f01529f7e22ef3`
- Preservation commit: `fbf7e9968671cffc18760e9a4aadaf165843d4cc`
- Review baseline: `f0f638f`

The artifact was removed after individual review and maintained replacement.
Its exact bytes remain recoverable from the preservation commit and Git blob
listed above.

## Review evidence

The notebook contains 16 cells: 6 markdown cells, 9 nonempty code cells, and 1
empty code cell. It stores 25 outputs, including 22 PNG images, 2 widget-view
records, and no saved error output.

A clean in-memory execution successfully completed the first seven nonempty code
cells and reproduced the static bands, state, phase, Fourier, and
cosine-potential figures. It then stopped at the first interactive cell because
`ipywidgets` is not an installed project dependency.

All 22 saved PNG displays were extracted and inspected before disposition. They
contain repeated folded-band plots, sampled-state amplitude/density/phase plots,
cell-grid Fourier weights, one cosine-potential band comparison, and duplicated
plots associated with the widget cells. Each valuable presentation has the
maintained destination listed below.

## Maintained destinations

| Prototype layer | Maintained destination |
| --- | --- |
| Exact folded free-electron branches | `FreeElectronBloch1D` |
| Endpoint-excluded finite-difference Hamiltonian | `periodic.pbc1d` and the free-electron finite-difference constructor |
| Grid-normalized sampled state | Maintained free-electron laboratory |
| Amplitude, density, and phase rendering | `BlochState1DVisualizer` |
| Cell-periodic Fourier weights | `BlochState1DPlotData` and `BlochState1DVisualizer` |
| Cosine potential and first-gap extension | Maintained cosine-potential model and laboratory |
| Interactive parameter choice | Explicit reduced-wave-number and band-index notebook inputs |
| Explanatory material | Free-electron lecture note and maintained laboratory |

## Adaptations and rejected behavior

The prototype overlays complex amplitude and probability density on one
ordinate even though a normalized physical amplitude has inverse-square-root
length units and its density has inverse-length units. The maintained visualizer
places them on distinct labeled axes.

The prototype applies a periodic discrete Fourier transform directly to the
nonzero-$k$ twisted Bloch state. The maintained adapter first removes the Bloch
carrier, $u_k=e^{-ikx}\psi_k$, and only then reports integer reciprocal-cell
weights. This prevents carrier leakage from being presented as physical
reciprocal-mode mixing.

Raw truncated SI constants, notebook-local Hamiltonian construction, dense
conversion as a reusable default, and the embedded cosine-potential duplicate
were not promoted. Individual eigenvectors at degeneracy remain explicitly
basis-dependent.

The two widget views are not package behavior. The maintained notebook keeps
parameter selection executable through ordinary explicit variables.

## Final disposition

**Removed after replacement.** Exact folded bands, finite-difference bands,
sampled states, phase, corrected cell-periodic Fourier weights, cosine-potential
bands, and explicit parameter selection all have maintained package, notebook,
and lecture-note destinations. Duplicated code, mixed-unit axes, carrier-leaked
Fourier weights, dense notebook-local defaults, and the optional widget runtime
were not migrated.
