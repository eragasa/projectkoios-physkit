# Oblique periodic-boundary prototype review

## Source identity

The reviewed artifact is
`notebooks/qm/piab2d/tise/numerical/oblique-periodic-boundary-prototype.ipynb`,
relocated byte-for-byte from `notebooks/solidstate/piab2d_pbc.ipynb`.

- Git blob: `5cbc204ca336a3601ecd7fbfa935c0ff521cab3c`
- SHA-256: `0874fadf3a29997ec9ad059087c8e5e74387b369310a1be2db44466d26ffd7aa`
- Preservation commit: `fbf7e9968671cffc18760e9a4aadaf165843d4cc`
- Review baseline: `4a70159`

The artifact remains present pending final disposition.

## Review evidence

The notebook contains seven code cells, six of which are nonempty. It stores six
outputs, including four PNG displays and no saved error output. A clean
execution stops in the second code cell with a `SyntaxError` at `N: int 50`.
Consequently, its saved eigenspectrum and probability-density displays are not
reproducible from the stored source.

## Maintained destinations

| Prototype layer | Maintained destination |
| --- | --- |
| Direct and reciprocal oblique geometry | `periodic.lattice` and `DirectLatticeMetric` |
| Endpoint-excluded fractional grid | `PeriodicFiniteDifferenceGrid` |
| Pure and mixed inverse-metric stencil | `BlochPeriodicLaplacianConstructor` |
| Periodic and Bloch-twisted seams | Boundary-twist records and quotient-seam construction |
| Sparse low-spectrum solve | Maintained oblique periodic free-particle laboratory |
| Jacobian-aware state normalization | Maintained laboratory and lecture note |
| Physical-parallelogram probability plot | Maintained laboratory |
| Numerical and scientific explanation | Periodic finite-difference architecture and new lecture note |

The maintained successor is
`notebooks/qm/piab2d/tise/numerical/oblique-periodic-free-particle.ipynb`.
It composes the package-owned geometric operator into the explicit
nondimensional Hamiltonian $H=-\Delta_h/2$.

## Rejected and corrected behavior

The prototype calls the model a particle in a box even though its stated intent
is periodic. Its active one-dimensional stencils omit periodic seam entries,
while an unused periodic helper omits spacing normalization. Fractional plotting
coordinates include both zero and one, duplicating the periodic endpoint. The
constructor ignores its intended resolution and spacing inputs in favor of
hard-coded values.

The prototype does not apply the primitive-cell Jacobian when normalizing or
interpreting eigenvectors. Experimental helpers also contain undefined names,
and saved outputs cannot repair the invalid stored source. These behaviors were
not migrated.

The successor uses unique endpoint-excluded points, quotient-seam phases on
pure and simultaneous mixed wraps, C-order state indexing, sparse Hermitian
storage, area weight $|\det A|/(N_1N_2)$, exact discrete plane-wave identities,
and an observed second-order refinement check.

## Current disposition

**Retain pending operator disposition.** The intended periodic free-particle
pedagogy now has a maintained destination, but the exploratory artifact has not
been removed. The separate generalized oblique-coordinate prototype still
requires its own individual review.
