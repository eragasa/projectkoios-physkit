# Package `projectkoios.physkit.dos3d`

`UniformBandDensityOfStates3D` normalizes unique equal-weight band samples over
a complete three-dimensional Brillouin zone. A physical
`brillouin_zone_measure` must have reciprocal-volume units. With degeneracy
$g$, every complete band integrates to $g/\Omega$ states per volume for a
primitive cell of volume $\Omega$. Its shared
`solve_chemical_potential(...)` façade assigns that same per-sample
reciprocal-volume weight to every exact sampled band energy.

- Code: `src/python/projectkoios/physkit/dos3d/model.py`
- Tests: `tests/projectkoios/physkit/dos3d/model/`
- [Shared DOS contract](../dos/index.md)
