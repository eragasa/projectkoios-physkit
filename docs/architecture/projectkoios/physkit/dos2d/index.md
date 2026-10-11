# Package `projectkoios.physkit.dos2d`

`UniformBandDensityOfStates2D` normalizes unique equal-weight band samples over
a complete two-dimensional Brillouin zone. A physical
`brillouin_zone_measure` must have reciprocal-area units. With degeneracy $g$,
every complete band integrates to $g/A$ states per area for a primitive cell of
area $A$. Its shared `solve_chemical_potential(...)` façade assigns that same
per-sample reciprocal-area weight to every exact sampled band energy.

- Code: `src/python/projectkoios/physkit/dos2d/model.py`
- Tests: `tests/projectkoios/physkit/dos2d/model/`
- [Shared DOS contract](../dos/index.md)
