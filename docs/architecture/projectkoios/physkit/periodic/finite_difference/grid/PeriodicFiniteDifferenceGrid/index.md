# Class `PeriodicFiniteDifferenceGrid`

## Contract

`PeriodicFiniteDifferenceGrid` is a frozen, slotted `DataObject` composing:

- one `DirectLatticeMetric` with dimension two or three;
- one dimension-matched `FinitePeriodicDomain` whose extents are all at least
  three; and
- one declared direct-basis unit, either `Unitless` or a physical length unit.

For extent `N_i`, the endpoint-excluded fractional sites are
`s_i=n_i/N_i`, `n_i=0,...,N_i-1`, with spacing `h_i=1/N_i`. The finite domain
fixes last-axis-fastest C-order flattening. The direct-basis magnitudes in the
metric are interpreted in the declared unit; the resulting Laplacian unit is
unitless or inverse-length squared.

## Limitations

The unit declaration does not convert the already represented basis
magnitudes. Callers must supply magnitudes in that unit. The record owns no
Cartesian wave vector, kinetic coefficient, or material metadata.

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/finite_difference/grid.py`
- Tests: `tests/projectkoios/physkit/periodic/finite_difference/grid/test__PeriodicFiniteDifferenceGrid__init.py`
