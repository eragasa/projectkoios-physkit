# Class `BlockHoppingPower`

## Current contract

`BlockHoppingPower` binds one exact lattice displacement to the raw squared
Frobenius norm $\sum_{ij}|T_{ij}|^2$ of its hopping block. The
`frobenius_power` value is a nonnegative `ScalarQuantity` carrying the squared
model energy unit, or `Unitless` for a dimensionless model.

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py::BlockHoppingPower`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__BlockHoppingPower__init.py`
