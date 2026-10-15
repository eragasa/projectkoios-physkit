# Class `BlockHoppingPowerResult`

## Current contract

The result retains a `BlockHoppingModel` and one canonical `BlockHoppingPower` per
ordered term. `total_power` returns the finite sum over all terms.
`omitted_power(retained_displacements)` returns the raw power outside a supplied
sorted unique model-subset inventory; the inventory may be empty and need not be
conjugate closed because this method reports a numerical fact rather than applying
a truncation.

No normalization, cutoff, locality label, or acceptance status is owned.

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py::BlockHoppingPowerResult`
- Invariant tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__BlockHoppingPowerResult__init.py`
- Omitted-power tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__BlockHoppingPowerResult__omitted_power.py`
