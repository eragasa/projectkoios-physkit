# Class `BlockHoppingCoefficientComparison`

## Current contract

The record names one source displacement and its negative target. When the target
exists, `frobenius_defect` is the nonnegative `ScalarQuantity`
$\lVert T[-R]-T[R]^\dagger\rVert_F$ in the model energy unit. It is `None` when
the target is absent. `target_is_present` distinguishes a missing comparison from
an exact zero defect.

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py::BlockHoppingCoefficientComparison`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__BlockHoppingCoefficientComparison__init.py`
