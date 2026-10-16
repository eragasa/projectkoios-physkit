# Class `BlockHoppingBandTolerance`

This page defines the numerical tolerance record used by block-band evaluation.

## Current contract

`BlockHoppingBandTolerance` records caller-supplied numerical Hermiticity
parameters without selecting scientific policy.

- `hermiticity_absolute` is a finite nonnegative `ScalarQuantity`.
- The absolute-tolerance unit equals the evaluated model energy unit.
- `hermiticity_relative` is an exact finite binary64 value in `[0, 1)`.
- The pointwise cutoff is $a+r\max_{ij}|A_{ij}|$.
- Construction establishes numerical validity rather than scientific acceptance.

## Mapping

The implementation and invariant tests are maintained together.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bands.py::BlockHoppingBandTolerance`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingBandTolerance__init.py`
