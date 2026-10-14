# Class `ScalarHoppingBlochEvaluation`

## Ownership

`ScalarHoppingBlochEvaluation` is the immutable correlated result of scalar
hopping Bloch evaluation. It retains:

- the exact `ScalarHoppingModel`;
- the unitless fractional reciprocal wave-vector matrix; and
- one complex Bloch value per wave vector in the model energy unit.

Construction validates model ownership, coordinate dimension, coordinate unit,
value count, and value unit. It does not claim Hermiticity, force values to be
real, choose a numerical tolerance, or establish scientific acceptance.

## Mapping

- Implementation: `src/python/projectkoios/physkit/solidstate/tb/bloch.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bloch/test__ScalarHoppingBlochEvaluation__init.py`
