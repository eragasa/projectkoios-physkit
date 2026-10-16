# Class `BlockHoppingRangeBandEvaluator`

This page defines accelerated evaluation of one source and an explicit ordered family
of retained hopping models.

## Current contract

The evaluator composes the policy-free truncator, shared-grid Bloch evaluation, and
batched complete-spectrum evaluation.

- Input includes one source `BlockHoppingModel`.
- Input selections contain explicit conjugate-closed inventories and distinct model
  identifiers.
- Duplicate inventories and identifiers reject before numerical work.
- Every `BlockHoppingTruncationResult` is constructed once and retained exactly.
- Source and retained models use one shared finite-Fourier channel operation.
- Every nonempty Bloch matrix is passed through one concatenated spectrum solve.
- Output order exactly matches selection order.
- Full source Bloch matrices remain available in `source_band_evaluation`.

## Empty and failure semantics

Empty coordinate grids retain structural evidence while skipping both numerical
backends.

- Source and retained Bloch matrix tuples are empty.
- Every energy matrix has shape `(0, basis_size)`.
- Point diagnostics are empty and aggregate diagnostics are `None`.
- Spectrum backend declaration strings remain present without an eigensolver call.
- Invalid selections raise `TypeError` or `ValueError` before Fourier execution.
- Numerical and eigensolver failures preserve the band evaluator's exception classes.

## Policy boundary

The evaluator does not interpret why inventories were supplied or compare their
outputs.

- No metric converts displacements into range.
- No cutoff sequence is constructed or sorted.
- No error, normalization, convergence, or acceptance calculation is present.
- No finite-cyclic alias or Nyquist treatment is inferred.

## Mapping

The implementation and focused construction/evaluation tests are maintained
together.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bands.py::BlockHoppingRangeBandEvaluator`
- Tests:
  - `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingRangeBandEvaluator__init.py`
  - `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingRangeBandEvaluator__evaluate.py`
