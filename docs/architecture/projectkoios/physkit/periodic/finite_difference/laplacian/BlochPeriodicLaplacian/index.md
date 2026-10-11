# Class `BlochPeriodicLaplacian`

The frozen, slotted `ResultsObject` retains the complete construction request,
the reduced `TwistFiber` tagged `QUOTIENT_SEAM`, and a canonical immutable
`ComplexSparseMatrixQuantity`. Matrix shape equals the square of the grid point
count, and its unit is the grid's unitless or inverse-length-squared Laplacian
unit.

The result establishes a represented geometric operator. It does not establish
kinetic-energy scaling, a Hamiltonian, or material-specific validity.

- Code: `src/python/projectkoios/physkit/periodic/finite_difference/laplacian.py`
- Tests: `tests/projectkoios/physkit/periodic/finite_difference/laplacian/test__BlochPeriodicLaplacianConstructor__action.py`
