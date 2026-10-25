# `SparseSpectralNormOptionsEnum`

## Responsibility

`SparseSpectralNormOptionsEnum` is the PhysKit-owned selector for the two explicit sparse
spectral representation paths. It is a pure `Enum` with opaque automatically assigned
values.

## Members

### `complete`

`complete` authorizes one explicit CSR-to-dense materialization followed by the complete
NumPy singular-value decomposition owner. Iterative controls are invalid for this path.

### `iterative`

`iterative` invokes the SciPy CSR spectral owner with one exact
`ScipyCsrSpectralNormControls` instance. The path prohibits dense materialization,
complete-SVD retry, and algorithm fallback.

## Selection contract

A sparse spectral evaluation requires one member. No member is inferred from matrix
shape, density, controls, or caller context. Dense evaluation accepts neither member.

The enum does not inherit from `str`. The API accepts no strings, compatibility aliases,
or default sparse choice. The public Optimization control field is
`krylov_dimension`; PhysKit defines no `ncv` name.
