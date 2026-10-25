# `SparseSpectralNormOptionsEnum`

`SparseSpectralNormOptionsEnum` is a pure PhysKit enum with lower-case members:

- `complete`
- `iterative`

The values are opaque and are not strings. The enum selects an explicit sparse
quantity-representation path; it does not redefine the canonical spectral norm semantic
identity.

- `complete` authorizes explicit CSR materialization followed by dense complete SVD.
- `iterative` selects the sparse Optimization evaluator and prohibits dense fallback.

There is no default sparse option, string shortcut, compatibility alias, or public `ncv`
name.
