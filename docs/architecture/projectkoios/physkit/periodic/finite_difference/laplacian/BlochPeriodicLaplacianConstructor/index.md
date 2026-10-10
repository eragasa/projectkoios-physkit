# Class `BlochPeriodicLaplacianConstructor`

## Mathematics

For `r=A s`, with inverse direct metric `g^{-1}=(A^T A)^{-1}`, the constructor
represents

```text
Delta_h = sum_i g^{ii} D_ii + 2 sum_{i<j} g^{ij} D_ij,
```

where

```text
D_ii f(n) = [f(n+e_i)-2f(n)+f(n-e_i)] / h_i^2,
D_ij f(n) = [f(n+e_i+e_j)-f(n+e_i-e_j)
             -f(n-e_i+e_j)+f(n-e_i-e_j)] / (4 h_i h_j).
```

If a raw stencil target reduces as `n+d = wrapped + N*q`, its contribution is
multiplied by

```text
exp(2 pi i kappa dot q) = product_a exp(2 pi i kappa_a q_a).
```

The constructor uses the reduced twist representative. Positive and negative
wraps therefore receive conjugate phases, and simultaneous wraps receive the
required phase product. COO contributions are canonicalized into immutable CSR.

## Conventions and limitations

- The Bloch convention is `psi(s+q)=exp(+2 pi i kappa dot q) psi(s)`.
- Primitive vectors are square intrinsic-coordinate bases; 2D `3×2` embeddings
  are outside this contract.
- The centered stencil is second-order accurate for smooth modes.
- No kinetic, mass-tensor, potential, or energy scaling is applied.

## Verification

Tests cover C-order seam indexing, face/corner phase products, Hermiticity, the
zero-twist constant null mode, orthogonal Kronecker-sum reduction, skew-2D and
fcc-primitive 3D plane-wave symbols, second-order convergence, physical-unit
conversion, and canonical immutable CSR storage. This is numerical
verification of the represented formula, not scientific validation of a
material model.

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/finite_difference/laplacian.py`
- Tests: `tests/projectkoios/physkit/periodic/finite_difference/laplacian/test__BlochPeriodicLaplacianConstructor__action.py`
