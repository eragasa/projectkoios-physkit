# Module `projectkoios.physkit.linalg.matrix.norm.maximum_absolute_entry`

## Responsibility

This module owns
[`MaximumAbsoluteEntryMatrixNorm`](MaximumAbsoluteEntryMatrixNorm/index.md), the
quantity-aware adapter for the largest absolute matrix entry. It delegates the reduction
to the exact dense or CSR Optimization owner.

## Numerical owners

- Dense real and complex quantities use `NumpyMaximumAbsoluteEntryNorm`.
- Sparse real and complex quantities use `ScipyCsrMaximumAbsoluteEntryNorm`.

The CSR owner includes implicit zeros in the mathematical representation without dense
materialization. The module contains no local reduction or fallback.

## Related pages

- [Class `MaximumAbsoluteEntryMatrixNorm`](MaximumAbsoluteEntryMatrixNorm/index.md)
- [Matrix-norm package](../index.md)
- [Base contracts](../base/index.md)
