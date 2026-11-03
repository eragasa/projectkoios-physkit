# Class `SparseBlockHoppingBlochEvaluator`

## Current contract

`SparseBlockHoppingBlochEvaluator` evaluates
$H(\mathbf q)=\sum_{\mathbf R}T(\mathbf R)e^{+2\pi i\mathbf q\cdot\mathbf R}$
from a `SparseBlockHoppingModel`. Each phase-weighted coefficient is accumulated in
CSR storage. The evaluator canonicalizes the sparse result and rejects nonfinite
arithmetic; it does not invoke the dense evaluator or silently densify.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/bloch.py::SparseBlockHoppingBlochEvaluator`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bloch/test__SparseBlockHoppingBlochEvaluator__evaluate.py`
