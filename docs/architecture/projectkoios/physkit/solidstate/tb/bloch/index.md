# Module `projectkoios.physkit.solidstate.tb.bloch`

## Current responsibility

The module owns physics-semantic bridges from ordered `ScalarHoppingModel` and
`BlockHoppingModel` records to scalar or matrix-valued Bloch Hamiltonians. For
unitless fractional reciprocal coordinates $\mathbf q$ and integer direct-lattice
displacements $\mathbf R$, it evaluates

$$
h(\mathbf q)=\sum_{\mathbf R} t_{\mathbf R}
\exp\!\left(+2\pi i\,\mathbf q\cdot\mathbf R\right)
$$

or the entrywise block analogue with $t_{\mathbf R}$ replaced by
$T_{\mathbf R}$.

The positive exponent, factor of $2\pi$, supplied hopping-term order, and absence
of normalization are fixed conventions. Models own coefficient order, energy
unit, energy-reference identity, and basis identity. Supplied wave vectors are
rows of a unitless `MatrixQuantity` with width equal to the model dimension.

Results remain complex because hopping models do not impose Hermiticity. The
block result retains one immutable energy-unit-bearing `ComplexMatrixQuantity`
per supplied reciprocal point. Diagonalization, diagnostics, tolerances,
acceptance, and interpretation remain separate responsibilities.

## Numerics boundary

The scalar adapter creates one `FiniteFourierSeriesEvaluationDefinition`; the
block adapter creates one definition per matrix entry using shared coordinates
and displacements. Each calls a nominal `FiniteFourierSeriesEvaluator` and wraps
the returned built-in complex values in the model energy unit. The dependency is
bound to exact Numerics owner commit
`cf551409872a621b4e47d2236ac4ffdc12c187c0`.

The upstream Numerics operation remains unit-free and does not receive lattice,
hopping, Bloch, basis, energy-reference, or acceptance meaning. PhysKit is not a
pass-through facade: these adapters restore physical semantics around the generic
finite-series calculation.

The default adapters construct the reviewed ordered NumPy backend through the
Numerics factory. A nominal evaluator may be injected for verification or a
separately reviewed implementation. Empty wave-vector matrices of shape
`(0, dimension)` produce an empty scalar vector or empty block matrix inventory.
Numerics validation and numerical failures propagate without broad translation.

## Public classes

- [`ScalarHoppingBlochEvaluator`](ScalarHoppingBlochEvaluator/index.md)
- [`ScalarHoppingBlochEvaluation`](ScalarHoppingBlochEvaluation/index.md)
- [`BlockHoppingBlochEvaluator`](BlockHoppingBlochEvaluator/index.md)
- [`BlockHoppingBlochEvaluation`](BlockHoppingBlochEvaluation/index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/bloch.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bloch/`
