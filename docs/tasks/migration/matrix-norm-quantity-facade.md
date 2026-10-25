# Matrix-quantity norm facade

**Status:** Implemented and qualified locally; uncommitted review pending.

## Purpose

PhysKit owns matrix quantities and restores their units around unit-free matrix-norm
algorithms from `projectkoios-optimization`. Raw NumPy and SciPy evaluation remains a
direct Optimization API and is not wrapped by PhysKit.

This design is independent of legacy KSD and PhysKit API shapes. Donor behavior supplies
provenance and verification context, not compatibility requirements.

## Package layout

```text
src/python/projectkoios/physkit/linalg/matrix_norms/
├── __init__.py
├── base.py
├── frobenius.py
├── maximum_absolute_entry.py
└── spectral.py
```

- `base.py` owns `BaseMatrixQuantityNorm` and the `MatrixNorm` facade. It imports the
  exact Optimization `MatrixNormEnum`; it does not redefine that semantic identity.
- `frobenius.py` owns `FrobeniusMatrixNorm`.
- `maximum_absolute_entry.py` owns `MaximumAbsoluteEntryMatrixNorm`.
- `spectral.py` owns `SpectralMatrixNorm` and `SparseSpectralNormOptionsEnum`.

No compatibility module, factory module, helper module, or utilities container is
introduced. Every helper behavior is an instance method of one of these cohesively owned
classes. The modules contain no behavioral module-level functions, static utility
methods, or classmethod helpers.

## Accepted quantity representations

Every PhysKit norm accepts only an existing nominal matrix quantity:

```python
MatrixQuantity | ComplexMatrixQuantity | SparseMatrixQuantity | ComplexSparseMatrixQuantity
```

No new union alias is introduced. A quantity may carry either `PhysicalUnit` or
`Unitless`; the exact unit object is preserved. Raw ndarray and CSR inputs are outside
the PhysKit facade and use Optimization directly.

## Public identities

`projectkoios.physkit.linalg.matrix_norms.__all__` deliberately exposes:

- `BaseMatrixQuantityNorm`
- `MatrixNorm`
- exact Optimization `MatrixNormEnum`
- `FrobeniusMatrixNorm`
- `MaximumAbsoluteEntryMatrixNorm`
- `SpectralMatrixNorm`
- `SparseSpectralNormOptionsEnum`
- exact Optimization `ScipyCsrSpectralNormControls`

`projectkoios.physkit.linalg.__all__` exposes the KSD-facing facade identities:

- `MatrixNorm`
- exact Optimization `MatrixNormEnum`
- `SparseSpectralNormOptionsEnum`
- exact Optimization `ScipyCsrSpectralNormControls`

Optimization's `MatrixNormABC` and concrete backend classes retain their defining
Optimization imports. PhysKit imports them for implementation and typing but does not
re-export them as convenience aliases.

The supported surface contains no string selector, mapping accessor, `ncv` alias,
compatibility alias, callable evaluator shortcut, or separate factory.

## Base contract

`BaseMatrixQuantityNorm` is the nominal PhysKit base for one quantity-aware norm. It
retains an exact matrix quantity, exposes the canonical `matrix_norm` identity, and
returns `ScalarQuantity` from `.result()`.

The base owns common instance behavior for:

- retaining the exact quantity;
- distinguishing dense from CSR quantity representations;
- extracting the exact dense magnitude;
- obtaining a fresh canonical CSR value from sparse quantities;
- delegating magnitude evaluation to a concrete class; and
- restoring the exact original unit object.

It does not own a numerical kernel, sparse algorithm selection, scientific threshold,
or exception translation.

Optimization's exact `MatrixNormABC` remains the separate unit-free numerical contract.
Its `.result()` returns a built-in `float`; PhysKit neither subclasses nor redefines that
identity.

Underscore attributes and methods are internal/non-public conventions, not claims of
Python-enforced privacy.

## Concrete quantity norms

### `FrobeniusMatrixNorm`

- identifies `MatrixNormEnum.frobenius`;
- dense quantities delegate to `NumpyFrobeniusNorm(matrix.magnitude).result()`;
- sparse quantities delegate to `ScipyCsrFrobeniusNorm(matrix.to_csr()).result()`;
- never densifies sparse storage; and
- returns the matrix-entry unit, not its square.

Frobenius Power remains a separate capability outside this facade.

### `MaximumAbsoluteEntryMatrixNorm`

- identifies `MatrixNormEnum.maximum_absolute_entry`;
- dense quantities delegate to `NumpyMaximumAbsoluteEntryNorm`;
- sparse quantities delegate to `ScipyCsrMaximumAbsoluteEntryNorm`;
- preserves implicit sparse-zero semantics; and
- never densifies sparse storage.

### `SpectralMatrixNorm`

The class accepts the retained quantity plus:

```python
sparse_spectral_option: SparseSpectralNormOptionsEnum | None
sparse_spectral_controls: ScipyCsrSpectralNormControls | None
```

Dense quantities require both sparse-only arguments to be absent and delegate to
`NumpySpectralNorm`.

Sparse quantities require one explicit option:

- `complete` forbids iterative controls, explicitly materializes the fresh canonical CSR
  value, checks materialization arithmetic, and delegates to `NumpySpectralNorm`;
- `iterative` requires exact Optimization controls and delegates to
  `ScipyCsrSpectralNorm(matrix, controls).result()` without dense conversion or fallback.

Optimization's sparse evaluator uses its deterministic nondensifying direct branch when
`min(shape) <= 2` and approximate ARPACK iteration when `min(shape) >= 3`. That branch
is part of the selected owner evaluator, not a third PhysKit option or a fallback.

## Sparse option identity

```python
class SparseSpectralNormOptionsEnum(Enum):
    complete = auto()
    iterative = auto()
```

The pure enum is PhysKit-owned because it selects a quantity-representation path. Its
values are opaque. It does not inherit from `str`, and no string-to-enum shortcut is
supported.

The exact Optimization control field is `krylov_dimension`. `ncv` appears only at
Optimization's internal SciPy call boundary and is never a PhysKit API name.

## `MatrixNorm` facade

```python
MatrixNorm(
    matrix: MatrixQuantity
    | ComplexMatrixQuantity
    | SparseMatrixQuantity
    | ComplexSparseMatrixQuantity,
    matrix_norm: MatrixNormEnum,
    *,
    sparse_spectral_option: SparseSpectralNormOptionsEnum | None = None,
    sparse_spectral_controls: ScipyCsrSpectralNormControls | None = None,
).result() -> ScalarQuantity
```

The immutable slotted facade owns immutable internal/non-public mappings from:

- canonical norm enum to stable diagnostic name; and
- canonical norm enum to exact quantity-aware implementation class.

The mappings use `MappingProxyType`, are not exports, and have no accessor.
`MatrixNormEnum.name` is the only supported textual semantic identity.

The facade validates the canonical enum and delegates to one quantity-aware norm class.
It performs no numerical calculation, unit conversion, dense conversion, retry,
fallback, result caching, or exception translation.

## Validation and failure ownership

- Wrong nominal quantity, enum, option, or controls identities raise `TypeError`.
- Missing, incompatible, or ignored sparse spectral options raise `ValueError`.
- Backend dtype, dimensionality, emptiness, and source-finiteness failures retain their
  Optimization classification.
- Nonfinite explicit dense materialization, backend output, or iterative residual
  arithmetic raises `FloatingPointError` at the owning operation.
- Dense or iterative decomposition nonconvergence remains `RuntimeError` from the
  numerical owner.
- Unrelated backend defects propagate unchanged.

No PhysKit layer catches a numerical-owner exception to retry, select another algorithm,
or provide compatibility behavior.

## Unit contract

All three norms have the same unit as each matrix entry:

```python
ScalarQuantity(owner_result, matrix.unit)
```

The exact input unit object is retained whether it is `PhysicalUnit` or `Unitless`. No
matrix values are converted to SI or another scale.

## Dependency boundary

Selected reviewed Optimization revision:

- commit `2a18e485109e307b2f7928bd0418bcf796aeb0d5`;
- tree `77262fed5d75e7ddf01334321f9fb24bf6d5961e`;
- remotely advertised on `origin/feature/scipy-csr-spectral-norm` at review time.

The implementation requests:

```text
projectkoios-optimization[numpy,scipy] @ ...@2a18e485109e307b2f7928bd0418bcf796aeb0d5
scipy>=1.18,<2
```

The exact pin is used directly; no unreviewed merge revision is substituted.

## Verification

```text
tests/projectkoios/physkit/linalg/matrix_norms/
├── test__matrix_norms_public_api.py
├── test__BaseMatrixQuantityNorm.py
├── test__MatrixNorm.py
├── test__FrobeniusMatrixNorm.py
├── test__MaximumAbsoluteEntryMatrixNorm.py
├── test__SpectralMatrixNorm__dense.py
└── test__SpectralMatrixNorm__sparse.py
```

Verification covers:

- exact Optimization enum and controls identity reuse;
- deliberate package exports and absence of backend re-exports;
- absence of mapping accessors, string selectors, raw-input paths, and `ncv`;
- nominal base relationships and immutable facade state;
- real/complex dense and sparse quantities with physical or unitless units;
- analytic values for every norm;
- exact matrix quantity, magnitude, unit, enum, option, and controls retention;
- direct dense delegation and canonical sparse delegation;
- sparse Frobenius/maximum and iterative spectral non-densification;
- complete spectral's sole explicit dense materialization boundary;
- iterative controls forwarding, deterministic small-shape behavior, nonconvergence, and
  no fallback;
- empty and unrepresentable-result owner boundaries, while nonfinite and noncanonical
  matrix construction remains owned by the quantity classes; and
- unchanged existing block-diagnostic behavior.

Qualification completed with 38 focused tests and the full suite at 1545 passed plus
20 passed subtests. Ruff checks and formatting, compileall, exact lock verification,
`git diff --check`, fresh-sdist wheel construction, metadata/payload inspection, and an
isolated installed-wheel dense/complete-sparse/iterative-sparse smoke test passed.
PhysKit has no configured strict static checker; that limitation is reported explicitly.

## Documentation and provenance completion

Current-state documentation mirrors each package, module, and class under:

```text
docs/architecture/projectkoios/physkit/linalg/matrix_norms/
```

Implementation also adds a KSD prior-art provenance page without promising API
compatibility and updates exact Optimization revision references.

## Explicit exclusions

The capability does not provide or preserve:

- raw ndarray or CSR input through PhysKit;
- manufactured `Unitless` wrappers for raw input;
- KSD `RepresentedMatrixNormAnalyzer` or aggregated result;
- reference normalization;
- compatibility aliases or KSD re-exports;
- generic factories;
- module-level helper functions;
- static or classmethod helpers;
- Helper/Utils containers;
- structural protocols;
- string selectors;
- implicit sparse spectral selection;
- public `ncv`;
- iterative-to-complete fallback;
- sparse densification outside explicit `complete`; or
- scientific acceptance, uncertainty, serialization, or workflow policy.

## Sequencing

The preceding Frobenius integration is committed locally as
`025e0d5f16899c85b6fa3bfb034972343e0e98dc` and has not been pushed.

1. Complete and review source-mirrored architecture documentation.
2. Implement the package exactly as documented.
3. Repin and regenerate the lock.
4. Add mapped tests and provenance.
5. Run complete qualification.
6. Stop with the facade diff uncommitted for review.

No push, second commit, release, tag, merge, or KSD mutation is authorized.
