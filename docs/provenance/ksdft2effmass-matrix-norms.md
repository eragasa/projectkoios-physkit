# KSD matrix-norm prior art and PhysKit adaptation

## Donor identity

- Repository: `https://github.com/eragasa/ksdft2effmass`
- Commit: `de61df02e88bf85a57f758ea8cf13c37b34111df`
- Tree: `5b550574a0ef5fc2b1d2e5bde5a9233b5b2d8f9a`
- Source path: `python/src/ksdft2effmass/operators/matrix_norms.py`
- Source blob: `dee104bb95d2ba1f8c4eb4c10b2c4f7553bce432`
- Source SHA-256:
  `a2e8606c459d49a03f85cee46be0ae76135bd33f19a6dcb56255caccc472f868`
- Test path:
  `python/tests/software_verification/ksdft2effmass/operators/test__RepresentedMatrixNormAnalyzer.py`
- Test blob: `59c9d551a36771739f36b9ddc89cd946fb9c985f`
- Test SHA-256:
  `5fb730af47807220bf8a5e08ca8bae273a63057e3846c48ce2f4b47cdf4f6dc0`
- License: Apache-2.0

## Numerical owner identity

Unit-free algorithms are consumed rather than copied from
`projectkoios-optimization`:

- Commit: `2a18e485109e307b2f7928bd0418bcf796aeb0d5`
- Tree: `77262fed5d75e7ddf01334321f9fb24bf6d5961e`

The consumed identities are the canonical matrix-norm enum and dense/CSR Frobenius,
maximum-absolute-entry, and spectral evaluators. Iterative controls and convergence
qualification remain Optimization-owned.

## PhysKit destination

```text
src/python/projectkoios/physkit/linalg/matrix/norm/
├── base.py
├── frobenius.py
├── maximum_absolute_entry.py
└── spectral.py
```

PhysKit adapts only reusable quantity semantics:

- nominal real/complex dense/sparse matrix quantities;
- exact unit restoration as `ScalarQuantity`;
- explicit complete or iterative sparse spectral selection; and
- an enum-dispatch facade over separately owned quantity-aware norms.

## Deliberate adaptations

The PhysKit API is designed independently rather than preserving the donor surface.

- Raw numerical matrices are excluded; callers use Optimization directly.
- Complex dense and sparse quantities are supported.
- Sparse complete and iterative spectral behavior is selected explicitly.
- Sparse complete evaluation is the only authorized dense materialization path.
- Iterative sparse evaluation cannot densify or fall back.
- The exact public iterative control name is `krylov_dimension`; no `ncv` alias exists.
- Each norm is implemented in its own module behind a common nominal quantity base.
- Numerical exceptions propagate from their owning operation without compatibility
  translation.

## Excluded donor behavior

PhysKit does not migrate or preserve:

- `RepresentedMatrixNormAnalyzer`;
- the aggregated `RepresentedMatrixNormResult` record;
- `normalized_by` reference ratios;
- legacy result-field names;
- campaign tolerances or acceptance;
- serialization or workflow behavior; or
- any KSD wrapper or re-export.

Those exclusions prevent donor API compatibility from constraining the reusable owner
contract. Passing software tests establishes implementation conformance only, not
scientific validation or human acceptance.
