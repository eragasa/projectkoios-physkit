# KSD model-system quantity extraction

## Frozen donor

The canonical PhysKit unit surface reconciles established PhysKit ownership with frozen
KSD quantity evidence.

- Donor repository: `https://github.com/eragasa/ksdft2effmass`
- Donor commit: `8aea8232a2d4bc33af095cac9a32bd3a564a7cee`
- Donor tree: `562a5d4bdd97681c23cdff37d4f15e2a162a3e95`
- Donor source path: `python/src/ksdft2effmass/operators/quantities.py`
- Donor source blob: `d98be473e19e17a59563fc7c3e00cca02fd822d3`
- Donor source SHA-256:
  `132ab2efb5dac7c56ddd5d199b75c8d22dab39941e28244ac88e106aca0db9b5`

## Donor verification evidence

The extraction was checked against the donor's focused unit and quantity tests.

- `test__PhysicalUnit.py`: blob `d987be8bb3de73198dac5ff6a715e2dc8474d5de`,
  SHA-256 `3fdbefe5084a317ca5bae239d6399ba2845d4c06c418ea6305f33f8e39dfc052`
- `test__PintUnitConverter.py`: blob `1875bdd1e764bb21b63fdf2fc700df0294a0fe42`,
  SHA-256 `953ee4bb805eb1046afed09efc2437a6272178afb266ffb2d8f03e8e5ffa17b0`
- `test__ScalarQuantity.py`: blob `7864cd1702ea54ace09851f9f0cbcb5309c20b57`,
  SHA-256 `286a4b9803374aa507fd2f7a19f2aa34e04f6b9f790abb58b2662f2be370b377`
- `test__VectorQuantity.py`: blob `5a0cf7f4569d88496f78dfca4987b4674987045a`,
  SHA-256 `fb74d74a75629139fc3ef62c94e6b1c2113bf60e76c9cd7bcdea748d3ce28868`
- `test__MatrixQuantity.py`: blob `de572616e0488f0e4d97404077911051c7461be0`,
  SHA-256 `675c8bd4072a1f17bd2c0a681b2e74cb3368d4f9d600a7c6ec3f251f052bd2e3`
- Operator complex/sparse evidence: blobs
  `2e81f0cc1189a7d818c17c25637be408ec286d1e`,
  `b23ed60ae2f58cd89e95b7aed542a21d1be4c3ba`, and
  `804b488fc235f42a897f3f6e70bcfdbfee3006d6`.

## Canonical owner contract

PhysKit retains and strengthens the reusable quantity semantics.

- `PhysicalUnit`, `Unitless`, nominal unit/quantity bases, immutable scalar and dense or
  sparse array records, the Pint registry, and conversion operations belong to PhysKit.
- Real and complex dense/sparse inputs retain binary64 or complex128 storage and explicit
  physical or unitless identity.
- CSR quantities retain canonical ordering, eliminate duplicate and explicit-zero
  entries, and do not densify implicitly.
- The fixed Pint registry is a public identity shared by every owner converter.
- `ModelSystemUnitConverter` is the canonical shared converter instance.
- Scalar and real-vector conversion supports affine units pointwise.
- Linear-algebra and sparse quantities accept only multiplicative conversions because a
  nonzero offset is not a linear map and would destroy sparse structure.
- `conversion_factor(...)` therefore rejects compatible conversions having nonzero
  affine offset rather than returning a misleading factor.

## Consumer audit

A read-only KSD audit bounds the eventual identity migration.

- Ninety source files and 133 test files import the local quantity family.
- Three maintained façades publicly re-export those identities.
- Seventy-five source files contain 144 exact-type and 138 nominal `isinstance` checks
  outside the quantity definitions.
- Twenty-eight source files make 130 converter calls; active production units are
  multiplicative, while the donor's factor-based affine conversion is incorrect.
- Three maintained serializers and 22 retained schema/artifact documents persist
  magnitude/unit representations.
- The separate KSD METAL calculator-conversion subsystem is not part of this owner
  extraction.

## Adaptation boundary

The destination does not preserve KSD as a second nominal owner.

- KSD quantity classes are donor and migration evidence, not runtime dependencies.
- Future KSD adoption must deliberately re-export or import the complete exact PhysKit
  family; wrappers, subclasses, structural substitutes, and piecemeal replacement remain
  unsupported.
- KSD owns its persisted schemas, result serializers, calculator-unit conventions,
  fixtures, rollback evidence, and campaign acceptance.
- A PhysKit checkpoint alone does not authorize KSD pin changes, schema changes, class
  deletion, or consumer migration.
