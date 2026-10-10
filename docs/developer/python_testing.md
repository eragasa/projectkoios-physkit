# Python testing conventions

## Mirrored module layout

Tests mirror the implementation module as a directory. A test for one class
method uses:

```text
tests/projectkoios/physkit/<subpackage>/<module>/test__<ClassName>__<method>.py
```

Examples:

```text
src/python/projectkoios/physkit/units/systems.py
tests/projectkoios/physkit/units/systems/test__UnitSystem__energy_unit.py

src/python/projectkoios/physkit/qm/piab1d/tise_fd/solver.py
tests/projectkoios/physkit/qm/piab1d/tise_fd/solver/
    test__Piab1dTiseFdSolver__solve.py
```

Use `init` for constructor tests and the property name for property tests. New
or materially changed top-level test functions use
`test__surface__facet__behavior`. Keep unrelated classes and methods in separate
files. Shared fixtures may live in the mirrored module directory's `conftest.py`.

All maintained PhysKit package tests, including their generated cache location,
belong under `tests/projectkoios/physkit/`. Repository-level conformance tests
remain under `tests/repository/`. Do not retain or recreate a parallel
`tests/physkit/` tree. The repository-layout test rejects that legacy test root
so test collection and bytecode caches cannot obscure the canonical mirrored
namespace.

## Test scope

Each test should state one observable requirement through its name and assertions.
Verify nominal type ownership, dimensional units, immutable storage, shape
correlation, numerical values, and failure behavior where applicable.

For migrated numerical code, retain the donor's exact represented formulas and
add selected-scale checks. In particular, physical quantum tests must verify that
matrices remain in the requested numerical scale rather than merely being
dimensionally convertible through SI.

## Execution

Run the narrowest relevant test directory while developing, followed by the
connected numerical layer. Examples:

```bash
python -m pytest -q tests/projectkoios/physkit/qm/piab1d
python -m pytest -q tests/projectkoios/physkit/operators/operators_1d/fd
python -m pytest -q tests/projectkoios/physkit/operators/qm/qm_1d
python -m pytest -q tests/projectkoios/physkit/numerics/eigenproblems/real_symmetric
```

Also compile changed Python modules and run `git diff --check`.

## Deprecated tests

Tests retained specifically for `projectkoios.physkit.legacy` use the registered
`@pytest.mark.deprecated` marker, either directly or through module-level
`pytestmark`. Maintained implementation tests must not carry this marker or
import deprecated modules.

## Evidence boundaries

Passing tests establishes only the requirements actually asserted. Keep software
verification, numerical verification, physical validation, pedagogical
validation, and uncertainty quantification distinct. Executing a notebook or
reproducing a finite matrix identity does not by itself validate a physical
model.
