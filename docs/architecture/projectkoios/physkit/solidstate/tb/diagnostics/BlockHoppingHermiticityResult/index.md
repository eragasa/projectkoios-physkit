# Class `BlockHoppingHermiticityResult`

## Current contract

The result retains one `BlockHoppingModel` and its exact canonical comparisons.
Each available unordered pair `{R,-R}` appears once; zero displacement compares
with itself; a one-sided term records one missing target. `missing_targets` exposes
coverage gaps and `maximum_frobenius_defect` returns the largest available
unit-aware defect or `None` when every target is missing.

No boolean status or tolerance is owned.

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py::BlockHoppingHermiticityResult`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__BlockHoppingHermiticityResult__init.py`
