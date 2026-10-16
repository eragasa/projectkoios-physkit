# Class `BlockHoppingRangeBandSelection`

This page defines one explicit retained inventory requested within an ordered band
family.

## Current contract

The record carries only caller-selected ordinary-lattice displacements and the new
retained-model identity.

- `retained_displacements` is nonempty, sorted, unique, and dimensionally uniform.
- The inventory is closed under ordinary integer negation.
- `truncated_identifier` is a nonempty string.
- Family evaluation requires inventories and identifiers to be unique across entries.
- Source membership and metadata compatibility are checked by the truncator.
- No range metric, distance, cutoff, cyclic equivalence, or Nyquist rule is inferred.

## Mapping

The implementation and invariant tests are maintained together.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bands.py::BlockHoppingRangeBandSelection`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingRangeBandSelection__init.py`
