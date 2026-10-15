# Class `ConjugateClosedBlockHoppingTruncator`

## Current contract

`execute(source, retained_displacements, truncated_identifier)` validates the
explicit inventory and distinct nonempty retained-model identity, canonicalizes
selection to source-owned displacement and term records, builds a model preserving
basis/unit/energy-reference/dimension metadata, and reports omitted block power.

The action enforces closure under displacement negation but does not choose the
inventory or test coefficient Hermiticity. Closure uses ordinary integer-lattice
negation, not cyclic finite-mesh equivalence. Selection criteria and any Nyquist
alias treatment remain with the consumer.

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/truncation.py::ConjugateClosedBlockHoppingTruncator`
- Tests: `tests/projectkoios/physkit/solidstate/tb/truncation/test__ConjugateClosedBlockHoppingTruncator__execute.py`
