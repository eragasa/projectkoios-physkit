# Module `projectkoios.physkit.solidstate.tb.truncation`

## Current responsibility

The module applies an exact caller-supplied block-hopping displacement selection.
The retained inventory must be nonempty, sorted, unique, contained in the source,
and closed under ordinary integer-lattice $R\mapsto-R$. The caller supplies a
nonempty retained-model identifier distinct from the source identifier. The action
reuses selected source terms; preserves dimension, basis, energy unit, and energy
reference; and reports source power outside the selection.

Conjugate closure is structural only. The module does not claim that paired
coefficients obey $T[-R]=T[R]^\dagger$; that remains a separate diagnostic.
There is no cyclic modulus or finite-DFT aliasing: an even-mesh Nyquist
representative is not self-paired. Consumers must exclude it from an ordinary
closed retained subset or keep cyclic selection separate. The module does not
choose a range, metric, cutoff, tolerance, normalized loss, campaign inventory,
or acceptance status.

## Navigation

- [`BlockHoppingTruncationResult`](BlockHoppingTruncationResult/index.md)
- [`ConjugateClosedBlockHoppingTruncator`](ConjugateClosedBlockHoppingTruncator/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/truncation.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/truncation/`

## Evidence boundary

Tests cover metadata and source-term identity preservation, canonical displacement
inventories, conjugate closure, omitted power, and invalid selections. They do not
establish that any chosen truncation is scientifically adequate.
