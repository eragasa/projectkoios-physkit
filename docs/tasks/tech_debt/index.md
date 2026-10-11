# PhysKit technical-debt tasks

This directory records bounded follow-up work that is known but is not part of
the current capability tranche. A task records the problem, intended contract,
migration constraints, and completion criteria. It does not itself establish
scientific validation or authorize implementation beyond an explicitly approved
scope.

## Open tasks

- [Backend-neutral linear-algebra solver boundary](backend-neutral-linear-algebra-solvers.md)
  — separate model code from NumPy dense and SciPy sparse solver APIs while
  retaining explicit storage and densification policy.
