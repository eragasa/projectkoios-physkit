# Linear Slater–Koster parameterization scientific contract

The module exposes the exact linear span of a previously admitted finite bond model. It
does not decide whether the basis, bond inventory, parameter sharing, or model class is
physically appropriate.

A zero coefficient is distinct from an absent channel: every group retains all seven
basis models, and a client supplies zero explicitly when required by its concrete model.
Parameter-group identifiers express sharing only.

The linear representation contains hopping terms only. Onsite composition, target-space
alignment, fitting objectives, weighting, train/validation policy, residuals,
identifiability decisions, uncertainty, evidence, and scientific acceptance remain
separate owners.

Software reconstruction tests establish represented linearity, not adequacy for any
material.

- [Index](index.md)
- [Numerics](numerics.md)
- [Schematics](schematics.md)
