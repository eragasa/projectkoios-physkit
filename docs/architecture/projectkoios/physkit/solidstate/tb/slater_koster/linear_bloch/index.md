# Module `projectkoios.physkit.solidstate.tb.slater_koster.linear_bloch`

This module evaluates every linear Slater–Koster channel model on one explicit ordered
reciprocal sample.

## Current contract

- `HomopolarSp3sStarLinearBlochEvaluator` delegates the complete compatible model family
  to PhysKit's dense block-hopping Bloch evaluator.
- `HomopolarSp3sStarLinearBlochEvaluation` retains matrices in parameter-group, channel,
  then sample order.
- Every matrix is nondimensional and acts on the exact parameterization basis.
- Exact lookup requires a declared group identifier, a channel enum member, and a
  built-in integer sample index.

## Class architecture

- [`HomopolarSp3sStarLinearBlochEvaluation`](HomopolarSp3sStarLinearBlochEvaluation/index.md)
- [`HomopolarSp3sStarLinearBlochEvaluator`](HomopolarSp3sStarLinearBlochEvaluator/index.md)

## Exclusions

The module does not flatten matrices, construct a regression design, define residuals,
select reciprocal samples, fit parameters, diagonalize matrices, or decide adequacy.

- [Numerics](numerics.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
- [Parent package](../index.md)
