# Package `projectkoios.physkit.numerics`

## Current responsibility

The package owns reusable unit-free numerical algorithms and the shared typing
of supported numerical representations. Model layers strip units before this
boundary and restore them after numerical execution.

## Public contract

Numerical algorithms accept explicit NumPy or SciPy representations according
to their defining-module contracts. Defining modules remain the preferred
import routes.

## Navigation

- [Package `eigenproblems`](eigenproblems/index.md)
- [Package `linear_algebra`](linear_algebra/index.md)
- [Package `typing`](typing/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/numerics/`
- Tests: `tests/projectkoios/physkit/numerics/`

## Evidence

Implementation conformance is supported by the mapped numerical tests. These
tests do not establish scientific validation for every consuming model.
