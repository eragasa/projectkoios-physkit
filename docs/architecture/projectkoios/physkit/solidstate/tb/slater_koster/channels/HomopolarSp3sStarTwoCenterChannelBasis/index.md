# Class `HomopolarSp3sStarTwoCenterChannelBasis`

## Current contract

Correlates one normalized bond direction with seven ordered immutable nondimensional coefficient matrices.

The public boundary accepts a `SlaterKosterBondDirection` and a seven-element tuple of `MatrixQuantity` values. Every matrix has shape `(5, 5)` and the exact nondimensional unit. Lookup uses the stable enum order.

The central represented relation is

$$
\mathcal B_{\mathrm{channel}}(d)=(B_1(d),\ldots,B_7(d)),\qquad B_\alpha\in\mathbb R^{5\times5}.
$$

## Ownership boundary

The result is an angular linear basis. It contains no energy coefficients and makes no claim about fit identifiability. It does not construct an energy-valued hopping matrix or choose parameter values.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/channels.py::HomopolarSp3sStarTwoCenterChannelBasis`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/channels/`
- [Owning module](../index.md)
