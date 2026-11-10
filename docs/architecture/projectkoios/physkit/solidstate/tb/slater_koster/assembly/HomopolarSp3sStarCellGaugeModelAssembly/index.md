# Class `HomopolarSp3sStarCellGaugeModelAssembly`

## Current contract

Correlates one exact structured model definition with its generated dense ordinary-cell block-hopping model.

The public boundary accepts a cell-gauge model definition and exact `DenseBlockHoppingModel`. Model identity, dimension, energy unit, energy reference, basis metadata, and basis size must agree with the retained definition.

The central represented relation is

$$
\mathcal A=(\mathcal D_{\mathrm{cell}},\{H_R\}_{R\in\mathbb Z^3}).
$$

## Ownership boundary

The retained definition remains necessary to interpret the structured basis; generic model metadata alone is not a substitute for exact basis compatibility. It does not verify scientific adequacy, Bloch spectra, or external provenance.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/assembly.py::HomopolarSp3sStarCellGaugeModelAssembly`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/assembly/`
- [Owning module](../index.md)
