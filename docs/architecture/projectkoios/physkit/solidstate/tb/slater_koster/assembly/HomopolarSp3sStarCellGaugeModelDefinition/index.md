# Class `HomopolarSp3sStarCellGaugeModelDefinition`

## Current contract

Correlates one exact structured basis, unique parameter groups, unique directed representatives, model identity, and energy reference.

The public boundary accepts a keyword-only model identifier, admitted basis, group tuple, bond tuple, and energy-reference string. Group identifiers and directed bonds are unique, explicit reverses are absent, endpoints belong to the basis, units agree, and every group is referenced.

The central represented relation is

$$
\mathcal D_{\mathrm{cell}}=(\mathcal B,\{V_g\},\{b\},E_{\mathrm{ref}}).
$$

## Ownership boundary

The definition records a selected finite model inventory. Its validity does not establish neighbor completeness, material fidelity, or fit acceptance. It does not generate geometry, onsite terms, parameters, residuals, or evidence.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/assembly.py::HomopolarSp3sStarCellGaugeModelDefinition`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/assembly/`
- [Owning module](../index.md)
