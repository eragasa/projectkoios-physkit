# Class `HomopolarSp3sStarCellGaugeModelAssembler`

## Current contract

Embeds directed five-orbital hopping representatives and generated Hermitian reverses into sorted dense ordinary-cell coefficients.

The public boundary accepts one keyword-only `HomopolarSp3sStarCellGaugeModelDefinition`. The action resolves groups exactly, accumulates into complex128 matrices, raises on nonfinite accumulation, and sorts displacement terms lexicographically.

The central represented relation is

$$
H_{R_b}\mathrel{+}=E_{s_b}T_bE_{t_b}^{\mathsf T},\qquad H_{-R_b}\mathrel{+}=E_{t_b}T_b^{\mathsf T}E_{s_b}^{\mathsf T}.
$$

## Ownership boundary

Assembly preserves declared representation conventions. It does not choose bonds, add onsite physics, fit coefficients, or establish model adequacy. It does not evaluate Bloch phases, diagonalize matrices, or construct linear fitting designs.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/assembly.py::HomopolarSp3sStarCellGaugeModelAssembler`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/assembly/test__HomopolarSp3sStarCellGaugeModelAssembler__execute.py`
- [Owning module](../index.md)
