# `HomopolarSp3sStarCellGaugeModelAssembler` scientific contract

Assembly preserves declared representation conventions. It does not choose bonds, add onsite physics, fit coefficients, or establish model adequacy.

The represented relation

$$
H_{R_b}\mathrel{+}=E_{s_b}T_bE_{t_b}^{\mathsf T},\qquad H_{-R_b}\mathrel{+}=E_{t_b}T_b^{\mathsf T}E_{s_b}^{\mathsf T}.
$$

is a mathematical or structural software contract. Passing construction and tests
establishes conformance to that contract only; it does not establish material accuracy,
parameter acceptance, convergence, independent verification, or experimental agreement.

It does not evaluate Bloch phases, diagonalize matrices, or construct linear fitting designs.

- [Index](index.md)
- [Schematics](schematics.md)
- [Numerics](numerics.md)
