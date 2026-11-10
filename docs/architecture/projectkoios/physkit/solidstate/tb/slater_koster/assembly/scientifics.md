# Slater–Koster assembly scientific contract

The module constructs a finite matrix representation from an already selected basis,
parameter inventory, and directed-bond inventory. It does not choose those scientific
inputs.

The ordinary-cell coefficients satisfy the generated Hermitian-support relation

$$
H_{-R}=H_R^{\dagger}
$$

when the supplied representative inventory is complete under the assembly rule. This is
a representation property, not evidence that the chosen physical model is adequate.

The integer displacement controls Bloch phase, while explicit center coordinates remain
separate. A full-center representation is obtained only through

$$
H_{\mathrm{full}}(q)=D(q)^{\dagger}H_{\mathrm{cell}}(q)D(q),
$$

where $D(q)$ is owned by `orbital_center_gauge`. Assembly never folds
$\tau_t-\tau_s$ into $R$.

Parameter groups state sharing only. Directed-bond uniqueness prevents double counting,
but it does not establish that a neighbor list is crystallographically complete or that
all physically relevant interactions are present.

Onsite operators, material composition, neighbor generation, distance scaling, fitting,
residuals, uncertainty, validation evidence, and acceptance remain separate owners.

- [Index](index.md)
- [Schematics](schematics.md)
- [Numerics](numerics.md)
