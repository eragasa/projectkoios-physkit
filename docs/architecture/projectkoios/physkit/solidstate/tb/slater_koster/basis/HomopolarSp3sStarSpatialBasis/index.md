# Class `HomopolarSp3sStarSpatialBasis`

## Current contract

Admits only exact ordered five-coordinate homopolar sp3s-star blocks for every center.

The public boundary accepts one keyword-only `SlaterKosterSpatialBasisDefinition`. The complete coordinate tuple must equal the generated expected tuple. Reordering, omission, insertion, or component substitution fails closed.

The central represented relation is

$$
\mathcal B_c=((c,s,0),(c,p_x,0),(c,p_y,0),(c,p_z,0),(c,s,1)),\qquad N_b=5N_c.
$$

## Ownership boundary

Admission establishes software-level basis structure only. It does not establish that a material is accurately represented by sp3s-star orbitals. It does not admit heteropolar semantics, onsite values, spin, overlap, or hopping parameters.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/basis.py::HomopolarSp3sStarSpatialBasis`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/basis/test__HomopolarSp3sStarSpatialBasis__init.py`
- [Owning module](../index.md)
