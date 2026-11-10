# Class `HomopolarSp3sStarTwoCenterHoppingConstructor`

## Current contract

Constructs one energy-valued two-center matrix by ordered linear combination of the seven angular channel matrices.

The public boundary accepts one keyword-only `HomopolarSp3sStarTwoCenterHoppingDefinition`. Binary64 accumulation is checked for nonfinite output. Overflow raises rather than clipping or silently regularizing values.

The central represented relation is

$$
T(d;V)=\sum_{\alpha=1}^{7}V_\alpha B_\alpha(d),\qquad T(-d;V)=T(d;V)^{\mathsf T}.
$$

## Ownership boundary

The action implements the declared homopolar sign convention only. It performs no fitting, shell inference, or validation against external data. It does not add onsite terms, overlap, spin, lattice support, or Bloch phases.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/sp3s_star.py::HomopolarSp3sStarTwoCenterHoppingConstructor`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/sp3s_star/test__HomopolarSp3sStarTwoCenterHoppingConstructor__execute.py`
- [Owning module](../index.md)
