# Class `HomopolarSp3sStarTwoCenterChannelBasisConstructor`

## Current contract

Constructs the seven nondimensional angular matrices in stable channel order from one normalized direction.

The public boundary accepts one keyword-only `SlaterKosterBondDirection`. The action uses exact channel dispatch and binary64 arrays; each branch writes only the entries owned by that channel.

The central represented relation is

$$
B_\alpha(d)=\frac{\partial T(d;V)}{\partial V_\alpha}.
$$

## Ownership boundary

The constructor implements a stated sign and orbital-order convention. It does not select a universal Slater-Koster convention for other models. It does not multiply by energies, fit coefficients, or assemble lattice blocks.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/channels.py::HomopolarSp3sStarTwoCenterChannelBasisConstructor`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/channels/test__HomopolarSp3sStarTwoCenterChannelBasisConstructor__execute.py`
- [Owning module](../index.md)
