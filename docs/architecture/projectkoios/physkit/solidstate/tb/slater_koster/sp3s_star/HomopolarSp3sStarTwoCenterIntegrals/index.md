# Class `HomopolarSp3sStarTwoCenterIntegrals`

## Current contract

Retains seven real physical energy quantities in stable homopolar sp3s-star channel order.

The public boundary accepts seven `ScalarQuantity` values and one explicit non-nondimensional `UnitSystem`. Every scalar must carry a physical energy unit, all units must match exactly, and that unit must equal the selected system energy unit.

The central represented relation is

$$
V=(V_{ss\sigma},V_{sp\sigma},V_{pp\sigma},V_{pp\pi},V_{s^*p\sigma},V_{ss^*\sigma},V_{s^*s^*\sigma}).
$$

## Ownership boundary

The values are represented model inputs. Valid storage does not establish fit quality, transferability, or scientific acceptance. It does not infer absent channels, convert units, or attach material or shell meaning.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/sp3s_star.py::HomopolarSp3sStarTwoCenterIntegrals`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/sp3s_star/`
- [Owning module](../index.md)
