# Class `HomopolarSp3sStarBondGeometry`

## Current contract

Binds a directed ordinary-cell bond to its nonzero physical Cartesian displacement and selected physical unit system.

The public boundary accepts a `SlaterKosterDirectedBond`, physical length `VectorQuantity`, and non-nondimensional `UnitSystem`. Length dimensionality, exact selected length unit, shape `(3,)`, and nonzero norm are validated before direction derivation.

The central represented relation is

$$
d=\frac{\Delta r}{\lVert\Delta r\rVert_2},\qquad [\Delta r]=U_{\mathrm{system}}^{\mathrm{length}}.
$$

## Ownership boundary

The geometry separates R from the physical vector. A client remains responsible for correlating R, fractional centers, and lattice vectors. It does not generate neighbors, verify crystallographic completeness, or apply distance scaling.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/bond.py::HomopolarSp3sStarBondGeometry`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/bond/test__HomopolarSp3sStarBondGeometry__init.py`
- [Owning module](../index.md)
