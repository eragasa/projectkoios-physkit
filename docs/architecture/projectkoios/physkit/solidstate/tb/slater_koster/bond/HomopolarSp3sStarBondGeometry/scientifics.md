# `HomopolarSp3sStarBondGeometry` scientific contract

The geometry separates R from the physical vector. A client remains responsible for correlating R, fractional centers, and lattice vectors.

The represented relation

$$
d=\frac{\Delta r}{\lVert\Delta r\rVert_2},\qquad [\Delta r]=U_{\mathrm{system}}^{\mathrm{length}}.
$$

is a mathematical or structural software contract. Passing construction and tests
establishes conformance to that contract only; it does not establish material accuracy,
parameter acceptance, convergence, independent verification, or experimental agreement.

It does not generate neighbors, verify crystallographic completeness, or apply distance scaling.

- [Index](index.md)
- [Schematics](schematics.md)
- [Numerics](numerics.md)
