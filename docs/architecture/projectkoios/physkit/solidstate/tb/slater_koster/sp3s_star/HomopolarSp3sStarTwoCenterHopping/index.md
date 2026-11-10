# Class `HomopolarSp3sStarTwoCenterHopping`

## Current contract

Correlates one exact two-center definition with its immutable real five-orbital energy matrix.

The public boundary accepts a hopping definition and `(5, 5)` energy-valued `MatrixQuantity`. The matrix shape and energy unit must match the definition integrals exactly. The fixed orbital order is exposed explicitly.

The central represented relation is

$$
T(d;V)=\sum_{\alpha=1}^{7}V_\alpha B_\alpha(d).
$$

## Ownership boundary

The result is a matrix representation for one directed bond. It is not a complete Hamiltonian or evidence of material adequacy. It does not retain lattice displacement, assemble reverse support, or evaluate Bloch spectra.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/sp3s_star.py::HomopolarSp3sStarTwoCenterHopping`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/sp3s_star/`
- [Owning module](../index.md)
