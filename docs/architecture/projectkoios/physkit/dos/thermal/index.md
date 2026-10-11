# Module `projectkoios.physkit.dos.thermal`

## Current responsibility

This module owns Fermi--Dirac transforms of an already normalized Gaussian or
histogram DOS. It returns occupations, occupied and unoccupied DOS, the thermal
window, thermally active DOS, and numerical occupied/unoccupied state-density
integrals.

The module does not solve band energies, choose state degeneracy, or infer a
particle-density target. Chemical-potential inversion uses exact sampled states
through `BaseUniformBandDensityOfStates.solve_chemical_potential(...)`; thermal
DOS transforms use the requested rendered DOS representation.

## Public contracts

- [`FermiDiracThermalDensityOfStatesEvaluator`](FermiDiracThermalDensityOfStatesEvaluator/index.md)
- [`FermiDiracThermalDensityOfStatesRequest`](FermiDiracThermalDensityOfStatesRequest/index.md)
- [`FermiDiracThermalDensityOfStates`](FermiDiracThermalDensityOfStates/index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/dos/thermal.py`
- Tests: `tests/projectkoios/physkit/dos/thermal/`
- Laboratory: `notebooks/solidstate/bloch1d/cosinepotential/cosine-potential-bands.ipynb`

## Navigation

- [DOS package](../index.md)
- [Numeric contract](../numeric.md)
- [Scientific contract](../scientific.md)
