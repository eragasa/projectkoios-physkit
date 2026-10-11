# Weighted Fermi--Dirac and thermal-DOS extraction lineage

## Donor artifact

The reviewed source is
`notebooks/solidstate/bloch1d/cosinepotential/exploratory-bands-dos-thermal-prototype.ipynb`,
relocated byte-for-byte from
`notebooks/scratch/solidstate/bloch_1d_cosine_potential_bands_dos.ipynb`.

- Git blob: `687879784df581a4308104614ff59bc796f5e55f`
- SHA-256: `a029505c63023d0a09477d71fce8f1bb11c8275016f5de57ef9d1d07d0fde674`
- Preservation commit: `fbf7e9968671cffc18760e9a4aadaf165843d4cc`
- Extraction baseline: `5e936a2`

The artifact was removed from the maintained notebook tree after extraction and
individual review. Its exact bytes remain recoverable from the preservation
commit and Git blob above.

## Review evidence

The notebook contains 25 cells: 9 markdown cells, 1 raw lineage cell, 13
nonempty code cells, and 2 empty code cells. It stores 17 outputs, including 12
PNG displays and no saved error output.

A clean linear execution fails at the first band-function definition that uses
`m_e` as a default before the constant is defined. The self-contained final
778-line calculation cell executed successfully in isolation and reproduced its
301-by-6 band inventory, normalized histogram, fixed-filling chemical potentials,
and thermal plots. Successful execution does not validate the prototype's state
normalization.

## Exact adaptations

| Donor layer | Maintained destination |
| --- | --- |
| Stable finite-temperature occupation | Existing `FermiDiracState` |
| Fixed-number root inversion | Existing `FermiDiracChemicalPotentialSolver`, extended by `solve_weighted(...)` |
| Full-zone state weights | `BaseUniformBandDensityOfStates.solve_chemical_potential(...)` |
| Occupied and unoccupied DOS | `projectkoios.physkit.dos.thermal` |
| Thermal window $-\partial f/\partial E$ | `FermiDiracThermalDensityOfStates` |
| Thermally active DOS | `FermiDiracThermalDensityOfStates` and the generic DOS visualizer |
| Pedagogical calculation | Maintained 1D cosine-potential laboratory and lecture note |

The weighted solver retains the existing stable occupation evaluation and
bisection controls. It adds explicit positive state weights, a target in the
same physical or unitless state-weight unit, occupied state weight, and a
unit-aware residual.

## Rejected and corrected behavior

The prototype includes both Brillouin-zone endpoints with equal full weight. It
then normalizes all six represented one-spin bands to one abstract state and
labels half of the resulting spin-degenerate total as an electron number. The
reported value is a filling fraction under that normalization, not a physical
one-dimensional electron density.

The prototype also solves the particle-number equation on 700 histogram-bin
centers. In the reviewed 300 K calculation, this produced
`3.425126995805 eV`, while applying the pre-existing finite-level solver to the
same 301-by-6 sampled energy inventory produced `3.427016887412 eV`, a binning
difference of about `1.89 meV`. The maintained full-zone façade therefore solves
on original endpoint-excluded sampled eigenvalues. Gaussian or histogram DOS is
used only for declared thermal visualization and approximate plotted integrals.

No zero-temperature Fermi-energy convention was extracted. An integer-filled
insulator admits any chemical potential inside its gap at zero temperature;
midgap or another convention requires a separate explicit contract.

## Disposition

**Removed after replacement.** The band model, normalized full-zone DOS,
weighted finite-temperature inversion, occupied and unoccupied DOS, thermal
window, thermally active DOS, visualization, and pedagogical explanations all
have named maintained destinations. Keeping the monolithic notebook live would
duplicate those contracts while retaining clean-execution, endpoint,
normalization, and histogram-inversion defects.
