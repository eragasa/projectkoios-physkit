# Notebook changes

This file records chronological notebook extraction, consolidation, relocation,
and retirement decisions. It preserves source paths, content identities,
presentation lineage, and disposition rationale. Current notebook review and
migration directions remain at the beginning of [`README.md`](README.md).

## First extraction: finite square-wave series

The byte-identical exploratory notebooks formerly at
`notebooks/math/fourier_series.ipynb` and
`notebooks/math/fourier/fourier_series.ipynb` were consolidated into
`notebooks/foundations/fourier-and-periodic-methods/square-wave-fourier-series.ipynb`.
The reusable odd-harmonic sum now belongs to
`projectkoios.physkit.math.fourier_series.OddSquareWaveFourierSeriesEvaluator`;
the notebook retains visualization, checks, interpretation, and exercises.

## Second extraction: Maxwell–Boltzmann speed density

The overlapping explorations formerly at
`notebooks/statmech/kinetic/mb_velocity_distribution.ipynb` and
`notebooks/gas/dev_maxwell_speed_pdf.ipynb` were consolidated into
`notebooks/thermal/statmech/kinetic/maxwell-boltzmann-speed-distribution.ipynb`.
Reusable SI state validation, characteristic speeds, and density evaluation now
belong to `projectkoios.physkit.thermal.statmech.kinetic_theory`. The distinct
root-level energy-distribution note is now the maintained
`notebooks/thermal/statmech/distributions/boltzmann-exponential-energy-distribution.ipynb`.
Its exact exponential formula was extracted to
`projectkoios.physkit.thermal.statmech.distributions.boltzmann_energy`, and its
constant-density-of-states assumption is explicit so it is not confused with
the three-dimensional translational kinetic-energy distribution.

## Third extraction: von Mises equivalent stress

The separate two- and three-dimensional exploratory notebooks under
`notebooks/materials-science/mechanical-properties/` were consolidated into
`notebooks/mechanics/continuum-mechanics/von-mises-equivalent-stress.ipynb`.
Reusable plane-stress and three-dimensional evaluation now belongs to
`projectkoios.physkit.mechanics.continuum_stress`, preserving pairwise
off-diagonal averaging and consistent-unit output.

## Fourth extraction: Townsend discharge current

The Townsend notebook formerly stored inside the Python package moved to
`notebooks/plasma-physics/gas-discharge/townsend-current.ipynb`. Reusable state,
result, closed-form, generational, and fixed-point behavior now belongs to
`projectkoios.physkit.plasmas.gas_discharge.townsend`. Empty package-local
`breakdown.ipynb` and `paschen.ipynb` placeholders were removed.

## Fifth consolidation: one-dimensional particle in a box

The byte-identical, incomplete notebooks formerly at
`notebooks/qm/qm.piab1d.comp.ipynb` and
`notebooks/qm/qm.piab1d.physkitlib.ipynb` were replaced by
`notebooks/qm/piab1d/piab1d.ipynb`.
The maintained laboratory uses the existing analytical, finite-difference, and
spectral wavefunction APIs under `projectkoios.physkit.qm.piab1d`; no notebook
implementation was promoted into a second library route. Both removed notebooks
had SHA-256
`80c05dda12c7fc2631a1bfec059a82db6de57e355b718687f7f814a8c232ddbf`.

## Sixth consolidation: one-dimensional tight binding

The byte-identical exploratory notebooks formerly at
`notebooks/solidstate/dev_tightbinding.ipynb` and
`notebooks/solidstate/solst01_tightbinding.ipynb` were replaced by
`notebooks/solidstate/tight-binding/one-dimensional-tight-binding.ipynb`.
The maintained laboratory uses the existing finite periodic domain, scalar
hopping model, boundary twist, and represented-operator constructor APIs. It
narrows the lesson to a fully executable reduced nearest-neighbor lattice model;
the duplicated continuum-well and LCAO experiments remain available through Git
history. Both removed notebooks had SHA-256
`e0b3f94ac89e4e1eca4a4318c6dd65cb953afc17241c5038526dcf3f22af597c`.

## Seventh cleanup: invalid notebook placeholders

Three files with no recoverable notebook payload were removed:

- `notebooks/scratch/qm_planck_blackbody.ipynb`, an empty file with SHA-256
  `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`;
- `notebooks/solidstate/lesson_01_directlattice.ipynb`, the same empty payload;
  and
- `notebooks/scratch/thermo_dulongpetit_data.ipynb`, a one-byte file containing
  `p` with SHA-256
  `148de9c5a7a44d19e56cd9ae1a554bf67847afb0c58f6e12fa29ac7ddfca9940`.

Repository history contains no valid notebook document at these paths. A
repository-layout test now parses every remaining `.ipynb` file as JSON.

## Eighth consolidation: direct and reciprocal lattices

The overlapping notebook-local lattice prototypes formerly at
`notebooks/scratch/ch1_direct_lattice.ipynb` and
`notebooks/solidstate/dev_DirectLattice.ipynb` were replaced by
`notebooks/solidstate/lattices/direct-and-reciprocal-lattices.ipynb`.
The maintained laboratory uses the existing physical `DirectLattice*` and
`ReciprocalLattice*` APIs and explicitly distinguishes them from integer
`FinitePeriodicDomain`. The removed notebooks had SHA-256 values
`81112bc7c3a5947c9601e8815a0be271494c9d1ffb87d9d4ba2f20719660bbc1`
and `dff7da097761fae7332dcb6f69f326e7d46ee53b4bf2571567ac33e7beb341b0`,
respectively.

## Ninth extraction: Antoine vapor pressure

The incomplete `projectkoios.physkit.thermal.thermo.vaporpressure` package was replaced
with explicit unit-aware `AntoineVaporPressureModel`, `AntoineTemperatureRange`,
typed evaluation request, internal ActionObject, and
`AntoineVaporPressureEvaluation` owners. The maintained synthetic-data
laboratory retains its semilogarithmic pressure plot at
`notebooks/thermal/thermo/vapor-pressure/antoine-vapor-pressure.ipynb`.
`notebooks/thinfilm/vapor_pressure.py` remains a migration adapter until the
legacy aluminum and aluminum–copper notebooks no longer depend on its former
curve-object interface.

## Tenth consolidation: remaining one-dimensional particle-in-a-box notebooks

The obsolete numerical narrative at `notebooks/qm/qm.piab1d.ipynb` and the
notebook-local symbolic model at `notebooks/qm/qm.piab1d.sym.ipynb` were folded
into the maintained particle-in-a-box lecture note and laboratory. The laboratory
now verifies the continuum eigenvalue equation and normalization directly with
bounded symbolic expressions while retaining the package-owned analytical and
finite-difference implementations. The removed notebooks had SHA-256 values
`e47775cda9a16a29a02b7b7c2b5d922497c106faa6f7b950779b0db967219dac`
and `4d9dacb1a62b2d731b2704be8064e0a9f52310ec73d0221861cfafaa57f71cd4`,
respectively.

## Eleventh extraction: two-dimensional particle in a box

The mixed and error-bearing `notebooks/qm/piab2d.ipynb` was replaced by
`notebooks/qm/piab2d/tise/analytical/piab2d.ipynb`.
The maintained laboratory uses the new `Piab2D`,
`Piab2DAnalyticalEvaluator`, and `Piab2DAnalyticalSolution` owners for the
rectangular homogeneous-Dirichlet continuum spectrum. Periodic reciprocal-space
experiments from the predecessor were not treated as particle-in-a-box behavior
and remain available only through Git history. The removed notebook had SHA-256
`90b4b7ffc0995119f4568bf88afca9c62a68cd3b0b8e1f7150fef41c1538c1dd`.

## Twelfth laboratory: one-dimensional spectral time evolution

`notebooks/qm/piab1d/piab1d__time_evolution.ipynb`
uses the existing analytical TISE and TDSE package APIs to propagate a normalized
two-mode superposition. It checks conserved spectral and coordinate-space norms,
conserved energy expectation, and probability-density recurrence after one beat
period. No notebook-local solver or time integrator was introduced.

## Thirteenth laboratory: finite-level fermion occupations

`notebooks/qm/piab1d/piab1d__fermions.ipynb`
extracts the Fermi–Dirac chemical-potential, discrete-level occupation, and
occupation-weighted spatial-density analysis behind the historical PIAB1D
visualizations. The maintained implementation lives in
`projectkoios.physkit.thermal.statmech.qm`; stationary
orbital evaluation remains intrinsic to `Piab1DAnalyticalResults`. The
laboratory declares a spinless independent-particle model and checks both the
occupation sum and spatial-density integral.

The historical source notebooks remain present during review so their additional
visualization content can be compared before final disposition.

## Fourteenth laboratory: sampled Gaussian basis projection

`notebooks/qm/piab1d/piab1d__basis_projection.ipynb`
extracts the Gaussian-state construction, explicit weighted normalization, and
energy-eigenbasis projection behind the historical PIAB1D expansion-probability
visualization. Reusable behavior lives in
`projectkoios.physkit.qm.sampled_states`; the laboratory owns its grid,
finite basis, plots, interpretation, and exercises.

The maintained PIAB laboratories now reside under the dimension-specific
hierarchy `qm/piab1d/`, `qm/piab2d/`, and `qm/piab3d/`.

## Fifteenth laboratory: analytical PIAB3D spectrum

`notebooks/qm/piab3d/piab3d.ipynb` evaluates finite positive quantum-number
triples for a homogeneous-Dirichlet rectangular cuboid. The laboratory checks
the analytical spectrum, deterministic ordering, cubic permutation degeneracy,
orthorhombic splitting, and the distinction between quantum-number inventory
and physical coordinate space. Reusable behavior lives under
`projectkoios.physkit.qm.piab3d`.

The corresponding lecture notes follow the same dimension-specific naming
standard. PIAB2D further groups material beneath `tise/analytical/` and reserves
`tdse/numerical/` for a separately specified future solver.

## Sixteenth laboratories: analytical multidimensional eigenfunctions

`notebooks/qm/piab2d/tise/analytical/piab2d__eigenfunctions.ipynb` evaluates
normalized product eigenfunctions on rectangular coordinate grids and visualizes
four probability-density heatmaps. `notebooks/qm/piab3d/piab3d__eigenfunctions.ipynb`
uses the general Cartesian plane-slice evaluator for `X`-, `Y`-, and
`Z`-normal slices while retaining the three-dimensional density interpretation. Both laboratories use tested package-owned evaluators and check
boundary values and represented normalization or slice weight.

## Seventeenth laboratory: thermal de Broglie wavelength

The former root-level `notebooks/thermal_debroglie_wavelength.ipynb` is now
`notebooks/thermal/statmech/quantum/thermal-de-broglie-wavelength.ipynb`.
Reusable unit-aware wavelength and degeneracy-parameter evaluation belongs to
`projectkoios.physkit.thermal.statmech.thermal_wavelength`. The laboratory
combines the maintained analytical derivation with an independent numerical
momentum-space Gaussian integral, checks $T^{-1/2}$ scaling, and retains a
logarithmic wavelength plot.

## Eighteenth consolidation: ideal-gas equation of state

Two obsolete root-level ideal-gas notebooks were removed after bounded source
inspection:

- `notebooks/ideal_gas_law.ipynb`, SHA-256
  `cfd858306e6028980b562f07fe232f9c62d2c3124448d143e62f0ed5d2822fa3`,
  contained only empty `EquationOfState` and `IdealGasEos` class placeholders;
- `notebooks/thermo/ideal_gas.ipynb`, SHA-256
  `3a88ab155a3cae0e443d88fbde5ef8cb8d04b4b2ddca82d1597f5aa50d729348`,
  contained a dimensionless local pressure calculation and a saved plot whose
  pressure and volume axis labels were reversed.

Their maintained successor is
`notebooks/thermal/thermo/eos/ideal_gas.ipynb`. Unit-aware reusable behavior,
including ordered isotherm construction, belongs to
`projectkoios.physkit.thermal.thermo.eos.ideal_gas`. Git history retains the two
removed exploratory payloads.

## Nineteenth laboratory: blackbody wavelength spectrum

The error-bearing root notebook `notebooks/wpd_blackbodyradiation.ipynb`, with
SHA-256
`b7104613efbf24989bdc2788de530a3e4411c2658d490262ad759e311db92a88`,
is now the maintained
`notebooks/thermal/statmech/radiation/blackbody-wavelength-spectrum.ipynb`.
The original mixed undefined names, duplicated constants, misspelled the Wien
approximation, and combined wavelength formulas with a frequency-domain local
class. Reusable unit-aware Planck, Wien, and Rayleigh--Jeans evaluation now
belongs to `projectkoios.physkit.thermal.statmech.radiation.blackbody`. The
laboratory checks each approximation in its dimensionless regime and retains
two plot panels.

## Twentieth laboratory: symmetric regular solution

The root notebook `notebooks/thermo/mixing_enthalpy.ipynb`, with SHA-256
`7c9be4451bd98032adf6fd203511c2fde1a7825977fcd0bc7637a3d7d38738a2`,
is now
`notebooks/thermal/thermo/mixtures/symmetric-regular-solution.ipynb`.
The original used dimensionless local parameters and labeled curves as phase
separation or complete miscibility without a coexistence analysis. Reusable
unit-aware mixing enthalpy, configurational entropy, and Gibbs free-energy
evaluation now belongs to
`projectkoios.physkit.thermal.thermo.mixtures.regular_solution`. The maintained
laboratory limits its interpretation to symmetry, endpoint values, and local
curvature; it retains the mixing-quantity plot.

## Twenty-first consolidation: Knudsen number

Three exploratory Knudsen notebooks were consolidated:

- `notebooks/gas/knudsen_number.ipynb`, SHA-256
  `ae4a9a8eaf451289d3d3b5cac8740482aecbc81cda0d86931aa9d7656e79540f`,
  contained one empty cell;
- `notebooks/thinfilm/knudsen_number.ipynb`, SHA-256
  `0337649447d51f7fe84a0cfd31ccbc1445e1be9074c1cbbd1eea20829a033f9a`,
  also contained one empty cell;
- `notebooks/vacuum/dev_knudsen.ipynb`, SHA-256
  `2ebca332d91cc65b3fce400c3180d5ccbf6e5d0b0625be1ca500e3f1a9d4b1de`,
  contained the ratio $\lambda/L$ and an inconsistent regime classifier that
  labeled both its smallest and largest represented ranges as molecular flow;
  and
- `notebooks/vacuum/knudsen_number.ipynb`, SHA-256
  `d6d35507d102a093e2cceaa96f1b6a11b7a15eba9e16536b16790be2b59ed6ff`,
  contained only an untyped one-line $\lambda/L$ function.

The maintained successor is
`notebooks/thermal/statmech/transport/knudsen-number.ipynb`. Reusable unit-aware
ratio evaluation belongs to
`projectkoios.physkit.thermal.statmech.transport.knudsen`. No regime thresholds
were inferred; the laboratory checks inverse-pressure scaling and retains a
logarithmic plot.

## Twenty-second consolidation: ideal occupation statistics

The oversized interactive workbook
`notebooks/thermo/statisticaldistributions/plot_fermi-dirac.ipynb`, SHA-256
`945385324fe8490a224331103d3cec234aac7eb29aab3e2d2cd09df590ca8253`,
and the malformed scratch file
`notebooks/scratch/fermi_distribution_function.ipynb`, SHA-256
`be46a93d7eb9617c4a03a6d44b34b031457293fc4d1dd51a37bbb3df6910728c`,
were consolidated into
`notebooks/thermal/statmech/distributions/occupations/ideal-occupation-statistics.ipynb`.
The maintained implementation separates concrete Fermi--Dirac, Bose--Einstein,
and Maxwell--Boltzmann ownership beneath
`projectkoios.physkit.thermal.statmech.distributions.occupations`. Public calls
accept explicit energy, chemical-potential, and temperature quantities; only
the resulting occupation factors are unitless. The laboratory checks the
dilute-limit convergence and retains a three-curve semilogarithmic plot.

## Twenty-third laboratory: planar three-dimensional Hertz–Knudsen flux

The exploratory `notebooks/thinfilm/hertz_knudsen.ipynb`, SHA-256
`79a33652d296a5a1fc7fc3079381c00c614e782da7b1c4fb7169610bed532f6e`,
is now
`notebooks/thermal/transport/phase-change/hertz-knudsen-net-flux.ipynb`.
Reusable signed interfacial flux evaluation belongs to
`projectkoios.physkit.thermal.transport.phase_change.hertz_knudsen3d`. The model
is explicitly a three-dimensional gas incident on a two-dimensional planar
interface, retains separate evaporation and condensation coefficients, and does
not clamp negative net flux. The laboratory numerically integrates the positive
normal-velocity Maxwell moment, composes a synthetic Antoine pressure, and
retains pressure and flux plots.

## Twenty-fourth consolidation: thin-film vapor-pressure adapters

`notebooks/thinfilm/Al_evaporation.ipynb`, SHA-256
`3f386e0eeade1cc815b48fb1ba50a57fd1f0eb72fc346a0bf072a9f4910d137c`,
is now the maintained
`notebooks/thermal/transport/phase-change/aluminum-evaporation-legacy-coefficients.ipynb`.
It uses the package Antoine and planar three-dimensional Hertz--Knudsen APIs and
labels the retained coefficient attribution as unverified legacy input.

The following incomplete or unsupported exploratory material was removed after
bounded source inspection and remains available through Git history:

- `notebooks/thinfilm/AlCu_evaporation.ipynb`, SHA-256
  `60b69461e8a86dba790786434e64469480550937ce970a19526ff6e169b0f93e`,
  combined incomplete activity calculations with an unverified secondary link
  for copper coefficients;
- `notebooks/thinfilm/evaporationrate_Al.ipynb`, SHA-256
  `155d5d38b6ea38be4b202cd3c25124d636fd7d19b462b9662e9296c1dddb4f33`,
  contained syntactically incomplete Antoine, Shomate, and data-holder classes;
- `notebooks/thinfilm/vapor_pressure.py`, SHA-256
  `d5f88e2df90e19e14ed35b6291907983fbffef09860fc0c4d23eb7dcf9ba2d55`,
  was superseded by the maintained unit-aware Antoine implementation;
- `notebooks/thinfilm/hertz_knudsen.py`, SHA-256
  `718de3391c317173aad4c6deb4f0eb1ed265f6fca698ef904381426b79963411`,
  clamped condensation-dominated net flux and used one coefficient for both
  directions;
- `notebooks/thinfilm/mixtures.py`, SHA-256
  `c0ee151f51719c777acf8f79a51d635586b038b5f14d77e5a26e7f8781f4a1ad`,
  depended on the removed vapor-pressure adapter and untyped dictionaries; and
- `notebooks/thinfilm/thinfilm.py`, SHA-256
  `a6f250a710e3071df627ba9c0f3b6a0a6ab28fe039a701415b12c7f0280d0e94`,
  contained an incomplete species record with undefined attributes.

## Twenty-fifth consolidation: planar deposition geometry

`notebooks/thinfilm/AlEvaporationLab.ipynb`, migrated working-copy SHA-256
`ccde18d413807097002950b24ac632585d188f2a6c926c68e7f5ce60e821d515`,
is now
`notebooks/deposition/geometry/planar-source-thickness-profiles.ipynb`.
The maintained laboratory retains the distinct point-source and finite circular
source geometry while removing notebook-local thermodynamics, phase-change,
geometry, and interactive image-click implementations. Reusable unit-aware
normalized-profile behavior belongs to
`projectkoios.physkit.deposition.geometry.plane_point` and
`projectkoios.physkit.deposition.surface_source.disk`. The laboratory explicitly
distinguishes the shrinking isotropic-disk limit from the Lambertian angular
law and adds quadrature-refinement checks.

The unverified aluminum coefficient example was already migrated separately in
the twenty-fourth consolidation. Placeholder source area, film density, and
apparatus dimensions were not promoted to maintained physical predictions.

## Twenty-sixth consolidation: molar and mass composition

Two exploratory mixture notebooks were consolidated:

- `notebooks/thinfilm/dev_molar_mixture.ipynb`, SHA-256
  `c738bd22c0b1f1fc19086e922401f089ee8717fb586dada19c70191e3cac1d67`,
  held mole- and mass-based dictionary prototypes; and
- `notebooks/mixtures/mixture_molefraction.ipynb`, migrated working-copy
  SHA-256
  `cd58c5952524b4dfd83cc5604a34dd1cfdea08f0a473e3ffa86ac918b30b5ff3`,
  represented mole fractions and an optional total with a duplicate local
  class.

Their maintained successor is
`notebooks/materials/mixtures/molar-and-mass-composition.ipynb`. Reusable
composition records and bidirectional conversion belong to
`projectkoios.physkit.materials.mixtures.molar`. The maintained API replaces
mutable dictionaries with deterministic species tuples and immutable
unit-bearing vectors. The laboratory uses synthetic component data, checks a
molar-to-mass-to-molar round trip, and compares mole and mass fractions in a
retained plot.

## Twenty-seventh disposition: unsupported thin-film explorations

`notebooks/thinfilm/janaf_tables.ipynb`, SHA-256
`a30ef0c71eb9b192ca285394913552e2753b08d2aedb365971a93c1ca2060bf1`,
was retired rather than promoted. Its notebook-local parser used mutable pandas
state, did not define a stable missing-value or phase-transition-row contract,
and loaded an empty liquid-state file as though it were a table. The five input
files were moved byte-for-byte to
`notebooks/materials/thermochemistry/resources/legacy-janaf-unverified/`, where
their individual hashes and missing provenance are recorded. They are not
maintained literature evidence.

`notebooks/thinfilm/ase_example.ipynb`, SHA-256
`9526f62de1ea470960a57587e1c0b0b6049c8f3ae2edf5a02aa90e80b90b2f2e`,
was retired. It attempted an in-notebook `ase` installation, contained
incomplete lattice prototypes with invalid NumPy construction and coordinate
checks, and produced no maintained calculation. ASE was not added as a project
dependency.

The unverified image `notebooks/thinfilm/img/ohring2ed_fig3.1.png`, SHA-256
`00a555b85b47f00c2c8868c4da3cc2b3cdb5c3593466182943ece2ddc16d74e5`,
was preserved byte-for-byte under
`notebooks/deposition/geometry/resources/legacy/`. It is not displayed because
its exact source and reproduction permission remain unverified.

With these dispositions, the exploratory `notebooks/thinfilm/` directory has
no remaining files. Its maintained successors are organized under `thermal`,
`deposition`, and `materials` according to model ownership.

## Twenty-eighth consolidation: sublattice site fractions

`notebooks/mixtures/sublattice.ipynb`, SHA-256
`a2ee5d27baf2811035af801d2b228a93180ea41403f2ed1e29730e5f6933e890`,
is now `notebooks/materials/mixtures/sublattice-site-fractions.ipynb`.
Notebook-local mutable dictionaries and inconsistent prototype type names were
replaced by immutable ordered records under
`projectkoios.physkit.materials.mixtures.sublattice`. The maintained model
normalizes each sublattice independently, applies explicit formula-unit
multiplicities, preserves deterministic species order, and records any vacancy
or pseudocomponent exclusion in the evaluation request. The laboratory uses
synthetic site counts and makes no equilibrium or real-material claim.

## Twenty-ninth consolidation: radial pair potentials

Two exploratory notebooks were consolidated:

- `notebooks/scratch/notebook-lennardjones.ipynb`, SHA-256
  `ad5327ee8a045bf38bc566724afa5e5c6eeadfc6dfc226f9cdf9ac73d244d1dc`,
  contained notebook-local Lennard--Jones energy and force functions; and
- `notebooks/solidstate/dev_InteratomicPotential.ipynb`, SHA-256
  `b83c56d675b8e5b563c713ee866cef06a02bde2ee669194db0c51ced77c7f1a1`,
  contained duplicate Morse prototypes and labeled $+dV/dr$ as the force
  $-dV/dr$.

Their maintained successor is
`notebooks/solidstate/interatomic-potentials/radial-pair-potentials.ipynb`.
Unit-aware shared contracts belong to
`projectkoios.physkit.solidstate.interatomic_potentials.pair.base`, with
concrete ownership under `pair.lennard_jones_12_6` and `pair.morse`. The
deterministic Morse sign correction restores positive repulsion, negative
attraction, and
$F_r=-dV/dr$. The laboratory matches illustrative well depth, equilibrium
distance, and local curvature, checks both forces by independent numerical
differentiation, and makes no material-fit claim.

## Thirtieth consolidation: dimensionless Cahn--Hilliard evolution

`notebooks/scratch/notebook_cahnhilliard_1d.ipynb`, SHA-256
`61eb02cd34db2ca27d577881fc7473e683cd4c7ea2d3284f19771a86f08fbe2e`,
is now
`notebooks/solidstate/phase-field/dimensionless-periodic-cahn-hilliard.ipynb`.
The original notebook contained both one- and two-dimensional periodic FFT
experiments despite its 1D filename. Both ideas are retained under
`projectkoios.physkit.solidstate.phase_field.cahn_hilliard` with shared reduced
parameters and dimension-specific immutable states, solve records, and spectral
solvers. The maintained laboratory uses deterministic synthetic fields and
checks zero-mode conservation. It does not claim dimensional parameter mapping,
convergence, energy stability, coarsening scaling, or material validation.

An unused exploratory chemical-potential helper had the wrong sign on its
spectral gradient contribution. The maintained solver follows the documented
semi-implicit update and omits that inconsistent unused helper.

## Thirty-first consolidation: generic complex plane waves

Three exploratory notebooks are consolidated into
`notebooks/waves/plane/complex-plane-waves.ipynb`:

- `notebooks/scratch/qm_planewave.ipynb`, SHA-256
  `f0f2451b501fb1e0995bea1de248e9c601cfa727192bd45576231eef8a5ea43d`;
- `notebooks/scratch/ch2_planewave_1d.ipynb`, SHA-256
  `19b6a73acd6affb2351219e5429f8e9a8279a3709914c044e24d6e78bf3349f2`;
  and
- `notebooks/scratch/ch2_planewaves.ipynb`, SHA-256
  `b4db965f62e4e59e61285f0b497b9dfdfa185f9ccd82bf5c7975167d622be5f0`.

Their useful ideas—angular frequency, wave number, complex one-dimensional
fields, two-dimensional wave vectors, and an FFT-space peak—are retained in
unit-aware one- and two-dimensional models under
`projectkoios.physkit.waves.plane`. The maintained laboratory treats the plane
wave as a generic unitless complex field. It does not silently promote it to a
quantum wavefunction or impose a dispersion relation.

## Thirty-second consolidation: remaining PIAB1D prototypes

Two remaining exploratory PIAB1D notebooks are retired as superseded:

- `notebooks/solidstate/dev_L3_piab1d.ipynb`, SHA-256
  `7d0e67c4f3fde727cb118d4ad07e94b8cc7de1de6dbb520379d996236e6de37a`;
  and
- `notebooks/scratch/dev_analytical_PIAB1d.ipynb`, SHA-256
  `7fd2266a3d474cb273f81ca3c9924c98f4c5f457b27d11d3268fad7aacaff01f`.

Their maintained successors separate the ideas by responsibility:

- `notebooks/qm/piab1d/piab1d.ipynb` compares continuum and represented
  finite-difference spectra and explicitly identifies discretization error;
- `notebooks/qm/piab1d/piab1d__basis_projection.ipynb` owns Gaussian-state
  normalization and finite-basis projection;
- `notebooks/qm/piab1d/piab1d__canonical_probabilities.ipynb` owns normalized
  single-particle canonical probabilities on a finite PIAB1D level inventory;
- `notebooks/qm/piab1d/piab1d__fermions.ipynb` owns finite-level Fermi--Dirac
  occupations and occupied spatial density; and
- `notebooks/qm/piab1d/piab1d__time_evolution.ipynb` owns spectral state
  evolution.

The package APIs used by these laboratories retain the analytical and
finite-difference PIAB solutions, continuously normalized sampled states,
projection coefficients, occupation solving, and immutable evolution results.
The graphical and instructional sequence is also retained deliberately:

- the full/interior/boundary grid presentation, absolute and scaled spectra,
  relative-error graph, ground-state wavefunction, probability density, and a
  visible Hamiltonian matrix block are in `piab1d.ipynb`;
- the expansion-probability bars and reconstructed Gaussian density are in
  `piab1d__basis_projection.ipynb`;
- the single-particle canonical-probability bars and their finite energy
  inventory are in `piab1d__canonical_probabilities.ipynb`;
- the continuous Fermi--Dirac curve, discrete-level samples, and occupied
  density are in `piab1d__fermions.ipynb`; and
- the time-dependent probability-density presentation is in
  `piab1d__time_evolution.ipynb`.

The original course, term, lecture purpose, author, and operator-to-matrix
presentation lineage are stated in the maintained primary laboratory. The
retired lecture prototype contained broken names and unreachable cells; those
failures are not retained as implementations.

## Thirty-third consolidation: direct-lattice vector prototypes

Two development notebooks are retired as superseded by the maintained direct-
and reciprocal-lattice laboratory:

- `notebooks/solidstate/dev_DirectLatttice2D.ipynb`, SHA-256
  `19976540132682ad2d6a9092db07094d20de4e16a914b6eb6998cdaf781ab701`;
  and
- `notebooks/solidstate/dev_unitcell.ipynb`, SHA-256
  `5ce25db16212a12c1b694400fed45b2e2866949be67246d6706dcf7d1ce79cd1`.

The first contained only an obsolete import. The second explored row and column
array shapes and an incomplete `DirectLattice` prototype; despite its filename,
it did not construct a basis-bearing unit-cell object. The maintained successor
is
`notebooks/solidstate/lattices/direct-and-reciprocal-lattices.ipynb`,
which imports the package-owned `DirectLattice1D`, `DirectLattice2D`, and
`DirectLattice3D` representations, checks direct--reciprocal duality, and plots
both direct and reciprocal points for an oblique two-dimensional lattice. The
original column-vector and basis-matrix presentation is retained in the linked
lecture note through $A=[\mathbf a_1\;\cdots\;\mathbf a_d]$ and
$\mathbf R_{\mathbf n}=A\mathbf n$. Neither retired prototype contained a
figure to preserve; the successor adds the maintained graphical presentation.
This consolidation does not merge physical direct-lattice ownership with
integer finite-periodic domains or introduce a new unit-cell dependency.

### Standalone-retention decision

| Retired notebook | Standalone value | Disposition |
|---|---|---|
| `dev_L3_piab1d.ipynb` | No. It mixed grid construction, finite differences, spectra, basis projection, canonical statistics, and broken unreachable cells. | Split its narrative and graphs among the primary, basis-projection, canonical-probability, and fermion PIAB1D laboratories. |
| `dev_analytical_PIAB1d.ipynb` | No. It was an incomplete duplicate with inconsistent names, but its continuous Fermi--Dirac overlay was useful. | Preserve that overlay in `piab1d__fermions.ipynb`. |
| `dev_DirectLatttice2D.ipynb` | No. It contained only one obsolete import and no narrative, calculation, or figure. | Retire without a standalone successor; use the maintained direct/reciprocal lattice laboratory. |
| `dev_unitcell.ipynb` | No. It explored array orientation and an incomplete direct-lattice prototype rather than a unit-cell model. | Preserve the column-basis convention in the direct/reciprocal lecture note and its geometry in the lattice laboratory. |

Thus the information and meaningful presentations are preserved, but none of
the four source notebooks merits continued standalone maintenance. This is an
artifact-ownership decision, not a claim of pedagogical validation.

## Thirty-fourth consolidation: finite-difference QHO1D

`notebooks/solidstate/dev_L3_qho1d.ipynb`, SHA-256
`08caf30410de1cb55347c0d6b757f89fbac4d36e92e72d54681f9710ebc47eaa`,
is now `notebooks/qm/qho1d/qho1d__finite_difference.ipynb`.

The original notebook had one useful and coherent numerical experiment: a
bounded second-order finite-difference representation of the dimensionless
one-dimensional quantum harmonic oscillator, with low eigenfunctions shifted by
their energies over the harmonic potential. That graph and presentation deserve
preservation, but the notebook-local Hamiltonian construction and normalization
do not. They are now owned by the tested
`projectkoios.physkit.qm.qho1d.tise.finite_difference` package.

The maintained laboratory adds explicit analytical-energy comparison,
continuous-normalization checks, a visible relative-error graph, prerequisites,
interpretation, and exercises. The finite interval is identified as numerical
truncation rather than a physical hard wall.

The repository-root `notebooks/qho-hermite-basis.ipynb` is not treated as a
duplicate. It uses a truncated Hermite function basis, so its basis-truncation
and reference-frequency questions remain distinct from coordinate-grid spacing
and finite-domain truncation. It remains a separate migration candidate.

### Standalone-retention decision

The exploratory path does not merit standalone maintenance because reusable
behavior was embedded locally and the presentation lacked assumptions and
error boundaries. Its finite-difference experiment does merit a dedicated
maintained successor rather than consolidation into the Hermite-basis notebook,
because the two representations have different state spaces and numerical
error controls.

## Thirty-fifth consolidation: incomplete potential object sketch

`notebooks/solidstate/dev_Potential.ipynb`, SHA-256
`77d7fb60602c8c9c309559dac9e2981e7d28d8f3f89d96c7c32c90c3aa393ae2`,
is retired without a standalone notebook successor.

The one-cell sketch contained five ideas with different dispositions:

| Prototype idea | Disposition |
|---|---|
| Generic interatomic-potential owner | Consolidated into `projectkoios.physkit.solidstate.interatomic_potentials`. |
| Radial pair-potential boundary | Consolidated into `pair/base.py` as `RadialPairPotential` and typed evaluation records. |
| Lennard--Jones parameter and potential owners | Consolidated into the maintained `lennard_jones_12_6` parameter, potential, and evaluator modules. |
| Chemical-species pair attached to parameters | Preserved as an explicit prospective force-field mapping concern, not added to the species-agnostic radial formula. |
| Mutable `Atom` and dictionary-backed `AtomList` | Preserved only as an atomistic-structure concept. The incomplete mutable implementation is not promoted, and this migration does not couple pair potentials to `periodic.unit_cell`. |

The notebook contained no calculation, output, narrative, or graph to preserve.
Its abstract classes returned the `NotImplemented` singleton instead of defining
an enforceable abstract contract; its constructors used erased dictionaries;
and the `super` calls were invalid. None of those implementations merits
maintenance. The useful boundary distinctions are recorded in the pair-potential
architecture and lecture note, while any future species-to-parameter registry or
atomistic structure requires its own explicit contract and task.

## Thirty-sixth retirement: empty intrinsic-carrier notebook

`notebooks/semphy/intrinsic_carrier_electron.ipynb`, SHA-256
`05e1090eac867e7fff7b38d58b203d90504e1ff80c574e2dc7f64056039ade0e`,
is retired. Both cells were empty: there was no equation, assumption,
implementation, output, graph, citation, or presentation to consolidate.

The filename alone is not treated as a scientific model or durable idea.
Intrinsic-carrier concentration remains future semiconductor-statistics work
that would require explicit band-edge, density-of-states, degeneracy,
temperature, and unit conventions. This empty artifact provides no authority
for choosing those contracts.

## Thirty-seventh disposition: incomplete effective-density-of-states sketch

`notebooks/semphy/effective_DOS_electron.ipynb`, SHA-256
`70e24ea3fd09ab25117fd34838cf81b6c6f0998f9c677fd072c0c94c140bfc3a`,
is retired without promoting its formula or classes.

The conceptual idea—effective conduction- and valence-band densities of states
for parabolic bands—deserves future treatment, but the notebook does not define
a usable scientific contract:

- the filename says electron, the title says holes, and the code mixes
  conduction- and valence-band owners;
- both class definitions are syntactically incomplete;
- the hole result returns an undefined name;
- temperature, mass, Planck constant, output volume, and unit conventions are
  unspecified;
- spin and valley degeneracies are not declared;
- the meaning of effective mass is not identified; and
- under the common three-dimensional parabolic-band convention, the written
  expression omits the reciprocal-space normalization factor associated with
  $2\pi$.

These are scientific-definition questions, not mechanical code corrections.
The one-cell notebook contained no output, graph, citation, or demonstrated
usage to preserve. A future maintained model must explicitly select spatial
dimensionality, band-edge and energy-zero conventions, density-of-states mass,
degeneracies, physical units, and a verified derivation before implementation.
No such choice is inferred from this artifact.

## Thirty-eighth consolidation: open-chain tight-binding level splitting

`notebooks/semphy/tb1d.level_splitting.ipynb`, SHA-256
`0a0b443fb15da453a20f1e7fed256e5862467bca298affb57819480f67ad47d6`,
is retired as a standalone notebook.

Its markdown-only idea is valuable: nearest-neighbor coupling splits an
$N$-fold isolated onsite level in a finite open chain, while a periodic chain
organizes the corresponding states by wave vector into
$E(k)=\varepsilon_0-2t\cos k$. The open-chain Hamiltonian, its exact finite
levels, and the illustrative separation-dependent hopping ansatz are now
consolidated into the one-dimensional tight-binding lecture note.

The artifact contained an empty code cell and no calculation, graph, output,
parameter provenance, or executable check. It therefore does not merit a
maintained laboratory. The existing maintained tight-binding laboratory remains
periodic and does not claim to implement physical separation-dependent hopping.
A future open-boundary computational laboratory would require package ownership
for that boundary condition and a clearly classified hopping law rather than
copying the unevaluated notebook sketch.

## Thirty-ninth consolidation: Kronig--Penney delta-comb bands

`notebooks/solidstate/kronig_penney.ipynb`, SHA-256
`452de58a6435320f0e86b00acbef04de8db1b9ab53cb4e82baf16b4af18075f9`,
is now
`notebooks/solidstate/electronic-structure/kronig-penney/delta-comb-bands.ipynb`.

The idealized band-structure graph deserves preservation, but the original
calculation does not. It assigned energy units to a delta-potential coefficient
that requires energy-times-length dimensions, used an arbitrary residual
tolerance to select points on a two-dimensional mesh, and supplied no derivation
or verification.

The maintained successor uses the explicit dimensionless convention

$$
\cos(ka)=\cos z+P\frac{\sin z}{z},
\qquad
E/E_a=z^2,
$$

with $P\geq0$. Package-owned typed sampling retains exactly those explicit
samples for which the dispersion right-hand side lies in $[-1,1]$. The
laboratory preserves the symmetric reduced-band presentation and adds an
allowed-value graph, executable residual checks, interpretation, and exercises.
It makes no material band-structure claim.

### Standalone-retention decision

The old path does not merit maintenance because its dimensional contract and
point-selection method were invalid. The corrected dimensionless model does
merit a dedicated maintained laboratory because the band/gap construction is a
coherent reusable solid-state concept not duplicated by the current
tight-binding laboratory.

## Fortieth relocation: PIAB2D periodic-boundary intent

`notebooks/solidstate/piab2d_pbc.ipynb`, SHA-256
`0874fadf3a29997ec9ad059087c8e5e74387b369310a1be2db44466d26ffd7aa`,
is preserved byte-for-byte at
`notebooks/qm/piab2d/tise/numerical/oblique-periodic-boundary-prototype.ipynb`.

The artifact is retained for its complete informational content, not merely its
current code. Its intended combination of periodic boundaries, oblique direct
and reciprocal geometry, metric-aware coordinate transformation, eigenproblem,
and physical-space probability-density visualization is a valuable prospective
PIAB2D model and presentation sequence. The saved notebook and outputs remain
available for joint review with `notebooks/scratch/generalized_piab2d.ipynb`.

The present implementation does not yet establish the intended periodic model:
the active solve path lacks seam couplings, endpoint and spacing conventions are
inconsistent, and the periodic helper is unused. Those findings label the
implementation maturity; they do not retire the concept, narrative intent, or
graphs. Consolidation is deferred until a maintained successor explicitly owns
and reproduces every valuable layer.

## Forty-first relocation: canonical solid-state notebook directory

The maintained notebook hierarchy formerly rooted at
`notebooks/solid-state-physics/` is now rooted at `notebooks/solidstate/`.
Solid-state laboratories, local README files, and their relative internal
structure moved together. Notebook contents and saved outputs were not rewritten
for this administrative relocation.

Live documentation links and the maintained-organization table now use
`notebooks/solidstate/`. Historical donor paths recorded in earlier entries
remain explicit where needed to preserve lineage. The Python package already
uses `projectkoios.physkit.solidstate`, so the notebook directory now matches
that domain spelling.

## Forty-second relocation: generalized oblique PIAB2D prototype

`notebooks/scratch/generalized_piab2d.ipynb`, SHA-256
`22c333ddea3d0eacf944124660f474bcd9e4428391d7d68b9c28661038e58322`,
is preserved byte-for-byte at
`notebooks/qm/piab2d/tise/numerical/generalized-oblique-coordinate-prototype.ipynb`.

The complete artifact remains intact, including 32 cells, 20 saved outputs, 10
saved graphical displays, and one saved error. Its value extends beyond current
implementation correctness: it records metric-tensor derivation, oblique
coordinate mapping, a mixed-derivative finite-difference route, reciprocal
basis construction, plane-wave energies, real-space reconstruction, and several
physical- and reciprocal-space presentations.

It is retained separately from
`oblique-periodic-boundary-prototype.ipynb`. The two artifacts overlap, but they
do not yet represent one identified mathematical problem: the generalized
notebook mixes bounded curvilinear and reciprocal-periodic routes, while the PBC
prototype records a more direct periodic-boundary intent. Consolidation is
deferred until the state spaces, boundary conditions, endpoint conventions,
normalizations, and numerical error sources are stated route by route and every
valuable narrative and graphical layer has a maintained destination.

## Forty-third extraction: periodic primitive-cell 2D methodological review

The periodic reciprocal-space layers identified in the two preserved PIAB2D
prototypes now have a maintained successor at
`notebooks/solidstate/semiconductors/effmass_2d/primitive-cell-methodological-review.ipynb`.
The successor is deliberately named as a periodic free-particle primitive-cell
model rather than a particle in a hard-wall box.

Reusable geometry, Cartesian effective mass, exact Bloch-mode energies, and
cell-normalized plane-wave sampling are owned by
`projectkoios.physkit.solidstate.semiconductors.effmass_2d`. Sibling
`effmass_1d` and `effmass_3d` packages reserve the dimension-specific naming
boundary but do not yet claim implementations. The maintained 2D laboratory
uses `UnitSystem.METAL`: direct lengths in angstrom, particle masses
in dalton, action in electron-volt picosecond, and energies in electron-volt.
The evaluator performs the required compound-unit conversion explicitly.

This 2D model is a methodological review before a separate 3D primitive-cell
contract. It verifies direct/reciprocal duality, reciprocal indexing, anisotropic
Cartesian mass tensors, rotation invariance, normalization, and Bloch
translation phases. Its basis and mass values are synthetic illustrative
inputs, not literature values or a material-specific band structure.

The two original prototypes remain preserved because the maintained successor
does not yet reproduce their bounded curvilinear finite-difference route or all
historical graphs and outputs. No DFT or Wannier calculation is implemented or
claimed.

## Forty-fourth extraction: three-dimensional primitive-cell effective-mass reference

The conventions reviewed by `effmass_2d` now have a maintained
three-dimensional implementation under
`projectkoios.physkit.solidstate.semiconductors.effmass_3d` and a laboratory at
`notebooks/solidstate/semiconductors/effmass_3d/primitive-cell.ipynb`.

The implementation owns a physical three-vector primitive basis, reciprocal
basis, cell volume, symmetric positive-definite Cartesian effective-mass tensor,
exact Bloch-mode kinetic energies, and volume-normalized plane-wave samples. It
retains the metal-unit convention and explicit compound-unit conversion reviewed
in 2D.

The laboratory uses an FCC-like nonorthogonal basis solely as illustrative
synthetic geometry. Its lattice scale and mass tensor are not identified as
silicon values or literature evidence. The calculated outputs verify the stated
idealized model; they do not reproduce DFT, fit a semiconductor, or establish
scientific validation.

`effmass_1d` remains reserved without an implementation. The exact historical
PIAB2D prototypes remain preserved because neither effective-mass successor
owns their bounded curvilinear finite-difference route or complete presentation
lineage.

## Forty-fifth extraction: `pbc1d` and free-electron `bloch1d`

`notebooks/semphy/freeelectron.bloch1d.ipynb`, SHA-256
`1dd413a740d143866141fd94af04b27d5c85dd85fdff423926f01529f7e22ef3`,
is preserved byte-for-byte at
`notebooks/solidstate/bloch1d/freeelectron/exploratory-finite-difference-prototype.ipynb`.
The artifact retains 16 cells, 25 saved outputs, 24 graphical displays, and its
interactive-widget presentation code without saved error output.

The maintained successor is
`notebooks/solidstate/bloch1d/freeelectron/free-electron-bloch-bands.ipynb`.
The physical folded-band model now belongs to
`projectkoios.physkit.solidstate.bloch1d.freeelectron`, not to
`semiconductors.effmass_1d`. It uses the bare electron mass and distinguishes the
free electron on the infinite line from each fixed-$k$ reference-cell Bloch
fiber.

Generic numerical boundary mechanics now belong independently to
`projectkoios.physkit.periodic.pbc1d`: an endpoint-excluded grid and centered
finite-difference Laplacian with conjugate Bloch seam phases. The
free-electron package combines that represented Laplacian with
$-\hbar^2/(2m_e)$ and compares its grid convergence with the exact folded
spectrum.

`semiconductors.effmass_1d` returns to reserved, unimplemented status. The
maintained notebook reproduces folded bands, finite-difference comparisons,
represented states, and Fourier content. The original prototype remains because
its interactive widgets do not yet have a maintained replacement.

## Forty-sixth extraction: cosine-potential `bloch1d`

Two related exploratory artifacts remain byte-identical under
`notebooks/solidstate/bloch1d/cosinepotential/`:

- `exploratory-band-gap-prototype.ipynb`, formerly
  `notebooks/semphy/cosine.bloch1d.ipynb`, SHA-256
  `e10afcc2b0fc2c73d77ba2934e570c4477c548f59ab9f4ecb8bf599679eba8cd`;
- `exploratory-bands-dos-thermal-prototype.ipynb`, formerly
  `notebooks/scratch/solidstate/bloch_1d_cosine_potential_bands_dos.ipynb`,
  SHA-256
  `a029505c63023d0a09477d71fce8f1bb11c8275016f5de57ef9d1d07d0fde674`.

The maintained band and first-gap successor is
`notebooks/solidstate/bloch1d/cosinepotential/cosine-potential-bands.ipynb`.
Reusable ownership is
`projectkoios.physkit.solidstate.bloch1d.cosinepotential`, which consumes the
independently maintained free-electron and `pbc1d` representations.

The maintained model verifies the represented free-electron limit,
Brillouin-zone endpoint agreement, weak-potential first-gap behavior, amplitude-
sign invariance, and grid refinement. Its potential amplitude is illustrative,
not a fitted semiconductor parameter. Reusable one-dimensional potential
rendering now belongs to `projectkoios.physkit.potentials.periodic1d`, while the
maintained band figure uses the dimension-independent scalar-path visualizer in
`projectkoios.physkit.solidstate.visualization` rather than notebook-local band
plotting behavior. The extended prototype remains preserved because its
finite-temperature chemical-potential, occupied/unoccupied DOS, thermal-window,
and interactive layers still require separate maintained contracts. Its
ordinary full-zone DOS is now superseded by the maintained dimensional DOS
contract described below.

## Forty-seventh relocation: canonical solid-state lecture-note root

The complete solid-state lecture-note hierarchy formerly rooted at
`docs/lecture-notes/solid-state-physics/` is now rooted at
`docs/lecturenotes/solidstate/`. The move preserves subject-relative structure
while aligning the canonical compact path requested for solid-state material.
Live documentation and notebook prerequisite links now use the new root.

## Forty-eighth extraction: nonorthogonal lattice metric geometry

The reviewed exploratory notebooks formerly staged at
`notebooks/solid-state-physics/crystal-structure/metric-tensor-2d.ipynb` and
`metric-tensor-3d.ipynb` were consolidated into
`notebooks/solidstate/lattices/nonorthogonal-metric-geometry.ipynb` and its
prerequisite lecture note.

Reusable direct-basis metric construction, fractional-to-Cartesian mapping, and
exact two- and three-dimensional nearest-image resolution now belong to
`projectkoios.physkit.periodic.lattice.metric`. The maintained laboratory uses
those package APIs to check direct--reciprocal duality, metric identities, a
skew-cell nearest-image counterexample, the two-dimensional mixed-derivative
coefficient, and equivalent Cartesian and metric Bloch factors.

The former bounded neighbor-shell search was not promoted as a general
algorithm. The maintained resolver derives a complete finite integer search
bound from the smallest singular value of the direct basis. Exact source hashes
and semantic adaptations are recorded in
`docs/provenance/nonorthogonal-lattice-geometry-notebook-migration.md`.

## Forty-ninth extraction: dimensional full-zone density of states

The exploratory one-dimensional histogram and maintained two-dimensional
checkerboard Gaussian DOS were rewritten behind shared normalization contracts
in `projectkoios.physkit.dos`. Concrete reciprocal-length,
reciprocal-area, and reciprocal-volume validation belongs to `dos1d`, `dos2d`,
and `dos3d` respectively. A dimension-independent DOS visualizer owns reusable
horizontal and vertical rendering.

The maintained 1D cosine-potential, 2D checkerboard, and 3D effective-mass
laboratories now use endpoint-excluded full-zone meshes and the dimensional DOS
models. Every represented complete band has an explicit state-count integral;
state degeneracy and Gaussian broadening remain declared inputs. The former 1D
prototype's duplicated endpoint and arbitrary-unit label were not retained.

Exact source identities and adaptations are recorded in
`docs/provenance/full-zone-density-of-states-library.md`.

## Fiftieth disposition: exploratory cosine band-gap prototype

`notebooks/solidstate/bloch1d/cosinepotential/exploratory-band-gap-prototype.ipynb`
was removed after individual review. All 12 nonempty code cells executed without
error, and its saved figures were inspected. Its cosine potential, Bloch bands,
first-gap checks, full-zone DOS, visualizers, and explanatory material now have
named maintained package, notebook, and lecture-note destinations.

The prototype's duplicated Brillouin-zone endpoint weighting and arbitrary-unit
DOS label were not retained. Exact source identity, evidence, destination
mapping, rejected behavior, and recovery information are recorded in
`docs/provenance/exploratory-band-gap-prototype-disposition.md`. The separate
extended thermal prototype remains present and requires its own review.

## Fifty-first extraction: weighted chemical potential and thermal DOS

The extended cosine-potential prototype was reviewed individually. Its
self-contained final calculation executes, but the notebook does not execute
linearly from a clean namespace. More importantly, it duplicates a periodic
Brillouin-zone endpoint, normalizes six one-spin bands to one abstract state,
and solves its particle-number equation on histogram-bin centers.

The existing `FermiDiracChemicalPotentialSolver` now supports explicit positive
state weights and a target in the same physical or unitless state-weight unit.
The full-zone DOS models map every original sampled eigenvalue to its exact
state-density weight before solving; neither a histogram nor Gaussian display
grid participates in the inversion.

`projectkoios.physkit.dos.thermal` owns occupied and unoccupied DOS, the
Fermi--Dirac thermal window, thermally active DOS, and numerical rendered-DOS
integrals. The maintained 1D cosine-potential laboratory demonstrates a
spin-degenerate fixed state density over temperature. It explicitly separates
the exact sampled-state solve from broadened thermal visualization.

Exact donor identity, execution evidence, the observed histogram-bin shift, and
adaptation decisions are recorded in
`docs/provenance/weighted-fermi-dirac-thermal-dos-library.md`. No
zero-temperature Fermi-energy convention was inferred or added.

## Fifty-second disposition: extended cosine thermal prototype

`notebooks/solidstate/bloch1d/cosinepotential/exploratory-bands-dos-thermal-prototype.ipynb`
was removed after the weighted chemical-potential and thermal-DOS extraction.
Every reviewed reusable and pedagogical layer has a maintained destination, and
the prototype's exact bytes remain recoverable from preservation commit
`fbf7e9968671cffc18760e9a4aadaf165843d4cc` and Git blob
`687879784df581a4308104614ff59bc796f5e55f`.

The removal does not promote the prototype's duplicated endpoint, abstract
six-band normalization, histogram-bin chemical-potential inversion, or broken
clean linear execution. Final disposition is recorded in
`docs/provenance/weighted-fermi-dirac-thermal-dos-library.md`.

## Fifty-third extraction: sampled one-dimensional Bloch-state visualization

The exploratory free-electron finite-difference notebook was reviewed
individually. Its seven static nonempty code cells execute from a clean
namespace; execution then stops because its saved widget cells import
`ipywidgets`, which is not a project dependency. The artifact remains present
pending final disposition.

`BlochState1DPlotData` and `BlochState1DVisualizer` now own unit-aware complex
amplitude, probability-density, phase, and cell-periodic Fourier rendering. The
maintained laboratory exposes reduced wave number and band index as ordinary
editable inputs without requiring a widget runtime.

The extraction corrects two presentation defects. Amplitude and probability
density no longer share one ordinate despite having different physical units,
and the discrete Fourier transform is applied to
$u_k=e^{-ikx}\psi_k$ after removing the Bloch carrier rather than directly to a
twisted nonzero-$k$ state. Exact source identity, execution evidence, destination
mapping, and remaining disposition are recorded in
`docs/provenance/free-electron-bloch-prototype-review.md`.
