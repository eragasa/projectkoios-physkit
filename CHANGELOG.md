# Changelog

## Unreleased

- Add raw scalar-hopping Hermiticity and supplied-operation diagnostics plus a
  unit-aware Bloch evaluator backed by the exact reviewed
  `projectkoios.numerics` finite Fourier-series owner.
- Advance the exact `projectkoios-optimization` dependency revision while
  preserving the finite Fourier-series implementation consumed by PhysKit.

## 0.1.2 - 2026-10-11

- Add immutable unit-aware stress tensors with explicit pressure-unit and sign-
  convention conversion and derived hydrostatic pressure.
- Add unit-aware, quotient-seam Bloch-periodic finite-difference Laplacians for
  intrinsic nonorthogonal 2D and full nonorthogonal 3D primitive cells, including
  inverse-metric mixed derivatives and scale-invariant 2D basis validation.
- Add an explicit ambient-3D/intrinsic-2D adapter for full-rank `3×2` direct
  bases, preserving frame, projector, dual-basis, and coordinate maps, together
  with a strict deterministic domain JSON codec for the defining basis and unit.
- Add a scale-invariant `BlochCell2D` contract for metal and explicitly
  nondimensional checkerboard models, canonical sparse plane-wave Hamiltonians,
  automatic sparse partial-spectrum solving, and explicit dense fallback for
  complete spectra.
- Add unit-aware periodic scalar-potential contracts, one-, two-, and
  three-dimensional cosine-product models, and explicit sampled-diagonal and
  reciprocal-Fourier sparse operator converters; route the maintained cosine and
  checkerboard Bloch models through the reusable potential boundary.
- Add a model-backed one-dimensional periodic-potential visualizer and a
  dimension-independent band-path visualizer with unit- and shape-checked
  optional reference bands.
- Add shared equal-weight full-zone DOS estimators with explicit one-, two-, and
  three-dimensional normalization, Gaussian and histogram representations,
  state-count results, and dimension-independent visualization; migrate the
  maintained 1D, 2D, and 3D laboratories to these APIs.
- Extend the Fermi--Dirac chemical-potential solver with explicit unit-aware
  state weights, add occupied, unoccupied, thermal-window, and thermally active
  DOS results, and migrate the maintained 1D laboratory to exact sampled-state
  inversion rather than histogram-bin inversion.
- Add unit-aware one-dimensional sampled Bloch-state visualization with
  separate amplitude, probability-density, and phase axes and Fourier weights
  computed after removing the Bloch carrier.
- Add a maintained nondimensional oblique periodic free-particle laboratory
  using the package-owned endpoint-excluded inverse-metric Laplacian, Bloch
  quotient seams, sparse eigensolving, plane-wave verification, and
  Jacobian-aware state normalization.
- Add a unit-aware bounded `ObliquePiab2D` model, sparse homogeneous-Dirichlet
  inverse-metric finite-difference solver, physical-area state normalization,
  and maintained hard-wall parallelogram laboratory.
- Remove the superseded generalized oblique-coordinate prototype after its
  bounded and periodic routes received separate maintained destinations.
- Remove the superseded oblique periodic-boundary prototype after separating
  its periodic intent from its bounded saved-output behavior.
- Remove the superseded free-electron Bloch prototype after migrating its bands,
  sampled-state visualization, and corrected cell-periodic Fourier analysis.

## 0.1.1 - 2026-10-10

- Extend `UnitCellJsonCodec` with the version-one `unit-cell` representation so
  an exact base `UnitCell` round-trips without primitive or conventional
  inference while preserving the existing representation bytes and decode types.
- Document the canonical Project Koios source and test namespaces, permit a
  path-mirrored transitional `src/physkit/` tree, and keep notebook organization
  independent under pedagogical ownership.
- Consolidate duplicate square-wave Fourier-series notebooks into one maintained
  computational laboratory and extract reusable odd-harmonic evaluation into
  `projectkoios.physkit.math.fourier_series`.
- Consolidate overlapping Maxwell–Boltzmann speed-distribution notebooks and
  extract explicit SI gas state, characteristic speeds, and speed-density
  evaluation into `projectkoios.physkit.thermal.statmech.kinetic_theory`.
- Consolidate separate two- and three-dimensional von Mises notebooks and
  extract equivalent-stress evaluators into
  `projectkoios.physkit.mechanics.continuum_stress`.
- Move the Townsend discharge notebook out of package source and extract its
  state, result, closed-form, generational, and fixed-point current behavior into
  `projectkoios.physkit.plasmas.gas_discharge.townsend`.
- Replace two byte-identical, incomplete particle-in-a-box notebooks with one
  maintained analytical and finite-difference comparison laboratory.
- Replace two byte-identical tight-binding explorations with one maintained
  nearest-neighbor finite-periodic laboratory using the represented lattice APIs.
- Remove three zero- or one-byte invalid notebook placeholders and add a
  repository check that every remaining notebook parses as JSON.
- Replace two notebook-local direct-lattice prototypes with one maintained
  direct-and-reciprocal-lattice laboratory using the physical lattice APIs.
- Replace the incomplete vapor-pressure package with explicit Antoine model,
  evaluator, and result owners plus a synthetic-data computational laboratory.
- Fold the remaining obsolete numerical and symbolic PIAB1D notebooks into the
  maintained lecture note and analytical/finite-difference laboratory.
- Add a unit-aware analytical two-dimensional rectangular particle-in-a-box
  model, ordered state result, evaluator, tests, note, and laboratory.
- Add a maintained analytical PIAB1D spectral time-evolution laboratory with
  norm, energy, and probability-density recurrence checks.
- Extract finite-level canonical and Fermi–Dirac occupations, chemical-potential
  solving, stationary-orbital evaluation, and occupied-orbital density analysis
  from historical PIAB1D visualizations into tested package APIs and a maintained
  laboratory.
- Extract Gaussian-state construction, weighted wavefunction normalization, and
  sampled basis projection into tested APIs and a maintained PIAB laboratory.
- Organize maintained particle-in-a-box laboratories and corresponding lecture
  notes by dimension under `qm/piab1d/`, `qm/piab2d/`, and `qm/piab3d/`.
- Add a unit-aware analytical three-dimensional rectangular particle-in-a-box
  model, ordered state result, evaluator, tests, lecture note, architecture
  documentation, and computational laboratory.
- Add tested analytical PIAB2D eigenfunction-grid evaluation and PIAB3D
  Cartesian-plane evaluation with probability-density visualization laboratories.
- Organize PIAB2D stationary analytical work beneath `tise/analytical/` and
  reserve `tdse/numerical/` without claiming an unimplemented solver.
- Rename the repository and distribution to `projectkoios-physkit`, move the
  Python package to the PEP 420 namespace `projectkoios.physkit`, and mirror
  source, tests, and current-state architecture documentation under Project
  Koios ownership paths.
- Add finite-periodic domains, boundary twists, scalar hopping and localized
  perturbation models, integral lattice operations, represented scalar operators,
  compatibility-gated sparse operator composition, centered uniform-link sparse
  operator construction, and site-diagonal twist-gauge bridge analysis.
- Consolidate reviewed two- and three-dimensional metric-tensor explorations into
  one package-backed nonorthogonal lattice lesson, add range-stable measure and
  scale-aware exact nearest-image metric geometry, and add unit-free sparse
  2D/3D Kronecker-sum assembly for floating and complex-floating later
  multidimensional represented operators.
- Relicense PhysKit from the MIT License to the Apache License 2.0 with the
  authorization of the sole copyright holder. Copies previously distributed
  under the MIT License remain available under those terms.
