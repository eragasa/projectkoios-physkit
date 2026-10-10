# Changelog

## Unreleased

- Add unit-aware, quotient-seam Bloch-periodic finite-difference Laplacians for
  intrinsic nonorthogonal 2D and full nonorthogonal 3D primitive cells, including
  inverse-metric mixed derivatives and scale-invariant 2D basis validation.
- Add an explicit ambient-3D/intrinsic-2D adapter for full-rank `3×2` direct
  bases, preserving frame, projector, dual-basis, and coordinate maps, together
  with a strict deterministic domain JSON codec for the defining basis and unit.

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
