# Finite-periodic lattice source mapping

This record binds the reusable finite-periodic lattice implementation to exact
`ksdft2effmass` donor source. The donor repository is
`https://github.com/eragasa/ksdft2effmass`, licensed under Apache-2.0. The
extraction preserves reusable represented mathematics while removing
semiconductor campaigns, retained calculation state, scientific acceptance
policy, and independent quotient-seam construction.

## Source identities

| Donor revision | Donor path | Git blob | SHA-256 | Destination module |
| --- | --- | --- | --- | --- |
| `7bd913151f7e61ed2bdba593df920be36573b502` | `python/src/ksdft2effmass/solid_state/geometry.py` | `726cbf49b49c4d5176f9cfdf39ff1ae9c1f7298c` | `80f1cfebc96ff3bb30c971a07dd07536ac55181598e958473af611730d691605` | `projectkoios.physkit.periodic.lattice.finite_domain` |
| `7bd913151f7e61ed2bdba593df920be36573b502` | `python/src/ksdft2effmass/solid_state/boundary_phases.py` | `ad14c46e8a08db30cd4941fb9ebb357f78132b49` | `1b08f34895ca007bf7a84bde71c97a7c12d65aab1f34512abaaccdb7911d8cee` | `projectkoios.physkit.periodic.lattice.boundary_phases` |
| `7bd913151f7e61ed2bdba593df920be36573b502` | `python/src/ksdft2effmass/solid_state/symmetry.py` | `2e62d76262094217bac767ca81c6acd2740591f5` | `938167c287cc604d8b92dd101b9a37ac561c074c1fb540a03c06207ef4a8320b` | `projectkoios.physkit.periodic.lattice.symmetry` |
| `7bd913151f7e61ed2bdba593df920be36573b502` | `python/src/ksdft2effmass/solid_state/lattice_models.py` | `4eeae126560345a1d8914670793fcfe2ea06faaf` | `387f05240eef2f695e3f520c0bcf907c411807ca81e301fcf56b3058f647ec28` | `projectkoios.physkit.periodic.lattice.hopping` |
| `7bd913151f7e61ed2bdba593df920be36573b502` | `python/src/ksdft2effmass/solid_state/represented_operators.py` | `d48c1f217df8a1085bcad6e1935fa118b236968b` | `535296e6cd952d7e416109cf9a834a17eb9958fe5ced606f69ea84ea36ce8050` | `projectkoios.physkit.periodic.lattice.represented_operators` |
| `2578d398b1d0aa3bfafb45f73b39c79ede5f5c42` | `python/src/ksdft2effmass/solid_state/operator_composition.py` | `481dcff3984f86b43d72daa93a0b4e3b37996c9c` | `2134c1ed4863a746b292cc0b398a32ed67b58ff9df71b42c5fbc9accb565bd3a` | `projectkoios.physkit.periodic.lattice.operator_composition` |
| `fcca7be70845aebe016dd6ba970f54d0bf8a879c` | `python/src/ksdft2effmass/solid_state/operator_construction.py` | `e85af86c298e7c488a53828e3517b44bb9a638b8` | `ce3f10f4249ddcea0310888bfb37885f16d678e0132c1ca4a960a026f2cc6d2c` | `projectkoios.physkit.periodic.lattice.operator_construction` |
| `a3a569064033cfa508d4f9ccbb0509f64745f8d6` | `python/src/ksdft2effmass/solid_state/gauge_bridges.py` | `a10a010b7457280e64d202ea8885c8d0cfe6b169` | `87a73884b7e222b2292ac17d00b82fe1d4942b21eae649a5b04b6fa0e7ccbe1f` | `projectkoios.physkit.periodic.lattice.gauge_bridges` |

## Destination-authored extensions

`projectkoios.physkit.solidstate.tb.diagnostics` was authored at the PhysKit
destination against the extracted `ScalarHoppingModel` and
`IntegralLatticeOperation` contracts. It does not claim exact implementation
lineage from a donor module. The module owns raw tolerance-free scalar and block
Hermiticity diagnostics plus displacement-correlated unit-aware block powers.

`BlockHoppingTerm` and `BlockHoppingModel` generalize the destination hopping
records to energy-unit-bearing square matrices with explicit ordered-basis and
energy-reference identity. `projectkoios.physkit.solidstate.tb.truncation`
selects a caller-supplied ordinary-integer-lattice conjugate-closed displacement
inventory with a distinct retained-model identity, without owning a cyclic
finite-mesh alias rule, range, cutoff, tolerance, normalization, or acceptance
policy. These contracts
were authored at the PhysKit destination in response to an internal KSD ownership
handoff and do not claim copied implementation lineage from a donor module.

`projectkoios.physkit.solidstate.tb.bloch` was also authored at the PhysKit
destination. It owns the physics-semantic mapping from fractional reciprocal
wave vectors and energy-unit-bearing scalar or block hopping coefficients to
correlated Bloch values or matrices. Its positive-exponent convention agrees
with the donor's hopping transform convention, but the adapters do not claim
copied implementation lineage from that module.

`projectkoios.physkit.solidstate.tb.bands` was authored at the PhysKit destination
as the unit-aware adapter from block Bloch matrices to complete ascending bands.
Its public model, coordinate, tolerance, diagnostic, energy-unit, and backend-
evidence records do not claim donor implementation lineage. Generic unit-free
Hermiticity validation, Hermitian-representative construction, diagnostics, and
complete spectra delegate directly to the exact pinned Numerics owner.

`projectkoios.physkit.solidstate.tb.finite_periodic` was authored at the PhysKit
destination as the physical boundary around generic finite-periodic matrix Fourier
operations. It correlates samples, complete or explicitly selected coefficients,
tensor-grid evaluations, Parseval evidence, and cyclic Hermiticity evidence with exact
energy units, energy references, and ordered-basis identities. Computed records are
created only by their owning PhysKit actions, and each diagnostic action invokes its
Numerics owner once. The package does not claim copied implementation lineage, choose
a support, or duplicate the generic transform, evaluation, pairing, or reduction
arithmetic owned by Numerics.

## Numerics owner binding

The generic numerical calculations are consumed from the
`projectkoios-optimization` distribution at exact commit
`d38e00ceef7f0fb1c137d5c3450d210c1ae4c99f`, tree
`6fe28b2fe52e80e5a1eb0199092da42770e9fb44`.

- `projectkoios.numerics.finite_periodic.finite_series` fixes
  $\sum_j c_j\exp(+2\pi i\,x\cdot m_j)$ without normalization.
- `projectkoios.numerics.finite_periodic.batched_finite_series` evaluates ordered
  coefficient channels on one shared coordinate/mode grid with point-major output.
- `projectkoios.numerics.finite_periodic.numpy.batched_finite_series` owns one dense
  shared-phase NumPy reduction and its backend evidence.
- `projectkoios.numerics.finite_periodic.meshes` owns centered odd/even tensor meshes,
  dimensionless offsets, origins, and parity-aware integer inventories.
- `projectkoios.numerics.finite_periodic.numpy.transforms` owns complete normalized
  centered-mesh matrix transforms for odd and even sizes.
- `projectkoios.numerics.finite_periodic.numpy.tensor_product_evaluation` owns bounded
  separable evaluation of complete or selected matrix series.
- `projectkoios.numerics.finite_periodic.numpy.transform_diagnostics` owns unit-free
  Parseval reductions and origin-aware cyclic Hermiticity pairing.
- `projectkoios.numerics.linalg.frobenius` returns deterministic raw power and its
  square root for finite built-in complex tuple matrices.
- `projectkoios.numerics.linalg.hermitian_spectrum` owns normalized batch
  definitions, tolerances, diagnostics, results, validation, and explicit
  Hermitian representatives.
- `projectkoios.numerics.linalg.numpy.hermitian_spectrum` owns the NumPy complete-
  spectrum implementation and its single nonempty batched `eigvalsh` call.
- PhysKit passes only unit-free built-in coordinates, integer modes, coefficients,
  matrices, and numerical tolerances to those owners.
- Lattice interpretation, model identity, basis identity, energy reference, and
  energy or squared-energy units remain in PhysKit adapters and results.

## Semantic adaptations

- `FiniteLatticeShape` became `FinitePeriodicDomain` because the record is an
  integer quotient domain rather than a physical direct lattice.
- Legacy ksdft compatibility names remain aliases in ksdft and are not duplicate
  nominal types in this distribution.
- Represented-operator finite geometry uses the field name `domain`.
- Compatibility issue code `SHAPE = "shape"` became
  `DOMAIN = "domain"`.
- ksdft quantity records were replaced by the corresponding native
  `projectkoios.physkit` quantity records before dependent behavior moved.
- Sparse construction, addition, and gauge comparison preserve CSR operation
  without implicit densification.
- The repository and distribution moved to Project Koios ownership under
  `projectkoios-physkit`; the import namespace is `projectkoios.physkit`.

## Exclusions

The extraction does not move:

- `QuotientSeamOperatorConstructor` or the independent seam route;
- cross-route reconciliation workflows;
- campaign tolerances or acceptance policy;
- retained calculations, provenance-bound execution scripts, or authorization
  records;
- semiconductor-specific physical interpretation; or
- scientific-validation or human-acceptance conclusions.

Those surfaces remain owned by `ksdft2effmass`.

## Verification mapping

The destination tests under
`tests/projectkoios/physkit/periodic/lattice/` cover finite-domain geometry,
boundary phases, symmetry operations, hopping records, represented operators,
compatibility-gated sparse addition, centered uniform-link construction, and
site-diagonal gauge bridges. Tests under
`tests/projectkoios/physkit/solidstate/tb/` cover destination-authored raw
scalar and block hopping diagnostics and the physics-semantic adapters to the exact
Numerics scalar, batched and finite-periodic matrix-series, Frobenius, and
Hermitian-spectrum owner revision. Finite-periodic adapter tests cover odd and even
inventories, explicit units and squared units, immutable correlation, reconstruction,
selected supports, and empty supports. The corresponding ksdft integration tests
exercise its compatibility
imports and independent seam-route reconciliation.

These tests establish only their documented software and bounded numerical
requirements. They do not establish physical-model adequacy, semiconductor
scientific validation, uncertainty quantification, publication authority, or
human acceptance.
