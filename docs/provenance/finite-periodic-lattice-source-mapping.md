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
lineage from a donor module. The module owns raw tolerance-free Hermiticity and
caller-supplied integral-operation comparisons.

`projectkoios.physkit.solidstate.tb.bloch` was also authored at the PhysKit
destination. It owns the physics-semantic mapping from fractional reciprocal
wave vectors and energy-unit-bearing scalar hopping coefficients to correlated
Bloch values. Its positive-exponent convention agrees with the donor's hopping
transform convention, but the adapter does not claim copied implementation
lineage from that module.

## Numerics owner binding

The generic ordered finite Fourier-series calculation is consumed from the
`projectkoios-optimization` distribution at exact commit
`23e1c88772f1bc9a1505cd25df5418cd945a76fa`, tree
`ef48b2f531d89a720f5bf5cdd9a64cf6dfb17a37`. The owning API is
`projectkoios.numerics.finite_periodic.finite_series` and fixes
$\sum_j c_j\exp(+2\pi i\,x\cdot m_j)$ without normalization.

PhysKit passes only unit-free built-in coordinate, integer-mode, and coefficient
tuples to that owner. Lattice interpretation, model identity, basis identity,
energy reference, and energy unit remain in the PhysKit adapter and result.

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
scalar hopping diagnostics and the physics-semantic adapter to the exact
Numerics finite-series owner revision. The corresponding ksdft integration
tests exercise its compatibility imports and independent seam-route
reconciliation.

These tests establish only their documented software and bounded numerical
requirements. They do not establish physical-model adequacy, semiconductor
scientific validation, uncertainty quantification, publication authority, or
human acceptance.
