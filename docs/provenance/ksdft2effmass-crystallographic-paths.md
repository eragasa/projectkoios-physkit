# KSD crystallographic-path extraction

## Source checkpoint

The first bounded path tranche adapts frozen KSD donor evidence into its PhysKit owner.

- Donor repository: `https://github.com/eragasa/ksdft2effmass`
- Donor base commit: `8aea8232a2d4bc33af095cac9a32bd3a564a7cee`
- Donor base tree: `562a5d4bdd97681c23cdff37d4f15e2a162a3e95`
- Donor source working-tree blob:
  `d7161bc90bc6f672fa4333d2b43e341099c75c05`
- Donor source SHA-256:
  `b656deb3bfd31758205248ddbdbe9274215a31e62e79d1e51dc54e857f467bb4`
- Donor source path:
  `python/src/ksdft2effmass/analysis/model_adequacy/band_equivalence_locality/strain.py`
- Donor test working-tree blob:
  `14126a86de062f45f31ce56c08cd0407e5839601`
- Donor test SHA-256:
  `05ceefbcb1e7ee47b248e096a859dfa2b58faccd84588f0f5868271e446d0e68`
- Donor test path:
  `python/tests/software_verification/ksdft2effmass/analysis/model_adequacy/band_equivalence_locality/test__BainPath.py`
- Donor bridge/decoration test working-tree blob:
  `49141506ea851d10d8e928198f8abe104aea8035`
- Donor bridge/decoration test SHA-256:
  `fcc1f8da6bbed91a5ac3c182477973edf2a13ecd2e036f0ab639545e8c6d3118`
- Donor bridge/decoration test path:
  `python/tests/software_verification/ksdft2effmass/analysis/model_adequacy/band_equivalence_locality/test__SyntheticStrainGaugeLocality.py`

## Extracted behavior

The destination preserves the donor's represented Bain geometry and numerical evidence.

- Axial ratio varies linearly from `sqrt(2)` to one on the closed parameter interval.
- Transverse scale is `(2 / axial_ratio)^(1 / 3)`.
- Longitudinal scale is `axial_ratio * transverse_scale`.
- The ordered body-centered-tetragonal primitive vectors are retained exactly.
- A 201-point scan retains the unit primitive-volume invariant at binary64 precision.
- Metric determinants remain correlated to squared primitive volume.

## Owner-native adaptations

PhysKit names and validates the convention as an independent owner contract.

- Donor `BainPath` maps to `UnitVolumeBainPath3D`.
- Donor `BainPathPoint` maps to `UnitVolumeBainPathPoint3D`.
- `strain_fraction` becomes the neutral `path_parameter`; `aspect_ratio` becomes the
  crystallographic `axial_ratio`.
- Wrong runtime parameter types raise `TypeError`; invalid values raise `ValueError`.
- The point validates complete scalar/lattice correlation and owns a direct-lattice
  copy.
- Public exports, architecture pages, and isolated-wheel evidence are destination-owned.

## Decorated-basis and coordinate-bridge extraction

A second bounded tranche adapts the frozen crystallographic constants into typed owner
records.

- `BAIN_FCC_UNIMODULAR_MAP` maps to `FccBainCoordinateBridge3D` with its exact inverse,
  integer/fractional actions, and finite-quotient evidence.
- `BAIN_DIAMOND_BASIS_OFFSETS` maps to the immutable, explicitly unitless basis owned by
  `DiamondDecoratedBainPathPoint3D`.
- `DiamondDecoratedBainPath3D` composes the decoration with the previously owned Bravais
  path without attaching species, scale, or electronic response.
- The area-preserving rhombic path remains a separate 2D owner tranche.
- Synthetic strain diagnostics, electronic response models, selected bonds, Fourier
  calculations, thresholds, artifacts, and interpretations remain KSD-owned.

## Evidence limits

Extraction establishes represented software behavior only.

- It does not establish a material transformation mechanism, activation barrier,
  kinetics, stress response, temperature dependence, or phase stability.
- Endpoint labels are crystallographic identities under the declared basis convention.
- Scientific use and acceptance require separate evidence and human review.
