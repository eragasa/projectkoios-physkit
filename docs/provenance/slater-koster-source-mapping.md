# Slater–Koster source mapping

This mapping records the reusable Slater–Koster migration tranches from the
`ksdft2effmass` calculation incubation into PhysKit, including the later reimplementation
of the linear parameterization and reciprocal-sample tranches that were not present in
the first published migration.

## Donor context

- Repository worktree: `ksdft2effmass-source-cleanup-audit`
- Initial revision context: `6fdaeea62207be4f74f44138d93b3b2298c1cbe2`
- Linear-donor reimplementation context: `784406e3` plus the uncommitted donor hashes
  listed below on branch `work/source-code-cleanup-audit`.
- Branch context: `work/source-code-cleanup-audit`
- Status: the donor files were uncommitted incubation files, so their SHA-256 hashes,
  rather than the revision alone, establish the exact source files.

## Exact file mapping

| Donor path | SHA-256 | PhysKit destination |
|---|---|---|
| `calculations/Si/bulk/solidstate/tb/basis_centers.py` | `f98f66c57c70a9b3a38ffb56755f524171188dbbe5d8fafc849bfc3eac749b5f` | `src/python/projectkoios/physkit/solidstate/tb/orbital_center_gauge.py` |
| `calculations/Si/bulk/solidstate/tb/slater_koster/basis.py` | `242bd936bad6441f55e467d2769a658496a5b3f93cca73757174c01b5be2515e` | `src/python/projectkoios/physkit/solidstate/tb/slater_koster/basis.py` |
| `calculations/Si/bulk/solidstate/tb/slater_koster/bond.py` | `d73241477f5d791a196f5166cb8f5e66a7fc7f5e32a4e4e29e6802acd9ac7111` | `src/python/projectkoios/physkit/solidstate/tb/slater_koster/bond.py` |
| `calculations/Si/bulk/solidstate/tb/slater_koster/channels.py` | `66f564cd3814ffd9c5c0e13955bb128bdab17d1617e558f84760ac6659582f24` | `src/python/projectkoios/physkit/solidstate/tb/slater_koster/channels.py` |
| `calculations/Si/bulk/solidstate/tb/slater_koster/parameters.py` | `bb4881a56f60e3b30e371333b6ae17b6f89d00f3ce8255ca8a401879567f047d` | `src/python/projectkoios/physkit/solidstate/tb/slater_koster/parameters.py` |
| `calculations/Si/bulk/solidstate/tb/slater_koster/sp3s_star.py` | `687fa6adca05db5b57ff9e33f740a36cb7d852d97436f72af6b614e4fd74942d` | `src/python/projectkoios/physkit/solidstate/tb/slater_koster/sp3s_star.py` |
| `calculations/Si/bulk/solidstate/tb/slater_koster/assembly.py` | `74370af68b4dc9af77dc6f90dcf3396e99798fc54a3135c5c65c55a2fb82c621` | `src/python/projectkoios/physkit/solidstate/tb/slater_koster/assembly.py` |
| `calculations/Si/bulk/solidstate/tb/slater_koster/linear_parameterization.py` | `61f21bb8014fb9d8cdc50230d71f22978a343908e2d635cc8d7921be9b7a700e` | `src/python/projectkoios/physkit/solidstate/tb/slater_koster/linear_parameterization.py` |
| `calculations/Si/bulk/solidstate/tb/slater_koster/linear_bloch.py` | `a26e64e08b24976065aac0ea740c36f702f8ed64b2262f3e0047158b27d26efa` | `src/python/projectkoios/physkit/solidstate/tb/slater_koster/linear_bloch.py` |
| `calculations/Si/bulk/solidstate/tb/slater_koster/mathematics.md` | `1a789ad1937d1868db2156c47aff647102e09e8ad190dea3dc048bd3190a9845` | `docs/architecture/projectkoios/physkit/solidstate/tb/slater_koster/` |
| `calculations/Si/bulk/solidstate/tb/slater_koster/README.md` | `dd3cd639c3b65b306eecd15cb5f01296586a94813cd2a492c43dac3ea19a43b1` | linear-module architecture and migration-boundary documentation |

## Deliberate destination adaptations

- The donor's direction-cosine input became the nominal immutable
  `SlaterKosterBondDirection`; its Cartesian components enter the formulas as
  direction cosines, and validation is not hidden behind an empty validator object.
- Computed channel-basis, hopping, and gauge values are correlated
  `ResultsObject` records.
- The generic orbital-center gauge transformation remains outside the
  Slater–Koster specialization.
- Assembly uses explicit `SlaterKosterDirectedBond` values rather than donor
  tuple-based bond keys, preserves exact basis and parameter-group correlation, and
  delegates the resulting coefficient inventory to PhysKit's existing dense
  block-hopping owner.
- The linear destination reuses PhysKit's internal block-accumulation and term-freezing
  operations rather than retaining the calculation donor's copied private assembly
  functions.
- Linear Bloch evaluation delegates the compatible channel-model family to PhysKit's
  existing dense block-hopping evaluator and restores group/channel/sample nesting; it
  does not flatten operator entries or introduce a fitting owner.
- Material selection, neighbor-shell selection, fitting, residual construction,
  campaign policy, evidence, and scientific disposition were not migrated.
- Physical lengths and energies use the existing `UnitSystem`; the default physical
  convention is `UnitSystem.ATOMIC`, while fractional coordinates and angular
  coefficient matrices use the existing nondimensional unit-system view.
- The two-center formulas, fixed `(s, p_x, p_y, p_z, s*)` order, seven-channel
  order, positive Bloch/gauge phase convention, and transpose-under-bond-reversal
  convention are preserved.
