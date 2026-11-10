# Slater–Koster source mapping

This mapping records the first reusable Slater–Koster migration tranche from the
`ksdft2effmass` calculation incubation into PhysKit.

## Donor context

- Repository worktree: `ksdft2effmass-source-cleanup-audit`
- Revision context: `6fdaeea62207be4f74f44138d93b3b2298c1cbe2`
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
| `calculations/Si/bulk/solidstate/tb/slater_koster/mathematics.md` | `1a789ad1937d1868db2156c47aff647102e09e8ad190dea3dc048bd3190a9845` | `docs/architecture/projectkoios/physkit/solidstate/tb/slater_koster/` |

## Deliberate destination adaptations

- The donor's direction-cosine input became the nominal immutable
  `SlaterKosterBondDirection`; its Cartesian components enter the formulas as
  direction cosines, and validation is not hidden behind an empty validator object.
- Computed channel-basis, hopping, and gauge values are correlated
  `ResultsObject` records.
- The generic orbital-center gauge transformation remains outside the
  Slater–Koster specialization.
- Material selection, neighbor-shell selection, fitting, residual construction,
  campaign policy, evidence, and scientific disposition were not migrated.
- Physical lengths and energies use the existing `UnitSystem`; the default physical
  convention is `UnitSystem.ATOMIC`, while fractional coordinates and angular
  coefficient matrices use the existing nondimensional unit-system view.
- The two-center formulas, fixed `(s, p_x, p_y, p_z, s*)` order, seven-channel
  order, positive Bloch/gauge phase convention, and transpose-under-bond-reversal
  convention are preserved.
