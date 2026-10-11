# Embedded-lattice JSON serialization adaptation

This record identifies the strict JSON wire-boundary implementation consulted
while designing PhysKit's domain-specific `EmbeddedDirectLattice2DJsonCodec`.
The donor repository remained read-only.

## Reviewed source identity

| Repository | Revision | Source path | Git blob | SHA-256 |
| --- | --- | --- | --- | --- |
| `ksdft2effmass` | `1c9a2f7719b64758a1e7ede4b867acb420a2b6ab` | `python/src/ksdft2effmass/operators/serialization.py` | `a095a5e072f0de7938a089ba567ce7945997c081` | `7d62165a5f7c8c8e74de36640d61bbe795c28504b9a932fa5276dc3a3c20fa93` |

## Retained strict-boundary patterns

The PhysKit codec adapts these reviewed mechanics:

- deterministic sorted JSON output with `allow_nan=False`;
- explicit schema version and representation discriminator;
- exact required object fields with unknown and missing fields rejected;
- duplicate-member rejection through `object_pairs_hook`;
- rejection of nonstandard `NaN` and infinity constants;
- rejection of booleans and numeric strings where numbers are required; and
- typed reconstruction through public domain constructors.

## PhysKit-specific adaptations

- Ownership remains adjacent to `projectkoios.physkit.periodic.lattice.embedding`;
  no runtime dependency on ksdft2effmass or a generic serializer framework was
  introduced.
- The fixed schema stores only two three-component ambient primitive-vector
  columns and an explicit physical/unitless unit discriminator.
- Input size and unit-expression length are bounded.
- Derived `G`, `G^{-1}`, `B`, `Q`, and `P` values are absent from the wire and
  reconstructed through `EmbeddedDirectLattice2DAdapter`.
- The public error and codec names follow PhysKit ownership rather than donor
  operator-record names.
- Existing `UnitCellJsonCodec` behavior and bytes were not changed.

## Excluded behavior

No legacy-compatible serializer that admits non-finite values was used. No
operator-record matrix schema, provenance model, package-global registry,
generic dataclass introspection, or cross-language conformance claim was
migrated.

## Evidence boundary

Golden, strict-failure, exact-type, and derived-reconstruction tests establish
the declared PhysKit schema behavior. They do not establish scientific
provenance or material validity for an encoded basis.
