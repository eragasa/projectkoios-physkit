# Package `projectkoios.physkit.units`

## Current responsibility

The package owns reusable physical-unit and immutable quantity semantics.

- `PhysicalUnit` and `Unitless` are the nominal unit identities for model-system data.
- Scalar, real/complex vector, real/complex matrix, and real/complex sparse matrix
  quantities retain explicit units and immutable numerical storage.
- `PintUnitConverter` owns compatibility checks, conversion factors, and represented
  conversions through one fixed registry.
- Domain-specific unit helpers and established unit-system records remain available
  through the package façade.
- KSD local quantity classes are legacy consumers/donor evidence, not a second owner.

## Nominal boundary

Runtime quantity ownership is explicit rather than structural.

- `ModelSystemUnit` is the nominal base for physical and explicitly unitless values.
- `ModelSystemQuantity` is the nominal base for immutable unit-bearing numerical values.
- Public numerical records use exact concrete quantity identities where correlation
  depends on ownership.
- Consumers may deliberately re-export exact PhysKit identities during migration.
- Wrappers, duck-typed bridges, and a PhysKit dependency on KSD are not owner contracts.

## Conversion boundary

Unit conversion mechanics remain separate from scientific normalization.

- Pint dimensional compatibility and physical conversion belong to PhysKit.
- Affine scalar/vector conversions preserve offset-unit behavior.
- Complex and matrix conversions use compatible multiplicative factors where affine
  semantics are not meaningful.
- Calculator-native scaling fields, external-file conventions, campaign normalization,
  and serialization compatibility remain with their consumer owners.
- Conversion success does not establish scientific comparability or acceptance.

## Migration boundary

KSD quantity adoption must remain a separately reviewed consumer migration.

- Existing KSD science may continue using its frozen local quantities until a bounded
  migration is authorized.
- Exact nominal incompatibility is an explicit integration blocker, not a reason to
  weaken PhysKit validation.
- Migration must audit exact-type checks, persisted schemas, converter behavior,
  fixtures, and retained numerical evidence.
- No crystallography or path tranche authorizes repository-wide quantity replacement.

## Local mapping

Maintained source and tests define the supported unit surface.

- Public façade: `src/python/projectkoios/physkit/units/__init__.py`
- Nominal quantities and converter: `src/python/projectkoios/physkit/units/quantities.py`
- Tests: `tests/projectkoios/physkit/units/`
- [Parent package](../index.md)
