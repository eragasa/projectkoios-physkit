# Module `projectkoios.physkit.periodic.unit_cell.serialization`

## Current responsibility

The module owns the version-one deterministic JSON representation of supported
periodic unit-cell types.

## Public contract

The schema records `schema_version`, caller-supplied `structure_id`, an explicit
representation discriminator, dimensionless direct-lattice vectors ordered as
columns, a physical lattice-parameter magnitude and unit, and the ordered
atomic basis in fractional coordinates.

The exact supported nominal types and discriminators are:

- `UnitCell` as `unit-cell`;
- `ConventionalUnitCell` as `conventional`; and
- `PrimitiveUnitCell` as `primitive`.

Encoding uses exact runtime types. Decoding selects the exact type solely from
the discriminator; it never infers primitive or conventional semantics from
geometry or atoms. Unsupported types, aliases, schema versions, identities,
vector conventions, units, and malformed values are rejected with
`UnitCellSerializationError`.

## Navigation

- [Class `UnitCellJsonCodec`](UnitCellJsonCodec/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/unit_cell/serialization.py`
- Tests: `tests/projectkoios/physkit/periodic/unit_cell/serialization/`

## Evidence

Focused tests verify byte-stable existing primitive and conventional payloads,
deterministic base-cell encoding, exact nominal decode types, vector and atom
ordering, units, strict discriminator handling, and malformed-input failures.
They establish the declared serialization behavior, not external compatibility
or physical validation.
