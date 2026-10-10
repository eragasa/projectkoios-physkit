# Class `projectkoios.physkit.periodic.unit_cell.serialization.UnitCellJsonCodec`

## Current contract

`UnitCellJsonCodec` encodes and decodes version-one unit-cell JSON. `dumps`
requires a nonempty structure identity and an exact `UnitCell`,
`ConventionalUnitCell`, or `PrimitiveUnitCell`; it returns sorted, indented JSON
with one trailing newline. `loads` requires the expected structure identity,
validates the schema and represented values, and returns the exact nominal type
named by `unit-cell`, `conventional`, or `primitive`.

## Invariants

- Existing `primitive` and `conventional` serialized bytes and decode types are
  unchanged.
- `unit-cell` decodes to exact type `UnitCell`, never to a semantic subclass.
- Direct-lattice vectors remain `a1`, `a2`, and `a3` in column order.
- The physical lattice-parameter unit and ordered atomic basis are retained.
- The schema version remains `1`; representation aliases are not accepted.
- Decode and validation failures cross the API boundary as
  `UnitCellSerializationError` while preserving applicable causes.

## Navigation

- [Parent module](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/unit_cell/serialization.py::UnitCellJsonCodec`
- Tests: `tests/projectkoios/physkit/periodic/unit_cell/serialization/test__UnitCellJsonCodec__*.py`

## Evidence

Focused positive and adversarial tests cover deterministic bytes, all three
exact nominal types, explicit non-inference for `unit-cell`, invalid aliases,
identity mismatch, column-vector declaration, invalid units, and invalid JSON.
This is implementation-conformance evidence only.
