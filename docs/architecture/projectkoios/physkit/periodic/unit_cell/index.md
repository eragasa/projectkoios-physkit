# Package `projectkoios.physkit.periodic.unit_cell`

## Current responsibility

The package owns immutable calculator-neutral atoms, ordered atomic bases, and
three-dimensional unit cells. It also owns the versioned deterministic JSON
boundary for those representations.

## Public contract

`UnitCell` stores an immutable dimensionless direct-lattice column basis, a
positive physical length parameter, and an ordered nonempty atomic basis.
`ConventionalUnitCell` and `PrimitiveUnitCell` add explicit nominal
representation semantics without changing that geometry.

The narrow package API exports these models, `UnitCellJsonCodec`, and
`UnitCellSerializationError`. It does not infer conventional or primitive
semantics from geometry and does not execute a calculator.

## Navigation

- [Module `projectkoios.physkit.periodic.unit_cell.serialization`](serialization/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/unit_cell/`
- Tests: `tests/projectkoios/physkit/periodic/unit_cell/`

## Evidence

Mapped tests cover nominal ownership, immutable geometry, units, atom ordering,
public exports, deterministic serialization, and invalid inputs. These checks
establish represented software behavior only.
