# Module `projectkoios.physkit.periodic.lattice.embedding.serialization`

## Current responsibility

This module owns the strict deterministic schema-version-1 JSON text boundary
for `EmbeddedDirectLattice2D`. It serializes only the defining ambient `3×2`
basis, column convention, and unit. Derived metric, inverse metric, dual basis,
frame, and projector values are reconstructed from the validated basis and are
never trusted from wire input.

## Wire contract

The top-level object has exactly:

- `schema_version`, exact JSON integer `1`;
- `representation`, exact string `embedded-direct-lattice-2d`; and
- `ambient_basis`, an exact object containing `a1`, `a2`, `vector_axis`, and
  `unit`.

Primitive vectors contain exactly three finite JSON numbers and
`vector_axis` is exactly `columns`. Physical and unitless values use explicit
discriminated unit objects.

Encoding uses sorted keys, compact separators, one trailing newline, and
`allow_nan=False`. Decoding rejects oversized input, malformed JSON, duplicate
members, unknown or missing fields, nonstandard constants, numeric strings,
booleans as numbers, invalid units, and invalid embedded geometry.

## Public contracts

- [`EmbeddedDirectLattice2DJsonCodec`](EmbeddedDirectLattice2DJsonCodec/index.md)
- [`EmbeddedDirectLattice2DSerializationError`](EmbeddedDirectLattice2DSerializationError/index.md),
  the domain-specific strict wire-boundary error.

## Ownership

The codec is PhysKit-owned and domain-specific. It has no dependency on
ksdft2effmass, no generic dataclass introspection, and no package-global codec
registry. Existing `UnitCellJsonCodec` bytes and behavior are unchanged.

## Local mapping

- [Adaptation provenance](../../../../../../../provenance/embedded-lattice-json-serialization.md)
- Code: `src/python/projectkoios/physkit/periodic/lattice/embedding/serialization.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/embedding/serialization/`
