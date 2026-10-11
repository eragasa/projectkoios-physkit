# Class `EmbeddedDirectLattice2DJsonCodec`

The stateless codec maps exactly between `EmbeddedDirectLattice2D` and strict
schema-version-1 JSON text. It emits deterministic compact sorted JSON with a
trailing newline and rejects non-finite output.

Decoding enforces exact fields at every object level, duplicate-member
rejection before dictionary collapse, standard finite numeric values, exact
representation and vector-axis discriminators, explicit units, a bounded input
size, and exact nominal domain reconstruction.

The wire carries no derived arrays. Callers use `EmbeddedDirectLattice2DAdapter`
after decoding to reconstruct `G`, `G^{-1}`, `B`, `Q`, and `P` from `A`.

- [Module index](../index.md)
- [Embedding numerical contract](../../numerical/index.md)
- Code: `src/python/projectkoios/physkit/periodic/lattice/embedding/serialization.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/embedding/serialization/test__EmbeddedDirectLattice2DJsonCodec__*.py`
