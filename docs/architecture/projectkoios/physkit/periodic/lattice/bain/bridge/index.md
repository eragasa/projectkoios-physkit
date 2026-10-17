# Module `projectkoios.physkit.periodic.lattice.bain.bridge`

## Current responsibility

The module owns exact coordinate maps between declared crystallographic bases.

- [`FccBainCoordinateBridge3D`](FccBainCoordinateBridge3D/index.md) maps reduced
  coordinates between one FCC primitive basis and the declared Bain basis.
- Integer and fractional maps preserve their distinct numerical representations.
- The package façade exports the bridge by its defining-module identity.
- No basis convention is inferred from arbitrary vectors.

## Ownership boundary

Exact crystallographic coordinate mechanics remain separate from campaign use.

- PhysKit owns the integral matrix, its inverse, mapping direction, and finite-quotient
  bijection.
- Consumers own selection of the bridge, path sampling, structure interpretation, and
  scientific acceptance.
- The bridge does not attach physical units, atomic species, or electronic state.
- It introduces no dependency on KSD.

## Local mapping

Source, tests, and related path documentation remain explicit.

- Code: `src/python/projectkoios/physkit/periodic/lattice/bain/bridge.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/bain/bridge/`
- [Bain deformation paths](../path/index.md)
- [Parent package](../index.md)
