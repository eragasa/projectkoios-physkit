# Class `projectkoios.physkit.periodic.basis.ReciprocalModeBasis3D`

## Current contract

The class represents one ordered finite inventory of 3D reciprocal-lattice modes.

- `reciprocal_lattice` is retained by exact canonical identity.
- `modes` has logical shape `(number_of_modes, 3)` and contains unique signed-64-bit integers.
- `vectors` contains the corresponding Cartesian vectors in the same order.
- Both arrays are independently owned, immutable NumPy storage.
- `size` is the represented finite basis dimension.

## Symmetric factory

`from_mode_limit(...)` constructs a deterministic cubic inventory.

- A nonnegative built-in integer cutoff `M` represents `[-M, M]^3`.
- The basis dimension is `(2 * M + 1)^3`.
- Iteration uses first-index outer, second-index middle, and third-index inner order.
- Cutoff zero retains exactly `(0, 0, 0)`.
- The factory does not select a physically adequate cutoff.

## Shifted vectors

`shifted_wavevectors(...)` restores one Cartesian carrier displacement.

- Input is one finite real vector with shape `(3,)`.
- The vector is added to every represented reciprocal vector.
- The operation does not reduce coordinates into a Brillouin zone.
- No route, neighbor, gauge, or sewing interpretation is assigned.

## Ownership boundary

The class owns physical representation while preserving numerical and campaign owners.

- `ReciprocalLattice3D.vector(...)` owns lattice-index conversion.
- NumPy owns dense array arithmetic used by the existing lattice implementation.
- Callers own cutoff selection, convergence, basis completeness, and scientific use.
- KSD source did not define this owner contract and is not a runtime dependency.

## Local mapping

Source, tests, and parent documentation remain adjacent.

- Code: `src/python/projectkoios/physkit/periodic/basis.py::ReciprocalModeBasis3D`
- Tests: `tests/projectkoios/physkit/periodic/basis/test__ReciprocalModeBasis3D__contract.py`
- Façade tests: `tests/projectkoios/physkit/periodic/test__periodic_reciprocal_basis_public_api.py`
- [Parent module](../index.md)
