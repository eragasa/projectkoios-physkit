# Package `projectkoios.physkit.solidstate.electronic`

## Current responsibility

The package owns reduced electronic models for solid-state systems.

- [`kronig_penney`](kronig_penney/index.md) owns one-dimensional periodic piecewise and
  singular-potential model families.
- The plane-wave 1D family composes a periodic potential, canonical reciprocal-mode
  basis, Bloch Hamiltonian, complete dense eigensolve, and represented band result.
- The plane-wave family consumes
  `projectkoios.physkit.periodic.basis.ReciprocalModeBasis1D`; it does not define a
  second reciprocal-lattice identity.
- [`frame_transport`](frame_transport/index.md) binds Optimization's unit-free
  closed-loop polar transport to exact route or one-dimensional periodic-mesh identity,
  ordered ambient basis, retained frame, and endpoint-sewing semantics.

## Numerical scope

The maintained plane-wave solver is a direct pedagogical dense implementation.

- Reciprocal vectors and their declared order come from the supplied basis.
- Kinetic and potential terms are assembled explicitly at each wave vector.
- Complete dense eigenpairs are ordered by NumPy's Hermitian eigensolver.
- Backend-neutral production spectrum operations remain owned by Optimization and are
  not implied by this older pedagogical family.
- Closed-loop frame transport delegates to Optimization's exact SciPy capability and
  adds no numerical threshold, retry, fallback, or scientific qualification policy.

## Evidence limits

Focused tests establish composition and represented numerical behavior only.

- Tests bind Hamiltonian assembly and solver correlation to the canonical reciprocal
  basis.
- Tests do not establish basis-cutoff convergence, material fidelity, or scientific
  validation.

## Local mapping

Maintained code and tests mirror the electronic package.

- Code: `src/python/projectkoios/physkit/solidstate/electronic/`
- Tests: `tests/projectkoios/physkit/solidstate/electronic/`
- [Reciprocal-basis architecture](../../periodic/basis/index.md)
