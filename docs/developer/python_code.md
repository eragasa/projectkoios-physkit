# Python code conventions

## Package ownership

Reusable physics, numerical methods, quantities, grids, operators, solvers, and
visualizers belong under `src/python/projectkoios/physkit/`. Keep
research-campaign orchestration, retained campaign records, manuscript workflows,
and generated evidence outside the package.

An optional legacy tree at `src/physkit/` may remain for workspace artifacts or
transitional source placement. Its presence alone is not a layout violation.
Every `.py` or `.pyi` file placed there must have a corresponding file at the
same relative path under `src/python/projectkoios/physkit/`. The canonical tree
remains the only packaging root and the owner of maintained code; path mirroring
does not create a supported top-level `physkit` import, compatibility facade, or
second implementation owner. Cache files and other non-source artifacts are not
part of the mirror check. Notebook paths are excluded: computational laboratories
follow their pedagogical organization and do not mirror package or test paths.

Organize dimensional model families explicitly. The clean particle-in-a-box
implementation uses:

```text
projectkoios/physkit/qm/piab1d/
├── __init__.py
├── base.py
├── tise_analytical/
├── tdse_analytical/
├── tise_fd/
└── tdse_fd/
```

Each solution package owns its solver and immutable results type. Analytical
TDSE uses exact phase evolution in a finite continuum eigenbasis. FD TDSE uses
exact spectral propagation of the complete finite-difference Hamiltonian; `fd`
identifies the spatial representation, not a finite-difference time integrator.

Do not add empty 2D or 3D placeholders. Add a dimensional package only when its
model and represented mathematics are implemented.

Generic one-dimensional finite-difference operators remain in
`projectkoios.physkit.operators.operators_1d.fd`. Schrödinger kinetic energy,
sampled quantum potentials, and Hamiltonian composition belong in
`projectkoios.physkit.operators.qm.qm_1d`.

Superseded quantum and numerical implementations are retained only under
`projectkoios.physkit.legacy.qm` and
`projectkoios.physkit.legacy.numerics`. Both namespaces are deprecated, warn
when imported, and must not be dependencies of maintained package code.

## Public imports

The defining module is the default import location. Do not recursively re-export
symbols through parent package `__init__.py` files merely for convenience.

A local `__init__.py` may aggregate a coherent public inheritance family when
that makes ownership clearer. Keep such aggregation at the narrowest owning
package and do not bubble it into broader namespaces without a concrete need.

## Runtime ownership

Base objects use an explicit `Base` prefix. The quantum-model hierarchy begins
with `BaseQuantumModel`; dimensional model bases inherit from it directly, for
example `BaseQuantumModel1D(BaseQuantumModel)`. A concrete one-dimensional model
such as `Piab1D` inherits `BaseQuantumModel1D` directly.

Use nominal base classes and explicit `isinstance` checks where runtime ownership
matters. Do not introduce structural protocols merely to avoid choosing an owner.
Keep immutable data and results in frozen, slotted dataclasses when practical.

## Unit-aware numerical calculations

`UnitSystem` selects a coherent numerical scale. Model and operator boundaries
retain units, but arrays passed into `projectkoios.physkit.numerics` are stripped
to plain numerical magnitudes. Represented operator or model layers wrap numerical
results with units again. Numerical matrices remain in the selected scale; they are
not converted wholesale to canonical SI.

For example, `UnitSystem.METAL` uses ångström, dalton, picosecond, electron-volt,
and electron-volt–picosecond. Unit conversion may supply a multiplicative
coefficient, but the represented Laplacian, Hamiltonian, and eigenvalues remain
well-scaled in the selected units.

Use `ScalarQuantity`, `VectorQuantity`, `ComplexVectorQuantity`,
`MatrixQuantity`, and `SparseMatrixQuantity` at package boundaries. Maintain
canonical immutable CSR storage for sparse matrices.

## Results objects

A computation returns one immutable object inheriting `ResultsObject` when its
outputs are correlated. Do not return unrelated arrays that require the caller
to reconstruct provenance or dimensional relationships.

A concrete results object retains the exact inputs and representations needed to
interpret it. `RealSymmetricEigendecomposition` correlates one unit-aware represented
operator with its eigenvalues and eigenvectors while the unit-free eigensolver result
remains an internal numerical-layer output. `Piab1dTiseFdSolution` composes that
represented decomposition with its model and interval and exposes the Hamiltonian,
energies, and eigenvectors as correlated accessors. Use `fd` for finite-difference
module and type names.

Results objects do not absorb campaign metadata, serialization policy, workflow
state, or manuscript evidence.

## External-source migrations

Treat external repositories as read-only donors unless changes are explicitly
authorized. Record the exact donor revision for migrated numerical behavior. Do
not add a runtime dependency on a local checkout.

Adapt names and ownership to Project Koios PhysKit while preserving represented mathematics,
units, numerical behavior, and relevant verification. Do not copy campaign
infrastructure when only reusable model or numerical machinery is needed.
