# Generalized oblique-coordinate prototype review

## Source identity

The reviewed artifact is
`notebooks/qm/piab2d/tise/numerical/generalized-oblique-coordinate-prototype.ipynb`,
relocated byte-for-byte from `notebooks/scratch/generalized_piab2d.ipynb`.

- Git blob: `2cf7a63a6722724537758c0e12cad60a3ca24445`
- SHA-256: `22c333ddea3d0eacf944124660f474bcd9e4428391d7d68b9c28661038e58322`
- Preservation commit: `fbf7e9968671cffc18760e9a4aadaf165843d4cc`
- Review baseline: `bf54f7d`

The artifact remains present pending final disposition.

## Inventory and execution evidence

The notebook contains 32 cells: 24 code cells, of which 21 are nonempty, and
eight Markdown cells. It stores 20 outputs, including ten PNG displays and one
saved `SyntaxError` from the prose entered as code in cell 30.

A clean sequential execution reaches cell 17 and stops with
`NameError: name 'B' is not defined`. The ten preceding nonempty code cells run,
but their eigenvalues use unit grid spacing rather than the declared fractional
interval spacing. Therefore the saved reciprocal-space results do not reproduce
from the stored execution sequence, and the preceding bounded-spectrum values
do not represent the stated unit cell.

## Route-by-route assessment

### Bounded curvilinear finite differences

The first route maps a fractional square through an oblique direct basis and
forms the correct structural combination

$$
g^{11}D_{11}+g^{22}D_{22}+2g^{12}D_{12}.
$$

It is a homogeneous-Dirichlet interior-stencil route, not a periodic route. The
code sets `h = 1` for 50 represented points while plotting them with
`linspace(0, 1, 50)`. It consequently neither uses the interior spacing
$1/(N+1)$ nor consistently represents boundary samples. The Markdown expansion
initially omits the factor two on the mixed term even though the later discrete
formula and implementation include it. Eigenvectors are not normalized with
the cell Jacobian.

This bounded oblique-parallelogram problem is scientifically distinct from the
periodic primitive-cell free particle. Its maintained destination is now the
separate `ObliquePiab2D` model and sparse finite-difference solver; the existing
rectangular `Piab2D` contract was not silently generalized.

### Reciprocal plane waves

The second route computes $B=2\pi A^{-\mathsf T}$ and the exact free-particle
energies $\hbar^2|B\mathbf m|^2/(2m)`. That conceptual layer is maintained by
`solidstate.semiconductors.effmass_2d.PeriodicFreeParticle2D` and
`notebooks/solidstate/semiconductors/effmass_2d/primitive-cell-methodological-review.ipynb`.
The nondimensional isotropic counterpart is also demonstrated by
`notebooks/qm/piab2d/tise/numerical/oblique-periodic-free-particle.ipynb`.

The prototype defines the missing reciprocal basis only in later duplicate
cells. Several subsequent cells mix arrays of 20 eigenvalues with 121 reciprocal
points, reassign a different direct basis, or reconstruct with coefficients
from an earlier basis. Its labeled “Fermi surface” uses a median single-particle
energy and a broad tolerance without a particle number, dimensional
normalization, or chemical-potential solve. It is at most an unverified
isoenergy-point sketch and was not migrated as a Fermi-surface claim.

## Maintained destinations

| Prototype layer | Destination or disposition |
| --- | --- |
| Direct metric and reciprocal duality | `projectkoios.physkit.periodic.lattice` |
| Periodic inverse-metric finite differences | `projectkoios.physkit.periodic.finite_difference` |
| Nondimensional periodic laboratory | `oblique-periodic-free-particle.ipynb` |
| Unit-aware reciprocal free-particle modes | `PeriodicFreeParticle2D` and its methodological-review notebook |
| Physical-parallelogram periodic density | Both maintained laboratories above |
| Bounded oblique homogeneous-Dirichlet operator | `ObliquePiab2D` and `qm.piab2d.tise.finite_difference` |
| Bounded physical-space state plots and refinement | `oblique-dirichlet-free-particle.ipynb` |
| Historical invalid or cross-route plots | Retain only in this recoverable prototype; do not migrate as evidence |

## Current disposition

**Retain pending final disposition.** The periodic, reciprocal, and bounded
homogeneous-Dirichlet layers now have separate maintained package and laboratory
destinations. The exploratory notebook has not been removed; deletion requires
an explicit authorization after the maintained replacement is validated.
Repairing the mixed notebook in place remains inappropriate because its two
routes use different state spaces and boundary conditions.
