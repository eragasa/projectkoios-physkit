# `BlochCell2D` scientific scope

## Represented scientific statement

`BlochCell2D` represents two ordered direct translations of a displayed
periodic cell in an intrinsic two-dimensional Cartesian plane. It retains the
basis convention required to interpret reciprocal vectors and plane-wave Bloch
wave vectors.

The record deliberately does not assert that the displayed cell is primitive.
For the checkerboard product potential, the supplied square conventional cell
has an additional half-diagonal translation symmetry. Calling that displayed
square primitive would overstate its crystallographic meaning.

## Supported conventions

### Metal convention

`UnitSystem.METAL` represents direct lengths in ångströms, reciprocal vectors in
inverse ångströms, mass in daltons, and energy in electron-volts. The
checkerboard model converts the physical bare-electron mass into that coherent
numerical convention.

### Nondimensional convention

`UnitSystem.NONDIMENSIONAL` uses first-class `Unitless` quantities throughout.
For the checkerboard Hamiltonian it declares the model convention

$$
\hbar=m_e=1.
$$

This is a nondimensional model definition, not a claim that physical constants
are dimensionless.

### Excluded SI convention

`UnitSystem.SI` is rejected at the `BlochCell2D` boundary. Supporting SI would
be a separate public-contract decision with its own selected-scale numerical
evidence; implicit acceptance is intentionally avoided.

## Assumptions

- Both basis vectors are finite and numerically independent in binary64.
- The basis is constant over the displayed cell.
- Coordinates are intrinsic two-dimensional Cartesian coordinates.
- A checkerboard consumer additionally requires orthogonal equal-length basis
  vectors.
- The reciprocal-mode cutoff is finite and separately chosen by the caller.
- Sparse storage is a numerical representation choice and does not alter the
  represented plane-wave Hamiltonian.

## Explicit separations

The cell does not own:

- a claim of primitive-cell minimality;
- a real-material crystal structure;
- an embedded `3×2` ambient orientation;
- a potential, particle filling, or temperature;
- basis-cutoff convergence or uncertainty estimates; or
- a mapping between the nondimensional convention and a physical material.

Ambient three-dimensional embeddings belong to the separate periodic-lattice
embedding capability. Semiconductor effective-mass primitive-cell records also
remain separate because they make a different ownership and scientific claim.

## Evidence classification

### Implementation conformance

Constructor, unit, geometry, and checkerboard tests establish the declared API
for tested cases.

### Numerical verification

Analytic direct/reciprocal identities, normalized shape checks, canonical sparse
storage, sparse-versus-dense spectral agreement, Hermiticity, free-particle
agreement, and the unitless kinetic diagonal provide bounded numerical evidence.

### Scientific validation

Not established. No comparison with measured material bands or an independent
electronic-structure calculation is included.

### Uncertainty quantification

Not established. The implementation supplies no uncertainty for cell
parameters, cutoff error, floating-point conditioning, or model form.

### Human acceptance

Not established by tests or documentation. A consuming project retains
responsibility for declared-use review and acceptance.

## Navigation

- [Module index](index.md)
- [Schematics](schematics.md)
- [Numeric page](numeric.md)
