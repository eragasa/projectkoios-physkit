# Slater–Koster scientific contract

The package represents a finite spin-independent orthogonal spatial basis and the
angular dependence of real homopolar two-center integrals. The basis definition is exact,
ordered, and center-resolved. The `s*` coordinate is a second radial channel of the
`s` component rather than a new universal angular-orbital kind.

Parameter-group identifiers state only explicit coefficient sharing. They do not
infer chemical species, material, neighbor shell, bond length law, or scientific
admissibility.

Ordinary-cell assembly accepts an explicit directed-representative inventory. The
integer displacement controls coefficient placement, while the physical Cartesian
vector controls only the bond direction. Assembly generates the Hermitian reverse and
rejects an explicitly duplicated reverse; it does not infer lattice geometry or
neighbor relationships.

The package constructs matrix representations only. It does not select a physical
model, add onsite terms, evaluate Bloch fields, fit parameters, align target spaces,
define residuals, qualify evidence, or judge agreement with first-principles
calculations.

- [Index](index.md)
- [Schematics](schematics.md)
- [Numerics](numerics.md)
