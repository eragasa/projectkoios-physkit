# Orbital-center gauge numerics

The constructor evaluates one binary64 dot product per ordered orbital center and
forms a diagonal complex matrix with the positive exponent
$\exp(+2\pi i q\cdot\tau_a)$. Inputs are nondimensional quantities of shapes `(basis_size, 3)`
and `(3,)`; no dimensional projection or reciprocal-axis reordering is performed.

Focused tests compare the matrix directly with an independently formed NumPy
expression and verify immutable output storage.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
