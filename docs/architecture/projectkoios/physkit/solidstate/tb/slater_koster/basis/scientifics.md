# Slater–Koster basis scientific contract

The module defines a finite represented spatial state space, not atomic wavefunctions or
a material model. For $N_c$ admitted centers, the represented vector space is

$$
\mathcal H_{\mathrm{spatial}}
\cong \bigoplus_{c\in\mathcal C}\mathbb C^5,
$$

with one exact ordered coordinate frame per center. The labels identify coordinates in
that frame; they do not establish completeness, orthogonality of underlying continuum
orbitals, chemical identity, or transferability of a parameter set.

The distinction

$$
(s,0)\neq(s,1)
$$

is radial-channel identity. The second coordinate is conventionally denoted $s^*$ in
the admitted model, but the module does not interpret the integer as a principal quantum
number.

Basis compatibility is exact and order-sensitive:

$$
\mathcal B_1\equiv\mathcal B_2
\quad\Longleftrightarrow\quad
(b^{(1)}_0,\ldots,b^{(1)}_{N_b-1})
=(b^{(2)}_0,\ldots,b^{(2)}_{N_b-1})
$$

under the same center definitions. Equal dimensions, matching strings, or equal orbital
counts are insufficient.

Spin, spin–orbit coupling, overlap matrices, generalized eigenproblems, orbital
normalization, and material admission remain outside this module.

- [Index](index.md)
- [Schematics](schematics.md)
- [Numerics](numerics.md)
