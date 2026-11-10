# Slater–Koster basis numerics

Let the ordered center inventory be

$$
\mathcal C=(c_0,c_1,\ldots,c_{N_c-1}),
\qquad \tau_c\in\mathbb R^3,
$$

where every $\tau_c$ is a nondimensional primitive fractional coordinate. A generic
ordered spatial coordinate is

$$
b_a=(c_a,\lambda_a,r_a),
$$

with center identifier $c_a$, real-cubic component $\lambda_a$, and nonnegative
built-in integer radial channel $r_a$.

For the admitted homopolar $sp^3s^*$ specialization, each center contributes

$$
\mathcal B_c=
\bigl(
(c,s,0),
(c,p_x,0),
(c,p_y,0),
(c,p_z,0),
(c,s,1)
\bigr).
$$

The complete basis is the exact concatenation

$$
\mathcal B=
\mathcal B_{c_0}\mathbin\Vert
\mathcal B_{c_1}\mathbin\Vert\cdots\mathbin\Vert
\mathcal B_{c_{N_c-1}},
\qquad N_b=5N_c.
$$

No permutation is admitted. The first orbital offset of center $c_j$ is therefore

$$
o(c_j)=5j.
$$

The ordered orbital-center matrix
$\boldsymbol{\tau}\in\mathbb R^{N_b\times3}$ is

$$
\boldsymbol{\tau}_{a,:}=\tau_{c_a}.
$$

Validation is exact and nominal: tuples remain tuples, identifiers are built-in strings,
radial channels are built-in integers, and fractional coordinates use
`UnitSystem.NONDIMENSIONAL.length_unit`. Array storage is inherited from immutable
PhysKit quantities.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
