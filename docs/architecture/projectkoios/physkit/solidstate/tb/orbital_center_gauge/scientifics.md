# Orbital-center gauge scientific contract

The transformation changes the coordinate representation of the same finite operator:

$$
H_{\mathrm{full}}(q)=D(q)^\dagger H_{\mathrm{cell}}(q)D(q).
$$

Because $D(q)$ is unitary,

$$
\operatorname{spec}H_{\mathrm{full}}(q)
=
\operatorname{spec}H_{\mathrm{cell}}(q),
$$

up to numerical eigensolver effects outside this module. Elementwise equality is not
expected unless the center phases are trivial or the matrix commutes with $D(q)$.

For an integer reciprocal-lattice shift $G\in\mathbb Z^3$,

$$
D(q+G)=D(q)S_G,
\qquad
(S_G)_{aa}=\exp(+2\pi iG\cdot\tau_a).
$$

If the ordinary-cell matrix is periodic,

$$
H_{\mathrm{cell}}(q+G)=H_{\mathrm{cell}}(q),
$$

then the full-center representation is covariant:

$$
H_{\mathrm{full}}(q+G)
=S_G^\dagger H_{\mathrm{full}}(q)S_G.
$$

Thus its spectrum is reciprocal-periodic even when its matrix entries require sewing.
The module constructs $D(q)$ only; it does not own a reciprocal-sewing policy or claim
that independently supplied center and basis objects are compatible.

The positive exponent is a fixed convention. Replacing it by a negative exponent would
reverse the declared transformation and is not an interchangeable formatting choice.

- [Index](index.md)
- [Schematics](schematics.md)
- [Numerics](numerics.md)
