# Orbital-center gauge numerics

Let the ordered fractional-center matrix be

$$
\boldsymbol{\tau}=
\begin{pmatrix}
\tau_0^{\mathsf T}\\
\tau_1^{\mathsf T}\\
\vdots\\
\tau_{N_b-1}^{\mathsf T}
\end{pmatrix}
\in\mathbb R^{N_b\times3},
$$

and let $q\in\mathbb R^3$ be a fractional reciprocal coordinate. The constructor first
computes the binary64 vector

$$
\theta=\boldsymbol{\tau}q,
\qquad \theta_a=q\cdot\tau_a,
$$

then evaluates

$$
d_a=\exp(+2\pi i\theta_a)
$$

and constructs

$$
D(q)=\operatorname{diag}(d_0,d_1,\ldots,d_{N_b-1})
\in\mathbb C^{N_b\times N_b}.
$$

For finite phases, every diagonal entry has unit modulus in exact arithmetic:

$$
|d_a|=1,
\qquad
D(q)^\dagger D(q)=I_{N_b}.
$$

The implementation uses one NumPy matrix-vector product, one vectorized complex
exponential, and `numpy.diag`. It preserves row order exactly and performs no modulo
reduction of $q$ or $\tau_a$. Inputs use `UnitSystem.NONDIMENSIONAL.length_unit`; this
unit records nondimensional coordinate ownership and is retained on the output
`ComplexMatrixQuantity`.

Focused tests compare against an independently formed NumPy expression, verify the
positive sign, and verify immutable output storage.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
