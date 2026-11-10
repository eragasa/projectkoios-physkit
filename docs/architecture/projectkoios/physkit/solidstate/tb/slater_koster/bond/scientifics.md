# Slater–Koster bond scientific contract

The ordinary-cell displacement and embedded-center displacement answer different
questions. The Bloch coefficient label is $R\in\mathbb Z^3$, while the physical bond
vector may satisfy a composing geometry relation such as

$$
\Delta r=A\bigl(R+\tau_t-\tau_s\bigr),
$$

for primitive-basis matrix $A$ and fractional centers $\tau_s,\tau_t$. This module does
not derive or verify that relation; the geometry owner supplies $\Delta r$ explicitly.

Consequently, the ordinary-cell phase is based on

$$
e^{+2\pi i q\cdot R},
$$

not on $R+\tau_t-\tau_s$. Fractional center offsets belong to an explicit gauge
transformation.

The parameter-group identifier in a parameterized bond states only

$$
g(b)=\text{declared coefficient-sharing group}.
$$

It does not mean first neighbor, second neighbor, chemical bond type, or material shell
unless a client model separately defines that interpretation.

The current geometry is three-dimensional, spin independent, and suitable for the
homopolar $sp^3s^*$ specialization. Heteropolar orientation rules, distance-dependent
scaling, lattice generation, and scientific admission remain outside the module.

- [Index](index.md)
- [Schematics](schematics.md)
- [Numerics](numerics.md)
