# Homopolar $sp^3s^*$ hopping scientific contract

The matrix represents a spin-independent orthogonal two-center operator between two
ordered five-orbital blocks. Rows belong to the source center and columns belong to the
target center:

$$
T_{ab}=\langle a,\mathrm{source}|\hat H|b,\mathrm{target}\rangle.
$$

The sign convention is explicit:

$$
\langle s|\hat H|p_i\rangle=d_iV_{sp\sigma},
\qquad
\langle p_i|\hat H|s\rangle=-d_iV_{sp\sigma},
$$

with the analogous convention for $s^*$ and $p_i$. This convention determines the
transpose-under-reversal relation and must not be mixed silently with another table
orientation.

The seven integrals are model inputs. Their presence in a valid record does not establish
that they were fitted, converged, independently verified, transferable, or scientifically
accepted. Explicit zero values mean zero represented contribution; they are not inferred
from a shell label.

The module assumes real homopolar integrals and an orthogonal spatial basis. It makes no
claim for heteropolar tables, overlap matrices $S\neq I$, generalized eigenproblems

$$
Hc=ESc,
$$

spin, spin–orbit coupling, magnetic interactions, many-body terms, or distance scaling.

- [Index](index.md)
- [Schematics](schematics.md)
- [Numerics](numerics.md)
