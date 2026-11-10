# `HomopolarSp3sStarLinearParameterizationAssembler` numerics

For each pair $(g,\alpha)$, the action embeds $B_\alpha(d_b)$ only for bonds explicitly
naming $g$, generates transpose reverses, accumulates `complex128` coefficients, and
sorts displacements lexicographically. Nonfinite accumulation fails closed.

- [Index](index.md)
