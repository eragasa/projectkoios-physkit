# Solid-state visualization scientific contract

A band structure evaluated along a reciprocal-space path can be rendered
against any ordered scalar parameter $s$:

$$
E_n(s_i),\qquad i=1,\ldots,N.
$$

The underlying wave vectors may belong to one, two, or three dimensions. The
visualizer receives only $s_i$ and the sampled energy columns; this makes the
rendering dimension-independent without treating different reciprocal-space
geometries as equivalent.

The caller owns the map from physical wave vectors to $s$, the order and
continuity of path segments, high-symmetry locations and labels, Brillouin-zone
conventions, and the meaning of any reference bands. The visualizer does not
infer degeneracies, gaps, band connectivity, occupations, density of states, or
material identity.

Drawing a curve establishes no scientific validation. Any claim about
convergence, a material band gap, or agreement with an analytical reference
must be supported by the owning model and its separate numerical or scientific
evidence.

## Navigation

- [Package index](index.md)
- [Schematics](schematics.md)
- [Numeric contract](numeric.md)
