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

A one-dimensional Bloch state factors as

$$
\psi_k(x)=e^{ikx}u_k(x),\qquad u_k(x+a)=u_k(x).
$$

Complex amplitude and probability density have different units and are not
placed on one ordinate. Integer reciprocal-cell Fourier labels describe the
periodic factor $u_k$, so the visualizer removes $e^{ikx}$ before applying the
discrete transform. Applying a periodic FFT directly to a nonzero-$k$ twisted
state would introduce carrier leakage unrelated to reciprocal-mode mixing.
Individual eigenvectors inside a degenerate numerical eigenspace are
basis-dependent; visualizing one does not establish a unique physical state.

Drawing a curve establishes no scientific validation. Any claim about
convergence, a material band gap, or agreement with an analytical reference
must be supported by the owning model and its separate numerical or scientific
evidence.

## Navigation

- [Package index](index.md)
- [Schematics](schematics.md)
- [Numeric contract](numeric.md)
