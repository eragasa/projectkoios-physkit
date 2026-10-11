# Package `projectkoios.physkit.potentials.periodic1d`

This package owns `CosinePotential1D`, the fractional-coordinate field

$$
V(u)=A\cos(2\pi h u+\phi).
$$

The model provides real-space evaluation, the exact support $m=\pm h$, and the
coefficients $V_{+h}=A e^{i\phi}/2$ and
$V_{-h}=A e^{-i\phi}/2$. It is consumed by the maintained one-dimensional
cosine-potential Bloch model through the sampled-diagonal operator converter.
`PeriodicPotential1DVisualizer` renders any nominal one-dimensional periodic
potential on supplied fractional coordinates without displaying or closing the
figure.

- Code:
  `src/python/projectkoios/physkit/potentials/periodic1d/cosine.py`
- Visualization:
  `src/python/projectkoios/physkit/potentials/periodic1d/visualization.py`
- Tests:
  `tests/projectkoios/physkit/potentials/periodic1d/`
- [Shared periodic-potential contracts](../periodic/index.md)
