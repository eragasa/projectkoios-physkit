# Package `projectkoios.physkit.potentials.periodic3d`

This package owns `CosineProductPotential3D`, the fractional-coordinate field

$$
V(u)=A\prod_{j=1}^{3}\cos(2\pi h_j u_j+\phi_j).
$$

The model has exact eight-mode Fourier support and can be converted to sampled
or reciprocal-mode multiplication operators. No three-dimensional Hamiltonian
or solver is implied by this field-level model.

- Code:
  `src/python/projectkoios/physkit/potentials/periodic3d/cosine_product.py`
- Tests:
  `tests/projectkoios/physkit/potentials/periodic3d/cosine_product/`
- [Shared periodic-potential contracts](../periodic/index.md)
