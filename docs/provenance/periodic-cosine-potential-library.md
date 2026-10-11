# Periodic cosine-potential library provenance

## Extraction baseline

The reusable periodic-potential contracts were introduced from behavior already
owned by this repository at commit
`e506a1d3e6397502dff20368c0d26bbc04973a4c`.

| Source | Git blob | Preserved behavior |
|---|---|---|
| `notebooks/solidstate/bloch1d/cosinepotential/exploratory-band-gap-prototype.ipynb` | `f4b38c1cba94e1ad8765eab60e44f54df98b12a0` | Real-space $V_0\cos(2\pi x/a)$ model, periodic sampling, and reciprocal coupling explanation |
| `src/python/projectkoios/physkit/solidstate/bloch1d/cosinepotential/model.py` | `e4809c10c33c95eb66225bc9ae93969cf0518080` | Unit-aware cosine amplitude and sampled finite-difference potential term |
| `src/python/projectkoios/physkit/solidstate/bloch2d/checkerboard/model.py` | `333ed5f154d88fd5f450808d7c5a0b1d089cb79c` | Checkerboard $V_0\cos(\mathbf b_1\cdot\mathbf r)\cos(\mathbf b_2\cdot\mathbf r)$ and four $V_0/4$ reciprocal couplings |

No runtime dependency on a notebook or donor checkout was introduced.

## Adaptation

The extracted models use dimensionless fractional coordinates rather than
embedding a cell length or direct basis. This separates scalar-field ownership
from the lattice mapping and from Hamiltonian construction. Physical reciprocal
vectors were replaced at the potential boundary by exact integer mode indices.

The one- and two-dimensional formulas retain their represented mathematics. The
three-dimensional cosine product is the direct dimensional extension of the
same documented separable identity and has eight exact Fourier-support modes.

Potential-to-operator converters make the previous implicit representation
steps explicit: sampled values become a diagonal CSR operator, while Fourier
coefficients become the convolution matrix $V_{m_r-m_c}$.

## Evidence boundary

Mapped tests verify the analytical sampled values and Fourier coefficients and
confirm that the existing one- and two-dimensional Bloch consumers retain their
represented results. This is implementation and numerical verification, not
experimental validation of a material potential.
