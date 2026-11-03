# Block-hopping scientific contract

A block term represents an ordinary-lattice coefficient
$T(\mathbf R)\in\mathbb C^{N\times N}$ in one explicitly identified ordered
primitive-cell basis. A model is the finite represented inventory
$\{T(\mathbf R)\}$ together with lattice dimension, energy unit, energy reference,
and basis identity.

Dense and sparse identify internal coefficient storage only; they do not identify
different physical models. Spatial locality is represented by the displacement
inventory, while sparse coefficient storage describes sparsity inside the retained
basis block.

The records do not infer $T(-\mathbf R)=T(\mathbf R)^\dagger$, cyclic equivalence,
finite-DFT representatives, material provenance, fitting status, or scientific
acceptance.

## Navigation

- [Module index](index.md)
- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
