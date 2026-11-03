# Block-hopping Bloch scientific contract

Bloch evaluation transforms an ordinary-lattice hopping representation into a
reciprocal-coordinate matrix field in the same retained basis, gauge, energy unit,
and energy reference. Fractional coordinates are paired with integer primitive-cell
displacements through the explicit $2\pi$ phase convention.

Dense and sparse evaluators preserve storage ownership; they do not define different
Hamiltonians. The returned matrices may be non-Hermitian because the source model is
not required to contain conjugate displacement partners.

Evaluation establishes represented Fourier assembly only. It does not establish
band accuracy, material validity, convergence, or physical adequacy.

## Navigation

- [Module index](index.md)
- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
