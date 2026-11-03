# `ScalarHoppingBlochEvaluator` scientific contract

For fractional reciprocal coordinate $\mathbf q$, the evaluator represents

$$
h(\mathbf q)=\sum_{\mathbf R}t(\mathbf R)
\exp(+2\pi i\,\mathbf q\cdot\mathbf R).
$$

The scalar output may be complex because the model does not require conjugate
hopping partners. Evaluation preserves the source model's energy convention and does
not infer Hermiticity, choose a path, truncate support, or establish model adequacy.

- [Index](index.md)
- [Schematics](schematics.md)
- [Numerics](numerics.md)
