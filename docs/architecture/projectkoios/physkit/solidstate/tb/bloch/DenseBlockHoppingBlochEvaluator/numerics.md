# `DenseBlockHoppingBlochEvaluator` numerics

The ordered route evaluates one scalar finite Fourier series per matrix entry. The
batched route forms a sorted displacement union, inserts exact complex zeros for
missing family terms, and evaluates channels ordered by model, row, and column.
Output reshaping preserves model and coordinate order.

The dense batched reduction is numerically equivalent within floating-point
roundoff but is not required to be bitwise identical to ordered scalar reduction.
Empty coordinate grids produce empty correlated result inventories.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
