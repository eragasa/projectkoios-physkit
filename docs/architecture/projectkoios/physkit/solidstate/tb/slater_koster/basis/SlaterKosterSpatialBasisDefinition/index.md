# Class `SlaterKosterSpatialBasisDefinition`

## Current contract

Retains exact ordered center and coordinate tuples for one finite spatial representation.

The public boundary accepts a display identifier, nonempty unique center tuple, and nonempty unique coordinate tuple. Every coordinate must reference a declared center. Equality and order are exact; dimensions and labels do not substitute for basis compatibility.

The central represented relation is

$$
\mathcal B=(b_0,b_1,\ldots,b_{N_b-1}),\qquad \boldsymbol{\tau}_{a,:}=\tau_{c_a}.
$$

## Ownership boundary

The definition identifies a represented coordinate frame. It does not establish continuum completeness, orthogonality, or material suitability. It does not admit a specific model specialization, construct operators, or infer missing coordinates.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/basis.py::SlaterKosterSpatialBasisDefinition`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/basis/`
- [Owning module](../index.md)
