# Project Koios PhysKit

`projectkoios-physkit` owns reusable computational-physics models, mathematical
representations, numerical methods, unit-aware quantities, and pedagogical
examples under the implicit namespace package `projectkoios.physkit`.

The component supports mathematically explicit teaching and bounded scientific
software development. Current implemented areas include units, discretization,
operators, quantum models, periodic lattices, solid-state model components,
thermodynamics, and visualization. Availability of an implementation or a
passing test does not by itself establish numerical verification, scientific
validation, pedagogical validation, uncertainty quantification, or human
acceptance.

## Ownership boundary

This repository owns reusable physics and numerical behavior that another
scientific project can use without adopting a material-specific research claim.
It does not own:

- research-campaign orchestration or retained campaign results;
- application workflow state, dispatch, retries, or execution authority;
- material-specific scientific acceptance criteria;
- electronic-structure or Wannier production execution; or
- optimization engines owned by `projectkoios-optimization`.

`ksdft2effmass` remains the owner of its semiconductor-specific models,
campaigns, provenance, tolerances, and scientific interpretations. Generic
finite-periodic domain and represented-operator mechanics adapted from that
repository are documented in
[`docs/provenance/finite-periodic-lattice-source-mapping.md`](docs/provenance/finite-periodic-lattice-source-mapping.md).
Quantity-aware matrix norms delegate their unit-free arithmetic to Optimization; their
KSD prior-art boundary is documented in
[`docs/provenance/ksdft2effmass-matrix-norms.md`](docs/provenance/ksdft2effmass-matrix-norms.md).

## Python package

Project Koios repositories share the PEP 420 implicit namespace `projectkoios`.
This distribution therefore provides `projectkoios.physkit` and does not provide
a top-level `projectkoios/__init__.py`.

```python
from projectkoios.physkit.periodic.lattice.finite_domain import (
    FinitePeriodicDomain,
    LatticeDimension,
)

finite_domain = FinitePeriodicDomain(LatticeDimension.ONE, (8,))
```

Python 3.14 or newer is required.

## Documentation

Current architecture documentation begins at
[`docs/architecture/projectkoios/physkit/index.md`](docs/architecture/projectkoios/physkit/index.md).
Developer conventions are under [`docs/developer/`](docs/developer/index.md), and
lecture-note ownership is described in
[`docs/lecture-notes/README.md`](docs/lecture-notes/README.md).

## Development

```bash
python3.14 -m venv .venv
.venv/bin/python -m pip install -e '.[dev]'
.venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check src/python tests
.venv/bin/python -m build --wheel
```

## License

Project Koios PhysKit is licensed under the
[Apache License 2.0](LICENSE). See [`NOTICE`](NOTICE) and the maintained
provenance records for adapted-source identities and limitations.
