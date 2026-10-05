# PhysKit lecture notes and computational laboratories

This directory owns durable lecture-note and manuscript material. The Python
package, lecture notes, notebooks, and computational studies have different
responsibilities and should not duplicate one another.

## Artifact responsibilities

### Maintained Python package

Reusable mathematics, numerical representations, solvers, data models, and
visualization behavior belong under `src/python/projectkoios/physkit/`. They
require ordinary
software tests. Do not place reusable implementation inside a manuscript figure
script or notebook merely to keep an example self-contained.

### Lecture notes

Lecture notes own the conceptual narrative:

- definitions and notation;
- derivations and assumptions;
- interpretation and limitations;
- worked examples;
- exercises; and
- references.

A chapter may use Markdown, MyST, or LaTeX according to its existing source
format. Do not reproduce the complete chapter narrative in its notebook.

### Notebooks

Notebooks are computational laboratories. They should contain:

- a link to the prerequisite lecture-note section;
- concise learning objectives and execution directions;
- executable code and visualizations;
- numerical or algebraic checks;
- exercises or questions; and
- brief interpretation of the calculated result.

Notebook execution and attractive output are not substitutes for package tests
or evidence that a lesson has been pedagogically validated.

### Computational studies

A figure derived from a retained computational study belongs with that study,
including its input data and reproduction script. A lecture note may reference
or reuse the published image without taking ownership of the study.

## Current subject notes

- [Finite-difference quantum harmonic oscillator](qm/qho1d/qho1d__finite_difference.md)
- [Complex plane waves](waves/plane/complex-plane-waves/index.md)
- [One-dimensional free-electron Bloch representation](../lecturenotes/solidstate/bloch1d/freeelectron/index.md)
- [One-dimensional cosine-potential Bloch bands](../lecturenotes/solidstate/bloch1d/cosinepotential/index.md)
- [Two-dimensional checkerboard-potential Bloch bands](../lecturenotes/solidstate/bloch2d/checkerboard/index.md)
- [Periodic primitive-cell free particle in two dimensions](../lecturenotes/solidstate/semiconductors/effmass_2d/index.md)
- [Three-dimensional periodic primitive-cell effective-mass reference](../lecturenotes/solidstate/semiconductors/effmass_3d/index.md)
- [Nonorthogonal lattice metric geometry](../lecturenotes/solidstate/lattices/nonorthogonal-metric-geometry/index.md)
- [Kronig--Penney delta-comb dispersion](../lecturenotes/solidstate/electronic-structure/kronig-penney/delta-comb/index.md)
- [Dimensionless periodic Cahn--Hilliard evolution](../lecturenotes/solidstate/phase-field/cahn-hilliard/index.md)
- [Lennard--Jones and Morse radial pair potentials](../lecturenotes/solidstate/interatomic-potentials/radial-pair-potentials/index.md)
- [Molar and mass composition of finite mixtures](materials/mixtures/molar-and-mass-composition/index.md)
- [Site-fraction composition on multiple sublattices](materials/mixtures/sublattice-site-fractions/index.md)
- [Planar deposition geometry from point and disk sources](deposition/planar-source-geometry/index.md)
- [Planar three-dimensional Hertz--Knudsen phase-change flux](thermal-transport/phase-change/hertz-knudsen/index.md)

## Chapter and figure layout

A chapter that owns generated figures should normally use:

```text
docs/lecture-notes/<subject>/<chapter>/
├── index.md                 # or main.tex for a LaTeX chapter
├── README.md                # local build and reproduction commands, when needed
└── figures/
    ├── make_figures.py
    ├── first-figure.png
    └── second-figure.png
```

Use one `make_figures.py` when the figures share a small coherent construction.
Split it into named scripts only when independent inputs or substantially
different reproduction steps make that clearer. Avoid a generic top-level
`scripts/` directory because it hides artifact ownership.

## Figure-script rules

A manuscript-specific figure script should:

1. import reusable behavior from `physkit`;
2. contain only the thin assembly needed for the owned figure;
3. use deterministic inputs and ordering;
4. write only declared outputs beneath its chapter;
5. avoid network access, repository discovery, and hidden working-directory
   assumptions;
6. provide a `main()` entry point and fail visibly on invalid inputs; and
7. document its reproduction command in the chapter README or source.

For example, a chapter heatmap should use:

```python
from projectkoios.physkit.visualize import HeatMapVisualizer
```

rather than reproduce color normalization and plotting behavior locally.

## Tracking and generated files

Track:

- manuscript and lecture-note sources;
- deterministic figure-generation scripts;
- small declared source data required by the figure; and
- rendered images required by the maintained documentation or manuscript.

Do not track:

- caches;
- temporary render directories;
- editor state;
- intermediate LaTeX files;
- exploratory images not referenced by maintained material; or
- duplicated computational-study outputs.

Notebook-only figures should normally remain notebook outputs generated by the
owning cells. If a notebook figure becomes a maintained manuscript figure,
create a thin chapter-owned script that calls the same package API rather than
scraping the notebook output.

## Reproduction

Run chapter scripts from the repository root unless the chapter explicitly
states otherwise:

```bash
python docs/lecture-notes/<subject>/<chapter>/figures/make_figures.py
```

The script and chapter README must identify any required optional dependencies.
Reusable behavior should be tested through the package test suite; figure scripts
need deterministic execution and output checks proportional to their role.
