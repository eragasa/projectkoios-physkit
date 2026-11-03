# PhysKit Capability Baseline

## Document status

**Status:** Proposed for human review.

**Proposed path:** `docs/harness/physkit.harness.01-capability-baseline.md`

This document has no accepted authority merely because it has been proposed or
written to a repository path. The path, scope, observations, and future use
remain subject to human review.

If accepted, its authority is limited to:

> An accepted, revision- and worktree-qualified inspection snapshot that may be
> used as planning evidence.

It must not become the active authority for:

- current project status;
- capability maturity;
- lifecycle state;
- canonical implementations;
- active Pi runtime state;
- current task status.

Later accepted capability records, lifecycle artifacts, task records, or other
authoritative documents may supersede the observations in this snapshot. The
historical snapshot should not then be edited to make it appear current. A later
artifact may cite, refine, or supersede it while preserving what was observed at
the recorded revision and working-tree state.

No implementation is selected as canonical by this document.

No capability maturity is assigned by this document. Every capability recorded
below has:

> **Classification: pending**

A future capability lifecycle policy may define maturity states and promotion
rules. This baseline does not design, anticipate, or apply those states.

## 1. Purpose

The purpose of this document is to record a reviewable baseline of what the
stated inspection found in the PhysKit repository.

It records:

- repository locations;
- source-declared objects and functions;
- associated tests;
- associated notebooks;
- associated documentation;
- observed syntax and import conditions;
- apparent competing or overlapping representations;
- observed evidence;
- repository claims that conflict with observed contents;
- unresolved human decisions.

This is an inspection and discrepancy snapshot. It is not:

- a capability lifecycle;
- a continuously updated status document;
- a canonical architecture;
- a public API specification;
- a declaration of support;
- an active task record;
- active Pi runtime state;
- a numerical-verification report;
- a physical-validation report;
- a pedagogical-validation report;
- an uncertainty-quantification report.

The presence of a class or function definition does not establish that it:

- parses;
- imports;
- executes;
- satisfies its intended contract;
- is part of the public API;
- is numerically verified;
- is physically validated;
- is pedagogically accepted.

Non-empty source definitions are therefore described as source declarations,
not as proof of a functional or supported capability.

## 2. Observation and decision labels

The following labels distinguish the status of statements in this document.

### Directly observed fact

`FACT` identifies a result obtained directly from repository files, Git state,
bounded static parsing, import probing, or file inventory.

A fact about a definition, test, import, or saved notebook output is not
automatically evidence that the underlying software contract, mathematics,
physics, or pedagogy is correct.

### Inference requiring confirmation

`INFERENCE` identifies an interpretation suggested by repository structure,
naming, or relationships but not established by an accepted contract.

### Proposed classification

`PROPOSED CLASSIFICATION` would identify a suggested lifecycle or governance
classification. No such classifications are made in this baseline. All
capability classifications remain pending.

### Protected human decision

`HUMAN DECISION` identifies a decision that this document does not make. These
include choices concerning physical models, mathematical formulations,
numerical representations, public APIs, pedagogical purpose, canonical
implementations, evidence acceptance, promotion, replacement, and deprecation.

## 3. Scope

The stated inspection covered:

- `README.md`;
- `pyproject.toml`;
- the current Git branch, revision, and working-tree state;
- Python source under `src/physkit/`;
- tests under `tests/`;
- notebooks under `notebooks/` and notebooks found under `src/`;
- tracked documentation under `docs/`;
- representative examples under `examples/`;
- representative modules and tests from units, grids, discretization, periodic
  systems, quantum mechanics, solid state, thermodynamics, deposition, and
  visualization.

Capability groups are used to keep material relationships and conflicts
visible. They are human-readable groupings, not accepted package boundaries,
public API boundaries, or maturity classifications.

## 4. Exclusions

The inspection did not:

- modify the inspected production source, tests, notebooks, examples, or
  documentation;
- install pytest or another tool;
- execute the test suite;
- execute notebooks;
- validate notebook results against independent references;
- establish numerical convergence;
- establish physical adequacy for an intended use;
- evaluate learning outcomes with students or instructors;
- quantify uncertainty;
- select public APIs;
- select canonical implementations;
- select canonical notebooks;
- classify capabilities by maturity;
- design a capability lifecycle;
- design Pi agents, chains, skills, schemas, or runtime-state conventions;
- determine whether existing magnetism notebook changes are complete,
  intentional, correct, or incorrect;
- create a deterministic complete per-file evidence inventory.

A complete per-file inventory may later be generated as a deterministic
evidence artifact if a human accepts its need, schema, location, generation
command, and retention policy. No such artifact is proposed or created by this
document.

## 5. Inspection provenance

### 5.1 Inspection identity

- **Inspection date:** 2026-08-04.
- **Inspection timezone:** PHT (UTC+08:00).
- **Inspection time:** not recorded; only the local inspection date was retained.
- **Inspected branch:** `main`.
- **Inspected revision:**
  `d307a18db2e2202e9c64863cceea4112db54ce11`.
- **Dirty working-tree paths at baseline collection:**
  - modified:
    `notebooks/magnetism/2d-ising-model-pedagogical.ipynb`;
  - modified:
    `notebooks/magnetism/qm-spin-magnetism-bloch-states.ipynb`;
  - untracked: `package-lock.json`.
- **Staged paths at baseline collection:** none reported.

The untracked proposed baseline file was written only after the initial
inspection, for review convenience. It was not part of the inspected corpus or
the dirty-path list above. This revised proposal does not reinterpret or modify
the pre-existing magnetism notebook changes.

> **Process deviation:** The planning instruction required the proposed content
> to be presented without creating a repository file. Because the complete
> artifact exceeded the available response buffer, it was instead written as
> this single untracked file by the bounded `1a` drafting step. No production
> source, test, notebook, example, or existing documentation file was modified.
> Human acceptance of this snapshot does not establish untracked repository
> files as the accepted transfer mechanism for future single-file reviews.

### 5.2 Committed versus working-tree contents

- Tracked-path counts used `git ls-files`; they therefore counted paths known to
  the Git index at the inspected revision and working tree.
- File contents used for empty-file detection, AST parsing, import probing,
  notebook inspection, documentation inspection, and relationship searches
  were read from the working tree.
- The two modified magnetism notebooks were therefore inspected in their
  working-tree form where included in bounded notebook searches.
- No conclusion was drawn from the difference between their committed and
  working-tree forms.
- The inspected source and test paths were not reported as modified by Git at
  baseline collection.

### 5.3 Python environment

System Python used for static scripts:

```text
Executable: /opt/homebrew/bin/python3
Version:    Python 3.14.6
```

Existing project virtual environment examined for imports and pytest
availability:

```text
Executable: /Users/eugene/repos/physkit/.venv/bin/python
Version:    Python 3.14.6
```

No environment was created or changed.

### 5.4 Reproduction commands and scripts

The following commands or equivalent inline scripts were used for the reported
counts and bounded static findings.

#### Repository identity and working-tree state

```bash
git branch --show-current
git rev-parse HEAD
git status --short --branch
```

#### Tracked-file counts

```bash
git ls-files 'src/**/*.py' | wc -l
git ls-files 'tests/**/*.py' | wc -l
git ls-files '*.ipynb' | wc -l
git ls-files 'docs/**' | wc -l
```

#### Empty-file detection and Python syntax parsing

```bash
python3 - <<'PY'
from pathlib import Path
import ast

for root in (Path("src"), Path("tests")):
    for path in sorted(root.rglob("*.py")):
        if not path.is_file():
            continue
        text = path.read_text(errors="replace")
        if not text.strip():
            print("EMPTY", path)
        try:
            ast.parse(text)
        except SyntaxError as error:
            print("SYNTAX", path, error.lineno, error.msg)
PY
```

#### Source-declared class and function discovery

```bash
python3 - <<'PY'
from pathlib import Path
import ast

for path in sorted(Path("src/physkit").rglob("*.py")):
    if not path.is_file():
        continue
    text = path.read_text(errors="replace")
    try:
        tree = ast.parse(text)
    except SyntaxError as error:
        print(path, "SYNTAX_ERROR", error.lineno, error.msg)
        continue
    names = [
        node.name
        for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef))
    ]
    if names:
        print(path, names)
PY
```

This discovery records declarations only. It does not establish that a declared
object imports, executes, or belongs to the public API.

#### Module discovery and import probing

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src .venv/bin/python - <<'PY'
import importlib
import pkgutil
import physkit

modules = sorted(
    module.name
    for module in pkgutil.walk_packages(
        physkit.__path__,
        physkit.__name__ + ".",
    )
)
print("modules discovered", len(modules))
for module in modules:
    try:
        importlib.import_module(module)
    except Exception as error:
        print("IMPORT_ERROR", module, type(error).__name__, str(error))
PY
```

`PYTHONDONTWRITEBYTECODE=1` was used to avoid generating bytecode during the
probe.

#### Notebook JSON parsing, code-cell references, and saved errors

```bash
python3 - <<'PY'
from pathlib import Path
import json

for path in sorted(Path("notebooks").rglob("*.ipynb")):
    if ".ipynb_checkpoints" in path.parts:
        continue
    try:
        notebook = json.loads(path.read_text(errors="replace"))
    except Exception as error:
        print("INVALID", path, type(error).__name__)
        continue

    cells = notebook.get("cells", [])
    code_cells = [
        cell for cell in cells
        if cell.get("cell_type") == "code"
    ]
    source = "\n".join(
        "".join(cell.get("source", []))
        for cell in code_cells
    )
    has_physkit_code = "physkit" in source
    has_saved_error = any(
        output.get("output_type") == "error"
        for cell in code_cells
        for output in cell.get("outputs", [])
    )
    print(
        path,
        len(cells),
        has_physkit_code,
        has_saved_error,
    )
PY
```

The `physkit` count was based on executable code-cell text, not arbitrary
notebook metadata or rendered output.

#### Duplicate-notebook detection

```bash
find notebooks -type f -name '*.ipynb' \
  ! -path '*/.ipynb_checkpoints/*' \
  -exec shasum {} + | sort
```

Files with identical complete-file SHA-1 output were treated as exact duplicate
content for this bounded check.

#### Test-import search

```bash
python3 - <<'PY'
from pathlib import Path
import ast
import collections

imports = collections.Counter()
for path in Path("tests").rglob("*.py"):
    if not path.is_file():
        continue
    tree = ast.parse(path.read_text(errors="replace"))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            if node.module and node.module.startswith("physkit"):
                imports[node.module] += 1
        elif isinstance(node, ast.Import):
            for name in node.names:
                if name.name.startswith("physkit"):
                    imports[name.name] += 1

for module, count in sorted(imports.items()):
    print(module, count)
PY
```

Where this document says that the stated inspection found no matching test
import, documentation page, or notebook relationship, that statement refers
only to these bounded repository searches. It is not an unrestricted claim that
no evidence exists elsewhere.

### 5.5 Unavailable tools

- `pytest` was not installed in the existing project virtual environment.
- `pytest` was also unavailable through the inspected system Python.
- No test dependency was installed.
- No accepted notebook-execution harness was identified or invoked.

### 5.6 Checks explicitly not performed

- pytest collection;
- unit or integration test execution;
- notebook execution;
- documentation builds;
- package build or installation;
- type checking;
- linting;
- code coverage;
- numerical convergence studies;
- comparison with independent physical data;
- pedagogical evaluation;
- uncertainty propagation;
- network, cluster, cloud, or HPC execution.

## 6. Repository-level findings

### 6.1 Stated project status conflicts with observed contents

- `FACT` — `README.md` describes PhysKit as “Pre-implementation.”
- `FACT` — `README.md` states “No functional modules implemented yet.”
- `FACT` — the stated inspection found non-empty source declarations for units,
  constants, grids, boundary conditions, periodic lattices, quantum models,
  solid-state models, thermodynamic utilities, deposition geometry, and other
  topics.
- `FACT` — associated tests and notebooks were found for some, but not all, of
  those source groups.
- `FACT` — non-empty source declarations do not establish that the declared
  objects parse, import, execute, satisfy a contract, or form supported
  capabilities.
- `INFERENCE` — a single project-wide status sentence may be insufficient to
  summarize the heterogeneous repository contents.
- `HUMAN DECISION` — determine whether future project status should be
  represented by capability-level records, with `README.md` containing only a
  derived summary.
- `HUMAN DECISION` — determine which accepted records and process may produce
  or update such a summary.

This snapshot must not itself become the source of current project status.

### 6.2 Source breadth exceeds the stated inspection’s test and documentation relationships

- `FACT` — the bounded AST test-import search found imports from constants,
  core grids, discretization, numeric helpers, and units.
- `FACT` — the stated inspection found no matching test imports for the sampled
  chemistry, deposition, materials, math-operator, mechanics, newer numerics,
  periodic, newer quantum, thermodynamic, or visualization source groups.
- `FACT` — solid-state lattice tests import test-local lattice classes rather
  than the current `physkit.periodic` package.
- `FACT` — five tracked documentation files were counted under `docs/`.
- `FACT` — the tracked documentation found by the inspection is concentrated on
  temperature, vapor pressure, and one semiconductor mobility topic.
- `INFERENCE` — repository breadth may reflect different histories or intended
  uses, but no accepted lifecycle exists with which to classify them.
- `HUMAN DECISION` — a future policy must decide what source, test, notebook,
  documentation, and review relationships are required for any lifecycle state.

### 6.3 Observed tests are not executed passing tests

- `FACT` — 20 tracked Python test paths were counted.
- `FACT` — three inspected test Python files were empty.
- `FACT` — pytest collection and execution were unavailable.
- `FACT` — test source contains assertions concerning selected software
  contracts, unit conversions, grid conventions, lattice relationships, and
  error behavior.
- `FACT` — these assertions were observed but not executed.
- `FACT` — even successful execution would not by itself establish numerical
  verification, physical validation, pedagogical validation, or UQ.

### 6.4 Notebook breadth and condition vary

For notebooks under `notebooks/`, excluding checkpoint copies:

- `FACT` — 127 notebook files were inspected.
- `FACT` — 124 parsed as JSON.
- `FACT` — three did not parse as JSON.
- `FACT` — 27 valid notebooks contained `physkit` in executable code-cell text.
- `FACT` — 15 valid notebooks retained at least one saved error output.
- `FACT` — 11 notebooks contained no more than one cell.
- `FACT` — exact duplicate notebook content was found in multiple pairs,
  including:
  - `notebooks/math/fourier_series.ipynb` and
    `notebooks/math/fourier/fourier_series.ipynb`;
  - `notebooks/solidstate/dev_tightbinding.ipynb` and
    `notebooks/solidstate/solst01_tightbinding.ipynb`;
  - `notebooks/qm/qm.piab1d.comp.ipynb` and
    `notebooks/qm/qm.piab1d.physkitlib.ipynb`.
- `INFERENCE` — duplicate notebooks may be intentional variants, historical
  copies, or accidental duplication.
- `HUMAN DECISION` — no interpretation or canonical selection is made here.

### 6.5 Repository organization mixes artifact types

- `FACT` — notebooks occur under `notebooks/`, `examples/`, `dev/`, and
  `src/physkit/plasmas/gas_discharge/`.
- `FACT` — standalone Python definitions also occur under `notebooks/`.
- `INFERENCE` — the current layout may mix reusable source, teaching material,
  experiments, course-specific material, and historical artifacts.
- `HUMAN DECISION` — establish future boundaries among production source,
  canonical notebooks, exploratory notebooks, examples, generated evidence,
  and historical artifacts.

### 6.6 Repository-wide validation evidence gaps

The stated inspection found no accepted repository artifact that establishes:

- a capability lifecycle;
- a canonical-notebook inventory;
- a complete source-to-test-to-notebook-to-documentation mapping;
- a repository-wide import gate;
- a notebook execution gate;
- a numerical-verification hierarchy;
- general tolerance ownership;
- physical-validation acceptance;
- pedagogical-validation acceptance;
- UQ applicability or acceptance;
- capability-level current status authority.

These are bounded findings from the stated repository inspection. They do not
claim that no evidence exists outside the inspected repository or that an
unrecorded capability lacks value.

## 7. Capability groups discovered

### 7.1 Numerical values, constants, and unit conversion

**Classification:** pending

**Repository locations**

- `src/physkit/numeric.py`
- `src/physkit/types.py`
- `src/physkit/constants.py`
- `src/physkit/units/`
- `src/physkit/chemistry/stoichiometry.py`
- `src/physkit/units/mass_molar.py`

**Source-declared objects and functions**

- `as_f64_array`
- `ConstantsSI`, `ConstantsGaussianCGS`
- `UnitQuantityProtocol`
- `Pressure`, `Temperature`, `Energy`, `Force`, `Length`, `Mass`, `Time`,
  `Charge`, `Velocity`, `Torque`, `Viscosity`, `Dipole`, `Density`, and
  `ElectricField`
- named unit-system containers
- molar-mass and particle-mass representations

These are source declarations only. This list does not declare them public,
correct, supported, or complete.

**Associated tests**

- `tests/physkit/constants/`
- `tests/physkit/units/test_energy.py`
- `tests/physkit/units/test_force.py`
- `tests/physkit/units/test_pressure.py`
- `tests/physkit/units/test_temperature.py`

**Associated notebooks**

- `notebooks/units/00_units_overview.ipynb`
- `notebooks/units/units_pressure.ipynb`
- `notebooks/units/units_temperature.ipynb`
- `notebooks/units.ipynb`
- additional physics notebooks containing unit or constant imports

**Associated documentation**

- `docs/physkit.units.temperature.md`
- links from `docs/physkit.00.md`

**Capability-specific condition**

- `FACT` — sampled units modules imported during the package sweep.
- `FACT` — the temperature source and temperature documentation use
  inconsistent API examples and names.

**Specific overlap**

- `FACT` — molar- and particle-mass declarations occur in both
  `chemistry/stoichiometry.py` and `units/mass_molar.py`.
- `INFERENCE` — these may be complementary or competing representations.

**Distinctive observed evidence**

- `FACT` — test source contains known conversion-factor, protocol,
  vectorization, and round-trip assertions.

**Unresolved human decisions**

- `HUMAN DECISION` — determine the intended public quantity API.
- `HUMAN DECISION` — determine whether unit-system containers and individual
  quantity classes belong to one coherent contract.
- `HUMAN DECISION` — determine the relationship between the two molar-mass
  representations.
- `HUMAN DECISION` — determine required baseline notebooks and
  numerical-verification cases.

### 7.2 Grids, state representations, boundaries, and discretization

**Classification:** pending

**Repository locations**

- `src/physkit/core/grids.py`
- `src/physkit/core/state.py`
- `src/physkit/core/boundaries.py`
- `src/physkit/discretization/grid_1d.py`
- `src/physkit/visualize/grids.py`
- `notebooks/basic/grid_1d.py`

**Source-declared objects and functions**

- `CartesianAxis`
- `CartesianGrid1D`, `CartesianGrid2D`, `CartesianGrid3D`
- `ActiveSet1D`
- `core.state.Grid1D`, `Wavefunction1D`
- `discretization.grid_1d.Grid1D`, `ActiveSetType1D`
- Dirichlet, Neumann, Robin, periodic, Bloch, asymptotic-decay, and radiation
  boundary-condition classes
- grid plotting functions

**Associated tests**

- `tests/physkit/core/test_CartesianAxis.py`
- `tests/physkit/core/test_CartesianGrid1D.py`
- `tests/physkit/core/test_CartesianGrid2D.py`
- `tests/physkit/core/test_CartesianGrid3D.py`
- `tests/physkit/discretization/test_ActiveSetType1D.py`
- `tests/physkit/discretization/test_Grid1D.py`

**Associated notebooks**

- `notebooks/basic/01_discretization_1d.ipynb`
- `notebooks/core/grids/cartesian-axis.ipynb`
- `notebooks/core/grids/cartesian-grid-1d.ipynb`
- `notebooks/scratch/discretization/grid_1d_conventions.ipynb`
- `notebooks/numerics/differentiation/finite-difference-laplacian-1d.ipynb`

**Associated documentation**

- The stated inspection found no dedicated tracked grid or discretization page
  under `docs/`.

**Capability-specific condition**

- `FACT` — `physkit.core.grids` and `physkit.discretization.grid_1d`
  imported during inspection.
- `FACT` — `CartesianGrid1D` tests expect names including `x_min`, `x_max`,
  `x_axis`, `Lx`, and `dx`.
- `FACT` — the inspected source uses different names including `x_lower`,
  `x_upper`, `x`, and `delta`.

**Specific overlap**

- `FACT` — at least three distinct one-dimensional grid representations were
  found:
  - an interior-state-oriented grid in `core.state`;
  - endpoint-configurable Cartesian grids in `core.grids`;
  - a closed reference grid with selectable active sets in
    `discretization.grid_1d`.
- `FACT` — another grid definition occurs under `notebooks/basic/`.
- `INFERENCE` — these may represent different mathematical purposes,
  historical alternatives, or unresolved API experiments.

**Distinctive observed evidence**

- `FACT` — test source asserts spacing, endpoints, active-index sets, shapes,
  dtypes, array immutability, and constructor validation.
- `FACT` — notebooks expose grid and discretization constructions.

**Unresolved human decisions**

- `HUMAN DECISION` — identify the physical and mathematical purpose of each
  grid representation.
- `HUMAN DECISION` — select or reject naming and endpoint conventions.
- `HUMAN DECISION` — decide whether active sets belong to grids,
  discretizations, boundary handling, or another object.
- `HUMAN DECISION` — define a baseline agreement case before selecting any
  canonical notebook.

### 7.3 Mathematical and finite-difference operators

**Classification:** pending

**Repository locations**

- `src/physkit/core/operator.py`
- `src/physkit/math/operators/`
- `src/physkit/numerics/finite_difference.py`
- `src/physkit/numerics/differentiation/laplacian.py`

**Source-declared objects and functions**

- two `LinearOperator1D` declarations
- symbolic `ContinuousOperator1D` and `PoissonOperator1D` declarations
- another operator module using the same class names
- `Laplacian1D`
- a module named `laplacian.py` that declares `Piab1D`

**Associated tests**

- The bounded test-import inspection found no direct matching imports.

**Associated notebooks**

- `notebooks/math_poisson1d.ipynb`
- `notebooks/fem/dev.01_discretization1d.ipynb`
- `notebooks/numerics/differentiation/finite-difference-laplacian-1d.ipynb`
- related scratch numerical notebooks

**Associated documentation**

- The stated inspection found no dedicated tracked operator or
  finite-difference page under `docs/`.

**Capability-specific condition**

- `FACT` — `physkit.math.operators` fails import because it requests
  `physkit.math.operators.state`.
- `FACT` — `physkit.numerics.finite_difference` requests
  `physkit.core.bc`.
- `FACT` — `physkit.numerics.differentiation` requests
  `physkit.core.constants`.
- `FACT` — those requested modules were not found at the referenced paths.
- `FACT` — static inspection also found unresolved or inconsistent names in
  sampled symbolic operator modules.

**Specific overlap**

- `FACT` — `core/operator.py` and `math/operators/base.py` contain near-duplicate
  `LinearOperator1D` declarations.
- `FACT` — continuous and discrete operator modules reuse class names for
  different apparent representations.
- `INFERENCE` — boundaries among mathematical operators, matrices,
  discretization methods, and physics models are not consistently expressed.

**Distinctive observed evidence**

- `FACT` — a finite-difference pedagogical notebook exists.
- `FACT` — the stated inspection found no matching convergence-test import or
  declared convergence requirement for this source group.

**Unresolved human decisions**

- `HUMAN DECISION` — define the distinction among mathematical operator,
  discrete representation, matrix construction, and solver.
- `HUMAN DECISION` — determine whether duplicate declarations are alternatives,
  layers, or historical remnants.
- `HUMAN DECISION` — determine their intended relationship to grid and boundary
  objects.

### 7.4 Periodic lattices, reciprocal space, modes, and visualization

**Classification:** pending

**Repository locations**

- `src/physkit/periodic/`
- `src/physkit/periodic/lattice.py`
- `src/physkit/periodic/lattice/`
- `src/physkit/periodic/visualization/`
- portions of `src/physkit/solidstate/`

**Source-declared objects and functions**

- abstract lattice, direct-lattice, reciprocal-lattice, Wigner–Seitz-cell, and
  first-Brillouin-zone classes
- one-, two-, and three-dimensional lattice families
- `ReciprocalModeBasis1D`
- `BlochPhase1D`, `BlochMode1D`
- one-, two-, and three-dimensional k-point path classes
- periodic-mode visualization
- additional standalone `DirectLattice1D` and `ReciprocalLattice1D`
  declarations

**Associated tests**

- lattice tests under `tests/physkit/solidstate/`
- those tests import `tests/physkit/solidstate/lattice.py`, a test-local
  implementation, rather than the current `physkit.periodic` package

**Associated notebooks**

- `notebooks/lattice/lattice-1d.ipynb`
- `notebooks/lattice/lattice-2d.ipynb`
- `notebooks/lattice/lattice-3d.ipynb`
- `notebooks/periodic/periodic-subpackage.ipynb`
- related solid-state and semiconductor notebooks

**Associated documentation**

- The stated inspection found no dedicated tracked periodic or lattice page
  under `docs/`.

**Capability-specific condition**

- `FACT` — the `physkit.periodic` package imported during inspection.
- `FACT` — both `periodic/lattice.py` and `periodic/lattice/` exist.
- `FACT` — the inspected package initializer exports objects from the package
  directory rather than from the standalone module.

**Specific overlap**

- `FACT` — direct- and reciprocal-lattice declarations occur in:
  - `periodic/lattice.py`;
  - `periodic/lattice/lattice1d.py`;
  - `periodic/reciprocal.py`;
  - `tests/physkit/solidstate/lattice.py`;
  - notebook-local code.
- `INFERENCE` — the standalone module may be shadowed by the package layout,
  but this snapshot does not classify it as obsolete or dead.

**Distinctive observed evidence**

- `FACT` — test source asserts basis storage, dimensions, cell measure,
  lattice-vector construction, reciprocal duality, and invalid-shape behavior.
- `FACT` — because those tests use test-local classes, they do not directly
  exercise the current `physkit.periodic` declarations.
- `FACT` — several notebooks demonstrate periodic concepts and some contain
  PhysKit imports.

**Unresolved human decisions**

- `HUMAN DECISION` — identify the intended public lattice representation.
- `HUMAN DECISION` — determine whether test-local lattice classes are fixtures,
  prototypes, or alternatives.
- `HUMAN DECISION` — establish matrix orientation, units, indexing, and
  reciprocal conventions.
- `HUMAN DECISION` — define an explicit-construction/library agreement case
  before selecting a canonical notebook.

### 7.5 Quantum particle-in-a-box and one-dimensional solvers

**Classification:** pending

**Repository locations**

- `src/physkit/qm/models/`
- `src/physkit/qm/solvers/piab1d/`
- `src/physkit/qm/solver1d.py`
- `src/physkit/qm/well1d.py`
- `src/physkit/solidstate/piab1d.py`

**Source-declared objects and functions**

- abstract quantum model, result, and solver classes
- generic PIAB model, result, and solver classes
- one-, two-, and three-dimensional PIAB model families
- analytical and numerical one-dimensional TISE solver classes
- `InfiniteSquareWell1D` and Hamiltonian/reconstruction helpers
- a separate solid-state `ParticleInABox1D` declaration

**Associated tests**

- `tests/physkit/solidstate/test_ParticleInABox1D.py`
- `FACT` — the inspected test file contains pass-only class declarations and
  does not contain substantive PIAB assertions.

**Associated notebooks**

- root-level and `notebooks/qm/` PIAB notebooks
- analytical, symbolic, numerical, comparison, and library-named variants
- additional PIAB notebooks under `notebooks/solidstate/` and
  `notebooks/scratch/`

**Associated documentation**

- The stated inspection found no dedicated tracked PIAB page under `docs/`.

**Capability-specific condition**

- `FACT` — `src/physkit/qm/models/piab3d.py` does not parse.
- `FACT` — multiple PIAB modules fail import because Python cannot construct the
  declared `ABC`/`Generic` method-resolution order.
- `FACT` — `src/physkit/qm/solver1d.py` imports a boundary name not supplied by
  the inspected `core.boundaries` module.
- `FACT` — sampled legacy solver code contains unresolved attribute and
  assignment discrepancies.
- `FACT` — at least one PIAB notebook requests `physkit.core.bc`.
- `FACT` — `examples/solidstate/smoke_particleinabox.py` contains a syntax
  error.

**Specific overlap**

- `FACT` — PIAB concepts are declared in several model, solver, well, and
  solid-state modules.
- `FACT` — several notebooks construct PIAB mathematics directly without an
  accepted relationship to one library representation.
- `INFERENCE` — these may represent multiple generations or experiments.

**Distinctive observed evidence**

- `FACT` — notebooks contain analytical and numerical constructions and saved
  outputs.
- `FACT` — some saved PIAB notebooks contain error outputs.
- `FACT` — the stated inspection found no accepted agreement case or tolerance
  connecting an explicit construction to one source-declared library path.

**Unresolved human decisions**

- `HUMAN DECISION` — select the physical model, boundary conditions, state
  space, units, and energy convention.
- `HUMAN DECISION` — select the intended mathematical and numerical
  formulations.
- `HUMAN DECISION` — decide whether any existing API should be retained.
- `HUMAN DECISION` — define baseline cases and agreement tolerances.
- `HUMAN DECISION` — decide whether PIAB is suitable as an initial vertical
  pedagogical capability.

### 7.6 Solid-state electronic and phonon models

**Classification:** pending

**Repository locations**

- `src/physkit/solidstate/electronic/`
- `src/physkit/solidstate/phonons/`
- `src/physkit/solidstate/visualization/`
- supporting declarations under `src/physkit/periodic/`

**Source-declared objects and functions**

- Fourier representations of periodic potentials
- a one-dimensional plane-wave Bloch model, Hamiltonian, solver, and band result
- a one-dimensional monatomic-chain model, dynamical matrix, solver, and phonon
  result
- electronic- and phonon-band visualizers

**Associated tests**

- The bounded test-import inspection found no imports of these current source
  modules.

**Associated notebooks**

- Bloch, tight-binding, Kronig–Penney, direct-lattice, phonon-adjacent, and
  electronic-band notebooks under `notebooks/solidstate/`,
  `notebooks/semphy/`, and `notebooks/scratch/solidstate/`

**Associated documentation**

- The stated inspection found no dedicated tracked documentation for these
  source groups.

**Capability-specific condition**

- `FACT` — sampled electronic and phonon modules imported during the package
  sweep.
- `FACT` — their numerical behavior was not executed.
- `FACT` — they depend on periodic declarations for which competing
  representations exist.

**Specific overlap**

- `FACT` — related Bloch and tight-binding constructions occur across source
  modules and several independent notebooks.
- `INFERENCE` — source declarations may have been extracted from notebook
  experiments, but their history and intended authority require confirmation.

**Distinctive observed evidence**

- `FACT` — saved notebook calculations and plotting outputs exist for related
  topics.
- `FACT` — the stated inspection found no direct test imports for the current
  electronic or phonon source modules.

**Unresolved human decisions**

- `HUMAN DECISION` — determine intended physical models, assumptions, and
  units.
- `HUMAN DECISION` — determine which calculations are illustrative and which
  may become reference cases.
- `HUMAN DECISION` — decide whether these topics are in scope for an early
  PhysKit capability.
- `HUMAN DECISION` — determine which, if any, notebook should later be proposed
  as canonical.

### 7.7 Thermodynamics, vapor pressure, mixtures, and transport-adjacent material

**Classification:** pending

**Repository locations**

- `src/physkit/thermo/vaporpressure/`
- `src/physkit/materials/mixtures/`
- related scripts under `notebooks/thinfilm/`
- transport-related notebooks and documentation

**Source-declared objects and functions**

- `VaporPressureBase`
- `VaporPressureCurve` protocol
- empty Antoine, data, fitting, and models modules
- `MolarMixture`
- standalone vapor-pressure, mixture, and thin-film definitions under
  `notebooks/`

**Associated tests**

- The bounded test-import inspection found no direct imports of these source
  modules.

**Associated notebooks**

- vapor-pressure, Hertz–Knudsen, evaporation, ideal-gas, mixing-enthalpy,
  mixture, and semiconductor transport notebooks

**Associated documentation**

- `docs/physkit.thermo.vaporpressure.md`
- `docs/physkit.thermo.vaporpressure.VaporPressureCurve.md`
- `docs/sempy/transport/physkit.semi.transport.mobility_lattice.md`
- links from `docs/physkit.00.md`

**Capability-specific condition**

- `FACT` — `physkit.thermo.vaporpressure` fails import because its initializer
  requests `VaporPressureCurveBase`, which is not declared by the inspected
  base module.
- `FACT` — `VaporPressureBase` refers to a validity-range field name different
  from its declared dataclass field.
- `FACT` — it calls `Temperature.convert` and
  `Temperature.check_in_range` with keyword signatures inconsistent with the
  current temperature source.
- `FACT` — several named vapor-pressure modules are empty.

**Specific overlap**

- `FACT` — reusable-looking vapor-pressure and mixture definitions occur both
  under `src/` and in notebook-adjacent Python files.
- `FACT` — documentation references module and object paths that differ from
  the inspected source layout.
- `INFERENCE` — some notebook-side definitions may predate or compete with the
  source package.

**Distinctive observed evidence**

- `FACT` — conceptual documentation and saved notebook calculations exist.
- `FACT` — validity-range metadata appears in source, but its provenance and
  acceptance were not assessed.

**Unresolved human decisions**

- `HUMAN DECISION` — define the intended vapor-pressure API and correlation
  ownership.
- `HUMAN DECISION` — establish parameter provenance, units, validity ranges,
  and extrapolation policy.
- `HUMAN DECISION` — determine which material is general PhysKit capability
  versus course- or application-specific material.
- `HUMAN DECISION` — determine physical-validation and UQ requirements for
  empirical correlations.

### 7.8 Deposition, thin-film geometry, and plasma notebooks

**Classification:** pending

**Repository locations**

- `src/physkit/deposition/`
- `src/physkit/plasmas/gas_discharge/`
- `notebooks/thinfilm/`
- `notebooks/vacuum/`

**Source-declared objects and functions**

- plane-point geometry
- disk-source deposition classes and quadrature
- radial kernel functions
- a minimal Townsend Python module
- thin-film helper definitions under `notebooks/`

**Associated tests**

- The bounded test-import inspection found no direct imports for this source
  group.

**Associated notebooks**

- deposition, evaporation, Knudsen-number, vapor-pressure, gas-discharge, and
  vacuum notebooks
- three notebooks located inside the source tree

**Associated documentation**

- The stated inspection found no dedicated tracked deposition or plasma page
  under `docs/`.

**Capability-specific condition**

- `FACT` — static inspection found unresolved deposition annotation and typing
  names.
- `FACT` — two source-tree plasma notebooks do not parse as JSON.
- `FACT` — `townsend.py` contains no class or function declaration.

**Specific overlap**

- `FACT` — related computations occur in both source modules and notebook-local
  Python files.
- `INFERENCE` — the intended source/notebook boundary for these topics has not
  been established.

**Distinctive observed evidence**

- `FACT` — saved notebook calculations and plots exist for related topics.

**Unresolved human decisions**

- `HUMAN DECISION` — determine the student-facing capability supported by this
  source group.
- `HUMAN DECISION` — decide whether empirical and geometric assumptions are
  sufficiently explicit.
- `HUMAN DECISION` — determine whether notebooks belong inside production
  source.

### 7.9 Chemistry, elements, and mixture utilities

**Classification:** pending

**Repository locations**

- `src/physkit/elements.py`
- `src/physkit/chemistry/`
- `src/physkit/materials/mixtures/`

**Source-declared objects and functions**

- element and isotope classes
- an element dictionary and SQLite-generation helpers
- molar- and particle-mass classes
- molar-mixture calculations

**Associated tests**

- The bounded test-import inspection found no direct imports for these source
  groups.

**Associated notebooks**

- `notebooks/basic/elements.ipynb`
- `notebooks/mixtures/mixture_molefraction.ipynb`
- `notebooks/mixtures/sublattice.ipynb`
- related material and thin-film notebooks

**Associated documentation**

- links from `docs/physkit.00.md`
- The stated inspection found no complete dedicated tracked capability page.

**Capability-specific condition**

- `FACT` — sampled modules imported during inspection.
- `FACT` — data provenance and runtime behavior were not evaluated.

**Specific overlap**

- `FACT` — mass representations occur under both chemistry and units.
- `INFERENCE` — element database generation may be infrastructure supporting
  multiple future pedagogical capabilities.

**Distinctive observed evidence**

- `FACT` — notebooks demonstrate element and mixture calculations.

**Unresolved human decisions**

- `HUMAN DECISION` — identify the intended student- or instructor-facing
  capabilities.
- `HUMAN DECISION` — establish data sources, units, and provenance.
- `HUMAN DECISION` — resolve the relationship among mass representations.

### 7.10 Mechanics and general numerical-method notebooks

**Classification:** pending

**Repository locations**

- `src/physkit/mechanics/`
- notebook groups under `notebooks/kinematics/`, `notebooks/math/`,
  `notebooks/fem/`, `notebooks/numerics/`, and repository-root notebook paths

**Source-declared objects and functions**

- `src/physkit/mechanics/oscillator.py` is empty.
- notebook-local definitions cover motion, integration, Fourier methods,
  Monte Carlo, differentiation, Poisson problems, and related topics.

**Associated tests**

- The bounded test-import inspection found no matching mechanics or
  general-method source imports.

**Associated notebooks**

- simple harmonic motion
- kinematics
- trapezoid rule
- Fourier series
- Monte Carlo Ising
- Poisson and finite-element experiments
- other numerical-method notebooks

**Associated documentation**

- The stated inspection found no dedicated tracked mechanics or general
  numerical-method page under `docs/`.

**Capability-specific condition**

- `FACT` — the inspected mechanics source module is empty.
- `FACT` — some notebook-local numerical definitions overlap with operator and
  finite-difference source modules that have import discrepancies.

**Specific overlap**

- `FACT` — algorithms are declared directly in multiple notebooks while related
  reusable source is absent, empty, or inconsistent.
- `INFERENCE` — these may be lesson prototypes rather than library
  capabilities.

**Distinctive observed evidence**

- `FACT` — narrative, equations, code, and saved outputs occur in the related
  notebooks.

**Unresolved human decisions**

- `HUMAN DECISION` — decide which notebook topics should be considered for
  reusable PhysKit abstractions.
- `HUMAN DECISION` — decide whether a notebook-only lesson may receive a
  capability record before a library abstraction exists.

### 7.11 Notebook-only or predominantly notebook-based physics domains

**Classification:** pending

**Repository locations**

Predominantly notebook-based material was observed for:

- electromagnetism;
- magnetism;
- optics and waves;
- statistical mechanics;
- special relativity;
- viscoelasticity;
- materials-science course material;
- vacuum and thin-film topics;
- semiconductor teaching material.

**Source-declared objects and functions**

- primarily notebook-local functions and direct NumPy, SciPy, or SymPy
  constructions
- no complete corresponding source-package group was identified for several of
  these domains by the stated inspection

**Associated tests**

- The bounded test-import inspection found no matching library imports for most
  of these domains.

**Associated notebooks**

- numerous domain-specific notebooks under corresponding root, domain,
  scratch, and course-specific directories

**Associated documentation**

- limited tracked documentation, including one lattice-scattering mobility page

**Capability-specific condition**

- `FACT` — notebook JSON and saved-output conditions vary.
- `FACT` — the two magnetism notebooks were already modified in the working
  tree at baseline collection.
- `FACT` — their changes were not interpreted, compared with `HEAD`, or
  modified.

**Specific overlap**

- `FACT` — similar topics occur in root notebook paths, domain directories,
  scratch directories, and course-specific directories.
- `INFERENCE` — directory naming may encode audience, history, or intended use,
  but no accepted convention was found.

**Distinctive observed evidence**

- `FACT` — narrative, equations, visible code, and saved outputs exist in many
  notebooks.
- `FACT` — saved outputs and absence of saved errors do not establish agreement
  with a library abstraction or attainment of a learning objective.

**Unresolved human decisions**

- `HUMAN DECISION` — determine whether and how notebook-only work enters future
  capability-level records.
- `HUMAN DECISION` — identify intended learners and learning objectives.
- `HUMAN DECISION` — classify notebook roles only after a lifecycle policy is
  accepted.
- `HUMAN DECISION` — preserve and separately review the existing magnetism
  notebook changes without interference from control-plane planning.

## 8. Source–test–notebook–documentation relationships

The following table summarizes bounded relationships found by the stated
inspection without assigning maturity.

| Capability group | Source declarations | Observed test relationship | Observed notebook relationship | Observed documentation relationship | Classification |
|---|---|---|---|---|---|
| Constants and units | Broad quantity and conversion declarations | Focused test source for constants and selected units | Several notebooks containing PhysKit unit use | Temperature page | pending |
| Grids and discretization | Several competing grid/state declarations | Substantial test source with API drift | Direct-construction and API-oriented notebooks | No dedicated page found by the stated search | pending |
| Mathematical/numerical operators | Overlapping declarations with import discrepancies | No direct test import found by the stated search | Several numerical notebooks | No dedicated page found by the stated search | pending |
| Periodic lattices | Multiple source and test-local declarations | Tests exercise test-local alternatives | Several periodic and lattice notebooks | No dedicated page found by the stated search | pending |
| PIAB quantum models | Several model, solver, and well declarations | Pass-only test surface | Many overlapping notebooks | No dedicated page found by the stated search | pending |
| Solid-state bands and phonons | Model/operator/solver/result declarations | No direct test import found by the stated search | Many related notebooks | No dedicated page found by the stated search | pending |
| Thermodynamics and vapor pressure | Partial declarations with API conflicts | No direct test import found by the stated search | Several related notebooks | Multiple partial pages | pending |
| Deposition and plasma | Partial source declarations and source-tree notebooks | No direct test import found by the stated search | Several application notebooks | No dedicated page found by the stated search | pending |
| Chemistry and mixtures | Element, mass, and mixture declarations | No direct test import found by the stated search | Elements and mixture notebooks | Index links only | pending |
| Mechanics/general methods | Empty mechanics source and notebook-local definitions | No matching import found by the stated search | Broad notebook collection | No dedicated page found by the stated search | pending |
| Other physics domains | Predominantly notebook-local definitions | No matching library import found for most sampled domains | Broad domain coverage | Minimal tracked coverage | pending |

This table is not a support matrix. An observed test relationship means test
source was found; it does not mean that the tests collect, pass, or establish a
validation claim.

## 9. Competing and overlapping representations

Material overlaps requiring later human review include:

1. multiple one-dimensional grid abstractions;
2. duplicate `LinearOperator1D` declarations;
3. repeated continuous/discrete operator class names;
4. both a `periodic/lattice.py` module and a `periodic/lattice/` package;
5. multiple direct- and reciprocal-lattice declarations;
6. test-local lattice classes separate from current package declarations;
7. multiple PIAB model, solver, well, and solid-state declarations;
8. reusable-looking algorithms duplicated between source and notebooks;
9. `visualize` and `viz` naming;
10. overlapping chemistry and units mass representations;
11. source-tree notebooks mixed with Python modules;
12. exact duplicate notebook files.

No overlap is resolved by this baseline. In particular:

- no older representation is declared obsolete;
- no newer representation is declared preferred;
- no notebook is declared canonical;
- no declaration is identified as the accepted public API;
- no package structure is declared canonical beyond the observed packaging
  configuration in `pyproject.toml`.

## 10. Syntax, import, and API discrepancies

### 10.1 Syntax and file-format discrepancies

- `src/physkit/qm/models/piab3d.py` did not parse as Python.
- `examples/solidstate/smoke_particleinabox.py` did not parse as Python.
- three notebooks under `notebooks/` did not parse as JSON.
- two notebooks under `src/physkit/plasmas/gas_discharge/` did not parse as
  JSON.

### 10.2 Import discrepancies

The import probe reported failures involving:

- missing `physkit.math.operators.state`;
- missing `physkit.core.constants`;
- missing `physkit.core.bc`;
- inconsistent `ABC` and `Generic` inheritance in PIAB model classes;
- a boundary name expected by a legacy solver but not supplied by the inspected
  boundary module;
- missing `VaporPressureCurveBase`;
- cascading PIAB solver import failures.

These are software-condition observations. They do not determine the correct
physical, mathematical, API, or architectural repair.

### 10.3 API discrepancies

The stated inspection observed:

- Cartesian-grid tests expecting names different from current source
  declarations;
- vapor-pressure source calling temperature conversion with incompatible
  keywords;
- vapor-pressure source referring to inconsistent validity-range field names;
- temperature documentation using names and signatures different from source;
- thermodynamics documentation referring to paths and classes not found at
  those paths;
- notebooks requesting source paths not found by the inspection;
- deposition annotations referring to unresolved names.

`HUMAN DECISION` — determine whether each discrepancy should eventually be
resolved in source, tests, notebooks, documentation, or an accepted
combination. This snapshot does not assume that the source declaration is
correct merely because it is under `src/`.

## 11. Evidence gaps

The repository-wide evidence gaps are recorded once in Section 6.6. Additional
cross-cutting gaps relevant to later planning are:

- no accepted mapping from source declarations to software requirements;
- no accepted mapping from explicit notebook constructions to library APIs;
- no declared owner for comparison tolerances;
- no deterministic complete per-file inventory;
- no accepted rule for when a bounded search is sufficient to support a
  capability decision;
- no accepted distinction in repository markers among software verification,
  numerical verification, physical validation, pedagogical validation, and UQ.

Absence of recorded evidence does not establish that a capability is invalid,
incorrect, or without pedagogical value. It limits the claims that may be made
from this inspection snapshot.

## 12. Unresolved human decisions

### 12.1 Snapshot and repository decisions

1. Accept, correct, or reject the provisional path of this snapshot.
2. Accept, correct, or reject the inspection findings and limitations.
3. Decide whether future current project status should be composed from
   capability-level records.
4. Decide whether `README.md` should contain only a derived project summary.
5. Define the accepted authority and update process for any derived summary.
6. Decide whether a complete deterministic per-file inventory is required.
7. If required later, separately approve its schema, path, generator, and
   retention policy.
8. Define future boundaries among source, tests, canonical notebooks,
   exploratory notebooks, examples, generated evidence, and historical
   material.

Acceptance of this snapshot would not resolve decisions 3–8.

### 12.2 Capability decisions

1. Define a capability lifecycle without retroactively applying unstated
   maturity states to this snapshot.
2. Select one bounded capability for subsequent contract planning.
3. Define that capability’s pedagogical purpose and intended learners.
4. Select its physical model and assumptions.
5. Select its mathematical formulation.
6. Select its numerical representation and conventions.
7. Select its public API.
8. Resolve only the competing representations needed by that capability.
9. Define software-verification requirements.
10. Define numerical-verification references, invariants, and tolerances.
11. Define physical-validation requirements for a declared intended use.
12. Define pedagogical-validation requirements.
13. Define UQ requirements or an explicit not-applicable rationale.
14. Decide who may accept each evidence class.
15. Decide when implementation, promotion, replacement, or deprecation may
    occur.

## 13. Limitations of the inspection

- Pytest was unavailable and was not installed.
- Tests were not collected or executed.
- Notebooks were not executed.
- Saved notebook outputs may be stale or environment-dependent.
- Static parsing and import probing do not establish behavioral correctness.
- Import success does not establish contract satisfaction or numerical
  correctness.
- Source declarations do not establish public API membership.
- File and directory names do not establish intended authority.
- Documentation was sampled for structural and API relationships, not reviewed
  as a complete physics or pedagogy specification.
- Notebook code was inventoried structurally and textually, not exhaustively
  reviewed cell by cell.
- Bounded searches may miss indirect, dynamically generated, externally stored,
  or differently named relationships.
- The working tree contained pre-existing human changes.
- Modified magnetism notebooks were not compared with their committed versions
  or interpreted.
- Capability groups are human-readable aggregations, not a complete per-file
  manifest.
- No conclusion in this document constitutes lifecycle classification.
- This snapshot cannot remain current automatically and must not be maintained
  as a competing current-state authority.

## 14. Relationship to future capability and lifecycle records

If accepted, this baseline may be cited by later planning or policy artifacts
as evidence of what the stated inspection observed.

A future lifecycle policy would need to define, in a separate human-reviewed
artifact:

- lifecycle states;
- entry and exit criteria;
- required capability contracts;
- protected decisions;
- required evidence classes;
- promotion and deprecation authority;
- treatment of exploratory and historical artifacts;
- relationship between capability records and a README summary.

This baseline does not design those rules.

Later accepted capability records or lifecycle artifacts may supersede this
snapshot’s observations. They should cite the newer evidence and applicable
authority rather than editing this historical snapshot to make it appear
current.

If a later repository condition differs from this snapshot, the difference is
new evidence. It does not make the historical observation false, provided the
recorded revision, working-tree state, methods, and limitations remain attached
to it.

## 15. Explicit non-decisions

This document explicitly records that:

1. **No implementation was selected as canonical.**
2. **No capability maturity was assigned.**
3. Every capability classification remains **pending**.
4. No notebook was selected as canonical.
5. No source declaration was declared supported or public.
6. No class or function definition was assumed to parse, import, execute, or
   satisfy its intended contract solely because it exists.
7. No test was declared passing.
8. No numerical result was declared verified.
9. No physical model was declared validated.
10. No pedagogical objective was declared achieved.
11. No uncertainty analysis was declared complete.
12. No repository-wide lifecycle state was assigned.
13. No current project status was established by this snapshot.
14. No active task or Pi runtime state was established.
15. No control plane was created or implemented.
