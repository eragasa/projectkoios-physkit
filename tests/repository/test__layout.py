"""Repository-layout conformance for the Project Koios component package."""

from __future__ import annotations

import json
from pathlib import Path


def test__repository_layout__uses_projectkoios_namespace() -> None:
    """Require canonical namespace packaging and legacy-source path mirrors.

    An optional ``src/physkit`` tree may remain for workspace artifacts or
    transitional sources, but every Python file there must have the same relative
    path in the canonical source tree. Tests must live only in the mirrored
    ``tests/projectkoios/physkit`` namespace.
    """
    root: Path = Path(__file__).resolve().parents[2]
    pyproject_text: str = (root / "pyproject.toml").read_text(encoding="utf-8")

    assert 'name = "projectkoios-physkit"' in pyproject_text
    assert 'package-dir = { "" = "src/python" }' in pyproject_text
    assert 'where = ["src/python"]' in pyproject_text
    assert "namespaces = true" in pyproject_text
    package_root = root / "src/python/projectkoios/physkit"
    assert (package_root / "__init__.py").is_file()
    assert not tuple(package_root.rglob("*.ipynb"))
    assert not (root / "src/python/projectkoios/__init__.py").exists()

    legacy_package_root = root / "src/physkit"
    legacy_source_paths = tuple(
        sorted(
            path
            for path in legacy_package_root.rglob("*")
            if path.is_file() and path.suffix in {".py", ".pyi"}
        )
    )
    missing_mirrors = tuple(
        str(path.relative_to(legacy_package_root))
        for path in legacy_source_paths
        if not (package_root / path.relative_to(legacy_package_root)).is_file()
    )
    assert not missing_mirrors, (
        "Python sources under src/physkit must be mirrored under "
        f"src/python/projectkoios/physkit: {missing_mirrors}"
    )

    test_root = root / "tests/projectkoios/physkit"
    assert test_root.is_dir()
    assert not (root / "tests/physkit").exists()


def test__repository_notebook_layout__preserves_pedagogical_laboratories() -> None:
    """Require named laboratories without mirroring package/test namespaces.

    Notebook paths are organized for pedagogical ownership. They may coexist in
    distinct teaching trees and are not required to mirror Python source paths.
    """
    root = Path(__file__).resolve().parents[2]

    assert (root / "notebooks/qm/piab1d").is_dir()
    assert (root / "notebooks/qm/piab2d").is_dir()
    assert (root / "notebooks/qm/piab3d").is_dir()
    assert (root / "notebooks/qm/piab1d/piab1d.ipynb").is_file()
    assert (root / "notebooks/qm/piab2d/tise/analytical/piab2d.ipynb").is_file()
    assert (root / "notebooks/qm/piab2d/tdse/numerical").is_dir()
    assert (root / "notebooks/qm/piab3d/piab3d.ipynb").is_file()


def test__repository_notebooks__are_valid_json_documents() -> None:
    """Require every repository notebook to contain a valid JSON document."""
    root = Path(__file__).resolve().parents[2]
    notebook_paths = tuple(sorted((root / "notebooks").rglob("*.ipynb")))

    assert notebook_paths
    for notebook_path in notebook_paths:
        json.loads(notebook_path.read_text(encoding="utf-8"))
