from __future__ import annotations

import re
from pathlib import Path


def test_publication_repository_files() -> None:
    root = Path(__file__).resolve().parents[1]
    required = {
        ".github/workflows/tests.yml",
        "CITATION.cff",
        "LICENSE",
        "README.md",
        "docs/FIGURE_MAP.md",
        "docs/METHODS.md",
        "docs/REPRODUCIBILITY.md",
        "docs/SOFTWARE_ENVIRONMENT.md",
        "environment.yml",
        "figures/figure_10_matched_constructions.png",
        "figures/figure_9_appendix_b.png",
        "pyproject.toml",
        "requirements-lock.txt",
        "requirements.txt",
        "results/appendix_b1_crossover.csv",
        "results/appendix_b2_estimator_summary.csv",
        "results/appendix_b3_pair_correlation.csv",
        "scripts/generate_appendix_b.py",
        "scripts/run_all.ps1",
        "src/correlated_abiogenesis/__init__.py",
        "src/correlated_abiogenesis/appendix_b.py",
    }
    assert all((root / relative).is_file() for relative in required)

    windows_users = r"[A-Za-z]:" + r"[\\/]+" + "Users" + r"[\\/]+" + r"[^\\/]+"
    posix_users = re.escape("/" + "Users" + "/") + r"[^/]+"
    machine_path = re.compile(f"(?:{windows_users}|{posix_users})")
    excluded_top_level = {"legacy", "local", "logs", "outputs"}
    for path in root.rglob("*"):
        relative_parts = path.relative_to(root).parts
        if (
            not path.is_file()
            or any(part.startswith(".") for part in relative_parts)
            or (relative_parts and relative_parts[0] in excluded_top_level)
        ):
            continue
        if path.suffix.lower() not in {".cff", ".csv", ".md", ".ps1", ".py", ".toml", ".txt", ".yml"}:
            continue
        text = path.read_text(encoding="utf-8")
        assert not machine_path.search(text), f"Machine-specific path in {path}"
