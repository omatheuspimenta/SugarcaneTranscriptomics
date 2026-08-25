"""Guards on the documentation itself.

These tests exist because documentation rots in specific, detectable ways:
source files stop parsing after a hand edit, a function is added without a
docstring, a reference is added to the code but not to the README (or the other
way round), or a rule is added to the configuration without being described.
Each of those is cheap to check and expensive to notice by eye.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

from rga import report
from rga.config import load_config

REPO = Path(__file__).resolve().parents[2]
CODE = REPO / "code"
DOCS = REPO / "docs" / "rga"
CONFIG = CODE / "config" / "rga_config.yaml"

#: The preserved original script is kept verbatim as a historical record and is
#: therefore exempt from the documentation rules applied to current code.
LEGACY = CODE / "rgas_prediction_legacy.py"


def _sources() -> list[Path]:
    """Every Python file that is part of the current pipeline."""
    return sorted(
        p
        for p in CODE.rglob("*.py")
        if p != LEGACY and ".venv" not in p.parts and "__pycache__" not in p.parts
    )


@pytest.mark.parametrize("path", _sources(), ids=lambda p: p.name)
def test_source_file_parses(path: Path) -> None:
    """Every source file is valid Python.

    A hand edit once stripped the ``#`` from ten section-banner comments, which
    left bare ``---------`` lines at module level and made the entry point
    unrunnable. Nothing else in the suite catches that, because an unimportable
    entry point cannot be imported to be tested.
    """
    ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


@pytest.mark.parametrize("path", _sources(), ids=lambda p: p.name)
def test_everything_has_a_docstring(path: Path) -> None:
    """Every module, class and function carries a docstring.

    Nested helper closures inside tests are exempt: they are local to one
    assertion and naming them well is enough.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    assert ast.get_docstring(tree), f"{path.name} has no module docstring"

    is_test = path.parent.name == "tests"
    undocumented = [
        f"{node.name}:{node.lineno}"
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
        and not ast.get_docstring(node)
        and not (is_test and node.col_offset > 0)
    ]
    assert not undocumented, f"{path.name}: undocumented {undocumented}"


def test_every_code_reference_is_in_the_readme() -> None:
    """The bibliography in ``report.py`` and the one in the README agree.

    The reports cite by DOI, the README by DOI link; a reference added to one
    and not the other means a published report cites something the methodology
    document does not explain.
    """
    readme = (DOCS / "README.md").read_text(encoding="utf-8")
    missing = [doi for _, doi in report.REFERENCES if doi not in readme]
    assert not missing, f"cited in code but absent from README: {missing}"


def test_every_readme_reference_is_in_the_code() -> None:
    """No DOI is listed in the README that the generated reports omit."""
    readme = (DOCS / "README.md").read_text(encoding="utf-8")
    listed = set(re.findall(r"doi:\[(10\.[^\]]+)\]", readme))
    cited = {doi for _, doi in report.REFERENCES}
    assert not listed - cited, f"in README but not cited by the reports: {listed - cited}"


def test_every_rule_is_documented() -> None:
    """Each configured rule id appears in the README rule table."""
    cfg = load_config(CONFIG)
    readme = (DOCS / "README.md").read_text(encoding="utf-8")
    missing = [rule.id for rule in cfg.rules if rule.id not in readme]
    assert not missing, f"rules absent from the README: {missing}"


def test_every_rule_has_a_description() -> None:
    """Each rule carries a human-readable description, reused in the report."""
    cfg = load_config(CONFIG)
    undescribed = [rule.id for rule in cfg.rules if not rule.description.strip()]
    assert not undescribed, f"rules without a description: {undescribed}"


def test_architecture_document_exists_and_is_linked() -> None:
    """``ARCHITECTURE.md`` exists and both READMEs point at it."""
    assert (DOCS / "ARCHITECTURE.md").is_file()
    assert "ARCHITECTURE.md" in (DOCS / "README.md").read_text(encoding="utf-8")
    assert "docs/rga/ARCHITECTURE.md" in (REPO / "README.md").read_text(encoding="utf-8")
