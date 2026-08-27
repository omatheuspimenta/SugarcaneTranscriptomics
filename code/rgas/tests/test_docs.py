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

REPO = Path(__file__).resolve().parents[3]
CODE = REPO / "code" / "rgas"
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


def test_the_rule_table_lists_exactly_the_configured_rules() -> None:
    """README section 6.2 must name every rule, at its priority, and no other.

    ``test_every_rule_is_documented`` only asks whether each rule id appears
    somewhere in the file, so it stays green when a rule is *removed* from the
    configuration and its row is left behind in the table. That happened once
    (the ``other-RLK`` row survived the switch to the Rody et al. RLK rule),
    and the table then documented a class the pipeline could not produce.
    """
    cfg = load_config(CONFIG)
    readme = (DOCS / "README.md").read_text(encoding="utf-8")
    table = readme.split("### 6.2 The rule table", 1)[1].split("### 6.3", 1)[0]
    documented = {
        m.group(2): int(m.group(1))
        for m in re.finditer(r"^\| (\d+) \| `([^`]+)` \|", table, re.M)
    }
    configured = {rule.id: rule.priority for rule in cfg.rules}
    assert documented == configured, (
        f"documented but not configured: "
        f"{sorted(set(documented) - set(configured))}; "
        f"configured but not in the table: "
        f"{sorted(set(configured) - set(documented))}; "
        f"priority mismatches: "
        f"{ {k: (documented[k], configured[k]) for k in documented.keys() & configured.keys() if documented[k] != configured[k]} }"
    )


def test_worked_example_counts_match_the_reference_run() -> None:
    """The per-class counts in README section 6.3 must match the shipped run.

    Skipped when the reference run is not present in the checkout. These counts
    are quoted in the text and were hand-maintained; they went stale the first
    time a rule was disabled, reporting an ``other-RLK`` class that no longer
    existed and an ``Other`` count 5,527 proteins short.
    """
    import csv

    counts = REPO / "results" / "rgas" / "SaccharumR570" / "rga_summary_counts.tsv"
    if not counts.is_file():
        pytest.skip("reference run not present in this checkout")

    with counts.open(encoding="utf-8", newline="") as handle:
        actual = {
            row["rga_subclass"]: int(row["n_proteins"])
            for row in csv.DictReader(handle, delimiter="\t")
            if row["level"] == "subclass"
        }

    readme = (DOCS / "README.md").read_text(encoding="utf-8")
    section = readme.split("### 6.3 One worked example per class", 1)[1]
    # bound at the next heading, so a later section's tables can never be
    # scraped into this one
    section = re.split(r"^#{2,3} ", section, maxsplit=1, flags=re.M)[0]
    documented = {}
    for line in section.splitlines():
        row = re.match(r"\|\s*`([^`]+)`\s*\|\s*([\d,]+|—)\s*\|", line)
        if row:
            value = row.group(2)
            documented[row.group(1)] = 0 if value == "—" else int(value.replace(",", ""))

    assert documented, "no class-count table found in README section 6.3"
    mismatched = {
        cls: (n, actual.get(cls, 0))
        for cls, n in documented.items()
        if n != actual.get(cls, 0)
    }
    assert not mismatched, f"README (documented, actual): {mismatched}"
    assert not set(actual) - set(documented), (
        f"classes in the run but not in the table: {sorted(set(actual) - set(documented))}"
    )


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


def test_every_prediction_column_is_documented():
    """The data dictionary must describe every column, and no column it lacks.

    The dictionary is what a reader consults instead of the source, so a column
    added without a row here is a column nobody can interpret -- and a row left
    behind after a rename documents something that no longer exists.
    """
    import re

    from rga.rules import PREDICTION_COLUMNS

    readme = (DOCS / "README.md").read_text(encoding="utf-8")
    start = readme.index("Data dictionary — `rga_predictions.tsv`")
    end = readme.index("`rga_domain_evidence_long.tsv`", start)
    documented = set(re.findall(r"^\| `([a-z_0-9]+)` \|", readme[start:end], re.M))
    assert documented == set(PREDICTION_COLUMNS), (
        f"undocumented: {sorted(set(PREDICTION_COLUMNS) - documented)}; "
        f"documented but absent: {sorted(documented - set(PREDICTION_COLUMNS))}"
    )


def test_readme_example_command_matches_the_recorded_reference_run():
    """The documented example and the reference run must be the same command.

    Skipped when the reference run is not present in the checkout (the results
    directory is large and may be absent); when it is present, the two must
    agree token for token, so the README cannot document a command that is not
    the one that produced the numbers beside it.
    """
    import json
    import shlex

    metadata = REPO / "results" / "rgas" / "SaccharumR570" / "run_metadata.json"
    if not metadata.is_file():
        pytest.skip("reference run not present in this checkout")
    recorded = shlex.split(json.loads(metadata.read_text())["command"])

    readme = (DOCS / "README.md").read_text(encoding="utf-8")
    # the last fenced block of section 3.1 is the run command; the first is the
    # one-time `uv sync`
    section = readme.split("### 3.1 The shortest path", 1)[1].split("### 3.2", 1)[0]
    block = section.split("```")[-2]
    documented = shlex.split(block.replace("\\\n", " ").replace("bash\n", "", 1))
    # drop the interpreter prefix the documented form carries
    assert documented[:3] == ["uv", "run", "python"]
    assert documented[3:] == recorded, (
        f"README example: {documented[3:]}\nreference run: {recorded}"
    )
