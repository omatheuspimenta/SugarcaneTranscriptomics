"""Tests for the progress-reporting contract.

The progress callback is easy to break silently: a refactor drops the
``on_progress`` argument, or a loop stops advancing, and the only symptom is a
bar that sits at 0 % for twenty minutes. These tests assert the contract
directly — a total is announced, and the advances add up to it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import rgas_prediction
from rga import evidence as evidence_mod
from rga import parsers, rules
from rga.progress import null_progress
from conftest import write_deepcoil, write_interproscan

IPS_ROWS = [
    ("P1", 400, "Pfam", "PF00931", "NB-ARC domain", 20, 200, "1e-40", "IPR002182"),
    ("P1", 400, "Pfam", "PF13855", "Leucine rich repeat", 250, 300, "1e-9", "IPR001611"),
    ("P2", 300, "Pfam", "PF00069", "Protein kinase domain", 40, 280, "1e-30", "IPR000719"),
]


class Recorder:
    """A :class:`~rga.progress.ProgressCallback` that remembers what it was told."""

    def __init__(self) -> None:
        """Start with no total and no advances."""
        self.total: int | None = None
        self.advances: list[int] = []

    def __call__(self, advance: int = 0, total: int | None = None) -> None:
        """Record one progress report."""
        if total is not None:
            self.total = total
        if advance:
            self.advances.append(advance)

    @property
    def completed(self) -> int:
        """Total units reported as done."""
        return sum(self.advances)


@pytest.fixture
def ips_tsv(tmp_path: Path) -> Path:
    """A tiny InterProScan TSV."""
    path = tmp_path / "ips.tsv"
    write_interproscan(path, IPS_ROWS)
    return path


def test_interproscan_reports_bytes(cfg, ips_tsv: Path) -> None:
    """The InterProScan bar is measured in bytes and reaches the file size."""
    recorder = Recorder()
    parsers.parse_interproscan(ips_tsv, cfg, recorder)
    assert recorder.total == ips_tsv.stat().st_size
    assert recorder.completed == ips_tsv.stat().st_size


def test_interproscan_result_is_identical_with_and_without_progress(cfg, ips_tsv) -> None:
    """Wrapping the file handle must not change what is parsed."""
    with_progress = parsers.parse_interproscan(ips_tsv, cfg, Recorder())
    without = parsers.parse_interproscan(ips_tsv, cfg)
    assert with_progress.n_rows == without.n_rows
    assert with_progress.protein_ids == without.protein_ids
    assert with_progress.hits.equals(without.hits)


def test_deepcoil_reports_one_unit_per_source(cfg, tmp_path: Path) -> None:
    """The DeepCoil2 bar counts sources, and finishes them all."""
    for name in ("A", "B"):
        directory = tmp_path / f"part_{name}"
        directory.mkdir()
        write_deepcoil(directory, f"PROT{name}1p", [0.9] * 30)
    recorder = Recorder()
    parsers.parse_deepcoil(tmp_path, cfg, workers=1, on_progress=recorder)
    assert recorder.total == 2
    assert recorder.completed == 2


def test_evidence_and_classification_report_every_protein(cfg, ips_tsv: Path) -> None:
    """Both per-protein loops announce the protein count and reach it."""
    ips = parsers.parse_interproscan(ips_tsv, cfg)
    options = rgas_prediction.resolve_options(
        rgas_prediction.build_parser().parse_args([]), cfg
    )
    protein_ids = sorted(ips.protein_ids)

    ev_recorder = Recorder()
    ev = evidence_mod.build_evidence(
        cfg, protein_ids, ips, None, None, None, None, None, {}, options,
        on_progress=ev_recorder,
    )
    assert ev_recorder.total == len(protein_ids)
    assert ev_recorder.completed == len(protein_ids)

    cls_recorder = Recorder()
    predictions = rules.classify_proteome(cfg, ev, options, on_progress=cls_recorder)
    assert cls_recorder.total == len(protein_ids)
    assert cls_recorder.completed == len(protein_ids)
    assert len(predictions) == len(protein_ids)


def test_stages_run_without_a_callback(cfg, ips_tsv: Path) -> None:
    """``null_progress`` is the default, so no stage requires a console."""
    ips = parsers.parse_interproscan(ips_tsv, cfg)
    assert ips.n_rows == len(IPS_ROWS)
    assert null_progress(1, total=2) is None


def test_tracked_yields_a_working_callback() -> None:
    """The CLI helper binds a real rich task and survives a total set late."""
    with rgas_prediction.tracked("test stage") as report:
        report(0, total=4)
        report(2)
        report(2)
