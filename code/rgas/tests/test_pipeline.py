"""End-to-end tests on a synthetic mini-proteome covering the documented edge cases."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

import rgas_prediction
from conftest import (
    write_deepcoil,
    write_deeploc,
    write_deeptmhmm,
    write_interproscan,
    write_phobius,
    write_signalp,
)

# Nine synthetic proteins, each exercising a different part of the pipeline.
IPS_ROWS = [
    # CNL with a canonical N-terminal coiled coil, LRR seen by three databases
    ("CNL1", 900, "Pfam", "PF00931", "NB-ARC domain", 178, 455, "1e-40", "IPR002182"),
    (
        "CNL1",
        900,
        "Pfam",
        "PF13855",
        "Leucine rich repeat",
        512,
        560,
        "1e-9",
        "IPR001611",
    ),
    ("CNL1", 900, "SMART", "SM00369", "LRR_typ_2", 515, 558, "1e-8", "IPR003591"),
    ("CNL1", 900, "Gene3D", "G3DSA:3.80.10.10", "-", 510, 562, "1e-7", "IPR032675"),
    ("CNL1", 900, "Coils", "Coil", "Coil", 21, 48, "-", "-"),
    # CNL-like protein whose only coiled coil sits C-terminal to the NB-ARC
    ("CNL2", 900, "Pfam", "PF00931", "NB-ARC domain", 50, 300, "1e-40", "IPR002182"),
    (
        "CNL2",
        900,
        "Pfam",
        "PF00560",
        "Leucine Rich Repeat",
        400,
        450,
        "1e-9",
        "IPR001611",
    ),
    # TNL
    ("TNL1", 800, "Pfam", "PF01582", "TIR domain", 10, 140, "1e-20", "IPR000157"),
    ("TNL1", 800, "Pfam", "PF00931", "NB-ARC domain", 180, 450, "1e-40", "IPR002182"),
    (
        "TNL1",
        800,
        "Pfam",
        "PF00560",
        "Leucine Rich Repeat",
        500,
        550,
        "1e-9",
        "IPR001611",
    ),
    # RNL (RPW8 through the PROSITE profile, the only domain-level evidence)
    (
        "RNL1",
        700,
        "ProSiteProfiles",
        "PS51153",
        "RPW8 domain profile.",
        5,
        150,
        "20.0",
        "IPR008808",
    ),
    ("RNL1", 700, "Pfam", "PF00931", "NB-ARC domain", 200, 460, "1e-40", "IPR002182"),
    (
        "RNL1",
        700,
        "Pfam",
        "PF00560",
        "Leucine Rich Repeat",
        520,
        570,
        "1e-9",
        "IPR001611",
    ),
    # LRR-RLK with a signal peptide and one internal helix
    (
        "RLK1",
        1000,
        "Pfam",
        "PF00560",
        "Leucine Rich Repeat",
        60,
        110,
        "1e-9",
        "IPR001611",
    ),
    (
        "RLK1",
        1000,
        "Pfam",
        "PF00069",
        "Protein kinase domain",
        600,
        850,
        "1e-50",
        "IPR000719",
    ),
    # LRR-RLP: same but no kinase
    (
        "RLP1",
        800,
        "Pfam",
        "PF00560",
        "Leucine Rich Repeat",
        60,
        110,
        "1e-9",
        "IPR001611",
    ),
    # TM-CC: coiled coil overlapping the transmembrane helix (the artefact case)
    (
        "TMCC1",
        300,
        "MobiDBLite",
        "mobidb-lite",
        "consensus disorder prediction",
        1,
        40,
        "-",
        "-",
    ),
    # Protein with only MobiDBLite -> Non-RGA
    (
        "NONE1",
        200,
        "MobiDBLite",
        "mobidb-lite",
        "consensus disorder prediction",
        5,
        60,
        "-",
        "-",
    ),
]

PHOBIUS_ROWS = [
    ("CNL1", 0, "0", "o"),
    ("CNL2", 0, "0", "o"),
    ("TNL1", 0, "0", "o"),
    ("RNL1", 0, "0", "o"),
    ("RLK1", 1, "Y", "n4-15c25/26o500-520i"),
    ("RLP1", 2, "Y", "n4-15c25/26o600-620i700-720o"),
    ("TMCC1", 1, "0", "o40-60i"),
    ("NONE1", 0, "0", "o"),
    ("ONLYIPS", 0, "0", "o"),
]

DEEPTMHMM_BLOCKS = [
    ("CNL1", 900, [("inside", 1, 900)]),
    ("CNL2", 900, [("inside", 1, 900)]),
    ("TNL1", 800, [("inside", 1, 800)]),
    ("RNL1", 700, [("inside", 1, 700)]),
    # DeepTMHMM mis-calls the signal peptide as a helix; the SP filter must drop it
    (
        "RLK1",
        1000,
        [
            ("TMhelix", 3, 24),
            ("outside", 25, 499),
            ("TMhelix", 500, 520),
            ("inside", 521, 1000),
        ],
    ),
    # conflict: Phobius sees two helices, DeepTMHMM none
    ("RLP1", 800, [("signal", 1, 25), ("outside", 26, 800)]),
    ("TMCC1", 300, [("TMhelix", 40, 60), ("inside", 61, 300)]),
    ("NONE1", 200, [("inside", 1, 200)]),
    ("ONLYIPS", 100, [("inside", 1, 100)]),
]

SIGNALP_ROWS = [
    ("CNL1", "OTHER", 0.0, ""),
    ("CNL2", "OTHER", 0.0, ""),
    ("TNL1", "OTHER", 0.0, ""),
    ("RNL1", "OTHER", 0.0, ""),
    ("RLK1", "SP", 0.99, "CS pos: 25-26. Pr: 0.95"),
    ("RLP1", "SP", 0.99, "CS pos: 25-26. Pr: 0.95"),
    ("TMCC1", "OTHER", 0.0, ""),
    ("NONE1", "OTHER", 0.0, ""),
    ("ONLYIPS", "OTHER", 0.0, ""),
]

# DeepLoc deliberately omits ONLYIPS and CNL2 to exercise partial coverage.
DEEPLOC_ROWS = [
    ("CNL1", "Cytoplasm", {"Cytoplasm": 0.87}),
    ("TNL1", "Extracellular", {"Extracellular": 0.8}),  # inconsistent for an NLR
    ("RNL1", "Cytoplasm", {"Cytoplasm": 0.9}),
    ("RLK1", "Cell membrane", {"Cell membrane": 0.93}),
    ("RLP1", "Cell membrane", {"Cell membrane": 0.88}),
    ("TMCC1", "Cell membrane", {"Cell membrane": 0.7}),
    ("NONE1", "Nucleus", {"Nucleus": 0.95}),
]

# DeepCoil scores: 0 outside a segment, a constant plateau inside one.
DEEPCOIL = {
    "CNL1": [0.0] * 20 + [0.8] * 30 + [0.0] * 850,  # N-terminal, also seen by Coils
    "CNL2": [0.0] * 500 + [0.8] * 30 + [0.0] * 370,  # C-terminal to the NB-ARC
    "TNL1": [0.0] * 800,
    "RNL1": [0.0] * 700,
    "RLK1": [0.0] * 1000,
    "RLP1": [0.0] * 800,
    "TMCC1": [0.0] * 39 + [0.9] * 25 + [0.0] * 236,  # overlaps the TM helix
    "NONE1": [0.0] * 200,
    "ONLYIPS": [0.0] * 100,
}


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory) -> Path:
    """Write a complete synthetic input tree in the layout the pipeline expects."""
    root = tmp_path_factory.mktemp("inputs")
    (root / "InterProScan").mkdir()
    (root / "phobius").mkdir()
    (root / "DeepTMHMM").mkdir()
    (root / "SignalP6").mkdir()
    (root / "DeepLoc2").mkdir()
    write_interproscan(root / "InterProScan" / "ips.tsv", IPS_ROWS)
    write_phobius(root / "phobius" / "r.phobius", PHOBIUS_ROWS)
    write_deeptmhmm(root / "DeepTMHMM" / "TMRs.gff3", DEEPTMHMM_BLOCKS)
    write_signalp(root / "SignalP6" / "prediction_results.txt", SIGNALP_ROWS)
    write_deeploc(root / "DeepLoc2" / "results_test.csv", DEEPLOC_ROWS)
    for protein_id, scores in DEEPCOIL.items():
        write_deepcoil(root / "DeepCoil" / "part_001", protein_id, scores)
    return root


def _run(inputs: Path, outdir: Path, extra: list[str] | None = None) -> pd.DataFrame:
    """Run the pipeline and return the prediction table."""
    argv = [
        "--input-dir",
        str(inputs),
        "--outdir",
        str(outdir),
        "--organism-name",
        "Synthetic",
        "--log-level",
        "WARNING",
        *(extra or []),
    ]
    assert rgas_prediction.main(argv) == 0
    # keep_default_na=False: the pipeline writes the literal string "NA" for
    # missing values, and reading it back as a float NaN would hide that.
    return pd.read_csv(
        outdir / "rga_predictions.tsv", sep="\t", dtype=str, keep_default_na=False
    )


@pytest.fixture(scope="module")
def predictions(synthetic, tmp_path_factory) -> pd.DataFrame:
    """The prediction table of the default run."""
    return _run(synthetic, tmp_path_factory.mktemp("out"))


def test_every_protein_appears_exactly_once(predictions):
    """Including the protein no tool but InterProScan ever saw."""
    assert len(predictions) == 9
    assert predictions["protein_id"].is_unique
    assert "ONLYIPS" in set(predictions["protein_id"])


def test_expected_classes(predictions):
    """Each synthetic protein lands in its intended class."""
    calls = dict(zip(predictions["protein_id"], predictions["rga_subclass"]))
    assert calls["CNL1"] == "CNL"
    assert calls["TNL1"] == "TNL"
    assert calls["RNL1"] == "RNL"
    assert calls["RLK1"] == "LRR-RLK"
    assert calls["RLP1"] == "LRR-RLP"
    assert calls["TMCC1"] == "TM-CC"
    assert calls["NONE1"] == "NA"
    assert calls["ONLYIPS"] == "NA"


def test_overlapping_lrr_hits_count_as_one_copy(predictions):
    """Three databases calling the same LRR region give n_lrr == 1."""
    row = predictions.set_index("protein_id").loc["CNL1"]
    assert int(row["n_lrr"]) == 1


def test_tm_helix_inside_the_signal_peptide_is_discarded(predictions):
    """DeepTMHMM's helix at 3-24 lies inside the 1-25 signal peptide."""
    row = predictions.set_index("protein_id").loc["RLK1"]
    assert int(row["n_tm_deeptmhmm"]) == 1
    assert int(row["n_tm_deeptmhmm_raw"]) == 2
    assert int(row["n_tm_dropped_in_sp"]) == 1


def test_conflicting_tm_calls_resolve_by_union(predictions):
    """Phobius sees helices where DeepTMHMM sees none; the union policy keeps them."""
    row = predictions.set_index("protein_id").loc["RLP1"]
    assert int(row["n_tm_phobius"]) == 2
    assert int(row["n_tm_deeptmhmm"]) == 0
    assert row["tm_consensus"] == "True"


def test_cc_called_by_deepcoil_but_not_by_coils(predictions):
    """CNL2's coiled coil is DeepCoil2-only, so the call is demoted."""
    row = predictions.set_index("protein_id").loc["CNL2"]
    assert row["cc_deepcoil"] == "True"
    assert row["cc_coils"] == "False"
    assert row["cc_source"] == "deepcoil_only"
    assert row["confidence"] in {"medium", "low"}


def test_cc_supported_by_both_tools_is_high_confidence(predictions):
    """CNL1 is called by DeepCoil2 and corroborated by InterProScan Coils."""
    row = predictions.set_index("protein_id").loc["CNL1"]
    assert row["cc_source"] == "deepcoil+coils"
    assert row["cc_is_n_terminal"] == "True"
    assert row["confidence"] == "high"


def test_cc_c_terminal_to_nbarc_is_flagged_not_reclassified(predictions):
    """A C-terminal coiled coil keeps the CNL call but raises a warning."""
    row = predictions.set_index("protein_id").loc["CNL2"]
    assert row["rga_subclass"] == "CNL"
    assert row["cc_is_n_terminal"] == "False"
    assert "C-terminal" in str(row["warnings"])


def test_cc_overlapping_a_tm_helix_is_flagged(predictions):
    """The TM-CC class is the one most exposed to the CC/TM artefact."""
    row = predictions.set_index("protein_id").loc["TMCC1"]
    assert row["cc_tm_ambiguous"] == "True"
    assert "overlaps a predicted TM helix" in str(row["warnings"])
    assert row["confidence"] == "low"


def test_nlr_predicted_extracellular_is_flagged(predictions):
    """DeepLoc never changes the class, but the inconsistency is reported."""
    row = predictions.set_index("protein_id").loc["TNL1"]
    assert row["rga_subclass"] == "TNL"
    assert "inconsistent" in str(row["warnings"])


def test_protein_missing_from_deeploc_still_appears(predictions):
    """A protein absent from DeepLoc keeps NA localisation and is not dropped."""
    row = predictions.set_index("protein_id").loc["ONLYIPS"]
    assert row["predicted_localization"] == "NA"


def test_reason_field_cites_provenance(predictions):
    """The justification quotes the accession and coordinates behind the call."""
    reason = predictions.set_index("protein_id").loc["CNL1", "reason"]
    assert "Rule CNL (priority 1)" in reason
    assert "PF00931 @ 178-455" in reason
    assert "Confidence: high" in reason


def test_domain_architecture_is_n_to_c(predictions):
    """Architecture strings are ordered by coordinate."""
    row = predictions.set_index("protein_id").loc["CNL1"]
    assert row["domain_architecture"] == "CC-NB-ARC-LRR"


def test_all_outputs_are_written(synthetic, tmp_path):
    """Every promised output file exists after a run."""
    outdir = tmp_path / "out"
    _run(synthetic, outdir)
    for name in (
        "rga_predictions.tsv",
        "rga_predictions_rga_only.tsv",
        "rga_domain_evidence_long.tsv",
        "rga_summary_counts.tsv",
        "unmatched_ids_report.tsv",
        "rga_predictions_by_locus.tsv",
        "report.html",
        "report.md",
        "run_metadata.json",
        "logs/run.log",
    ):
        assert (outdir / name).is_file(), name


def test_counts_agree_across_outputs(synthetic, tmp_path):
    """The TSV, the JSON and the report must never disagree on a count."""
    outdir = tmp_path / "out"
    predictions = _run(synthetic, outdir)
    counts = pd.read_csv(outdir / "rga_summary_counts.tsv", sep="\t")
    metadata = json.loads((outdir / "run_metadata.json").read_text())
    subclass_total = int(counts.loc[counts["level"] == "subclass", "n_proteins"].sum())
    assert subclass_total == len(predictions) == metadata["counts"]["n_proteins"]
    report = (outdir / "report.md").read_text()
    assert f"{len(predictions):,} proteins were examined" in report


def test_html_report_is_self_contained(synthetic, tmp_path):
    """No external stylesheet, script or image may be referenced."""
    outdir = tmp_path / "out"
    _run(synthetic, outdir)
    html = (outdir / "report.html").read_text()
    assert "<script" not in html
    assert "http://" not in html
    assert "<img" not in html
    assert "cdn" not in html.lower()


def test_pipeline_runs_without_deepcoil(synthetic, tmp_path):
    """Missing DeepCoil2 falls back to Coils and demotes every CC-dependent call."""
    outdir = tmp_path / "out_nocc"
    predictions = _run(
        synthetic,
        outdir,
        ["--deepcoil", str(tmp_path / "missing")]
        if (tmp_path / "missing").exists()
        else [],
    )
    # Re-run with the DeepCoil directory renamed away.
    stripped = tmp_path / "inputs_nocc"
    stripped.mkdir()
    for child in synthetic.iterdir():
        if child.name != "DeepCoil":
            (stripped / child.name).symlink_to(child)
    predictions = _run(stripped, tmp_path / "out_nocc2")
    row = predictions.set_index("protein_id").loc["CNL1"]
    assert row["cc_source"] == "coils_only"
    assert row["confidence"] == "low"
    assert "unavailable" in str(row["warnings"])


def test_pipeline_runs_without_deeploc(synthetic, tmp_path):
    """The pipeline degrades gracefully when an optional tool is absent."""
    stripped = tmp_path / "inputs_noloc"
    stripped.mkdir()
    for child in synthetic.iterdir():
        if child.name != "DeepLoc2":
            (stripped / child.name).symlink_to(child)
    predictions = _run(stripped, tmp_path / "out_noloc")
    assert (predictions["predicted_localization"] == "NA").all()
    assert "deeploc" in str(predictions.iloc[0]["warnings"])


def test_duplicated_protein_ids_do_not_duplicate_output_rows(synthetic, tmp_path):
    """A protein reported twice by a tool still yields exactly one output row."""
    duplicated = tmp_path / "inputs_dup"
    duplicated.mkdir()
    for child in synthetic.iterdir():
        if child.name != "phobius":
            (duplicated / child.name).symlink_to(child)
    (duplicated / "phobius").mkdir()
    write_phobius(
        duplicated / "phobius" / "r.phobius", PHOBIUS_ROWS + [PHOBIUS_ROWS[0]]
    )
    predictions = _run(duplicated, tmp_path / "out_dup")
    assert predictions["protein_id"].is_unique


def test_rga_only_output(synthetic, tmp_path):
    """``--rga-only`` filters the main table but keeps the invariants."""
    outdir = tmp_path / "out_rgaonly"
    predictions = _run(synthetic, outdir, ["--rga-only"])
    assert (predictions["is_rga"] == "True").all()
    assert len(predictions) < 9


def test_reason_cites_the_channel_that_called_each_feature(predictions):
    """CC provenance must name DeepCoil2, not the InterProScan Coils interval."""
    reason = predictions.set_index("protein_id").loc["CNL1", "reason"]
    assert "CC [deepcoil+coils @" in reason
    assert reason.endswith("Confidence: high.")


def test_non_rga_reason_names_the_features_that_were_found(predictions):
    """A Non-RGA protein carrying a lone coiled coil says so explicitly."""
    reason = predictions.set_index("protein_id").loc["NONE1", "reason"]
    assert "no protein feature of interest was detected." in reason


def test_repeat_level_lrr_copy_number_is_reported(predictions):
    """`n_lrr_repeats` counts Pfam/SMART repeats, `n_lrr` merges every source."""
    row = predictions.set_index("protein_id").loc["CNL1"]
    assert int(row["n_lrr"]) == 1
    assert int(row["n_lrr_repeats"]) == 1


def test_counts_are_written_as_integers(predictions):
    """A count column must never be written as ``447.0``."""
    for column in ("sequence_length", "n_lrr", "n_tm_phobius", "rule_priority"):
        values = [v for v in predictions[column] if v != "NA"]
        assert all("." not in v for v in values), column


def test_canonical_nlr_subdomains_are_not_integrated_domains(cfg):
    """The structural sub-domains of an NLR must not be reported as integrations.

    Without the exclusion list, the winged-helix domain of the NB-ARC module and
    the Rx-type N-terminal coiled coil are flagged on nearly every NLR and the
    integrated-domain column stops meaning anything.
    """
    excluded = {str(a) for a in cfg.raw["integrated_domain_exclusions"]}
    assert {"PF23559", "PF18052", "PF25019"} <= excluded
    # ... while the textbook integrated domains stay reportable
    assert not ({"PF03106", "PF02892", "PF00403", "PF00069"} & excluded)
