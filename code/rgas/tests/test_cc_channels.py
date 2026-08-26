"""The three coiled-coil channels: policies, provenance and confidence.

Config v1.1.0 promoted the Rx N-terminal domain to a CC evidence channel of its
own. These tests pin the three properties that change makes load-bearing:

* the domain channel can call a CC on its own, with no predictor support;
* a call it backs is *not* demoted, because a curated profile HMM is the same
  grade of evidence as the ``PF00931`` hit that made the protein an NLR at all;
* every channel is still reported separately, so a reader can see which one
  fired rather than having to trust the consensus.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

import rgas_prediction
from conftest import (
    write_deepcoil,
    write_deeploc,
    write_deeptmhmm,
    write_interproscan,
    write_phobius,
    write_signalp,
)
from rga.config import ConfigError, load_config
from rga.evidence import apply_cc_policy

# Three NB-ARC + LRR proteins that differ only in which CC channel supports them.
# Under the old two-predictor configuration RXONLY would have been called `NL`:
# an NLR reported as carrying no coiled coil while holding the very domain that
# is the coiled coil of an Rx-type CNL.
IPS_ROWS = [
    # CC from the domain model alone
    ("RXONLY", 900, "Pfam", "PF00931", "NB-ARC domain", 200, 460, "1e-40", "IPR002182"),
    ("RXONLY", 900, "Pfam", "PF00560", "Leucine Rich Repeat", 520, 570, "1e-9",
     "IPR001611"),
    ("RXONLY", 900, "Pfam", "PF18052", "Rx N-terminal domain", 20, 95, "1e-12",
     "IPR041118"),
    # CC from DeepCoil2 alone
    ("DCONLY", 900, "Pfam", "PF00931", "NB-ARC domain", 200, 460, "1e-40", "IPR002182"),
    ("DCONLY", 900, "Pfam", "PF00560", "Leucine Rich Repeat", 520, 570, "1e-9",
     "IPR001611"),
    # CC from InterProScan Coils alone
    ("COONLY", 900, "Pfam", "PF00931", "NB-ARC domain", 200, 460, "1e-40", "IPR002182"),
    ("COONLY", 900, "Pfam", "PF00560", "Leucine Rich Repeat", 520, 570, "1e-9",
     "IPR001611"),
    ("COONLY", 900, "Coils", "Coil", "Coil", 20, 60, "-", "-"),
    # every channel agreeing
    ("ALLTHREE", 900, "Pfam", "PF00931", "NB-ARC", 200, 460, "1e-40", "IPR002182"),
    ("ALLTHREE", 900, "Pfam", "PF00560", "LRR", 520, 570, "1e-9", "IPR001611"),
    ("ALLTHREE", 900, "Pfam", "PF18052", "Rx N-terminal domain", 20, 95, "1e-12",
     "IPR041118"),
    ("ALLTHREE", 900, "Coils", "Coil", "Coil", 20, 60, "-", "-"),
    # no CC at all
    ("NOCC", 900, "Pfam", "PF00931", "NB-ARC domain", 200, 460, "1e-40", "IPR002182"),
    ("NOCC", 900, "Pfam", "PF00560", "Leucine Rich Repeat", 520, 570, "1e-9",
     "IPR001611"),
]

DEEPCOIL = {
    "RXONLY": [0.0] * 900,
    "DCONLY": [0.0] * 20 + [0.8] * 40 + [0.0] * 840,
    "COONLY": [0.0] * 900,
    "ALLTHREE": [0.0] * 20 + [0.8] * 40 + [0.0] * 840,
    "NOCC": [0.0] * 900,
}
PROTEINS = list(DEEPCOIL)


@pytest.fixture(scope="module")
def inputs(tmp_path_factory) -> Path:
    """A five-protein proteome differing only in CC channel support."""
    root = tmp_path_factory.mktemp("cc_inputs")
    for name in ("InterProScan", "phobius", "DeepTMHMM", "SignalP6", "DeepLoc2"):
        (root / name).mkdir()
    write_interproscan(root / "InterProScan" / "ips.tsv", IPS_ROWS)
    write_phobius(
        root / "phobius" / "r.phobius", [(p, 0, "0", "o") for p in PROTEINS]
    )
    write_deeptmhmm(root / "DeepTMHMM" / "TMRs.gff3", [(p, 900, []) for p in PROTEINS])
    write_signalp(
        root / "SignalP6" / "prediction_results.txt",
        [(p, "OTHER", 1.0, "") for p in PROTEINS],
    )
    write_deeploc(
        root / "DeepLoc2" / "results_test.csv",
        [(p, "Nucleus", {"Nucleus": 0.9}) for p in PROTEINS],
    )
    for protein_id, scores in DEEPCOIL.items():
        write_deepcoil(root / "DeepCoil" / "part_001", protein_id, scores)
    return root


def _run(inputs: Path, outdir: Path, extra: list[str] | None = None) -> pd.DataFrame:
    """Run the pipeline and return the prediction table indexed by protein."""
    argv = [
        "--input-dir", str(inputs),
        "--outdir", str(outdir),
        "--organism-name", "CCChannels",
        "--log-level", "WARNING",
        *(extra or []),
    ]
    assert rgas_prediction.main(argv) == 0
    frame = pd.read_csv(
        outdir / "rga_predictions.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    return frame.set_index("protein_id")


@pytest.fixture(scope="module")
def predictions(inputs, tmp_path_factory) -> pd.DataFrame:
    """Predictions under the shipped default policy (``union``)."""
    return _run(inputs, tmp_path_factory.mktemp("cc_out"))


# --------------------------------------------------------------------------
# the policy function itself
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("policy", "channels", "expected"),
    [
        ("union", {"rx_domain": True, "deepcoil": False, "coils": False}, True),
        ("union", {"rx_domain": False, "deepcoil": False, "coils": False}, False),
        ("intersection", {"rx_domain": True, "deepcoil": True, "coils": True}, True),
        ("intersection", {"rx_domain": True, "deepcoil": True, "coils": False}, False),
        ("rx_domain", {"rx_domain": True, "deepcoil": False, "coils": False}, True),
        ("rx_domain", {"rx_domain": False, "deepcoil": True, "coils": True}, False),
        ("deepcoil", {"rx_domain": True, "deepcoil": False, "coils": True}, False),
        ("coils", {"rx_domain": True, "deepcoil": True, "coils": False}, False),
    ],
)
def test_apply_cc_policy(policy, channels, expected):
    """Each policy selects exactly the channels its name promises."""
    assert apply_cc_policy(policy, channels) is expected


def test_missing_channel_never_vetoes_an_intersection():
    """An unavailable tool is skipped, not counted as a negative."""
    channels = {"rx_domain": True, "deepcoil": None, "coils": True}
    assert apply_cc_policy("intersection", channels) is True


def test_naming_an_unavailable_channel_falls_back_to_the_union():
    """``--cc-policy deepcoil`` without DeepCoil2 must not silently mean "no CC"."""
    channels = {"rx_domain": True, "deepcoil": None, "coils": False}
    assert apply_cc_policy("deepcoil", channels) is True


def test_no_channel_available_is_no_evidence():
    """All-``None`` is the honest negative, not an exception."""
    assert apply_cc_policy("union", {"rx_domain": None, "deepcoil": None}) is False


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------


def test_domain_channel_alone_calls_a_cnl(predictions):
    """The whole point of the change: PF18052 with no predictor support is a CNL."""
    row = predictions.loc["RXONLY"]
    assert row["cc_rx_domain"] == "True"
    assert row["cc_deepcoil"] == "False"
    assert row["cc_coils"] == "False"
    assert row["cc_consensus"] == "True"
    assert row["rga_subclass"] == "CNL"


def test_a_domain_backed_cc_is_not_demoted(predictions):
    """A curated profile HMM is not second-class evidence."""
    row = predictions.loc["RXONLY"]
    assert row["confidence"] == "high"
    assert row["confidence_demotions"] == "NA"


def test_predictor_only_calls_are_still_demoted(predictions):
    """Without a domain model behind it, a propensity score is still hedged."""
    deepcoil_only = predictions.loc["DCONLY"]
    assert deepcoil_only["cc_source"] == "deepcoil_only"
    assert "cc_deepcoil_only" in deepcoil_only["confidence_demotions"]
    assert deepcoil_only["confidence"] == "medium"

    coils_only = predictions.loc["COONLY"]
    assert coils_only["cc_source"] == "coils_only"
    assert "cc_coils_only" in coils_only["confidence_demotions"]
    assert coils_only["confidence"] == "low"


def test_cc_source_names_every_contributing_channel(predictions):
    """Provenance is reported, not summarised away."""
    assert predictions.loc["ALLTHREE"]["cc_source"] == "rx_domain+deepcoil+coils"
    assert predictions.loc["RXONLY"]["cc_source"] == "rx_domain_only"
    assert predictions.loc["NOCC"]["cc_source"] == "NA"


def test_reason_states_the_domain_channel(predictions):
    """A biologist reading one row can see which channel carried the call."""
    assert "Rx domain yes" in predictions.loc["RXONLY"]["reason"]
    assert "Rx domain no" in predictions.loc["DCONLY"]["reason"]


def test_domain_channel_coordinates_are_the_domain_hit(predictions):
    """Coordinates come from the channel that fired, not from a default."""
    assert predictions.loc["RXONLY"]["cc_coords"] == "20-95"


def test_deepcoil_coordinates_win_when_both_fired(predictions):
    """DeepCoil2 resolves a segment per residue, so it is the precise channel."""
    assert predictions.loc["ALLTHREE"]["cc_coords"] == "21-60"


def test_protein_without_any_cc_channel_is_nl(predictions):
    """The negative control: no CC channel means no CC."""
    row = predictions.loc["NOCC"]
    assert row["cc_consensus"] == "False"
    assert row["rga_subclass"] == "NL"


def test_domain_hits_reach_the_long_evidence_table(inputs, tmp_path):
    """The audit trail carries the domain hit that made the call."""
    _run(inputs, tmp_path / "long_out")
    long = pd.read_csv(
        tmp_path / "long_out" / "rga_domain_evidence_long.tsv",
        sep="\t",
        keep_default_na=False,
    )
    rx = long[(long["protein_id"] == "RXONLY") & (long["accession"] == "PF18052")]
    assert len(rx) == 1
    assert int(rx.iloc[0]["start"]) == 20


def test_deepcoil_rows_are_not_attributed_a_domain_interval(inputs, tmp_path):
    """A CC called by the domain model must not be exported as a DeepCoil2 hit."""
    _run(inputs, tmp_path / "prov_out")
    long = pd.read_csv(
        tmp_path / "prov_out" / "rga_domain_evidence_long.tsv",
        sep="\t",
        keep_default_na=False,
    )
    deepcoil_rows = long[long["tool"] == "DeepCoil2"]
    assert set(deepcoil_rows["protein_id"]) == {"DCONLY", "ALLTHREE"}


def test_policy_flag_selects_a_single_channel(inputs, tmp_path):
    """``--cc-policy rx_domain`` ignores both predictors."""
    predictions = _run(inputs, tmp_path / "rx_out", ["--cc-policy", "rx_domain"])
    assert predictions.loc["RXONLY"]["rga_subclass"] == "CNL"
    assert predictions.loc["DCONLY"]["rga_subclass"] == "NL"
    assert predictions.loc["COONLY"]["rga_subclass"] == "NL"


def test_policy_sensitivity_covers_all_five_policies(inputs, tmp_path):
    """The sensitivity table must not silently keep reporting only four."""
    _run(inputs, tmp_path / "sens_out")
    table = pd.read_csv(
        tmp_path / "sens_out" / "cc_policy_sensitivity.tsv", sep="\t"
    )
    assert set(table.columns) == {
        "rga_subclass",
        "rx_domain",
        "deepcoil",
        "coils",
        "union",
        "intersection",
    }
    counts = table.set_index("rga_subclass")
    assert counts.loc["CNL", "union"] == 4
    assert counts.loc["CNL", "intersection"] == 1
    assert counts.loc["CNL", "rx_domain"] == 2


def test_contingency_reports_the_domain_channel(inputs, tmp_path):
    """``run_metadata.json`` records how often the domain model fired alone."""
    import json

    _run(inputs, tmp_path / "meta_out")
    meta = json.loads((tmp_path / "meta_out" / "run_metadata.json").read_text())
    contingency = meta["cc_contingency_deepcoil_vs_coils"]
    assert contingency["rx_domain"] == 2
    assert contingency["rx_domain_only"] == 1


# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------


def test_domain_accessions_may_not_double_as_feature_evidence(config_path, tmp_path):
    """An accession feeding two CC channels would make ``intersection`` a lie."""
    raw = yaml.safe_load(config_path.read_text())
    raw["cc_domain_accessions"] = ["Coil"]
    broken = tmp_path / "broken.yaml"
    broken.write_text(yaml.safe_dump(raw))
    with pytest.raises(ConfigError, match="may not repeat an accession"):
        load_config(broken)


def test_empty_domain_list_restores_the_two_predictor_behaviour(config_path, tmp_path):
    """The channel is switchable off from configuration alone, as documented."""
    raw = yaml.safe_load(config_path.read_text())
    raw["cc_domain_accessions"] = []
    path = tmp_path / "no_domain.yaml"
    path.write_text(yaml.safe_dump(raw))
    cfg = load_config(path)
    assert cfg.cc_domain_accessions() == ()
    assert "PF18052" not in cfg.accession_to_features()


def test_shipped_config_enables_the_domain_channel(cfg):
    """The reference run is the one documented in docs/rga/README.md."""
    assert cfg.cc_domain_accessions() == ("PF18052", "IPR041118")
    assert cfg.policies["cc"] == "union"


def test_domain_accessions_appear_in_the_accession_audit(inputs, tmp_path):
    """The channel driving the CNL count must not be the one missing from the audit."""
    _run(inputs, tmp_path / "audit_out")
    audit = pd.read_csv(
        tmp_path / "audit_out" / "accession_audit.tsv", sep="\t", keep_default_na=False
    ).set_index("accession")
    assert audit.loc["PF18052", "feature"] == "CC (domain channel)"
    assert int(audit.loc["PF18052", "n_hits"]) == 2
    assert audit.loc["IPR041118", "feature"] == "CC (domain channel)"


def test_feature_coords_and_accessions_are_machine_readable(predictions):
    """The provenance in `reason` is prose; these two columns are parseable."""
    row = predictions.loc["RXONLY"]
    coords = dict(part.split(":", 1) for part in row["feature_coords"].split(";"))
    assert coords["NB-ARC"] == "200-460"
    assert coords["CC"] == "20-95"
    accessions = dict(
        part.split(":", 1) for part in row["feature_accessions"].split(";")
    )
    # Signature accessions, as in `reason`: the InterPro accession in column 12
    # contributes to feature matching but a hit is only ever recorded once,
    # under the signature that produced it.
    assert accessions["NB-ARC"] == "PF00931"
    assert accessions["CC"] == "PF18052"


def test_feature_columns_list_only_features_the_protein_has(predictions):
    """A feature absent from `features_found` must not appear in either column."""
    row = predictions.loc["NOCC"]
    assert "CC:" not in row["feature_coords"]
    assert "CC:" not in row["feature_accessions"]
    for column in ("feature_coords", "feature_accessions"):
        listed = {part.split(":", 1)[0] for part in row[column].split(";")}
        assert listed == set(row["features_found"].split(";"))


def test_feature_accessions_survive_accessions_containing_a_colon():
    """Gene3D accessions look like `G3DSA:3.80.10.10`.

    The column uses `FEATURE:value` pairs, so it must be split on the *first*
    colon only. Pinned because the reference proteome is full of Gene3D hits and
    a naive `split(":")` silently truncates every one of them.
    """
    from rga.rules import _format_feature_map

    rendered = _format_feature_map(
        {"LRR": ["G3DSA:3.80.10.10", "PF00560"]},
        ["LRR"],
        lambda values: ",".join(values) if values else "",
    )
    assert rendered == "LRR:G3DSA:3.80.10.10,PF00560"
    feature, value = rendered.split(":", 1)
    assert feature == "LRR"
    assert value.split(",") == ["G3DSA:3.80.10.10", "PF00560"]


def test_report_prints_a_runnable_command(inputs, tmp_path):
    """The 'Reproduce this run' block must be pasteable, not merely indicative.

    `" ".join(sys.argv)` loses the quoting around an organism name with spaces,
    so the recorded command would not have been the command that ran.
    """
    import shlex

    outdir = tmp_path / "cmd_out"
    argv = [
        "--input-dir", str(inputs),
        "--outdir", str(outdir),
        "--organism-name", "Genus species cv. Two Words",
        "--log-level", "WARNING",
    ]
    assert rgas_prediction.main(argv) == 0
    recorded = json.loads((outdir / "run_metadata.json").read_text())["command"]
    assert shlex.split(recorded) == ["code/rgas/rgas_prediction.py", *argv]
    report = (outdir / "report.md").read_text()
    assert "### Reproduce this run" in report
    assert "'Genus species cv. Two Words'" in report


def test_report_is_light_only(inputs, tmp_path):
    """The report commits to one palette rather than following the reader's theme.

    It is a printable record that gets screenshotted and pasted into documents,
    so it must render the same for everyone.
    """
    outdir = tmp_path / "light_out"
    _run(inputs, outdir)
    page = (outdir / "report.html").read_text()
    assert "prefers-color-scheme" not in page
    assert "color-scheme: light" in page
