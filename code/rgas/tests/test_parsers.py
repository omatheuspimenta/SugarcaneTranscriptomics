"""Parser tests against synthetic files that reproduce each tool's real format."""

from __future__ import annotations

import pytest

from conftest import (
    make_archive,
    write_deepcoil,
    write_deeploc,
    write_deeptmhmm,
    write_interproscan,
    write_phobius,
    write_signalp,
)
from rga import parsers


def test_interproscan_matches_by_accession_not_description(cfg, tmp_path):
    """A description mentioning "coil" must not create a CC feature."""
    path = tmp_path / "ips.tsv"
    write_interproscan(
        path,
        [
            (
                "P1",
                500,
                "Pfam",
                "PF00931",
                "NB-ARC domain",
                178,
                455,
                "1e-30",
                "IPR002182",
            ),
            (
                "P1",
                500,
                "ProSiteProfiles",
                "PS51808",
                "Coiled coil-helix-coiled coil-helix (CHCH) domain profile.",
                10,
                60,
                "12.0",
                "-",
            ),
        ],
    )
    result = parsers.parse_interproscan(path, cfg)
    assert set(result.hits["feature"]) == {"NB-ARC"}


def test_interproscan_excludes_mobidblite(cfg, tmp_path):
    """MobiDBLite is never allowed to contribute evidence."""
    path = tmp_path / "ips.tsv"
    write_interproscan(
        path,
        [
            (
                "P1",
                200,
                "MobiDBLite",
                "mobidb-lite",
                "consensus disorder prediction",
                1,
                40,
                "-",
                "-",
            )
        ],
    )
    result = parsers.parse_interproscan(path, cfg)
    assert result.hits.empty
    assert result.protein_ids == {"P1"}


def test_interproscan_matches_through_the_interpro_column(cfg, tmp_path):
    """A hit is recognised through its integrated InterPro accession too."""
    path = tmp_path / "ips.tsv"
    write_interproscan(
        path,
        [("P1", 400, "Gene3D", "G3DSA:9.99.99.99", "-", 10, 90, "1e-5", "IPR001611")],
    )
    result = parsers.parse_interproscan(path, cfg)
    assert list(result.hits["feature"]) == ["LRR"]


def test_interproscan_hit_matching_twice_is_not_duplicated(cfg, tmp_path):
    """A row whose signature and InterPro accessions both map to LRR yields one row."""
    path = tmp_path / "ips.tsv"
    write_interproscan(
        path,
        [
            (
                "P1",
                400,
                "Pfam",
                "PF00560",
                "Leucine Rich Repeat",
                10,
                40,
                "1e-5",
                "IPR001611",
            )
        ],
    )
    result = parsers.parse_interproscan(path, cfg)
    assert len(result.hits) == 1


def test_phobius_separates_signal_peptide_from_helices(cfg, tmp_path):
    """The ``n…c…/…`` block must never be read as a transmembrane helix."""
    path = tmp_path / "r.phobius"
    write_phobius(path, [("P1", 1, "Y", "n4-15c20/21o396-419i")])
    frame = parsers.parse_phobius(path, cfg)
    assert frame.loc[0, "tm_intervals_phobius"] == [(396, 419)]
    assert frame.loc[0, "sp_end_phobius"] == 20
    assert bool(frame.loc[0, "sp_phobius"]) is True


def test_phobius_parses_multiple_helices(cfg, tmp_path):
    """All helices of a topology string are recovered, in order."""
    path = tmp_path / "r.phobius"
    write_phobius(path, [("P1", 2, "0", "o128-143i164-190o")])
    assert parsers.parse_phobius(path, cfg).loc[0, "tm_intervals_phobius"] == [
        (128, 143),
        (164, 190),
    ]


def test_deeptmhmm_keeps_proteins_without_helices(cfg, tmp_path):
    """TM-free proteins must be recorded with zero helices, not dropped."""
    path = tmp_path / "TMRs.gff3"
    write_deeptmhmm(
        path,
        [
            ("P1", 460, [("inside", 1, 460)]),
            (
                "P2",
                780,
                [
                    ("signal", 1, 20),
                    ("outside", 21, 395),
                    ("TMhelix", 396, 416),
                    ("inside", 417, 780),
                ],
            ),
        ],
    )
    frame = parsers.parse_deeptmhmm(path, cfg).set_index("protein_id")
    assert int(frame.loc["P1", "n_tm_deeptmhmm"]) == 0
    assert frame.loc["P2", "tm_intervals_deeptmhmm"] == [(396, 416)]
    assert frame.loc["P2", "sp_end_deeptmhmm"] == 20


def test_deeptmhmm_reads_the_final_block(cfg, tmp_path):
    """The legacy parser lost the last block; this one must not."""
    path = tmp_path / "TMRs.gff3"
    text = (
        "##gff-version 3\n"
        "# P1 Length: 10\n# P1 Number of predicted TMRs: 0\nP1\tinside\t1\t10\t\t\t\t\n//\n"
        "# P2 Length: 12\n# P2 Number of predicted TMRs: 1\nP2\tTMhelix\t1\t12\t\t\t\t\n"
    )
    path.write_text(text, encoding="utf-8")
    assert set(parsers.parse_deeptmhmm(path, cfg)["protein_id"]) == {"P1", "P2"}


def test_signalp_strips_the_fasta_header(cfg, tmp_path):
    """Column 1 is a whole FASTA header; only the first token is the ID."""
    path = tmp_path / "prediction_results.txt"
    write_signalp(path, [("P1", "SP", 0.999286, "CS pos: 30-31. Pr: 0.9323")])
    frame = parsers.parse_signalp(path, cfg)
    assert frame.loc[0, "protein_id"] == "P1"
    assert frame.loc[0, "sp_end_signalp"] == 30
    assert frame.loc[0, "cleavage_site"] == "30-31"
    assert bool(frame.loc[0, "sp_signalp"]) is True


def test_signalp_other_is_not_a_signal_peptide(cfg, tmp_path):
    """The OTHER class must not be counted as a signal peptide."""
    path = tmp_path / "prediction_results.txt"
    write_signalp(path, [("P1", "OTHER", 0.0, "")])
    frame = parsers.parse_signalp(path, cfg)
    assert bool(frame.loc[0, "sp_signalp"]) is False


def test_deeploc_multilabel_localizations(cfg, tmp_path):
    """``Cytoplasm|Nucleus`` reports Cytoplasm as primary with its own probability."""
    path = tmp_path / "results.csv"
    write_deeploc(
        path, [("P1", "Cytoplasm|Nucleus", {"Cytoplasm": 0.71, "Nucleus": 0.64})]
    )
    frame = parsers.parse_deeploc(path, cfg)
    assert frame.loc[0, "predicted_localization"] == "Cytoplasm"
    assert frame.loc[0, "localization_prob"] == pytest.approx(0.71)
    assert frame.loc[0, "all_localizations"] == "Cytoplasm|Nucleus"


def test_deepcoil_segments_split_on_a_change_of_plateau(cfg, tmp_path):
    """Adjacent segments with different scores are two segments, not one."""
    scores = [0.0] * 5 + [0.43] * 10 + [0.66] * 10 + [0.0] * 5
    directory = tmp_path / "part_001"
    write_deepcoil(directory, "P.1.p", scores)
    frame = parsers.parse_deepcoil(directory, cfg)
    assert list(zip(frame["start"], frame["end"], frame["cc"])) == [
        (6, 15, 0.43),
        (16, 25, 0.66),
    ]


def test_deepcoil_reads_from_an_archive(cfg, tmp_path):
    """Segments can be streamed straight out of a ``.tar.xz`` without extracting."""
    directory = tmp_path / "part_001"
    write_deepcoil(directory, "P.1.p", [0.0, 0.7, 0.7, 0.0])
    archive = make_archive(directory, tmp_path / "part_001.tar.xz")
    frame = parsers.parse_deepcoil(archive, cfg)
    assert list(frame["deepcoil_id"]) == ["P1p"]
    assert list(frame["start"]) == [2]


def test_deepcoil_positions_are_one_based():
    """Residue 1 is the first row after the header."""
    lines = [
        "aa\tcc\traw_cc\tprob_a\tprob_d",
        "M\t0.900\t0.9\t0\t0",
        "A\t0.000\t0\t0\t0",
    ]
    assert parsers._segments_from_lines(lines) == [(1, 1, 0.9)]
