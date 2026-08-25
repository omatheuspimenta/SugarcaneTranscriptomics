"""Shared fixtures: the real configuration and a synthetic mini-proteome."""

from __future__ import annotations

import sys
import tarfile
from pathlib import Path

import pytest

CODE_DIR = Path(__file__).resolve().parents[1]
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from rga.config import load_config  # noqa: E402

CONFIG_PATH = CODE_DIR / "config" / "rga_config.yaml"


@pytest.fixture(scope="session")
def cfg():
    """The production configuration, as shipped."""
    return load_config(CONFIG_PATH)


@pytest.fixture
def config_path() -> Path:
    """Path to the production configuration."""
    return CONFIG_PATH


def write_interproscan(path: Path, rows: list[tuple]) -> None:
    """Write a synthetic 15-column InterProScan TSV.

    Each row is ``(protein_id, length, analysis, accession, description,
    start, end, score, interpro_accession)``.
    """
    lines = []
    for (
        protein_id,
        length,
        analysis,
        accession,
        description,
        start,
        end,
        score,
        ipr,
    ) in rows:
        lines.append(
            "\t".join(
                [
                    protein_id,
                    "0" * 32,
                    str(length),
                    analysis,
                    accession,
                    description,
                    str(start),
                    str(end),
                    str(score),
                    "T",
                    "09-08-2026",
                    ipr,
                    "-",
                    "-",
                    "-",
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_phobius(path: Path, rows: list[tuple[str, int, str, str]]) -> None:
    """Write a synthetic Phobius short-format file."""
    lines = ["SEQENCE ID                     TM SP PREDICTION"]
    lines.extend(f"{pid}  {tm}  {sp} {topology}" for pid, tm, sp, topology in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_deeptmhmm(
    path: Path, blocks: list[tuple[str, int, list[tuple[str, int, int]]]]
) -> None:
    """Write a synthetic DeepTMHMM ``TMRs.gff3``."""
    lines = ["##gff-version 3"]
    for protein_id, length, regions in blocks:
        n_tm = sum(1 for region, _, _ in regions if region == "TMhelix")
        lines.append(f"# {protein_id} Length: {length}")
        lines.append(f"# {protein_id} Number of predicted TMRs: {n_tm}")
        lines.extend(
            f"{protein_id}\t{region}\t{start}\t{end}\t\t\t\t"
            for region, start, end in regions
        )
        lines.append("//")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_signalp(path: Path, rows: list[tuple[str, str, float, str]]) -> None:
    """Write a synthetic SignalP 6.0 ``prediction_results.txt``."""
    lines = [
        "# SignalP-6.0\tOrganism: Other\tTimestamp: 20260811172837",
        "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tLIPO(Sec/SPII)\tTAT(Tat/SPI)\t"
        "TATLIPO(Tat/SPII)\tPILIN(Sec/SPIII)\tCS Position",
    ]
    for protein_id, prediction, probability, cs in rows:
        other = 1.0 - probability
        lines.append(
            "\t".join(
                [
                    f"{protein_id} pacid=1 transcript=x locus=y",
                    prediction,
                    f"{other:.6f}",
                    f"{probability:.6f}",
                    "0.000000",
                    "0.000000",
                    "0.000000",
                    "0.000000",
                    cs,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_deeploc(path: Path, rows: list[tuple[str, str, dict[str, float]]]) -> None:
    """Write a synthetic DeepLoc 2.0 results CSV."""
    classes = [
        "Cytoplasm",
        "Nucleus",
        "Extracellular",
        "Cell membrane",
        "Mitochondrion",
        "Plastid",
        "Endoplasmic reticulum",
        "Lysosome/Vacuole",
        "Golgi apparatus",
        "Peroxisome",
    ]
    header = [
        "Protein_ID",
        "Localizations",
        "Signals",
        "Membrane types",
        *classes,
        "Peripheral",
        "Transmembrane",
        "Lipid anchor",
        "Soluble",
    ]
    lines = [",".join(header)]
    for protein_id, localization, probabilities in rows:
        values = [str(probabilities.get(name, 0.0)) for name in classes]
        lines.append(
            ",".join(
                [protein_id, localization, "", "Soluble", *values, "0", "0", "0", "1"]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_deepcoil(directory: Path, protein_id: str, scores: list[float]) -> Path:
    """Write one synthetic DeepCoil2 ``.out`` file, named the DeepCoil way."""
    directory.mkdir(parents=True, exist_ok=True)
    stem = protein_id.replace(".", "")
    path = directory / f"{stem}.out"
    lines = ["aa\tcc\traw_cc\tprob_a\tprob_d"]
    lines.extend(f"A\t{value:.3f}\t{value:.3f}\t0.000\t0.000" for value in scores)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def make_archive(directory: Path, archive: Path) -> Path:
    """Pack a directory of ``.out`` files into a ``.tar.xz`` like DeepCoil2 ships."""
    with tarfile.open(archive, "w:xz") as tar:
        tar.add(directory, arcname=directory.name)
    return archive
