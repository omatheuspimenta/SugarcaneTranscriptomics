"""Configuration loading, validation and identifier normalisation."""

from __future__ import annotations


import pytest
import yaml

from rga.config import ConfigError, load_config, normalize_id


def test_shipped_config_loads(cfg):
    """The production configuration is valid as shipped."""
    assert cfg.raw["config_version"]
    assert cfg.rules[0].priority == 1
    assert [rule.priority for rule in cfg.rules] == sorted(
        rule.priority for rule in cfg.rules
    )


def test_accession_map_is_inverted_without_loss(cfg):
    """Every configured accession appears in the inverted lookup."""
    inverted = cfg.accession_to_features()
    for feature, accessions in cfg.raw["interproscan_features"].items():
        for accession in accessions:
            assert feature in inverted[str(accession)]


def test_mobidblite_is_excluded(cfg):
    """MobiDBLite must never be usable as evidence."""
    assert "MobiDBLite" in cfg.raw["excluded_analyses"]


def test_the_commented_out_other_rlk_rule_still_works(cfg, config_path, tmp_path):
    """Uncommenting the priority-13 `other-RLK` block must produce a valid rule set.

    The block is the documented opt-in for RGAugury RLK scope (README section
    6.4), and it is commented out, so nothing else in the suite would notice it
    rotting -- a renamed feature or a dropped key would only surface for the
    first person who tried to enable it.
    """
    import re

    from rga import rules

    text = config_path.read_text(encoding="utf-8")
    block = re.search(
        r"^  # - id: other-RLK\n(?:  #.*\n)+", text, re.M
    )
    assert block, "the commented-out other-RLK block is no longer in the config"
    uncommented = re.sub(r"^  # ?", "  ", block.group(0), flags=re.M)

    path = tmp_path / "with_other_rlk.yaml"
    path.write_text(text.replace(block.group(0), uncommented, 1), encoding="utf-8")
    enabled = load_config(path)

    assert "other-RLK" in {rule.id for rule in enabled.rules}
    assert len(enabled.rules) == len(cfg.rules) + 1
    assert rules.find_overlapping_rules(enabled) == []

    core = frozenset(enabled.core_immune_features)
    for anchor in ("TM", "SP"):
        call = rules.classify_features(enabled.rules, frozenset({"STTK", anchor}), core)
        assert (call.rule_id, call.family) == ("other-RLK", "RLK"), anchor

    # and the id the confidence section already carries is now a live class
    assert "other-RLK" in enabled.raw["confidence"]["classes_using_tm_sp"]


def test_ectodomain_token_is_expanded(cfg):
    """The ECTODOMAIN_FEATURES placeholder resolves to the configured list."""
    rule = next(r for r in cfg.rules if r.id == "other-RLP")
    assert tuple(cfg.ectodomain_features) in rule.any_of


def _write(tmp_path, mutate):
    """Write a mutated copy of the production configuration."""
    from conftest import CONFIG_PATH

    raw = yaml.safe_load(CONFIG_PATH.read_text())
    mutate(raw)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(raw))
    return path


def test_unknown_policy_is_rejected(tmp_path):
    """An invalid consensus policy must fail at load time, not at run time."""
    path = _write(tmp_path, lambda raw: raw["policies"].update({"cc": "magic"}))
    with pytest.raises(ConfigError, match="invalid cc policy"):
        load_config(path)


def test_rule_referencing_an_unknown_feature_is_rejected(tmp_path):
    """A typo in a rule must be caught by validation."""

    def mutate(raw):
        raw["rules"][0]["all_of"] = ["NB-ARK"]

    with pytest.raises(ConfigError, match="unknown feature"):
        load_config(_write(tmp_path, mutate))


def test_duplicate_rule_priorities_are_rejected(tmp_path):
    """Two rules may not share a priority: evaluation order would be ambiguous."""

    def mutate(raw):
        raw["rules"][1]["priority"] = raw["rules"][0]["priority"]

    with pytest.raises(ConfigError, match="duplicate rule priorities"):
        load_config(_write(tmp_path, mutate))


def test_missing_top_level_key_is_rejected(tmp_path):
    """A truncated configuration fails loudly."""

    def mutate(raw):
        del raw["confidence"]

    with pytest.raises(ConfigError, match="missing configuration keys"):
        load_config(_write(tmp_path, mutate))


def test_unknown_id_operation_is_rejected(tmp_path):
    """Only implemented ID normalisations may be configured."""

    def mutate(raw):
        raw["ids"]["per_tool"]["phobius"] = ["strip_everything"]

    with pytest.raises(ConfigError, match="unknown ID normalisation"):
        load_config(_write(tmp_path, mutate))


@pytest.mark.parametrize(
    "raw_id,ops,expected",
    [
        ("PROT.1.p pacid=1 locus=x", ["strip_after_whitespace"], "PROT.1.p"),
        ("SoffiXsponR570.01Bg000200.1.p", ["strip_dots"], "SoffiXsponR57001Bg0002001p"),
        ("sp|Q9XYZ1|NAME", ["strip_after_pipe"], "sp"),
        ("PROT.1.p", ["strip_after_whitespace", "strip_dots"], "PROT1p"),
    ],
)
def test_normalize_id(raw_id, ops, expected):
    """Each documented normalisation behaves as advertised."""
    assert normalize_id(raw_id, ops) == expected


def test_rstrip_suffixes():
    """TransDecoder-style suffixes can be removed when configured."""
    assert normalize_id("PROT.1.p1", ["rstrip_suffixes"], [".p1"]) == "PROT.1"


def test_no_accession_is_both_evidence_and_watch_only(cfg):
    """`watch_accessions` means "counted, never used"; using one as evidence too
    would make the audit table contradict itself."""
    evidence = {
        str(accession)
        for accessions in cfg.raw["interproscan_features"].values()
        for accession in accessions
    }
    evidence |= set(cfg.cc_domain_accessions())
    clash = sorted(evidence & {str(a) for a in cfg.raw.get("watch_accessions", {})})
    assert not clash, f"accessions used as evidence and also watch-only: {clash}"


def test_nacht_is_not_nb_arc_evidence(cfg):
    """NACHT (PF05729/IPR007111) is a distinct NTPase domain from NB-ARC.

    Treating it as NB-ARC would report NACHT proteins as NLRs with an NB-ARC
    domain. Harmless in R570 (0 hits) but wrong for any fungal or animal
    proteome, so it is pinned rather than left to a comment.
    """
    nb_arc = {str(a) for a in cfg.raw["interproscan_features"]["NB-ARC"]}
    assert not nb_arc & {"PF05729", "IPR007111"}
    assert "PF05729" in cfg.raw["watch_accessions"]
