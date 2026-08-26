"""Rule-engine tests: mutual exclusivity, one fixture per rule, architecture."""

from __future__ import annotations

import pytest

from rga import rules

#: One synthetic feature set per rule of the shipped configuration.
FIXTURES: dict[str, set[str]] = {
    "CNL": {"NB-ARC", "CC", "LRR"},
    "TNL": {"NB-ARC", "TIR", "LRR"},
    "RNL": {"NB-ARC", "RPW8", "LRR"},
    "NL": {"NB-ARC", "LRR"},
    "CN": {"NB-ARC", "CC"},
    "TN": {"NB-ARC", "TIR"},
    "RN": {"NB-ARC", "RPW8"},
    "N": {"NB-ARC"},
    "TX": {"TIR"},
    "RX": {"RPW8"},
    "LRR-RLK": {"STTK", "LRR", "TM"},
    "LysM-RLK": {"STTK", "LysM", "SP"},
    "other-RLK": {"STTK", "TM"},
    "LRR-RLP": {"LRR", "TM"},
    "LysM-RLP": {"LysM", "SP"},
    "TM-CC": {"TM", "CC"},
    "Other": {"STTK"},
    "Non-RGA": {"TM", "SP"},
}


def test_rules_are_mutually_exclusive(cfg):
    """No two non-fallback rules may fire on the same feature combination."""
    assert rules.find_overlapping_rules(cfg) == []


def test_every_feature_combination_matches_exactly_one_rule(cfg):
    """Ordered evaluation must assign exactly one class to every feature set."""
    import itertools

    core = frozenset(cfg.core_immune_features)
    for size in range(len(cfg.features) + 1):
        for combination in itertools.combinations(cfg.features, size):
            call = rules.classify_features(cfg.rules, frozenset(combination), core)
            assert call.rule_id


@pytest.mark.parametrize("rule_id", sorted(FIXTURES))
def test_one_fixture_per_rule(cfg, rule_id):
    """Each rule fires on its canonical synthetic protein."""
    core = frozenset(cfg.core_immune_features)
    call = rules.classify_features(cfg.rules, frozenset(FIXTURES[rule_id]), core)
    assert call.rule_id == rule_id


def test_other_rlp_is_unreachable_by_design(cfg):
    """`other-RLP` cannot fire while LRR and LysM are the only ectodomains."""
    core = frozenset(cfg.core_immune_features)
    fired = {
        rules.classify_features(cfg.rules, frozenset(c), core).rule_id
        for c in _all_combinations(cfg.features)
    }
    assert "other-RLP" not in fired


def test_protein_without_any_feature_is_non_rga(cfg):
    """A protein with no hits at all must still appear, as Non-RGA."""
    call = rules.classify_features(
        cfg.rules, frozenset(), frozenset(cfg.core_immune_features)
    )
    assert (call.family, call.subclass) == ("Non-RGA", "NA")


def test_lone_coiled_coil_is_not_an_rga(cfg):
    """CC alone is not core immune evidence under the shipped configuration."""
    call = rules.classify_features(
        cfg.rules, frozenset({"CC"}), frozenset(cfg.core_immune_features)
    )
    assert call.subclass == "NA"


def test_tir_kinase_stays_in_tx(cfg):
    """A TIR protein carrying a kinase is TX, not an RLK, under this rule order."""
    core = frozenset(cfg.core_immune_features)
    call = rules.classify_features(cfg.rules, frozenset({"TIR", "STTK", "TM"}), core)
    assert call.rule_id == "TX"


def test_domain_architecture_is_ordered_n_to_c():
    """Architecture strings follow coordinate order and collapse repeats."""
    intervals = {
        "NB-ARC": [(178, 455)],
        "LRR": [(512, 700), (701, 889)],
        "CC": [(21, 48)],
        "TM": [],
        "SP": [],
    }
    architecture = rules.domain_architecture(intervals, ["CC", "NB-ARC", "LRR"])
    assert architecture == "CC-NB-ARC-LRR"


def test_domain_architecture_keeps_features_without_coordinates():
    """A called feature with no interval is still reported."""
    architecture = rules.domain_architecture({"LRR": [(10, 40)]}, ["LRR", "TM"])
    assert architecture == "LRR-TM"


def _all_combinations(features):
    """Yield every subset of the feature vocabulary."""
    import itertools

    for size in range(len(features) + 1):
        yield from itertools.combinations(features, size)
