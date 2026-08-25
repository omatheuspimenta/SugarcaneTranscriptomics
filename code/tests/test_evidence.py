"""Evidence-layer tests: interval merging, CC segment calling, consensus policies."""

from __future__ import annotations

import pytest

from rga import evidence


# ---------------------------------------------------------------------------
# interval merging
# ---------------------------------------------------------------------------


def test_merge_intervals_collapses_overlaps():
    """Overlapping intervals from redundant databases collapse into one."""
    assert evidence.merge_intervals([(10, 20), (18, 30), (50, 60)]) == [
        (10, 30),
        (50, 60),
    ]


def test_overlapping_lrr_hits_from_three_databases_count_once():
    """Pfam, SMART and Gene3D calling the same LRR must yield a single copy."""
    hits = [(512, 560), (515, 558), (510, 562)]
    assert len(evidence.merge_intervals(hits)) == 1


def test_adjacent_intervals_are_not_merged():
    """Intervals that only touch end-to-start stay separate."""
    assert evidence.merge_intervals([(1, 10), (11, 20)]) == [(1, 10), (11, 20)]


def test_overlap_fraction_is_relative_to_the_first_interval():
    """The fraction is covered residues over the length of the query interval."""
    assert evidence.overlap_fraction((1, 10), [(1, 5)]) == pytest.approx(0.5)
    assert evidence.overlap_fraction((1, 10), []) == 0.0


# ---------------------------------------------------------------------------
# coiled-coil segment calling
# ---------------------------------------------------------------------------


def test_segment_just_below_min_length_is_dropped():
    """A 20-residue segment does not survive a 21-residue minimum."""
    called = evidence.call_cc_segments(
        [(1, 20, 0.8)], threshold=0.5, min_length=21, max_gap=2
    )
    assert called.n == 0


def test_segment_just_above_min_length_is_kept():
    """A 21-residue segment is exactly at the limit and is kept."""
    called = evidence.call_cc_segments(
        [(1, 21, 0.8)], threshold=0.5, min_length=21, max_gap=2
    )
    assert called.n == 1
    assert called.total_length == 21


def test_segment_below_threshold_is_dropped():
    """Plateau scores under the threshold never produce a call."""
    called = evidence.call_cc_segments(
        [(1, 40, 0.49)], threshold=0.5, min_length=21, max_gap=2
    )
    assert called.n == 0


def test_two_segments_at_the_gap_boundary_merge():
    """A gap of exactly ``max_gap`` residues is merged."""
    called = evidence.call_cc_segments(
        [(1, 15, 0.8), (18, 30, 0.7)], threshold=0.5, min_length=21, max_gap=2
    )
    assert called.n == 1
    assert called.intervals == [(1, 30)]


def test_two_segments_past_the_gap_boundary_stay_separate():
    """A gap of ``max_gap + 1`` residues is not merged, so both fail the length filter."""
    called = evidence.call_cc_segments(
        [(1, 15, 0.8), (19, 30, 0.7)], threshold=0.5, min_length=21, max_gap=2
    )
    assert called.n == 0


def test_merged_segment_keeps_the_best_score():
    """Merging two segments reports the stronger plateau score."""
    called = evidence.call_cc_segments(
        [(1, 15, 0.6), (17, 40, 0.9)], threshold=0.5, min_length=21, max_gap=2
    )
    assert called.max_prob == pytest.approx(0.9)


def test_mean_probability_is_length_weighted():
    """The reported mean weights each segment by its length."""
    called = evidence.call_cc_segments(
        [(1, 30, 0.6), (100, 159, 0.9)], threshold=0.5, min_length=21, max_gap=2
    )
    expected = (30 * 0.6 + 60 * 0.9) / 90
    assert called.mean_prob == pytest.approx(expected)


def test_sensitivity_grid_reports_one_row_per_combination():
    """The sensitivity table covers the full parameter grid."""
    grid = evidence.cc_sensitivity(
        {"p1": [(1, 25, 0.6)], "p2": [(1, 15, 0.9)]},
        thresholds=[0.2, 0.5],
        min_lengths=[14, 21, 28],
        max_gap=2,
    )
    assert len(grid) == 6
    row = grid[(grid["threshold"] == 0.5) & (grid["min_length"] == 14)]
    assert int(row["n_proteins_with_cc"].iloc[0]) == 2


# ---------------------------------------------------------------------------
# consensus and the SP/TM artefact
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "policy,first,second,expected",
    [
        ("union", True, False, True),
        ("union", False, False, False),
        ("intersection", True, False, False),
        ("intersection", True, True, True),
        ("first", False, True, False),
        ("second", True, False, False),
        ("intersection", None, True, True),
        ("union", None, None, False),
    ],
)
def test_apply_policy(policy, first, second, expected):
    """Consensus policies behave as documented, including for missing channels."""
    assert evidence.apply_policy(policy, first, second) is expected


def test_helix_fully_inside_the_signal_peptide_is_discarded():
    """A helix covered by the signal peptide is a classic false positive."""
    kept, dropped = evidence.filter_helices_in_signal(
        [(3, 22)], sp_end=30, fraction=0.5
    )
    assert kept == [] and dropped == 1


def test_helix_outside_the_signal_peptide_is_kept():
    """A genuine internal helix survives the filter."""
    kept, dropped = evidence.filter_helices_in_signal(
        [(396, 419)], sp_end=30, fraction=0.5
    )
    assert kept == [(396, 419)] and dropped == 0


def test_helix_partially_overlapping_the_signal_peptide_is_kept_below_the_fraction():
    """Only helices mostly inside the signal peptide are removed."""
    kept, _ = evidence.filter_helices_in_signal([(25, 44)], sp_end=30, fraction=0.5)
    assert kept == [(25, 44)]


def test_no_signal_peptide_means_no_filtering():
    """Without a signal peptide every helix is kept."""
    kept, dropped = evidence.filter_helices_in_signal(
        [(1, 20)], sp_end=None, fraction=0.5
    )
    assert kept == [(1, 20)] and dropped == 0
