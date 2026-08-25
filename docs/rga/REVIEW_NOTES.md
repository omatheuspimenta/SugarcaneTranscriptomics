# Second-pass review — RGA prediction pipeline

Written after the implementation worked and the reference run completed, by re-reading
every module as a reviewer looking for reasons to reject it. Each section states what was
checked, what was changed as a result, and what risk remains.

---

## 1. Rule mutual exclusivity and ordering

**Checked.** Whether two rules can fire on the same protein, and whether the specification's
priority order is load-bearing.

**What I found.** The specification's rule list, taken literally, overlaps in several
places: `TIR + CC + NB-ARC + LRR` matches both `CNL` and `TNL`; `RPW8` without NB-ARC
matches both `RX` and (via LRR + TM) `LRR-RLP`; a TIR protein with a kinase matches both
`TX` and `other-RLK`. Relying on `break`-at-first-match would have hidden all of this.

**What I changed.** Every rule from priority 1 to 17 carries an explicit `none_of` set that
makes it disjoint from every other rule **independently of order**. Disjointness is then
*proved*, not asserted: `rules.find_overlapping_rules()` enumerates all 2⁹ = 512 subsets of
the feature vocabulary and reports every pair of non-fallback rules that fire on the same
one. `rgas_prediction.run()` calls `assert_mutually_exclusive()` before reading a single
input file, and `test_rules.py` runs the same check plus a second test asserting that
ordered evaluation assigns exactly one class to each of the 512 combinations.

Two of these disjointness decisions are real biological choices and are documented in
`README.md` §5.2 rather than buried in YAML: **TIR outranks RPW8 outranks CC** for NLR
N-terminal domains, and **TIR outranks kinase** so a TIR-kinase without NB-ARC is `TX`,
not an RLK.

**Residual risk.** The proof covers the feature vocabulary, not the *config*: someone who
adds a rule without the necessary `none_of` entries will get a loud `AssertionError` at
startup — which is the intended behaviour, but they must read it rather than delete the
check.

## 2. The `if`/`elif` leakage bug in the legacy code

**Checked.** The legacy failure mode where a protein matches an early `if`, is classified,
and then falls into a later independent `if` block that overwrites the call.

**Status.** Structurally impossible here. Classification is a single `for` loop over an
ordered rule list with `break` on the first match (`rules.classify_features`), and the
returned `Call` is immutable (`@dataclass(frozen=True)`). No later code path assigns
`rga_subclass`. The legacy script never reached classification at all — it stopped at the
feature matrix — so there was no working behaviour to preserve; it is kept verbatim as
`code/rgas_prediction_legacy.py`.

## 3. Interval coordinates and off-by-one errors

**Checked.** Every coordinate convention, empirically, against the real files, before any
interval was merged.

| Source | Convention | Verified by |
|---|---|---|
| InterProScan TSV cols 7/8 | 1-based inclusive | inspection of real rows |
| DeepTMHMM GFF3 cols 3/4 | 1-based inclusive; regions tile the protein exactly (`1..n`) | `inside 1 460` for a 460-residue protein |
| Phobius topology | 1-based inclusive | `o396-419i` matches DeepTMHMM's `TMhelix 396 416` |
| SignalP `CS pos: 30-31` | signal peptide is residues **1–30** | DeepTMHMM independently emits `signal 1 30` for the same protein |
| DeepCoil2 `.out` | row *i* after the header is residue *i* | `test_deepcoil_positions_are_one_based` |

**What I changed.** `_overlap()` is written as `min(ends) - max(starts) + 1`, so two
intervals sharing exactly one residue have overlap 1, and `merge_intervals` therefore does
**not** merge `(1,10)` with `(11,20)`. Both behaviours are pinned by tests
(`test_merge_intervals_collapses_overlaps`, `test_adjacent_intervals_are_not_merged`).
Segment lengths are consistently `end - start + 1`, tested at the `min_length` boundary in
both directions (20 residues rejected, 21 accepted).

**Residual risk.** If a future tool emits 0-based or half-open coordinates it will be
silently off by one. The convention is documented at the top of `parsers.py` and a new
parser must state which convention it converts from.

## 4. pandas correctness

**Checked.** Merge cardinality, silent NaN coercion, and `itertuples` behaviour.

- **`merge` cardinality.** The pipeline contains **no** `pd.merge` call. Joining six
  300k-row tables on a string key was rejected in favour of building one dictionary per
  protein: it is cheaper, and it makes "one output row per input protein" a property of
  the loop rather than something a join has to be trusted to preserve. There is therefore
  no `validate=` to add. The equivalent guarantee is enforced directly by
  `assert_invariants()`: `len(output) == len(unique input IDs)`, no duplicated
  `protein_id`, and `set(output ids) == set(input ids)`.
- **Duplicate keys.** `_lookup_frame()` originally used
  `frame.set_index(key).to_dict(orient="index")`, which raises on a duplicated protein ID.
  A tool reporting the same protein twice is a realistic input, so it now drops duplicates
  (keeping the first) and logs a warning; `test_duplicated_protein_ids_do_not_duplicate_output_rows`
  covers it end to end.
- **`itertuples` renames invalid identifiers.** *This was a real bug.* The evidence table
  has columns such as `feat_NB-ARC`, which is not a valid Python identifier, so
  `itertuples` silently renamed it to a positional `_23`. The result was that `NB-ARC` was
  missing from every `reason` string and from the confidence calculation, even though the
  classification itself was right. The wide DataFrame is gone; the rule engine iterates the
  list of dictionaries the evidence layer already builds. This also halved peak memory,
  which mattered on a 15 GB machine with ~2 GB free.
- **Silent NaN coercion.** `read_csv` is called with `dtype=str, na_filter=False` for the
  InterProScan TSV, so the literal `-` in the score column is never turned into `NaN` and
  a protein named `NA` would survive. Every `int()` conversion of a possibly-missing value
  is guarded by `pd.isna`. The output table is passed through `_tidy_dtypes()`, which casts
  counts to nullable `Int64` — otherwise a column holding `None` for a handful of proteins
  is promoted to float and a length of 447 is written as `447.0` — and rounds probabilities
  (`localization_prob` was being written as `0.6050000190734863`, a float32→float64
  artefact from DeepLoc). `test_counts_are_written_as_integers` pins this.
- **Missing values.** Written as the literal string `NA` everywhere, never blank, never
  `NaN`. The tests read the outputs back with `keep_default_na=False` so they would notice
  if that regressed.

## 5. Bugs found and fixed during the review

| # | Bug | Effect | Fix |
|---|---|---|---|
| 1 | `itertuples` renamed `feat_NB-ARC` | NB-ARC absent from all reasons and from confidence grading | iterate dictionaries, drop the wide frame |
| 2 | Integrated-domain scan treated **every** configured accession as canonical | a kinase or LysM fused to an NLR — a textbook integrated domain — was never flagged | new `integrated_domain_canonical_features: [NB-ARC, TIR, RPW8, CC, LRR]`; STTK and LysM deliberately excluded |
| 3 | `_lookup_frame` raised on duplicated protein IDs | crash on a realistic input | de-duplicate with a warning |
| 4 | `n_lrr` merged region-level and repeat-level signatures together | `--min-lrr-copies > 1` was a near-useless filter, since Gene3D spans the whole LRR region and collapses everything to one interval | added `n_lrr_repeats` from repeat-level analyses only; documented the caveat in the config and README |
| 5 | Reason strings cited the InterProScan Coils interval as provenance for a CC that DeepCoil2 had actually called | misleading trace (`CC [Coil @ 39-59]` while the call was `70-91`) | `_feature_label()` now cites the channel that made the call, for CC, TM and SP alike |
| 6 | Counts written as floats | `447.0`, `0.6050000190734863` | `_tidy_dtypes()` |
| 7 | Non-RGA reason read "no positive domain evidence" for a protein carrying a coiled coil | confusing | the fallback rule now names the features that *were* found |
| 8 | `locus_summary` used `groupby(...).apply(...)` over ~195k groups | **the process was killed with no traceback** on the real proteome, silently truncating the run after `accession_audit.tsv`: no locus table, no report, no `run_metadata.json` | rewritten as a single linear pass over Python tuples; peak RSS for that step fell to 426 MB and the run completes |
| 9 | The integrated-domain scan treated every non-feature Pfam domain as integrated | the structural sub-domains of an ordinary NLR — the NB-ARC winged-helix domain (`PF23559`, 3,407 NLR hits), the Rx-type N-terminal coiled coil (`PF18052`, 2,611) and the family-specific LRR models — were flagged on 3,716 of 4,023 NLRs, making the column meaningless | new `integrated_domain_exclusions` list, audited accession by accession against the R570 data; the textbook integrated domains (WRKY, BED, HMA, kinase, jacalin) are deliberately not excluded |

The confidence formula originally contained a `single_database_domain` demotion. It was
removed after the review: with the shipped accession list NB-ARC is only ever reported by
Pfam, so the rule fired on **every** NLR and carried no information. The count is still
reported as the informational column `defining_domain_databases`.

## 6. Determinism and reproducibility

**Checked.** Whether two runs on the same inputs can differ.

Every sort is `kind="stable"`; sets are always materialised through `sorted()` before they
reach an output; `iter_deepcoil_sources()` yields directories and archives in sorted order;
the parallel DeepCoil reader collects results and then sorts the segment table by
`(deepcoil_id, start, end)`, so the worker count cannot change the output. Dictionary
iteration order is never relied on for anything that reaches a file.

`run_metadata.json` records the fully resolved configuration, the SHA-256 of every input,
the environment and every count.

**Caught here.** Bug 8 above was found by this pass and not by the tests: the synthetic
fixture has 9 proteins over 8 loci, so the pandas `groupby.apply` path never came close to
the memory ceiling that killed it on 194,593 real loci. The lesson is recorded rather than
papered over — a green test suite said nothing about whether the run would finish, and only
reading the real run's log next to the list of files it actually produced revealed that the
last four outputs were missing. Output-file existence is now asserted by
`test_all_outputs_are_written`, but that test would still have passed on synthetic data;
the real guard is checking the real run.

**Residual risk.** The DeepCoil cache is keyed by path, not by checksum. Pointing
`--input-dir` at a *different* DeepCoil dataset while reusing an existing `--outdir` would
silently reuse the old cache. Mitigation for now is documentation plus
`--refresh-deepcoil-cache`; a checksum-keyed cache would be better and is left as a known
limitation rather than a silent hazard.

## 7. Cross-output consistency

**Checked.** Whether `report.html`, `report.md`, `rga_summary_counts.tsv` and
`run_metadata.json` can disagree.

Every count comes from a single `Summary` object built once by `report.summarize()`; no
output recomputes anything. `assert_report_consistency()` then verifies that the TSV totals,
the JSON `counts` block and the `Summary` agree before anything is written, and
`test_counts_agree_across_outputs` checks the same across all four files, including the
prose sentence in the Markdown report.

## 8. Graceful degradation

**Checked.** Each optional tool removed in turn.

The pipeline runs with any subset of the optional tools; only InterProScan is required. A
missing channel produces a prominent warning, `NA` values, an entry in
`evidence_tools_available`, a per-protein `warnings` entry, and a two-level confidence
demotion for any class that depended on it. Two end-to-end tests cover DeepCoil2 absent
(CC falls back to Coils, calls become `low`) and DeepLoc absent.

`apply_policy()` treats a missing channel as "no evidence" for `union` and skips it for
`intersection`, so an absent tool can never silently veto a call.

## 9. Style and structure

Type hints everywhere, numpydoc docstrings on every public function, PEP 8, no accession
or threshold in the Python source (`test_config.py` covers the validation that keeps the
configuration authoritative).

Four functions still exceed 60 lines of code once docstrings are discounted:
`_prediction_record` (66), `_build_context` (66), `build_parser` (62) and `_cc_fields` (61).
All four are flat literal construction — a dict of output columns, an argparse
declaration — with no branching. Splitting them would add indirection without reducing
complexity, so they are left as they are deliberately.

## 10. Remaining risks, stated plainly

1. **The default CC channel materially understates CNLs.** `PF18052`/`IPR041118`, the Rx
   N-terminal coiled-coil domain, is carried by 2,611 of the 4,023 R570 NLRs but is not
   used as CC evidence, because the specified channels are DeepCoil2 and InterProScan
   Coils. The effect is quantified in `README.md` §4.4 (CNL 396 → 1,790 → 2,648 depending
   on policy and whether the Rx domain is enabled). This is the single most consequential
   open decision in the pipeline and it is the user's to make.
2. **`TM-CC` is a screening bucket, not a class.** Its highest-confidence R570 member is a
   VAMP/synaptobrevin R-SNARE. Documented in the limitations section.
3. **The DeepCoil cache is path-keyed** (§6).
4. **The InterProScan release is not recoverable** from the TSV — only the run date. Left
   as an explicit `TODO` in the README rather than guessed.
5. **`other-RLP` is unreachable** with the default ectodomain list. Intentional, tested,
   and reported with a count of 0.
6. **Logging starts after the configuration is loaded**, so a malformed configuration
   raises to stderr without reaching `logs/run.log`. Harmless, but worth knowing when
   debugging a config error.
7. **The SP/TM overlap filter uses the most permissive `sp_end`** across all tools. If one
   tool over-calls a long signal peptide, a genuine N-terminal helix can be discarded. The
   raw counts and `n_tm_dropped_in_sp` are reported so the effect is always visible; in
   R570 the filter removed helices from a small minority of proteins.

## 11. Verification performed

- 146 pytest tests, all passing: one synthetic fixture per rule (all 19), the exhaustive
  512-combination exclusivity proof, every documented edge case (no hits, MobiDBLite only,
  TM inside the SP, conflicting Phobius/DeepTMHMM, ID present in InterProScan but missing
  from DeepLoc, duplicated IDs, overlapping LRR hits from three databases, CC by DeepCoil2
  but not Coils and vice versa, CC just below and just above `min_length`, two segments at
  the `max_gap` boundary in both directions, CC overlapping a TM helix, CC C-terminal to
  the NB-ARC, DeepCoil2 absent entirely), and cross-output count agreement.
- A full run on the real R570 data: 299,731 proteins, 0 unmatched identifiers, all
  invariants holding.
- Five predictions traced by hand from the output back to the raw tool files (one NLR, one
  RLK, one RLP, one TM-CC, one Non-RGA with partial evidence); all five reconcile exactly.
- All 16 DOIs resolved against the Crossref API, and the bioRxiv preprint against the
  bioRxiv API, on 2026-08-24; author, year, journal, volume and pages checked field by
  field against the Crossref record (see §12).

---

## 12. Third pass — documentation review (2026-08-24)

**Checked.** The state of the tree after the interactive edits that added `rich` output,
plus every reference, re-verified rather than trusted.

| # | Finding | Effect | Fix |
|---|---|---|---|
| 10 | `code/rgas_prediction.py` did not parse: ten section-banner comments had lost their `# ` prefix, leaving bare `---------` lines at module level | `SyntaxError` — the pipeline could not run at all | banners restored; `ast.parse` now checked for every file in `code/` |
| 11 | `code/rgas_prediction_legacy.py` was absent | the promise to preserve the original script was broken | restored byte-for-byte (414 lines, 16,545 bytes) from the session transcript and re-verified with `ast.parse` |
| 12 | `setup_logging()` called `root.handlers.clear()` and installed a plain `StreamHandler`, discarding the `RichHandler` installed at import | the `rich` console output was dead code at run time | the console handler is now a `RichHandler` bound to a shared `Console`; `run.log` keeps the plain timestamped formatter and is verified free of ANSI escapes |
| 13 | `RichHandler(markup=True)` | log messages carry bracketed accession lists (`['PF00931', ...]`), which rich parses as style tags and **deletes** | `markup=False`, with the reason recorded next to the setting |
| 14 | `_make_progress()` was never called | a progress bar that never appeared during the ~20-minute parse | every long stage now reports through a `ProgressCallback` (`rga/progress.py`) that the CLI binds to a `_make_progress()` bar with the `tracked()` context manager: InterProScan by bytes consumed, DeepCoil2 by archive, evidence and classification by protein, writing by file |
| 15 | The bibliography cited neither InterPro nor Pfam | every feature call rests on an accession issued by those resources; the tool papers do not cover the databases | added Blum et al. 2025 (`10.1093/nar/gkae1082`) and Paysan-Lafosse et al. 2025 (`10.1093/nar/gkae997`); Resistify's now-published volume and Rody's full title corrected |
| 16 | The byte-counting wrapper behind the InterProScan bar counted only `read` | pandas wraps a binary handle in `io.TextIOWrapper`, which calls **`read1`** — so the bar would have sat at 0 % for the entire parse | `read1` counted too; `test_progress.py` asserts the reported bytes equal the file size, which is what caught it |

**Determinism re-verified.** The pipeline was re-run from scratch after these changes and
every data table compared by MD5 against the previous run: `rga_predictions.tsv`,
`rga_summary_counts.tsv`, `accession_audit.tsv`, `cc_policy_sensitivity.tsv`,
`rga_predictions_by_locus.tsv` and `rga_domain_evidence_long.tsv` are **byte-identical**.
`report.md` differs only in the generation timestamp, the organism name and the recorded
command line. Counts unchanged: 299,731 proteins, 29,151 RGA candidates.

**What was checked and found correct.** Docstring coverage is complete — every module,
class and function in `code/rga/` and `code/rgas_prediction.py` carries a docstring
(checked programmatically with `ast`, not by eye). The 512-combination exclusivity proof,
all 103 pre-existing tests and every invariant still pass. Thirty-three further tests
were added in `code/tests/test_docs.py` — every source file parses, everything carries a
docstring, the two bibliographies agree in both directions, and every configured rule is
described in the README — plus `code/tests/test_progress.py`, which asserts that every
long stage announces a total and that its advances add up to it. The suite stands at
**146**.
