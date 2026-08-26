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
`README.md` §6.2 rather than buried in YAML: **TIR outranks RPW8 outranks CC** for NLR
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
`code/rgas/rgas_prediction_legacy.py`.

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

1. ~~**The default CC channel materially understates CNLs.**~~ **RETIRED in config
   v0.0.1 — see §13.** The paragraph below is kept as the record of the open question as
   it stood after the second pass; the decision it asked for has since been made.
   **The default CC channel materially understates CNLs.** `PF18052`/`IPR041118`, the Rx
   N-terminal coiled-coil domain, is carried by 2,611 of the 4,023 R570 NLRs but is not
   used as CC evidence, because the specified channels are DeepCoil2 and InterProScan
   Coils. The effect is quantified in `README.md` §4.4 of the two-channel configuration (CNL 396 → 1,790 → 2,648 depending
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

- 185 pytest tests, all passing: one synthetic fixture per rule (all 19), the exhaustive
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
| 10 | `code/rgas/rgas_prediction.py` did not parse: ten section-banner comments had lost their `# ` prefix, leaving bare `---------` lines at module level | `SyntaxError` — the pipeline could not run at all | banners restored; `ast.parse` now checked for every file in `code/` |
| 11 | `code/rgas/rgas_prediction_legacy.py` was absent | the promise to preserve the original script was broken | restored byte-for-byte (414 lines, 16,545 bytes) from the session transcript and re-verified with `ast.parse` |
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
class and function in `code/rgas/rga/` and `code/rgas/rgas_prediction.py` carries a docstring
(checked programmatically with `ast`, not by eye). The 512-combination exclusivity proof,
all 103 pre-existing tests and every invariant still pass. Thirty-three further tests
were added in `code/rgas/tests/test_docs.py` — every source file parses, everything carries a
docstring, the two bibliographies agree in both directions, and every configured rule is
described in the README — plus `code/rgas/tests/test_progress.py`, which asserts that every
long stage announces a total and that its advances add up to it. The suite stood at **146** after this pass.


---

## 13. Fourth pass — the coiled-coil evidence channel (2026-08-26)

**Checked.** The one open question §10 left on the table: whether the CC channel should
rest on a propensity predictor at all. The answer came out of the sources rather than out
of preference, so the sources are recorded first.

### 13.1 What the literature actually says

| Question | What was found | Where |
|---|---|---|
| What threshold does DeepCoil2 recommend? | **None.** The documentation defines `cc` only as "sharpened coiled coil propensity". The paper reports AUC/ROC and F1 — threshold-free metrics — and names a cut-off exactly once: *"a very strict cut-off of 0.9"*, used to mine the **human genome** for high-confidence novel coiled coils | <https://github.com/labstructbioinf/DeepCoil>; Ludwiczak et al. 2019 |
| Is there a defensible minimum length? | Yes. Coiled coils under three heptads are generally unstable, and the PDB-wide benchmark applied exactly two cut-offs, **14 and 21 residues**, reporting better performance at 21 — with *"the increase in sensitivity comes to the cost of precision, which decreases by 5–17%"* | Simm et al. 2021 |
| How good are the predictors? | *"the MCC indicates random prediction in case of NCOILS (MCC of 0.02) and close to random prediction for all other tools (MCC of 0.22 for MultiCoil2 being the highest value)"*, with a 30-fold spread in how many coiled coils each calls | Simm et al. 2021 |
| Was DeepCoil in that benchmark? | **No** — excluded for a 500-residue input cap and a mis-prediction at the first IQ motif of myosin-X. So it has no independent structural evaluation either way | Simm et al. 2021 |

So `min_length: 21` was already grounded and is now cited as such; `threshold: 0.5` had no
support and is now labelled in the configuration as a deliberate midpoint choice rather
than left to read as a default.

### 13.2 The measurement that decided it

Computed from the cached raw DeepCoil2 segments and the reference predictions, within the
4,023 R570 NLRs:

| | NLRs |
|---|---|
| Rx N-terminal domain (`PF18052`) | 2,611 |
| InterProScan Coils | 1,836 |
| DeepCoil2 (0.5 / 21) | 479 |
| Rx **and** DeepCoil2 | 377 |
| DeepCoil2 with neither Rx nor Coils | 33 |

**79 % of DeepCoil2's NLR calls were proteins that already carried the domain model**, at
roughly 14 % recall — and **2,003 of the 3,038 `NL` calls (66 %) carried `PF18052`**, i.e.
two thirds of the largest NLR class were proteins holding the coiled coil of an Rx-type
CNL, reported as having none.

Scored against the domain model as an independent reference, tightening the threshold does
not improve the calls — precision is flat within two points from 0.4 to 0.7 (82.8 % → 84.2 %)
while recall falls from 18.6 % to 2.7 %. Because classification is a **total partition**,
each `CNL` lost this way becomes a positive `NL` rather than an abstention: the error moves
instead of shrinking. At the authors' own 0.9, R570 yields 476 CC-positive proteins in the
entire proteome and `CNL` = 0 — the maximum plateau score anywhere in R570 is 0.922.

### 13.3 What changed

| # | Change | Effect |
|---|---|---|
| 17 | `PF18052`/`IPR041118` promoted from `watch_accessions` to a CC evidence channel of their own (`cc_domain_accessions`), entering as the reserved pseudo-feature `CC_domain` so the nine-feature vocabulary and the 2⁹ proof are untouched | `CNL` 396 → **2,648**; `NL` 3,038 → 786; `CN` 79 → 358; `N` 470 → 191 |
| 18 | `policies.cc` default `deepcoil` → `union`; the policy vocabulary gains `rx_domain` and `union`/`intersection` now range over three channels | `cc_policy_sensitivity.tsv` gains a fifth column. **The two-channel and three-channel `union` columns are not comparable** — hence the config version in `run_metadata.json` |
| 19 | New `cc_rx_domain` output column; `cc_source` now names every contributing channel (`rx_domain+deepcoil+coils`), keeping the `_only` spelling for a single channel | Provenance is visible per protein instead of inferred from the consensus |
| 20 | Confidence demotions rewritten as boolean channel predicates carrying `cc_rx_domain: false`; `_demotion_applies` gained an explicit `_BOOLEAN_WHEN_KEYS` whitelist so a typo in the config raises instead of matching everything | A domain-backed CC is never demoted. `CNL`: 1,617 `high` / 782 `medium` / 249 `low` |
| 21 | **Provenance bug found while implementing this.** `cc_intervals` holds whichever channel called the protein, and the long evidence table exported it as a DeepCoil2 hit. With three channels that would have attributed domain-model spans to DeepCoil2 — the same class of error as finding 5. A separate `deepcoil_intervals` column now feeds the long table | `rga_domain_evidence_long.tsv` 548,086 → 516,723 rows (518,059 after finding 24 added two SMART LRR models): the ~31.4 k removed rows are Coils spans that had been mislabelled as DeepCoil2 output since the two-channel configuration, partly offset by the new domain hits |

| 22 | The domain-level accessions were absent from `accession_audit.tsv`, which iterates `interproscan_features` and `watch_accessions` only | The accessions now driving the `CNL` count would have been the sole evidence accessions missing from the audit table the documentation calls exhaustive | `accession_audit()` gained a pass over `cc_domain_accessions`, reported under the feature label `CC (domain channel)` (2,827 hits each) |

### 13.4 The cost, stated rather than buried

Making `union` the default promotes InterProScan Coils — the channel Simm et al. benchmark
at MCC 0.02 — from cross-check to full channel. `TM-CC`, defined by the two least specific
features in the vocabulary, absorbs it: **3,960 → 8,105**, against 11 under `rx_domain`
alone. The confidence machinery flags this correctly without being asked to: **5,934 of the
8,105 (73 %) are `low`**, 420 `high`. `CNL` is barely affected by the same choice
(2,328 under `rx_domain`, 2,648 under `union`), so this is in practice a `TM-CC` trade-off,
and setting `policies.cc: rx_domain` reverses it in one line.

### 13.5 Verification

- **176 tests pass** (146 + 30). The new `code/rgas/tests/test_cc_channels.py` covers the policy
  function over all five policies including the missing-channel fallbacks, a domain-only CNL
  end to end, the absence of a demotion behind it, the predictor-only demotions still
  firing, `cc_source` provenance, coordinate precedence, the long-table attribution of
  finding 21, the five-column sensitivity table, and the two configuration guards (an
  accession may not feed two CC channels; an empty list restores the DeepCoil2-only behaviour).
- The bibliography-consistency test caught the missing Simm et al. entry in
  `report.REFERENCES` before it reached a report — which is what that test is for.
- The implementation **reproduced the projection made in the two-channel documentation
  exactly**: §4.4 of that README predicted `CNL` 2,648 / `NL` 786 / `CN` 358 / `N` 191 /
  `TM-CC` 8,105 for "union + Rx N-terminal domain", and the run produced those five numbers.
  That projection had been computed by a separate ad-hoc route, so the agreement is a real
  cross-check on the new channel rather than a restatement.
- Invariants hold on the full R570 run: 299,731 proteins, 0 unmatched identifiers,
  4,023 NLRs (unchanged, as expected — the CC channel partitions the NLRs, it does not
  create them), 33,296 RGA candidates.
- The DOI for Simm et al. 2021 was verified against Crossref on 2026-08-26 (title, four
  authors, Scientific Reports, volume 11, article 12439, 2021-06-14).

### 13.6 Risk retired, risks remaining

Risk 1 of §10 — "the default CC channel materially understates CNLs" — is **retired**: it
was the single most consequential open decision in the pipeline and it has now been made,
with the evidence recorded above. What remains in its place is narrower and is documented
in README §9: `PF18052` models Rx/Gpa2-type CNLs specifically, so the 1,006 R570 NLRs
positive for no CC channel stay CC-negative and `CNL` remains a floor. Risks 2–7 are
unchanged, except that risk 2 (`TM-CC` as a screening bucket) is now larger in absolute
terms and better instrumented.


---

## 14. Fifth pass — accession audit and output review (2026-08-26)

**Checked.** Every configured accession against the data and against InterPro; the
prediction table against what a reader actually needs from it; the HTML report against
what it claims to be.

### 14.1 Accessions — what was verified, and how

All 58 configured accessions were cross-checked against the signature description
InterProScan itself reported for them in the R570 run (the authoritative record of what
the run saw), and the five that never fire were resolved against the InterPro REST API.
Two defects and two non-defects came out of it.

| # | Finding | Effect | Fix |
|---|---|---|---|
| 23 | **`PF05729` (NACHT) was listed as `NB-ARC` evidence.** NACHT is a distinct NTPase domain — InterPro keeps `IPR007111` (NACHT) separate from `IPR002182` (NB-ARC) — and this file already said so, in a `nacht_accessions` key that no code read while `PF05729` sat in the NB-ARC list | 0 hits in R570, so nothing in the reference run; but in a fungal or animal proteome, where NACHT NLRs are the norm, NACHT proteins would have been reported as carrying NB-ARC. A portability defect, and P1 makes portability a promise | moved to `watch_accessions` with the reasoning; the dead `nacht_accessions` key deleted; pinned by `test_nacht_is_not_nb_arc_evidence` |
| 24 | **`SM00370` is a retired accession.** SMART deleted it on 2011-11-22 (the InterPro API answers HTTP 410, `deletion_date` 2011-11-22); it could never have fired. Meanwhile its two live siblings present in the data, `SM00364` (LRR, bacterial type) and `SM00368` (LRR, ribonuclease-inhibitor type), were not configured | none on classification — all 212 proteins carrying `SM00364`/`SM00368` already had LRR evidence from another signature, which is why this was invisible in the counts | `SM00370` removed, `SM00364` and `SM00368` added; `n_lrr_repeats` is now correct for those 212 proteins |

**Two candidate gaps were measured and rejected**, recorded so nobody re-derives them:

- `IPR038005` ("Virus X resistance protein-like, coiled-coil domain", 1,711 proteins)
  looked like a fourth CC channel. It adds **2 proteins** beyond `PF18052`/`IPR041118`,
  moving exactly one protein from `NL` to `CNL`. Not worth a channel.
- `PF25019` and `PF23622` are LRR domain models that sit in
  `integrated_domain_exclusions` but not in the `LRR` feature list, which looked
  asymmetric. They would rescue **0** NLRs: every protein carrying them already has LRR
  evidence, and all 582 NLRs currently without LRR (`CN`, `N`, `TN`) carry neither.

**Verified correct and left alone.** `SM00367` is `LRR_CC_2` — a *cysteine*-containing
LRR, correctly assigned to `LRR`, and the very example §4 of this file uses to explain why
description matching fails. `G3DSA:3.80.10.10` is described as "Ribonuclease Inhibitor",
the archetypal LRR fold, and is correctly treated as a region-level LRR signature.
`IPR000727` ("Target SNARE coiled-coil homology domain", 570 hits) is correctly *not* CC
evidence — it is the domain behind the VAMP/synaptobrevin false positive in §10 risk 2.
`PF05659` and `IPR008808` are the RPW8 models: InterPro types them `family` and `domain`
respectively, and `IPR008808` ("Powdery mildew resistance protein, RPW8 domain") is
domain-level, so the §12.4 argument that RPW8 evidence is domain-level holds.

### 14.2 The prediction table

`rga_predictions.tsv` was already the most detailed output in the pipeline, but everything
about *where* a domain sits and *which signature* found it lived only in the prose of
`reason`. Three columns make that machine-readable, and one closes a gap in the counts:

| # | Column | Why |
|---|---|---|
| 25 | `feature_coords` | `NB-ARC:191-361;LRR:595-1086,615-638` — the coordinates behind the call, filterable |
| 26 | `feature_accessions` | `NB-ARC:PF00931;LRR:G3DSA:3.80.10.10,PF00560` — the provenance behind the call |
| 27 | `integrated_domain_descriptions` | `Protein kinase domain` instead of `PF00069` |
| 28 | `n_tm_consensus` | the helix count behind the `tm_consensus` boolean, which had per-tool counts but no consensus count |

**A bug caught while adding them.** The first implementation put `feature_coords` and
`feature_accessions` into `_INTEGER_COLUMNS` — the tuple lives directly below
`PREDICTION_COLUMNS` and shares its entries, so a mechanical edit hit both. `_tidy_dtypes`
then cast both string columns to nullable `Int64` and every value came out `NA`. Caught by
the new tests before the reference run, not after.

Both columns use `FEATURE:value;FEATURE:value`, and Gene3D accessions contain colons
(`G3DSA:3.80.10.10`), so they must be split on the **first** colon only.
`test_feature_accessions_survive_accessions_containing_a_colon` pins it.

### 14.3 The HTML report

The page was already self-contained, responsive and script-free. What it lacked was
navigation and any way to see that a large class can be mostly untrustworthy.

| # | Change | Why |
|---|---|---|
| 29 | Table of contents, anchored `<h2>`s with self links, generated from one `_SECTIONS` tuple | Eleven sections with no way to jump between them. One source for the headings and the contents list, so they cannot drift |
| 30 | New section 6, **"How much of each class is trustworthy"**: a subclass × confidence cross-table in both HTML and Markdown | The most important thing this run has to say about `TM-CC` is that 5,934 of its 8,105 members are `low`. That was computable from the TSV and invisible in the report |
| 31 | A caveat callout below the first heading | The per-protein/per-locus distinction, the `--cc-policy` dependence and the `TM-CC` caveat now appear *before* the numbers, not in a limitations section of a different document |
| 32 | Sticky table headers, scroll-capped long tables, print stylesheet | Ordinary usability. The report is read on screen and printed |
| 33 | `_CC_NOTE` and `_METHODS` rewritten | **Both still described the two-channel design and called DeepCoil2 "the primary coiled-coil channel".** The generated report was contradicting the configuration it was generated from — the worst kind of stale documentation, because it ships with the result |

No JavaScript was added: the report loads no external stylesheet, script, font or image,
and that property is now worth more than a sortable table.

### 14.4 Two documented figures were wrong

Re-measuring rather than re-quoting turned up two stale performance claims, both of which
a reader would have acted on:

| # | Claim | Measured | Note |
|---|---|---|---|
| 34 | "about twenty minutes on a laptop" | **72 s** warm, 87 s with a cold DeepCoil2 cache | the 20-minute figure predates the cache and the chunked parse |
| 35 | "peak RSS ≈ 426 MB" | **2.84 GB** | 426 MB was real but misattributed: it was the peak of the *locus-summary step* after bug 8, not of the run. Identical at `--workers 1`, so it is the parent process, not the pool |

The 426 MB number was being used to size a machine, so it is corrected in place with a
note saying what it actually measured, rather than silently replaced.

`--workers 1` and `--workers 6` produce byte-identical outputs, which re-confirms the
determinism claim from a direction the earlier passes had not tested.

### 14.5 Sixth pass — report presentation (2026-08-26)

| # | Change | Why |
|---|---|---|
| 36 | **The report is light-only.** The `prefers-color-scheme: dark` block added in finding 32 was removed and `color-scheme: light` declared instead | The report is a printable record: it gets screenshotted into slides, pasted into theses and printed on paper. A page that changes appearance with the reader's desktop theme is the wrong thing for that job, and the dark variant was never going to be the one that reached a printer. `color-scheme: light` also stops the browser auto-darkening the scrollbars and controls around it |
| 37 | **A "Reproduce this run" block**, in both the HTML and the Markdown report, wrapping the invocation one flag per line | The person reading the report is the person who needs to re-run it, and the command was previously visible only as an unquoted one-liner in the metadata card |
| 38 | `run_metadata.json` records `shlex.join(...)`, not `" ".join(sys.argv)` | The old form dropped the quoting: `--organism-name Saccharum officinarum x spontaneum cv. R570` is five arguments, not one. The recorded command was **not the command that ran**, and pasting it back would have produced a differently-named organism or an argparse error |
| 39 | `main(argv)` records the argv it was given, and the script names itself from a constant rather than `sys.argv[0]` | `sys.argv` is the *host* program's command line whenever `main` is called programmatically — from a test, a notebook, or the Nextflow wrapper this is headed for. That is exactly the situation in which a misleading recorded command does damage |

Verified by running the command the report prints, verbatim, into a scratch directory:
it completes and reproduces `rga_predictions.tsv` and `rga_summary_counts.tsv`
byte-identically. Two tests were added — one asserting an organism name containing
spaces round-trips through `shlex.split` and reaches the Markdown report quoted, one
asserting the page carries no `prefers-color-scheme` — plus a third in `test_docs.py`
pinning the example command in README §3.1 against the command recorded in
`run_metadata.json`, so the documented invocation cannot drift from the one that produced
the numbers printed beside it.

### 14.6 Verification

- **185 tests pass** (176 + 9). New: the data dictionary is pinned to
  `PREDICTION_COLUMNS` in both directions, so a column can no longer be added without
  being documented or documented after being renamed away; the colon-splitting contract;
  the NACHT assignment; and the rule that no accession is both evidence and watch-only.
- The reference run was regenerated. **Counts are unchanged** — 299,731 proteins,
  33,296 RGAs, 4,023 NLRs, `CNL` 2,648 — which is the expected result: findings 23 and 24
  were both verified to be no-ops on R570 before they were applied.
- `accession_audit.tsv` grew from 56 to 60 rows: `PF18052`/`IPR041118` under
  `CC (domain channel)` (finding 22), `PF05729`/`IPR007111` under `(watch only)`, and
  `SM00370` replaced by `SM00364`/`SM00368`.
- All internal `§` references and anchor links in `README.md` were re-resolved
  programmatically after the new §3 shifted every downstream section by one, and the
  cross-file references in `ARCHITECTURE.md`, `REVIEW_NOTES.md`, the root `README.md` and
  `rga_config.yaml` were repointed with them.
