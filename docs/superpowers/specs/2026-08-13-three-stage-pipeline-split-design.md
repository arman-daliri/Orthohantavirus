# Three-Stage Pipeline Split — Design Spec

**Date:** 2026-08-13
**Branch:** `main`
**Status:** approved design, pending implementation plan

## 1. Purpose

The repo currently has two entry points that don't read as a pipeline: `all_hanta.py`, whose
`--phylogenetics` flag conditionally runs a ~730-line MAFFT+IQ-TREE analysis inline, and
`reassortment_detector.py`, a fully separate tool for scoring a brand-new candidate genome. A
user cloning the repo has to read code to learn that these are related, in what order to run
them, or that they're independent of each other.

This spec splits the phylogenetics analysis out of `all_hanta.py` into its own script and
documents all three as an explicit, numbered pipeline:

1. **`all_hanta.py`** — ingest the NCBI package, build accession/isolate/species tables, the ML
   feature matrix, and `ortho_ml_sequence_table.csv`.
2. **`phylogeny.py`** *(new)* — real MAFFT+IQ-TREE reassortment analysis over the reference
   isolate set (what `--phylogenetics` does today).
3. **`reassortment_detector.py`** — score a new, not-yet-published candidate genome against the
   5-gate framework (unchanged).

No new capability is added. This is a structural refactor plus documentation, verified to be
behavior-preserving by keeping every existing test passing throughout.

## 2. Architecture

Each stage keeps its own CLI (`argparse`) and can be run independently against the same raw NCBI
package root (`--input`). None of the three stages reads another stage's output files as a
required input — this was evaluated and rejected (§4) because it would add a false dependency
without saving real work, given what each stage's candidate-selection logic actually needs.

```
all_hanta.py --input <pkg> --output <out>          # Stage 1, always run first for the tables
phylogeny.py --input <pkg> --output <out>           # Stage 2, independent, needs mafft/iqtree
reassortment_detector.py --input <pkg> --output <out> --candidate ...   # Stage 3, independent
```

`environment.yml` (already pinned, see the reproducibility work from 2026-08-13) supplies
`mafft`/`iqtree` for Stages 2 and 3; Stage 1 has no phylogenetics dependency at all after this
split — a change from today, where `all_hanta.py` itself declares the `--phylogenetics` argument
group even though most runs never use it.

## 3. Component changes

### 3.1 `all_hanta.py` (Stage 1)

Loses the following, moving to `phylogeny.py`. This list was produced by grepping every
occurrence of each symbol in the current file, not by re-reading the phylogenetics section and
guessing at its boundaries — two of the entries below sit well outside the main contiguous block
and would be easy to miss without that check.

**The main block — `MIN_CDS_LENGTH` (line 2329) through `_phylo_call_by_isolate` (line 3093),
contiguous, 24 definitions, verified zero gaps:**

`MIN_CDS_LENGTH`, `is_human_affecting_or_likely`, `parse_cds_fasta`, `select_primary_cds`,
`collect_phylo_candidates`, `_phylogroup_for_isolate`, `_phylogroup_for_candidate`,
`_assignment_basis`, `_candidate_group_key`, `group_candidates`, `group_candidates_by_species`,
`resolve_binary`, `write_segment_fasta`, `run_mafft`, `run_iqtree`,
`_patristic_distance_matrix`, `tree_divergence_for_species`, `_nearest_supported_descendant`,
`bootstrap_support_for_pair`, `classify_reassortment_call`, `_phylo_row`,
`_resolve_outgroup_label`, `_analyze_group`, `_run_phylogenetics_scope`,
`run_phylogenetics_stage`, `_phylo_call_by_isolate`.

**Two outliers, defined outside that range, moving anyway because their only call site is inside
it** (each confirmed by a full-file grep — neither has any other caller):

- `PhyloCandidate` (line 733, a `@dataclass`) — used throughout the block above.
- `HOST_ORDER_PHYLOGROUP_FALLBACK` (line 451, a module-level dict) — used once, by
  `_phylogroup_for_isolate` (line 2520).
- `_segment_rank_metrics` (line 2028) — used once, by `tree_divergence_for_species` (line 2691).

**One shared dependency that does NOT move — `_rank_positions` (line 1997).** It's called both
by Stage 1's own k-mer divergence path (line 2153, `compute_segment_divergence`) and by the
moving `tree_divergence_for_species` (line 2691). It stays in `all_hanta.py`; `phylogeny.py`
calls `ah._rank_positions`. This is the only cross-module call `phylogeny.py` makes back into
`all_hanta.py` beyond the standard ingest primitives every stage already uses
(`resolve_inputs`, `parse_fasta`, `iter_jsonl`, `build_isolate_records`, `IsolateRecord`, etc. —
unchanged, still live in `all_hanta.py`, still used by Stage 2 and Stage 3 exactly as before).

**CLI surface:** `all_hanta.py` loses the `--phylogenetics`, `--mafft-bin`, `--iqtree-bin`,
`--phylo-min-taxa`, `--bootstrap`, `--min-bootstrap-support`, `--phylo-seed`, `--phylo-threads`,
`--phylo-scope`, `--phylo-outgroup` arguments and the `phylo_call_by_key` / try-except block in
`run_pipeline()` that invokes the phylogenetics stage. `ortho_risk_table.csv` is now always built
phylo-blind — this is byte-identical to today's behavior when `--phylogenetics` is *not* passed,
which is already the common case, so it's a no-op for most existing users. `all_hanta.py` also
stops writing `ortho_segment_divergence.csv` and `ortho_ml_segments.csv` (§5).

### 3.2 `phylogeny.py` (new, Stage 2)

**Correction from plan-writing (2026-08-13):** the original version of this section said Stage
3's import fix touched five symbols. Re-deriving the check found a sixth:
`detector_reference.py` line 40 calls `ah.parse_cds_fasta`, which is also moving. The corrected,
complete list is six symbols: `resolve_binary`, `PhyloCandidate`, `select_primary_cds`,
`MIN_CDS_LENGTH`, `_analyze_group`, and `parse_cds_fasta`. §3.3 below is corrected accordingly.

Receives everything listed in §3.1. Gets its own `parse_args`/`main`/pipeline-runner, taking the
`--phylogenetics`-era CLI flags (minus the `--phylogenetics` toggle itself, since running this
script *is* the toggle) plus `--input`/`--output`/`--verbose` to match the other two stages'
conventions. Re-derives isolates, the genomic/protein FASTA index, and the CDS index from
`--input` independently — the same ~10 lines of ingest-calling code `run_pipeline()` already runs
before touching accessions. This small duplication was a deliberate choice (§4), not an oversight.

Writes `ortho_reassortment_table.csv` and `phylogenetics/<species>/<segment>/` — same names,
same schemas, same content as today's `--phylogenetics` output. The withdrawn-k-mer-classifier
warning (`"Reassortment is NOT called by this pipeline..."`, logged at run_pipeline's
divergence-computation step, line 4015) stays in `all_hanta.py` unchanged — it's warning about
Stage 1's own k-mer diagnostic fields, which are still present in `ortho_isolate_table.csv` and
`ortho_ml_feature_matrix.csv` post-split, not about anything moving. `phylogeny.py` gets its own
new top-level "stage starting"/"stage skipped" log messages in its `main()`, replacing the
`if args.phylogenetics:` branch messages (lines 4067/4080) that don't carry over verbatim since
there's no longer a conditional branch to log around — running the script *is* the condition.

### 3.3 `reassortment_detector.py` / `detector_gates.py` / `detector_reference.py` (Stage 3)

No logic changes. Import-surface fix only: these three files currently do `import all_hanta as
ah` and call `ah.resolve_binary` (`reassortment_detector.py`), `ah.PhyloCandidate`,
`ah.select_primary_cds`, `ah.MIN_CDS_LENGTH`, `ah._analyze_group` (`detector_gates.py`), and
`ah.parse_cds_fasta` (`detector_reference.py`, line 40) — all six are moving. Every other `ah.`
reference in these files (`IsolateRecord`, `FastaRecord`, `kmer_profile`,
`cosine_distance_matrix`, `profiles_to_matrix`, `accession_base`, `parse_fasta`, `iter_jsonl`,
`build_annotation_index`, `build_biosample_index`, `transform_one_record`,
`build_isolate_records`, `resolve_inputs`, `configure_logging`, `SEGMENT_LENGTH_RANGE`) is not
moving and needs no change.

Fix: add `import phylogeny as ph` to all three files; repoint the six call sites above from
`ah.` to `ph.`.

## 4. Rejected alternative: stage-to-stage file dependencies

Considered having `phylogeny.py` and `reassortment_detector.py` consume Stage 1's output CSVs
(`ortho_ml_sequence_table.csv` or `ortho_isolate_table.csv`) instead of re-parsing the raw
package. Rejected because it doesn't actually reduce work: Stage 2's candidate selection
(`human_affecting_or_likely` / `divergence_candidate` / `complete_triplet` tags) runs over the
**full** reconstructed isolate population, not the 206 rows `ortho_ml_sequence_table.csv` is
deliberately filtered down to (complete-triplet, non-lab, non-vaccine only), and neither of
Stage 1's tables carries raw FASTA sequences for the full population. Stage 3's Gate 4 comparator
building has the identical requirement. Consuming Stage 1's output would add a hard run-order
dependency while Stage 2/3 still had to open the raw FASTA regardless — no benefit, real cost.

## 5. Output file changes

- **`ortho_segment_divergence.csv` — removed.** Every column (`sl_rank_separation`,
  `m_rank_offset`, `reassortment_status`, k-mer nearest-neighbor) already exists in
  `ortho_isolate_table.csv` and `ortho_ml_feature_matrix.csv`. It was residue from when the k-mer
  divergence diagnostic was the only reassortment-adjacent signal available; now that Stage 2 is
  a real, independently-documented pipeline stage, the standalone file adds nothing.
- **`ortho_ml_segments.csv` — removed.** Column-subset of `ortho_accession_table.csv`
  (same fields, narrower row set), fully derivable from it.
- **`ortho_risk_table.csv` — unconditionally phylo-blind.** `phylo_reassortment_call` and
  `priority_reassortment_term` are always empty/`0.0` now (previously only when
  `--phylogenetics` wasn't passed). Combining with Stage 2's evidence is a join on `isolate_key`
  against `ortho_reassortment_table.csv` — already how the README documents the fallback case
  today; now the only case. No schema change to `ortho_risk_table.csv` itself.
- Everything else is unchanged.

`output/` is git-tracked (checked, not assumed — 5 prior commits regenerate it deliberately,
e.g. "Regenerate default-pipeline outputs for the lab-derived ML flag change"), and the current
committed copy was generated *without* `--phylogenetics` (verified: no `ortho_reassortment_table.csv`
or `phylogenetics/` present in it today). So Stage 1's regenerated output after this split should
be identical to the current commit's, minus the two removed files — no other content should
change. The implementation plan re-runs `all_hanta.py` and re-commits `output/` as its own step,
following the repo's existing convention, with a diff review confirming only the two files
disappear and nothing else shifts.

## 6. Testing strategy

Phylogenetics-related test classes move from `test_all_hanta.py` to a new `test_phylogeny.py`,
following the repo's existing one-test-file-per-module convention
(`test_detector_gates.py` next to `detector_gates.py`). Verified by class-name inspection of
`test_all_hanta.py`, moving: `ClassifyReassortmentCallTest`, `BootstrapSupportRootFallbackTest`,
`CollectPhyloCandidatesScopeTest`, `GroupCandidatesTest`, `PhylogeneticsScopeEndToEndTest`,
`PhyloCallByIsolateTest`. **`BuildRiskRowsPhyloIntegrationTest` is deleted, not moved** — reading
its body during plan-writing showed it tests `build_risk_rows`'s `phylo_call_by_key` parameter
directly (`ah.build_risk_rows(isolates, {}, {"iso1": "..."})`), which is being removed from
`build_risk_rows` entirely per this section's "unconditionally phylo-blind" decision. The
behavior it tests no longer exists in either module, so it isn't a move candidate. Staying in
`test_all_hanta.py`:
`SegmentRankMetricsBackwardCompatTest` (tests `_rank_positions`/`_sl_m_rank_metrics`/
`_segment_rank_metrics` against each other — `_rank_positions` and `_sl_m_rank_metrics` stay in
`all_hanta.py`; this test needs `ah._segment_rank_metrics` repointed to `ph._segment_rank_metrics`
but otherwise stays put since it's fundamentally a Stage-1 regression test), `SpeciesEffortStatsTest`,
`SurveillancePriorityReassortmentTermTest` (tests that `build_risk_rows` handles an empty
`phylo_call_by_key` correctly — still a Stage-1 concern post-split).

`test_detector_gates.py` / `test_detector_reference.py` / `test_reassortment_detector.py` need no
assertion changes — checked for `mock.patch()`/string-based patch targets against any moving
symbol and found none (every dependency is a live attribute access through real code paths, not a
mock), so nothing there references the `ah.`/`ph.` split by string path.

**Automated verification, not manual review:** after the move, a one-time check asserts none of
the 26 moved names remain accessible as `all_hanta.<name>` and all 26 are accessible as
`phylogeny.<name>`. Combined with the full existing suite passing (currently 96 tests, unchanged
count expected — this is a move, not new coverage), a missed import-surface fix cannot pass
silently: every moved symbol is exercised by a real (unmocked) code path already in the suite, so
a stale `ah.X` reference fails immediately with `AttributeError` the first time that test runs.

## 7. README changes

Replace the current `### Optional: --phylogenetics` subsection under Stage 1's docs with a new
top-level "Pipeline stages" section listing Stage 1/2/3, what each needs (`--input`, plus
`mafft`/`iqtree` on PATH for Stages 2 and 3), and what each produces. The `Outputs` table drops
the two removed files and gets a one-line note on `ortho_risk_table.csv`'s row that phylogenetic
evidence requires a manual join with Stage 2's `ortho_reassortment_table.csv` on `isolate_key`.

## 8. Out of scope

- No new shared module or package restructuring beyond what's listed above — `all_hanta.py`
  keeps its ingest primitives, `phylogeny.py` keeps the tree-analysis primitives, nothing is
  further split into utility modules. If the plan phase finds the ~10-line ingest duplication
  between `all_hanta.py` and `phylogeny.py` (§3.2) is larger than expected in practice, revisit
  then — not speculatively now.
- No change to any output file's column schema, only which files exist.
- No change to Stage 3's behavior, tests, or CLI at all beyond the `ah.`/`ph.` import fix.
- No CI/automation pipeline wiring the three stages together into one command — each stays an
  independent, separately-invoked script, per the approved input-contract decision (§4).
