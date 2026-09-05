# Overnight full-pipeline rerun, reconciliation, and cleanup survey — design

**Written 2026-09-06. Branch `05-module-consolidation`.**

## 1. Purpose

Rerun the full pipeline (`all_hanta.py` → `phylogeny.py` → `rpi_core.py` +
`Algorithm.py`) from scratch into an isolated output folder, without letting the
existing `output/` bias the rerun, then reconcile the two and promote the
reconciled result. Separately: illustrate the phylogeny trees, survey (but do not
delete) redundant/junk artifacts, and tighten the report and docs. This also
absorbs last session's stalled gate4-confirmation run (§7f of
`docs/PROJECT-STATE.md`), which never completed (4/10 winners checkpointed, no
final JSON, driver process no longer running).

## 2. Non-goals

- No recalibration. `data/rpi_calibration.json` is reused as pinned — this run
  tests whether the pipeline reproduces itself deterministically downstream of
  calibration, not whether recalibrating changes anything.
- No autonomous code simplification. Simplification proposals for
  `all_hanta.py`/`phylogeny.py`/`rpi_core.py`/`Algorithm.py` are written up as
  proposals only, never applied without review.
- No deletion. The cleanup step produces a candidate list with reasons; nothing
  is removed tonight.
- No push. Everything lands as local commits only.

## 3. Ordered steps

1. **Baseline.** Working tree confirmed clean apart from 4 untracked stale files
   from the stalled gate4 run (`gate4_confirmation.log{,.stdout,.stderr}`,
   `gate4_confirmation_checkpoint.jsonl` under
   `output/ga_search_corrected_pool_long/`) — these stay untracked, go on the
   step-6 redundancy list, and are not committed as-is.

2. **Pre-run code review.** Review `all_hanta.py`, `phylogeny.py`, `rpi_core.py`,
   `Algorithm.py` for correctness/bias issues that would affect tonight's run.
   Score-neutral fixes (guard-tested) are applied and committed individually,
   following the existing code-review-fix-boundary rule. Anything that would
   change a reported number is written up as a proposal in the final report,
   not applied — this step is governed by the standing rule, not tonight's
   reconciliation override (see §4).

3. **Full rerun into `output_rerun_20260906/`** (isolated from `output/`, nothing
   promoted yet):
   - `all_hanta.py` ingest → accession/isolate/species/host tables, ML matrices.
   - `phylogeny.py` phylogenetics stage (both the plain and phylogroup scopes,
     matching `output/phylogenetics` and `output/phylogenetics_phylogroup`).
   - `rpi_core.py run` (reuses pinned calibration; calibrate → discriminate →
     stats → report).
   - `rpi_core.py validate --with-phylogeny --json .../validation_report.json`
     — all 5 gates including gate4 on the validation panel.
   - `Algorithm.py` GA search, both configs matching the existing runs
     (`ga_search_corrected_pool_default`, `ga_search_corrected_pool_long`).
   - `confirm_winners(..., with_phylogeny=True, checkpoint_path=...)` on the new
     run's 10 long-pool winners — the expensive step (~25+ min per
     hantanense-heavy winner observed last session; this is what makes the run
     overnight-scale). Checkpointed so a crash loses nothing already confirmed.
   - `figures.py` before/after rendering against the new run's trees.

4. **Compare** `output_rerun_20260906/` against `output/`: checksum diff first
   (same method as `ortho_output_checksums.txt`), then a semantic diff (values,
   not just bytes) on anything not byte-identical, since IQ-TREE/GA search
   could in principle differ run-to-run even under a fixed seed.

5. **Reconcile every real difference.** Root-cause it. Fix in code/config with
   full autonomy for *this reconciliation only* (guard-tested, one commit per
   fix) — this is the one place tonight's run overrides the standing
   score-neutral-only boundary. Log root cause + fix in the report. Where a fix
   changes a number, the corrected rerun's numbers are canonical.

6. **Promote.** Identical case: trivial promotion of
   `output_rerun_20260906/` over `output/`. Conflict case: promote the
   reconciled (corrected) result — same action, because step 5 already made it
   correct. Delete the now-redundant `output_rerun_20260906/` staging folder
   after promotion (it is intentionally short-lived; this is the one deletion
   tonight, since it is scaffolding, not project output).

7. **Redundancy survey.** Scan `output/`, repo root, `docs/` for redundant,
   unnecessary, or junk files/folders/code not needed even as evidence. Produce
   a list with reasons in the report. Includes at minimum: the 4 orphaned
   stale gate4 files from step 1, superseded scratch/driver scripts, stray
   `__pycache__`. No deletions.

8. **Report.** Written to close out the run: the 10 gate4 verdicts for the
   long-pool winners, one line each, verbatim (`no_conflict` and
   `insufficient_comparators` included per the binding reporting rule already
   set in `docs/PROJECT-STATE.md` §7f — 9/10 winners are within-hantanense
   combinations, so `no_conflict` there is gate4 doing its job, not a
   limitation to explain away); anything the stalled confirmation run covered
   that the standard rerun does not (stated explicitly, never dropped
   silently); every step-5 reconciliation with problem/fix; the step-7
   redundancy list; step-2 proposals not applied. `docs/PROJECT-STATE.md` §7f
   gets its verdict writeup from this report the same morning, not during the
   unattended run itself.

9. **Docs tightening.** README and `docs/PROJECT-STATE.md` tightened for the
   new state — de-duplicated, no content cut that documents a real limitation
   or decision.

10. **Commit, no push.** Everything committed locally in the ordered-step
    commits above plus a final commit for the report/docs. `git push` is not
    run.

## 4. Execution mechanics

Runs as a detached background process so it survives the session across the
overnight window, using checkpoint files for crash recovery (matching last
session's approach for the long IQ-TREE step). If a shell tool gets blocked at
any point, the run stops and the user is told, rather than working around it.
If something materially bigger than a numeric drift turns up mid-run (a real
correctness bug, not just an expected recalculation), the run pauses there and
flags it rather than pushing through.
