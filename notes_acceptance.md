# WS4 pipeline acceptance — final step

## SUMMARY

- Commit at start of session: `7ddd3bf` (Add pipeline.py, unpushed, branch
  `04-rpi-fitness-validation`).
- Fix commit (this session): `2cb0a70` — adds `"groups": scored_groups` to
  `pipeline._stage_report`'s payload and updates `StageReportTest` to assert on it.
- Verdict: **PASS**. Per-row scores were already byte-identical between
  `pipeline._stage_discriminate` and `rpi_discrimination.py`'s own panel build (step 1
  below). The wrapper gap was that `pipeline.py run`'s written report only ever included
  `disc.summarise_groups(...)`, never the per-row rows themselves — fixed and confirmed
  by a second seeded run (step 2 below).

---

## 1. Diff: `output/pipeline_acceptance/raw_scored_groups.json` vs `output/discrimination_scratch/discrimination_panel.json`

Compared the `"groups"` section of both files (the discrimination_panel.json has an
extra `"summary"` key the raw capture doesn't; that is expected and not part of this
diff). Same seed (20260829) in both.

- Group keys match: `a_random_fakes`, `b1_within_species_shuffle`,
  `b2_cross_species_sympatric`, `c_real_positive`.
- Group sizes match: 20 / 30 / 36 / 11.
- Every row, re-serialised with `json.dumps(row, indent=2, default=str)` (key order
  preserved from the original JSON load, so this catches ordering differences, not just
  value equality) is **byte-identical** between the two files, row for row, in order.
  **0 diffs across 97 rows.**

**Verdict: PASS — byte-identical.**

## 2. Fix applied: `pipeline.py run` must emit full per-row scores

`pipeline._stage_report` (pipeline.py:99-111) built its written report from
`disc.summarise_groups(scored_groups)` only — four small per-group aggregates (n,
n_scored, mean/min/max of `rpi_generation`) — and never wrote `scored_groups` itself.
`output/pipeline_acceptance/run_20260829.json` (the prior real run) confirms this: it
has `tau`, `group_summary`, `pass_criterion`, no per-row data at all.

**Fix:** added `"groups": scored_groups` to the payload dict in `_stage_report`, and
updated `StageReportTest.test_writes_valid_json_with_tau_summary_and_pass_criterion` to
assert `written["groups"] == scored_groups`. No other code path changed;
`_stage_discriminate` (which already produces correct per-row scores, confirmed in step
1) was untouched.

`python -m unittest pipeline` — 12 tests, all pass, run after the fix.

## 3. Rerun to confirm

Ran `python pipeline.py run --seed 20260829 --output output/pipeline_acceptance/run_confirm_20260830.json`
end to end (calibrate [skipped, artifact already present] -> discriminate -> stats ->
report), same seed as before.

- The written report now has top-level keys `tau`, `group_summary`, `groups`,
  `pass_criterion`.
- `report["groups"]` re-diffed against `output/discrimination_scratch/discrimination_panel.json`'s
  `groups` the same way as step 1: group sizes match (20/30/36/11), **0 diffs across all
  rows**.
- Console summary matched the prior run's group sizes and `overall_ranking_supported =
  False` (WS4 Task 13/14's known, unresolved-interpretation result — untouched by this
  fix, not in scope here).

**Verdict: PASS — confirmed.**

## 4. Committed

`git add pipeline.py notes_acceptance.md` and committed with message describing the
wrapper-gap fix. `output/pipeline_acceptance/*` and `notes_diagnosis.md` are left
uncommitted scratch/notes per the task instructions (not asked to commit those).
