# Single-File Consolidation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Finish branch `05-module-consolidation` by removing the six Phase-2
shims, the two generator scripts, and the scripts outside the six-module set
that still reference old filenames — leaving `rpi_core.py` as the one real
file for the RPI toolchain, with `all_hanta.py` and `phylogeny.py` untouched.

**Architecture:** `rpi_core.py` already contains every line of logic from the
six merged modules plus `_pipeline_main`'s subcommand dispatcher. This plan
(a) gives `rpi_core.py` its own `__main__` guard, (b) adds two subcommands
the dispatcher is currently missing (`calibrate`, `score`), (c) inlines the
two permanent diagnostics that currently live in separate files reached only
via lazy imports, (d) verifies nothing broke, then (e) deletes the eleven
now-redundant files and updates the docs that reference them.

**Tech Stack:** Python 3.14 stdlib (`argparse`, `unittest`, `unittest.mock`).
No new dependencies.

**Spec:** `docs/superpowers/specs/2026-08-31-single-file-consolidation-design.md`

## Global Constraints

- Zero change to any RPI scoring logic (τ, gates, weights, factors) — this is
  a file-layout and CLI-surface change only.
- `all_hanta.py` and `phylogeny.py` are out of scope; do not touch them.
- `python -m unittest rpi_core` must report exactly 322 tests, green, at
  every checkpoint from Task 3 onward.
- Docs (`CLAUDE.md`, `README.md`, `notes_consolidation.md`) are updated only
  in Task 5, after Task 3's gate passes and Task 4's deletions are committed.
- Every file deletion happens in one commit (Task 4); every doc edit happens
  in a separate, later commit (Task 5). Do not interleave them.

---

### Task 1: Add `__main__` guard, `calibrate`, and `score` subcommands to `rpi_core.py`

**Files:**
- Modify: `rpi_core.py:9315` (`COMMANDS` tuple)
- Modify: `rpi_core.py:9420-9438` (`_pipeline_main` dispatch body)
- Modify: `rpi_core.py` end of file (after line 9587, the last line of
  `MainDispatchTest`) — add `__main__` guard
- Modify: `rpi_core.py:9555-9587` (`MainDispatchTest`) — two new test methods

**Interfaces:**
- Consumes: `_calibration_main(argv=None) -> int` (already defined at
  `rpi_core.py:7154`), `_detector_main(argv=None) -> int` (already defined at
  `rpi_core.py:3248`), `_pipeline_main(argv=None) -> int` (already defined,
  being modified here).
- Produces: `python rpi_core.py calibrate ...` and `python rpi_core.py score
  ...` become valid invocations; `python rpi_core.py` (no args) still
  defaults to `run` (unchanged, `COMMANDS` membership check falls through to
  the `else` branch).

- [ ] **Step 1: Write the two failing dispatch tests**

Add to `MainDispatchTest` in `rpi_core.py`, immediately after
`test_validate_command_passes_through_remaining_args` (currently ending at
line 9572):

```python
    def test_calibrate_command_passes_through_remaining_args(self):
        with mock.patch(__name__ + "._calibration_main", return_value=0) as fake_calibrate:
            exit_code = _pipeline_main(["calibrate", "--stability"])
        fake_calibrate.assert_called_once_with(["--stability"])
        self.assertEqual(exit_code, 0)

    def test_score_command_passes_through_remaining_args(self):
        with mock.patch(__name__ + "._detector_main", return_value=0) as fake_detector:
            exit_code = _pipeline_main(["score", "--reference-input", ".", "--candidate", "x"])
        fake_detector.assert_called_once_with(["--reference-input", ".", "--candidate", "x"])
        self.assertEqual(exit_code, 0)
```

- [ ] **Step 2: Run the new tests to verify they fail**

Run: `python -m unittest rpi_core.MainDispatchTest -v`
Expected: `test_calibrate_command_passes_through_remaining_args` and
`test_score_command_passes_through_remaining_args` FAIL (dispatcher doesn't
recognize `"calibrate"`/`"score"` yet, so they fall into the `run` branch and
`_calibration_main`/`_detector_main` are never called).

- [ ] **Step 3: Update `COMMANDS` and the dispatch body**

At `rpi_core.py:9315`, change:

```python
COMMANDS = ("run", "validate", "diagnose-e", "diagnose-gate0")
```

to:

```python
COMMANDS = ("run", "validate", "calibrate", "score", "diagnose-e", "diagnose-gate0")
```

At `rpi_core.py:9420-9438`, change the body of `_pipeline_main` from:

```python
def _pipeline_main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in COMMANDS:
        command, rest = argv[0], argv[1:]
    else:
        command, rest = "run", argv

    if command == "validate":
        return _validate_main(rest)
    if command == "diagnose-e":
        import diagnose_e_stratification as diag_e
        diag_e.main()
        return 0
    if command == "diagnose-gate0":
        import diagnose_gate0_junk as diag_gate0
        diag_gate0.main()
        return 0
    return run(_parse_run_args(rest))
```

to:

```python
def _pipeline_main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in COMMANDS:
        command, rest = argv[0], argv[1:]
    else:
        command, rest = "run", argv

    if command == "validate":
        return _validate_main(rest)
    if command == "calibrate":
        return _calibration_main(rest)
    if command == "score":
        return _detector_main(rest)
    if command == "diagnose-e":
        import diagnose_e_stratification as diag_e
        diag_e.main()
        return 0
    if command == "diagnose-gate0":
        import diagnose_gate0_junk as diag_gate0
        diag_gate0.main()
        return 0
    return run(_parse_run_args(rest))
```

(The `diagnose-e`/`diagnose-gate0` branches are rewritten in Task 2, not
here — keep this step scoped to `calibrate`/`score` only, so Step 4's test
run isolates one change at a time.)

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m unittest rpi_core.MainDispatchTest -v`
Expected: all 8 `MainDispatchTest` methods (6 existing + 2 new) PASS.

- [ ] **Step 5: Add the `__main__` guard**

At the very end of `rpi_core.py` (after the last line of
`MainDispatchTest`), add:

```python


if __name__ == "__main__":
    sys.exit(_pipeline_main())
```

- [ ] **Step 6: Smoke-test the guard directly**

Run: `python rpi_core.py calibrate --stability` (expect the stability report
to print, same as `python rpi_calibration.py --stability` today — do not
compare output yet, that's Task 3; just confirm it runs without a traceback)
Run: `python rpi_core.py --help`
Expected: no traceback either way. (`--help` at the top level currently
falls through to `run`'s arg parser via `_parse_run_args`, since `--help`
isn't a recognized `COMMANDS` value — this is pre-existing `pipeline.py`
behavior, not something this task changes; confirm `python pipeline.py
--help` behaves identically before and after by running both.)

- [ ] **Step 7: Run full suite**

Run: `python -m unittest rpi_core`
Expected: 322 tests, green.

- [ ] **Step 8: Commit**

```bash
git add rpi_core.py
git commit -m "rpi_core.py: wire calibrate and score into the CLI dispatcher, add __main__ guard"
```

---

### Task 2: Inline `diagnose_e_stratification.py` and `diagnose_gate0_junk.py` into `rpi_core.py`

**Files:**
- Modify: `rpi_core.py` — add `_diagnose_e_main` and `_diagnose_gate0_main`
  near the end of the pipeline section, before `_pipeline_main`
- Modify: `rpi_core.py:9420-9438` region (post-Task-1) — `diagnose-e` /
  `diagnose-gate0` dispatch branches
- Modify: `rpi_core.py` `MainDispatchTest` — rewrite
  `test_diagnose_e_command_dispatches_with_no_args` and
  `test_diagnose_gate0_command_dispatches_with_no_args`
- Read only (source being moved, not edited): `diagnose_e_stratification.py`,
  `diagnose_gate0_junk.py`

**Interfaces:**
- Consumes: `stats.compare_groups` (already imported/aliased in `rpi_core.py`
  via the `stats = sys.modules[__name__]` self-alias), `rd.load_reference_data`,
  `rd._build_precomputed`, `rd.CandidateInput`, `rd.gate0_check`, `cal.genomic_sequence`,
  `cal.segment_accession` (all already defined in `rpi_core.py` under their
  original names — this task adds no new dependency, only new call sites).
- Produces: `_diagnose_e_main(argv: Optional[Sequence[str]] = None) -> None`,
  `_diagnose_gate0_main() -> None` — module-level functions in `rpi_core.py`.
  Module-level constants `E_THRESHOLD` (float, `0.5`) and `ISOLATE_KEYS`
  (list of 2 strings) also move in under those exact names (grep-confirmed:
  no existing `rpi_core.py` name collides with either).

- [ ] **Step 1: Write the failing tests (rewrite the two existing ones)**

In `MainDispatchTest`, replace:

```python
    def test_diagnose_e_command_dispatches_with_no_args(self):
        import diagnose_e_stratification as diag_e
        with mock.patch.object(diag_e, "main") as fake_diag:
            exit_code = _pipeline_main(["diagnose-e"])
        fake_diag.assert_called_once_with()
        self.assertEqual(exit_code, 0)

    def test_diagnose_gate0_command_dispatches_with_no_args(self):
        import diagnose_gate0_junk as diag_gate0
        with mock.patch.object(diag_gate0, "main") as fake_diag:
            exit_code = _pipeline_main(["diagnose-gate0"])
        fake_diag.assert_called_once_with()
        self.assertEqual(exit_code, 0)
```

with:

```python
    def test_diagnose_e_command_dispatches_with_remaining_args(self):
        with mock.patch(__name__ + "._diagnose_e_main") as fake_diag:
            exit_code = _pipeline_main(["diagnose-e", "some/panel.json"])
        fake_diag.assert_called_once_with(["some/panel.json"])
        self.assertEqual(exit_code, 0)

    def test_diagnose_gate0_command_dispatches_with_no_args(self):
        with mock.patch(__name__ + "._diagnose_gate0_main") as fake_diag:
            exit_code = _pipeline_main(["diagnose-gate0"])
        fake_diag.assert_called_once_with()
        self.assertEqual(exit_code, 0)
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m unittest rpi_core.MainDispatchTest -v`
Expected: both rewritten tests FAIL — `_diagnose_e_main`/`_diagnose_gate0_main`
don't exist yet as `rpi_core` module attributes, so `mock.patch(__name__ +
"._diagnose_e_main")` raises `AttributeError`.

- [ ] **Step 3: Add the two inlined functions to `rpi_core.py`**

Insert immediately before `def _pipeline_main(` (post-Task-1 location, still
in the "pipeline.py -- single orchestrating entry point" section):

```python
E_THRESHOLD = 0.5  # chosen after inspecting the distribution below -- see printed histogram


def _diagnose_e_main(argv: Optional[Sequence[str]] = None) -> None:
    """Standalone diagnostic (WS4 Task 13/14, PROJECT-STATE.md section 7c).

    Loads an already-generated discrimination panel JSON (does NOT regenerate it),
    stratifies groups.c_real_positive by rpi_e_total, and re-runs
    stats.compare_groups for reals(E~1) vs b1 only -- to check whether
    the reals-vs-b1 tie reported by evaluate_pass_criterion is driven by the
    published inter-species lab crosses (which carry a near-zero Factor E
    ecological prior) diluting an otherwise-intra-species real-positive panel.

    Pure read-and-recompute: does not modify any scoring logic, only calls
    the public compare_groups().
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    path = argv[0] if argv else "output/discrimination_scratch/discrimination_panel.json"
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    groups = data["groups"]
    c_real = groups["c_real_positive"]
    b1 = groups["b1_within_species_shuffle"]

    scored_real = [row for row in c_real if row.get("rpi_generation") is not None]
    unscored_real = len(c_real) - len(scored_real)

    e_counts = Counter(row["rpi_e_total"] for row in scored_real)
    print(f"Loaded: {path}")
    print(f"c_real_positive: n={len(c_real)}, scored={len(scored_real)}, "
          f"unscored/junk_sequence={unscored_real}")
    print()
    print("rpi_e_total distribution across scored c_real_positive rows:")
    for value, count in sorted(e_counts.items(), key=lambda kv: kv[0]):
        print(f"  rpi_e_total = {value!r}: {count} case(s)")
    print()
    print(f"Chosen threshold: rpi_e_total >= {E_THRESHOLD} "
          f"(distribution is bimodal -- 1.0 vs 0.01 -- this threshold falls cleanly "
          f"in the gap and does not need tuning)")
    print()

    reals_e1 = [row for row in scored_real if row["rpi_e_total"] >= E_THRESHOLD]
    reals_e_low = [row for row in scored_real if row["rpi_e_total"] < E_THRESHOLD]
    print(f"reals(E~1) stratum: n={len(reals_e1)} -- case_ids: "
          f"{[row['case_id'] for row in reals_e1]}")
    print(f"reals(E-low) stratum (excluded from this comparison): n={len(reals_e_low)} -- "
          f"case_ids: {[row['case_id'] for row in reals_e_low]}")
    print()

    b1_values = [row["rpi_generation"] for row in b1 if row.get("rpi_generation") is not None]
    reals_e1_values = [row["rpi_generation"] for row in reals_e1]

    result = stats.compare_groups(
        reals_e1_values, b1_values,
        "c_real_positive(E>=0.5)", "b1_within_species_shuffle",
        expect="reals(E~1) ~ b1: same D7 tie-expectation as the unstratified "
               "reals_vs_b1 comparison, restricted to the intra-species subset of "
               "reals to check whether Factor E composition (not M_pack/M_pol) drove "
               "the original result.")
    result["plausible_tie_not_proven_equivalent"] = (
        not result["significant_at_alpha_0_05"]
        or result["cliffs_delta_magnitude"] == "negligible")

    print("compare_groups(reals(E>=0.5), b1) result:")
    print(json.dumps(result, indent=2))


ISOLATE_KEYS = [
    "t2|Orthohantavirus sinnombreense|202003197",
    "t2|unclassified Orthohantavirus|msb89863",
]


def _diagnose_gate0_main() -> None:
    """Read-only diagnostic (roadmap Step 1): why REAL_ORTHO_sinnombreense_005
    (isolate_key t2|Orthohantavirus sinnombreense|202003197) and
    REAL_ORTHO_Orthohantavirus_006 (isolate_key t2|unclassified Orthohantavirus|msb89863)
    hit gate0's junk_sequence verdict in the c_real_positive group.

    Does not modify any scoring logic -- calls the public gate0_check /
    _build_precomputed / load_reference_data only.
    """
    reference = rd.load_reference_data(Path("all_hanta"))
    precomputed = rd._build_precomputed(reference)
    by_key = {i.isolate_key: i for i in reference.isolates}

    for key in ISOLATE_KEYS:
        isolate = by_key.get(key)
        print("=" * 100)
        print(f"isolate_key: {key}")
        if isolate is None:
            print("  NOT FOUND in reference.isolates")
            continue
        print(f"  species_name: {isolate.species_name}")
        print(f"  s_accession/m/l: {isolate.s_accession} / {isolate.m_accession} / {isolate.l_accession}")
        print(f"  s/m/l_completeness: {isolate.s_completeness} / {isolate.m_completeness} / {isolate.l_completeness}")
        print(f"  s/m/l_nt_length: {isolate.s_nt_length} / {isolate.m_nt_length} / {isolate.l_nt_length}")
        print(f"  triplet_status: {isolate.triplet_status}, has_complete_triplet: {isolate.has_complete_triplet}")
        print(f"  warnings: {isolate.warnings}")

        genomic = {
            seg: cal.genomic_sequence(reference, cal.segment_accession(isolate, seg))
            for seg in ("S", "M", "L")
        }
        for seg, seq in genomic.items():
            print(f"  {seg} genomic length actually loaded: {len(seq) if seq else 'EMPTY/None'}")

        candidate = rd.CandidateInput(
            candidate_id=key, s_genomic=genomic["S"], m_genomic=genomic["M"], l_genomic=genomic["L"])
        gate0 = rd.gate0_check(
            candidate, precomputed["protein_len_ranges"], precomputed["panhandle_consensus"],
            precomputed["protein_profiles"])
        print(f"  gate0.passed_structure: {gate0.passed_structure}")
        print(f"  gate0.passed_identity_floor: {gate0.passed_identity_floor}")
        print(f"  gate0.best_kmer_identity: {gate0.best_kmer_identity}")
        print("  gate0.failures:")
        for f in gate0.failures:
            print(f"    - {f}")
        print("  gate0.warnings:")
        for w in gate0.warnings:
            print(f"    - {w}")
        for seg, protein in gate0.orf_proteins.items():
            print(f"  {seg} longest_orf_protein length: {len(protein)}")
```

Confirm `Counter` is already imported at the top of `rpi_core.py` (it is —
`from collections import Counter` is part of the deduplicated import block
from Phase 1; if `grep -n "^from collections import Counter" rpi_core.py`
returns nothing, add it to the existing `collections` import line instead of
a new line).

- [ ] **Step 4: Update the dispatch branches**

Change (post-Task-1 state):

```python
    if command == "diagnose-e":
        import diagnose_e_stratification as diag_e
        diag_e.main()
        return 0
    if command == "diagnose-gate0":
        import diagnose_gate0_junk as diag_gate0
        diag_gate0.main()
        return 0
```

to:

```python
    if command == "diagnose-e":
        _diagnose_e_main(rest)
        return 0
    if command == "diagnose-gate0":
        _diagnose_gate0_main()
        return 0
```

- [ ] **Step 5: Run the dispatch tests to verify they pass**

Run: `python -m unittest rpi_core.MainDispatchTest -v`
Expected: all 8 methods PASS.

- [ ] **Step 6: Run full suite**

Run: `python -m unittest rpi_core`
Expected: 322 tests, green (no count change — no test was added or removed,
two were rewritten in place).

- [ ] **Step 7: Commit**

```bash
git add rpi_core.py
git commit -m "rpi_core.py: inline diagnose-e and diagnose-gate0, drop lazy cross-file imports"
```

---

### Task 3: Verification gate — full suite, seeded panel, CLI smoke tests

No code changes in this task. This is the checkpoint the spec requires
before any file deletion (Task 4) or doc update (Task 5) may proceed.

**Files:** none modified.

- [ ] **Step 1: Full suite**

Run: `python -m unittest rpi_core`
Expected: `Ran 322 tests ... OK`. If this fails, stop — do not proceed to
Task 4 until it's green. (Gate 4's tests inside this run are slow, ~13.5 min
total — this is expected, matches the documented baseline.)

- [ ] **Step 2: Seeded discrimination panel, bit-identical to the golden baseline**

Run:
```bash
python rpi_core.py score --help
python -c "
import json, subprocess, sys
subprocess.run([sys.executable, 'rpi_core.py', 'run', '--seed', '20260829',
                '--output', 'output/consolidation_baseline/discrimination_panel_phase5.json'],
               check=True)
"
```

If `output/consolidation_baseline/discrimination_panel_golden.json` from
Phase 0 stores the same shape as `run`'s `--output` report (check its
top-level keys first — `python -c "import json; print(json.load(open('output/consolidation_baseline/discrimination_panel_golden.json')).keys())"`),
diff the `groups` key of the new output against it:

```bash
python -c "
import json
old = json.load(open('output/consolidation_baseline/discrimination_panel_golden.json'))
new = json.load(open('output/consolidation_baseline/discrimination_panel_phase5.json'))
assert old.get('groups', old) == new.get('groups', new), 'MISMATCH -- do not proceed'
print('bit-identical: OK')
"
```

Expected: `bit-identical: OK`. If this prints `MISMATCH`, stop and
investigate before Task 4 — this is the reproducibility guarantee the whole
plan exists to protect.

- [ ] **Step 3: CLI help smoke test for every subcommand**

Run each of these and confirm none tracebacks:
```bash
python rpi_core.py --help
python rpi_core.py run --help
python rpi_core.py calibrate --help
python rpi_core.py validate --help
python rpi_core.py score --help
```

- [ ] **Step 4: Diagnostics smoke test**

Run:
```bash
python rpi_core.py diagnose-gate0
python rpi_core.py diagnose-e output/consolidation_baseline/discrimination_panel_golden.json
```

Expected: both print their report to stdout with no traceback. Compare the
printed `case_ids`/`rpi_e_total` numbers by eye against a run of the
original standalone scripts on the same input (`python
diagnose_gate0_junk.py`, `python diagnose_e_stratification.py
output/consolidation_baseline/discrimination_panel_golden.json` — these
still exist as separate files until Task 4 deletes them) to confirm the
inlined versions produce the same output.

- [ ] **Step 5: Record the gate result**

No commit for this task (nothing changed) — just confirm all four steps
passed before starting Task 4. If any step failed and required a code fix,
that fix belongs in Task 1 or Task 2 (amend those commits' follow-up, not a
new ad hoc commit here) and this task must be re-run from Step 1.

---

### Task 4: Delete the eleven redundant files

**Files:**
- Delete: `reassortment_detector.py`, `rpi_calibration.py`,
  `rpi_discrimination.py`, `rpi_discrimination_stats.py`,
  `validate_reassortment_detector.py`, `pipeline.py`
- Delete: `build_rpi_core.py`, `build_shims.py`
- Delete: `capture_raw_scored_groups.py`, `diagnose_e_stratification.py`,
  `diagnose_gate0_junk.py`

**Interfaces:** none — this task only removes files whose entire public
surface is now reachable through `rpi_core.py` (verified by Task 3).

- [ ] **Step 1: Confirm nothing else in the repo still imports the files being deleted**

Run:
```bash
python -c "
import re, pathlib
targets = ['reassortment_detector', 'rpi_calibration', 'rpi_discrimination',
           'rpi_discrimination_stats', 'validate_reassortment_detector',
           'pipeline', 'build_rpi_core', 'build_shims',
           'capture_raw_scored_groups', 'diagnose_e_stratification',
           'diagnose_gate0_junk']
keep = {'rpi_core.py'}
for f in pathlib.Path('.').glob('*.py'):
    if f.name in keep or f.stem in targets:
        continue
    text = f.read_text(encoding='utf-8')
    for t in targets:
        if re.search(rf'\b{re.escape(t)}\b', text):
            print(f'{f}: references {t}')
"
```

Expected: no output. If anything prints, that file needs updating (or is a
reason to reconsider deleting that target) before proceeding — do not delete
a file something else still imports.

- [ ] **Step 2: Delete the files**

```bash
git rm reassortment_detector.py rpi_calibration.py rpi_discrimination.py \
       rpi_discrimination_stats.py validate_reassortment_detector.py pipeline.py \
       build_rpi_core.py build_shims.py \
       capture_raw_scored_groups.py diagnose_e_stratification.py diagnose_gate0_junk.py
```

- [ ] **Step 3: Re-run the full suite and the seeded-panel check from Task 3, against `rpi_core.py` alone**

Run: `python -m unittest rpi_core`
Expected: 322 tests, green (this run no longer has the shims/generators on
disk at all — proves `rpi_core.py` has zero runtime dependency on any of
them).

Run: `python rpi_core.py diagnose-gate0` and `python rpi_core.py diagnose-e
output/consolidation_baseline/discrimination_panel_golden.json`
Expected: both still work (Task 3 verified them against the shim-era repo;
this confirms they don't secretly still need the now-deleted files).

- [ ] **Step 4: Commit**

```bash
git commit -m "$(cat <<'EOF'
Delete Phase-2 shims, generator scripts, and scratch/diagnostic scripts

rpi_core.py is now the only file for the RPI toolchain -- everything
these eleven files provided is reachable through it (calibrate/score
subcommands added in prior commits, diagnostics inlined, capture
script's one-off output already on disk). Repo root is now three
stage files: all_hanta.py, phylogeny.py, rpi_core.py.
EOF
)"
```

---

### Task 5: Update `CLAUDE.md`, `README.md`, `notes_consolidation.md`

**Files:**
- Modify: `CLAUDE.md:134-143`
- Modify: `README.md:12-51`, `README.md:155-167`, `README.md:904-911`,
  `README.md:976`
- Modify: `notes_consolidation.md` (append new section)

**Interfaces:** none — documentation only.

- [ ] **Step 1: Rewrite `CLAUDE.md`'s "Running things" section**

Replace `CLAUDE.md:134-143`:

```markdown
## Running things

```bash
python -m unittest rpi_calibration                # 39 tests, <1 s
python rpi_calibration.py                         # rebuilds data/rpi_calibration.json (~12 s)
python rpi_calibration.py --stability             # branch across seeds x windows (~2 min)
python -m unittest reassortment_detector          # 194 tests, ~13.5 min (gate 4 is slow)
python validate_reassortment_detector.py --json output/validation/validation_report.json
python reassortment_detector.py --rpi-report ...  # writes rpi_validation_report.json
```
```

with:

```markdown
## Running things

```bash
python -m unittest rpi_core                        # 322 tests, ~13.5 min (gate 4 is slow)
python rpi_core.py                                 # same as `run`: calibrate -> discriminate -> stats -> report
python rpi_core.py run --recalibrate               # force-rebuild data/rpi_calibration.json
python rpi_core.py calibrate                       # rebuild data/rpi_calibration.json directly (~12 s)
python rpi_core.py calibrate --stability           # branch across seeds x windows (~2 min)
python rpi_core.py validate --json output/validation/validation_report.json
python rpi_core.py score --reference-input . --candidate ... --rpi-report
```
```

- [ ] **Step 2: Rewrite README's "Module layout" paragraph and command block**

Replace `README.md:12-51` (the `## RPI fitness-validation toolchain
(\`pipeline.py\`)` section header through the paragraph ending "`rpi_core.py`
itself is a library module, not a CLI entry point.") with:

```markdown
## RPI fitness-validation toolchain (`rpi_core.py`)

**Purpose.** `rpi_core.py` is the single module and CLI entry point over the
Reassortment Plausibility Index (RPI) toolchain built on top of the base
pipeline above: it calibrates `M_pack`'s Boltzmann rescaling against a
held-out isolate panel, scores a four-group discrimination panel (random
fakes, within-species shuffles, cross-species sympatric shuffles, and real
published/supported reassortants), and runs the `reals ~ b1 > b2 > random`
pass criterion (Mann-Whitney U + Cliff's delta) over the result — see
["Reassortment Plausibility Index — limitations"](#reassortment-plausibility-index--limitations)
below for what the index does and does not claim.

**How to run:**

```bash
python rpi_core.py                       # same as `run`: calibrate -> discriminate -> stats -> report
python rpi_core.py run --recalibrate     # force-rebuild data/rpi_calibration.json
python rpi_core.py run --seed 20260829   # override the discrimination panel's seed (default shown)
python rpi_core.py calibrate             # rebuild data/rpi_calibration.json directly
python rpi_core.py calibrate --stability # branch the calibration across seeds x windows
python rpi_core.py validate --json output/validation/validation_report.json
python rpi_core.py score --reference-input . --candidate id:s.fasta,m.fasta,l.fasta
python rpi_core.py diagnose-e            # one-off: is the reals~b1 tie E-driven?
python rpi_core.py diagnose-gate0        # one-off: gate0 junk_sequence false positives
python -m unittest rpi_core              # full suite, 322 tests
```

`run` reuses the pinned calibration artifact (`data/rpi_calibration.json`) unless it is
missing or `--recalibrate` is passed — calibration is a distinct, auditable step, not
something a routine run silently redoes. A stage that raises halts the run there; each
stage's input is the previous stage's output, so nothing downstream runs on a half-built
result.

**Module layout.** `rpi_core.py` is the only file for this toolchain — the
six modules this consolidated from (`reassortment_detector.py`,
`rpi_calibration.py`, `rpi_discrimination.py`, `rpi_discrimination_stats.py`,
`validate_reassortment_detector.py`, `pipeline.py`) and the shims that
briefly stood in for them during the consolidation are gone; see
`notes_consolidation.md` for the full mechanical record of how the merge was
done and verified bit-identical at each phase.
```

- [ ] **Step 3: Update the three-stage table**

Replace `README.md:162`:

```markdown
| 3 | `reassortment_detector.py` | `mafft`/`iqtree` on PATH, Biopython, `dendropy` | Score a new, not-yet-published candidate genome against the 5-gate framework |
```

with:

```markdown
| 3 | `rpi_core.py score` | `mafft`/`iqtree` on PATH, Biopython, `dendropy` | Score a new, not-yet-published candidate genome against the 5-gate framework |
```

- [ ] **Step 4: Update the validation section**

Replace `README.md:904-911`:

```markdown
## Reassortment detector — validation against published reassortants

`validate_reassortment_detector.py` scores a panel whose ground truth is known
independently of the detector, and reports a confusion matrix. Run it with:

```bash
python validate_reassortment_detector.py --json output/validation/validation_report.json
```
```

with:

```markdown
## Reassortment detector — validation against published reassortants

`rpi_core.py validate` scores a panel whose ground truth is known
independently of the detector, and reports a confusion matrix. Run it with:

```bash
python rpi_core.py validate --json output/validation/validation_report.json
```
```

- [ ] **Step 5: Update the `--rpi-report` line**

Replace `README.md:976`:

```markdown
Run `reassortment_detector.py --rpi-report` to write `rpi_validation_report.json` alongside
```

with:

```markdown
Run `rpi_core.py score --rpi-report` to write `rpi_validation_report.json` alongside
```

- [ ] **Step 6: Append a closing section to `notes_consolidation.md`**

Append at the end of the file:

```markdown

## Phase 4 — single-file finish (2026-08-31)

The six Phase-2 shims, the two generator scripts (`build_rpi_core.py`,
`build_shims.py`), and the scripts outside the six-module set
(`capture_raw_scored_groups.py`, `diagnose_e_stratification.py`,
`diagnose_gate0_junk.py`) are deleted. `rpi_core.py` gained a `__main__`
guard and two subcommands the dispatcher was missing: `calibrate`
(→ `_calibration_main`, previously reachable only via `python
rpi_calibration.py` directly) and `score` (→ `_detector_main`, previously
reachable only via `python reassortment_detector.py` directly — this one was
a real gap in the original consolidation, not just a doc nit: nothing in
`pipeline.py`'s `COMMANDS` ever routed to it). The two permanent diagnostics
(`diagnose-e`, `diagnose-gate0`) moved from lazy cross-file imports to
functions defined directly in `rpi_core.py`.

Repo root is now three stage files: `all_hanta.py`, `phylogeny.py`,
`rpi_core.py`. Verification: `python -m unittest rpi_core` — 322 tests,
green; seeded panel (`--seed 20260829`) bit-identical to Phase 0's golden
baseline; all five CLI subcommands' `--help` sane; both diagnostics
smoke-tested and cross-checked against the pre-deletion standalone scripts.
See `docs/superpowers/specs/2026-08-31-single-file-consolidation-design.md`
for the full design.
```

- [ ] **Step 7: Commit**

```bash
git add CLAUDE.md README.md notes_consolidation.md
git commit -m "docs: update CLAUDE.md/README/notes_consolidation.md for the single-file rpi_core.py"
```

---

## Self-Review Notes

- **Spec coverage:** §1 (`__main__` guard) → Task 1 Step 5. §2 (`calibrate`
  command) → Task 1. §2a (`score` command, the mid-plan discovery) → Task 1.
  §3 (inline diagnostics) → Task 2. §4 (file deletions) → Task 4. §5 (test
  updates) → Task 1 Step 1/Task 2 Step 1. §6 (docs last) → Task 5, gated
  behind Task 3/4. Closing gate's four checks → Task 3 Steps 1-4.
- **Placeholder scan:** no TBD/TODO; every step shows the actual diff or
  command, not a description of one.
- **Type consistency:** `_diagnose_e_main(argv: Optional[Sequence[str]] =
  None) -> None` and `_diagnose_gate0_main() -> None` match their call sites
  in `_pipeline_main` (Task 2 Step 4) and their test mocks (Task 2 Step 1)
  exactly — `_diagnose_e_main` takes `rest` (a list), `_diagnose_gate0_main`
  takes nothing, matching the pre-existing external scripts' own signatures
  (`diag_e.main()` had no args in the old dispatch, but the file's `main()`
  itself read `sys.argv[1]`; this plan's inlined version explicitly threads
  `rest` through instead of falling back to `sys.argv`, which is the more
  correct behavior when called as a subcommand — noted here since it's a
  small, deliberate behavior improvement over the literal pre-image, not an
  oversight).
