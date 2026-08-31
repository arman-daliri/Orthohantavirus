# Single-file consolidation — finishing branch `05-module-consolidation`

## Context

`05-module-consolidation` merged six modules' logic into `rpi_core.py` (Phase
1) but kept the six original filenames as thin re-export shims (Phase 2), plus
two generator scripts (`build_rpi_core.py`, `build_shims.py`) and three
outside-the-six scripts that still import the old names
(`capture_raw_scored_groups.py`, `diagnose_e_stratification.py`,
`diagnose_gate0_junk.py`). The result is 12 Python files at repo root, which
the user found confusing when trying to identify the one file that actually
does the work.

This phase removes the shims and generator scripts entirely and folds the two
permanent diagnostics into `rpi_core.py` as CLI subcommands, leaving the repo
root at three stage files: `all_hanta.py`, `phylogeny.py`, `rpi_core.py`.
`capture_raw_scored_groups.py` is deleted outright — confirmed (by reading it)
to be a one-off scratch script whose output is already captured on disk.

**Out of scope:** `all_hanta.py`, `phylogeny.py`, and anything about the
accuracy of the RPI scoring itself (τ, gates, factors). This is a file-layout
and CLI-surface change only.

## What changes

### 1. `rpi_core.py` becomes directly executable

It currently has no `if __name__ == "__main__":` guard (each source file's own
guard was stripped by `build_rpi_core.py`). Add one at the very end:

```python
if __name__ == "__main__":
    sys.exit(_pipeline_main())
```

This makes `python rpi_core.py <command> ...` the one and only CLI entry
point, exactly reproducing what `python pipeline.py ...` does today (the
shim's own `__main__` guard calls the same `_pipeline_main`).

### 2. `calibrate` becomes a top-level command, not just an internal stage

`COMMANDS` (currently `("run", "validate", "diagnose-e", "diagnose-gate0")`)
gains `"calibrate"`, dispatching to the existing `_calibration_main(rest)` —
today only reachable via `python rpi_calibration.py` directly, or implicitly
via `run`'s `_stage_calibrate`. No change to `_calibration_main` itself.

New `COMMANDS`: `("run", "validate", "calibrate", "diagnose-e",
"diagnose-gate0")`.

### 3. The two permanent diagnostics move into `rpi_core.py`

`diagnose_e_stratification.py`'s `main()` (reads a panel JSON path from
`sys.argv[1]`, stratifies by `rpi_e_total`, calls `stats.compare_groups`) and
`diagnose_gate0_junk.py`'s `main()` (hardcoded `ISOLATE_KEYS`, calls
`rd.load_reference_data` / `rd._build_precomputed` / `rd.gate0_check`) move in
**verbatim**, renamed `_diagnose_e_main(argv)` / `_diagnose_gate0_main()` to
avoid colliding with the five existing renamed `_..._main` functions. Their
module-level constants (`E_THRESHOLD`, `ISOLATE_KEYS`) move in unrenamed — no
collision exists in `rpi_core.py` today (checked by grep). No argparse is
added to either; behavior stays exactly what it is now (confirmed: the user
chose "as described, no new argparse" over adding argparse to
`diagnose-gate0`).

`_pipeline_main`'s dispatch table changes from lazy `import
diagnose_e_stratification as diag_e; diag_e.main()` (and the gate0 twin) to
direct calls: `_diagnose_e_main(rest)` / `_diagnose_gate0_main()`. The lazy
imports existed only to break the Phase 2 circular-import (shim importing
`rpi_core` importing a shim); once the diagnostics live inside `rpi_core.py`
itself there is nothing left to import and no cycle.

### 4. Files deleted

- `reassortment_detector.py`, `rpi_calibration.py`, `rpi_discrimination.py`,
  `rpi_discrimination_stats.py`, `validate_reassortment_detector.py`,
  `pipeline.py` — the six shims. Nothing outside `rpi_core.py` imports these
  names once step 5 lands.
- `build_rpi_core.py`, `build_shims.py` — generators with nothing left to
  generate. History (this commit and the Phase 1/2 commits) preserves them if
  ever needed again.
- `capture_raw_scored_groups.py` — one-off scratch script, output already on
  disk at `output/pipeline_acceptance/raw_scored_groups.json`.
- `diagnose_e_stratification.py`, `diagnose_gate0_junk.py` — folded into
  `rpi_core.py` per step 3.

### 5. Test updates inside `rpi_core.py`

`MainDispatchTest.test_diagnose_e_command_dispatches_with_no_args` and
`test_diagnose_gate0_command_dispatches_with_no_args` currently do `import
diagnose_e_stratification as diag_e` / `import diagnose_gate0_junk as
diag_gate0` and mock `.main`. These become:

```python
with mock.patch(__name__ + "._diagnose_e_main", return_value=None) as fake_diag:
    exit_code = _pipeline_main(["diagnose-e"])
fake_diag.assert_called_once_with([])
```

(and the gate0 twin, `mock.patch(__name__ + "._diagnose_gate0_main")`, called
with no arguments). No other test changes anticipated — `StageCalibrateTest`,
`RunHaltsOnStageFailureTest`, `StageReportTest`'s existing assertions target
`_stage_calibrate`/`run`/`_stage_report`, none of which move or rename.

### 6. Docs updated last, only after the gate passes

`CLAUDE.md`'s "Running things" section and the README's module map get
rewritten to describe one file and one test command:

```bash
python -m unittest rpi_core                        # all tests
python rpi_core.py run [--recalibrate] [--seed N]
python rpi_core.py calibrate [--stability]
python rpi_core.py validate --json ...
python rpi_core.py diagnose-e [panel.json]
python rpi_core.py diagnose-gate0
```

The `notes_consolidation.md` Phase 1–3 record is left as-is (it's an accurate
history of what already happened); a new short section is appended recording
this phase's deletions, not rewriting the earlier phases.

## Closing gate (must pass before any doc change lands)

1. `python -m unittest rpi_core` — 322 tests green (same count; no tests
   added or removed, two existing tests rewritten in place).
2. Seeded discrimination panel (`--seed 20260829`) bit-identical to
   `output/consolidation_baseline/discrimination_panel_golden.json`.
3. `python rpi_core.py --help`, `run --help`, `calibrate --help`, `validate
   --help` produce sane output (byte-identical `--help` text is not required
   here, unlike Phase 1–2's shim work, since the CLI surface is intentionally
   changing — `pipeline.py` becomes `rpi_core.py`, and `calibrate` is newly
   reachable as a subcommand rather than a separate script).
4. `diagnose-e` and `diagnose-gate0` smoke-tested against a real panel file /
   the reference dataset, output compared by eye against a pre-change run of
   the original standalone scripts for the same input.

Only after all four pass: delete the eleven files listed in "Files deleted"
in one commit, update `CLAUDE.md`/README in a second commit.
