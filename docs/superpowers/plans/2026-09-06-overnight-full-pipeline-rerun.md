# Overnight Full-Pipeline Rerun Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. This plan runs mostly unattended overnight in a background process; there is no user available between tasks, so execution proceeds straight through unless a task's own verification fails or a materially-bigger-than-expected issue is found (see Global Constraints).

**Goal:** Rerun the full pipeline (`all_hanta.py` → `phylogeny.py` → `rpi_core.py`/`Algorithm.py`) from scratch into an isolated folder, reconcile it against the existing `output/`, promote the reconciled result, illustrate the phylogeny trees, survey (without deleting) redundant artifacts, and produce a morning report — then commit locally (no push).

**Architecture:** A single PowerShell/Python driver script runs the pipeline stages in sequence into `output_rerun_20260906/`, launched detached so it survives the session overnight, with checkpointing on the expensive gate4-confirmation step. Comparison and reconciliation happen after the driver completes, in-session, since they require judgment. No new product code architecture — this plan orchestrates existing entry points (`all_hanta.py`, `phylogeny.py`, `rpi_core.py`, `Algorithm.py`, `figures.py`).

**Tech Stack:** Python (existing repo scripts), PowerShell (detached process launch, sleep-block), git.

**Spec:** `docs/superpowers/specs/2026-09-06-overnight-full-pipeline-rerun-design.md`

## Global Constraints

- Reuse pinned `data/rpi_calibration.json` — do not recalibrate.
- No autonomous code simplification — proposals only, written into the report.
- No deletions except the short-lived `output_rerun_20260906/` staging folder after promotion (step 6). Everything else identified as redundant goes on a list only.
- No `git push`.
- Pre-run code review (Task 1) follows the standing score-neutral-fix-only boundary. Only the reconciliation step (Task 7) has full autonomy to change results.
- If a shell tool gets blocked, stop and tell the user — do not work around it.
- If a materially-bigger-than-expected issue turns up (a real correctness bug, not an expected recalculation), stop there and flag it in the report rather than pushing through.
- Gate4 verdicts for all 10 long-pool winners must be reported verbatim, one line each, `no_conflict`/`insufficient_comparators` included — no re-running for a "better" result.

---

### Task 1: Pre-run code review of the four pipeline scripts

**Files:**
- Read: `all_hanta.py`, `phylogeny.py`, `rpi_core.py`, `Algorithm.py`
- Modify (only if score-neutral fixes found, one file/fix per commit): any of the above
- Test: existing suites (`python -m unittest rpi_core`, `python -m unittest phylogeny`, `python -m unittest figures` if present)

**Interfaces:**
- Produces: a findings list (bugs/bias risks that could affect tonight's run vs. proposals for later) — carried forward into Task 10's report, not a separate file.

- [ ] **Step 1: Run the full existing test suite as a pre-flight check**

Run: `python -m unittest rpi_core -v 2>&1 | tail -20` and `python -m unittest phylogeny -v 2>&1 | tail -20`
Expected: all green (matches the 328-test count recorded in `docs/PROJECT-STATE.md` §1 — if the count differs, note it, don't assume it's wrong).

- [ ] **Step 2: Review `all_hanta.py` for correctness/bias issues affecting ingest**

Focus: anything that would make tonight's rerun non-reproducible or biased toward the existing `output/` (hardcoded paths into `output/`, any code that reads its own prior output as input, non-deterministic ordering/iteration over sets or dicts that affects written CSV row order).

- [ ] **Step 3: Review `phylogeny.py` for the same, plus the two known-open items**

Check whether `PhylogeneticsScopeEndToEndTest`'s Windows-portability fix (commit `6a75f48`) and the phylogroup-scope grouping logic are still consistent with how Task 3 (below) will invoke both scopes.

- [ ] **Step 4: Review `rpi_core.py` for the same, cross-checking open §7e findings 1, 3–11 in `docs/PROJECT-STATE.md`**

For each open finding, decide: does it affect tonight's run's correctness (fix now, score-neutral + guard test + commit) or not (leave as an open proposal, carry into Task 10's report unchanged)? Do not fix anything that would change a reported number here — that's Task 7's job, after the diff shows it's real.

- [ ] **Step 5: Review `Algorithm.py`, paying particular attention to open Findings 9–11** (duplicate `candidate_id="ga-confirm"` across winners, silent systemic-failure-as-flat-landscape, `--elitism 0` accepted)

Finding 9 (duplicate candidate_id colliding gate4 artifact directories) is likely to affect tonight's run directly if it's still unfixed — verify, and if so this is score-neutral (attribution bug, doesn't change any winner's fitness score) so fix it now with a guard test before Task 3 launches gate4 confirmation on the new run.

- [ ] **Step 6: Run the full test suite again to confirm any Task 1 fixes are green**

Run: `python -m unittest rpi_core -v 2>&1 | tail -20`
Expected: PASS, count noted for the report.

- [ ] **Step 7: Commit each Task 1 fix individually** (skip if no fixes were needed)

```bash
git add <file>
git commit -m "<script>.py: <fix description> (pre-rerun review)"
```

---

### Task 2: Write the overnight driver script

**Files:**
- Create: `<scratchpad>/overnight_rerun_driver.py` (ad hoc orchestration, not product code — matches last session's precedent for the gate4-confirmation driver; not committed)

**Interfaces:**
- Consumes: `all_hanta.py`, `phylogeny.py`, `rpi_core.py`, `Algorithm.py`, `figures.py` — all invoked as subprocesses or via their existing `main(argv)` entry points, writing into `output_rerun_20260906/` instead of `output/`.
- Produces: `output_rerun_20260906/` populated with the same file layout as `output/` (tables, `phylogenetics/`, `phylogenetics_phylogroup/`, `validation/validation_report.json`, `ga_search_corrected_pool_default/`, `ga_search_corrected_pool_long/` incl. `gate4_confirmation.json`/checkpoint, `figures/before/`, `figures/after/`), plus a top-level `overnight_rerun_driver.log`.

- [ ] **Step 1: Write the driver script**

```python
"""Ad hoc overnight driver: reruns the full pipeline into output_rerun_20260906/.
Not product code -- one-off orchestration, matching the precedent set by
last session's gate4-confirmation driver. Each stage's failure halts the
run there (checkpoint files preserve partial gate4 progress)."""
import logging
import subprocess
import sys
from pathlib import Path

OUT = Path("output_rerun_20260906")
OUT.mkdir(exist_ok=True)

logging.basicConfig(
    filename=str(OUT / "overnight_rerun_driver.log"),
    level=logging.INFO,
    format="%(asctime)s %(message)s",
)
log = logging.getLogger("driver")


def run_stage(name, argv):
    log.info("=== stage start: %s ===", name)
    log.info("command: %s", " ".join(argv))
    result = subprocess.run(argv, capture_output=True, text=True)
    log.info("exit code: %s", result.returncode)
    log.info("stdout tail:\n%s", result.stdout[-4000:])
    log.info("stderr tail:\n%s", result.stderr[-4000:])
    if result.returncode != 0:
        log.error("=== stage FAILED: %s ===", name)
        raise SystemExit(f"stage {name} failed, exit {result.returncode}")
    log.info("=== stage done: %s ===", name)


def main():
    run_stage("all_hanta ingest", [
        sys.executable, "all_hanta.py", "--output", str(OUT),
    ])
    run_stage("phylogeny (plain scope)", [
        sys.executable, "phylogeny.py", "--output", str(OUT),
        "--input-dir", str(OUT),
    ])
    run_stage("phylogeny (phylogroup scope)", [
        sys.executable, "phylogeny.py", "--output", str(OUT),
        "--input-dir", str(OUT), "--scope", "phylogroup",
    ])
    run_stage("rpi_core run", [
        sys.executable, "rpi_core.py", "run", "--reference-input", str(OUT),
        "--output", str(OUT),
    ])
    run_stage("rpi_core validate (with phylogeny, all 5 gates)", [
        sys.executable, "rpi_core.py", "validate", "--with-phylogeny",
        "--reference-input", str(OUT),
        "--json", str(OUT / "validation" / "validation_report.json"),
    ])
    run_stage("GA search (default pool)", [
        sys.executable, "Algorithm.py", "--reference-input", str(OUT),
        "--output", str(OUT / "ga_search_corrected_pool_default"),
    ])
    run_stage("GA search (long pool)", [
        sys.executable, "Algorithm.py", "--reference-input", str(OUT),
        "--output", str(OUT / "ga_search_corrected_pool_long"),
        "--generations", "300",
    ])
    run_stage("gate4 confirmation of long-pool winners", [
        sys.executable, "Algorithm.py", "confirm-winners",
        "--reference-input", str(OUT),
        "--winners", str(OUT / "ga_search_corrected_pool_long" / "seed_20260904.json"),
        "--with-phylogeny",
        "--checkpoint", str(OUT / "ga_search_corrected_pool_long" / "gate4_confirmation_checkpoint.jsonl"),
        "--output", str(OUT / "ga_search_corrected_pool_long" / "gate4_confirmation.json"),
    ])
    run_stage("figures (before/after)", [
        sys.executable, "figures.py", "--input-dir", str(OUT), "--output", str(OUT / "figures"),
    ])
    log.info("=== ALL STAGES COMPLETE ===")


if __name__ == "__main__":
    main()
```

**Note for the implementer:** the exact CLI flags above (`--reference-input`, `--output`, `--input-dir`, `--scope`, the `confirm-winners` subcommand name, `--generations` for the long pool) must be verified against each script's actual `parse_args`/`argparse` setup before this is run — do not trust this plan's flag names blindly, `argparse.ArgumentParser` in each of the four scripts is the ground truth. Fix the driver script's argv lists to match. This is expected to take a few minutes of cross-checking, not a redesign.

- [ ] **Step 2: Dry-run the first stage only to confirm the driver script and flags work**

Run: `python <scratchpad>/overnight_rerun_driver.py` and interrupt after the first "stage done" log line appears (or run with a temporary early `raise SystemExit` after stage 1 for the dry run), then check `output_rerun_20260906/overnight_rerun_driver.log` shows stage 1 succeeded and `output_rerun_20260906/` contains the expected ingest tables.

- [ ] **Step 3: Remove the temporary early-exit if one was added, restoring the full stage list**

---

### Task 3: Launch the overnight run detached

**Files:**
- None created; this is a process-launch step.

- [ ] **Step 1: Launch the driver detached so it survives the session, blocking sleep for its lifetime**

```powershell
$psi = Start-Process -FilePath "python" -ArgumentList "<scratchpad>/overnight_rerun_driver.py" -WindowStyle Hidden -PassThru
Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
public class Sleep {
    [DllImport("kernel32.dll")]
    public static extern uint SetThreadExecutionState(uint esFlags);
}
'@
[Sleep]::SetThreadExecutionState(0x80000003) # ES_CONTINUOUS | ES_SYSTEM_REQUIRED | ES_AWAYMODE_REQUIRED
$psi.Id
```

Expected: a process ID is printed; `output_rerun_20260906/overnight_rerun_driver.log` starts growing.

- [ ] **Step 2: Confirm the process is actually running (not exited immediately)**

Run: `Get-Process -Id <pid> -ErrorAction SilentlyContinue`
Expected: process present. If absent, check the log tail for an immediate failure before proceeding — do not re-launch in a retry loop; diagnose first.

---

### Task 4: Wait for completion

- [ ] **Step 1: Poll for completion via the log's terminal line, at intervals matched to the expected multi-hour runtime**

Use `ScheduleWakeup`/`Monitor` per session tooling rather than a tight sleep loop: check for `=== ALL STAGES COMPLETE ===` in `output_rerun_20260906/overnight_rerun_driver.log`, or a `SystemExit`/stage-FAILED line. Do not poll more often than every 15–20 minutes given the expected multi-hour runtime dominated by gate4 confirmation.

- [ ] **Step 2: On completion (success or failure), read the full driver log**

If a stage failed: diagnose from the log, decide whether it's a Task-1-missed bug (fix, guard test, commit, resume the driver from the failed stage) or an environment issue (report it, do not silently retry more than once).

---

### Task 5: Checksum diff

**Files:**
- Create: `<scratchpad>/checksum_diff.py` (ad hoc, not committed) or reuse whatever generated `output/ortho_output_checksums.txt`.

- [ ] **Step 1: Generate a checksums file for `output_rerun_20260906/`** using the same method as `output/ortho_output_checksums.txt` (check how that file was produced — likely `sha256sum`/`hashlib.sha256` over each output file, relative path as the label).

- [ ] **Step 2: Diff the two checksums files, listing every path present in one but not the other, and every path with a different hash**

```bash
diff <(sort output/ortho_output_checksums.txt) <(sort output_rerun_20260906/ortho_output_checksums.txt) > <scratchpad>/checksum_diff.txt
```

- [ ] **Step 3: Read `checksum_diff.txt` and classify each line: new/removed file, or changed content**

---

### Task 6: Semantic diff on changed files

- [ ] **Step 1: For each changed CSV/JSON from Task 5, load both versions and diff at the value level** (not just text diff, since column order or float formatting could differ trivially) — write findings inline into the report draft (Task 10), not a separate artifact.

- [ ] **Step 2: For each changed tree file (`.treefile`), compare topology (e.g. via `dendropy.calculate.treecompare.symmetric_difference`, already a dependency per `phylogeny.py`) rather than raw text, since Newick string formatting can differ without the topology differing.**

---

### Task 7: Reconcile real differences

- [ ] **Step 1: For each genuine difference (not just formatting/ordering noise) found in Task 6, root-cause it** — check whether it traces to a Task-1-missed bug, a genuinely non-deterministic step (unseeded randomness, dict/set iteration order, thread-timing-dependent IQ-TREE behavior), or an intentional prior-session workaround that's no longer being applied.

- [ ] **Step 2: Fix the root cause in code if it's a bug** (guard-tested, one commit per fix — this step has full autonomy per Global Constraints, unlike Task 1).

- [ ] **Step 3: If it's genuine non-determinism with no bug to fix, document which output(s) are affected and why, rather than "fixing" by picking one arbitrarily.**

- [ ] **Step 4: Re-run only the affected stage(s) if a fix was applied, to confirm the fix resolves the diff** (do not re-run the entire overnight pipeline for a single-stage fix).

---

### Task 8: Promote the reconciled output

- [ ] **Step 1: Confirm all Task 7 reconciliations are resolved or explicitly documented as non-deterministic-and-acceptable.**

- [ ] **Step 2: Replace `output/` with the reconciled `output_rerun_20260906/` contents**

```bash
rm -rf output.bak && mv output output.bak && mv output_rerun_20260906 output
```

- [ ] **Step 3: Spot-check the promoted `output/` looks complete** (same top-level entries as `output.bak` had, per Task 5's checksum file listing).

- [ ] **Step 4: Remove the `output.bak` staging copy** (the one deletion this plan performs outside the redundancy-survey list, since it's transient scaffolding from the promotion mechanics, not project output).

```bash
rm -rf output.bak
```

---

### Task 9: Redundancy survey

- [ ] **Step 1: List candidates for removal with reasons**, at minimum:
  - The 4 orphaned stale gate4 files from Task 1's baseline (now superseded by the fresh run's own `gate4_confirmation.json`).
  - Any superseded scratch/driver scripts left in the repo (there shouldn't be any — the driver lives in the scratchpad — but check for stray ones from prior sessions per `git log --all --diff-filter=A -- '*.py'` against what's actually referenced from README/docs).
  - Stray `__pycache__` directories.
  - Any `output/ga_search_*` run directories superseded by tonight's rerun (check `docs/PROJECT-STATE.md` §7d/§7e for which are already known-superseded vs. load-bearing evidence).

- [ ] **Step 2: For each candidate, confirm it is not referenced from README, `docs/PROJECT-STATE.md`, or any test**, so the list doesn't recommend removing something still load-bearing (grep for the filename across `*.py`, `*.md` before listing it).

- [ ] **Step 3: Write the list into the report draft (Task 10) — no deletions.**

---

### Task 10: Report and docs tightening

**Files:**
- Modify: `docs/PROJECT-STATE.md` (§7f verdict writeup, redundancy list, reconciliation log)
- Modify: `README.md` (if any headline numbers changed as a result of Task 7 reconciliation)

- [ ] **Step 1: Write the 10 gate4 verdicts for the long-pool winners into `docs/PROJECT-STATE.md` §7f, one line each, verbatim** (`no_conflict`/`insufficient_comparators` included, per the binding reporting rule already recorded there).

- [ ] **Step 2: State explicitly whether the stalled prior-session gate4 run covered anything the fresh rerun does not** (it shouldn't, since Task 2's driver reruns `confirm_winners` fresh on the new run's winners — but confirm the winner set is the same triplets before asserting equivalence, since a fresh GA search could in principle land on different tied triplets per the known structural-ceiling-flatness finding).

- [ ] **Step 3: Write Task 7's reconciliation log** (problem found → root cause → fix applied, or non-determinism documented) into `docs/PROJECT-STATE.md`.

- [ ] **Step 4: Write Task 9's redundancy list** into `docs/PROJECT-STATE.md` or a new short section, clearly marked "not yet actioned — awaiting review."

- [ ] **Step 5: Write Task 1's un-applied findings (score-changing, not fixed) as proposals**, same place.

- [ ] **Step 6: Update README's headline numbers/test count if Task 7 changed any of them**, and tighten any now-duplicated prose between README and `docs/PROJECT-STATE.md` without cutting content that documents a real limitation or decision.

- [ ] **Step 7: Update `docs/PROJECT-STATE.md` §1's snapshot table** (branch, HEAD description, test count) to reflect the new state.

---

### Task 11: Illustrate the phylogeny trees

**Files:**
- Verify: `output/figures/before/`, `output/figures/after/` (produced by Task 2's driver via `figures.py` against the promoted output)

- [ ] **Step 1: Confirm `figures.py`'s before/after PNGs were regenerated against the promoted `output/` and are non-empty/valid images** (open one file's header bytes or check file size > a trivial threshold, don't just check existence).

- [ ] **Step 2: If `figures.py` needs an update to point at the new run's paths (rather than hardcoded `output/` — check whether Task 2's driver already parameterized this via `--input-dir`/`--output`), fix and re-run just that stage, not the whole pipeline.**

---

### Task 12: Final commit (no push)

- [ ] **Step 1: Review `git status` and `git diff` for everything touched by Tasks 7–11.**

- [ ] **Step 2: Stage and commit the promoted output plus docs/report updates**

```bash
git add output/ docs/PROJECT-STATE.md README.md
git commit -m "$(cat <<'EOF'
Overnight full-pipeline rerun: promote reconciled output, gate4
verdicts for all 10 long-pool winners, redundancy survey

Co-Authored-By: Claude Haiku 4.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_014UNgjonny5S2VVFfKJ7WPy
EOF
)"
```

- [ ] **Step 3: Confirm no push was run** (`git log --oneline -1 origin/05-module-consolidation` vs local HEAD should still differ, or `git status` shows "ahead of origin").

- [ ] **Step 4: Report back to the user with a summary**: what ran, what matched, what was reconciled and why, the 10 gate4 verdicts, the redundancy list awaiting their review, and the un-applied Task-1/Task-7 proposals.
