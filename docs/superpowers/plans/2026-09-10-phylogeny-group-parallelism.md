# phylogeny.py group-level parallelism Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Parallelize `phylogeny.py`'s per-group (species/phylogroup) alignment+tree-building loop across a bounded worker pool, without changing what any individual group computes.

**Architecture:** Replace the serial `for group_key, members in sorted(by_group.items())` loop in `_run_phylogenetics_scope` with a `ThreadPoolExecutor`, sized by a new `--phylo-workers` flag (default 4). Groups are submitted largest-taxa-first (minimizes makespan given this dataset's extreme group-size skew) but results are always sorted by `group_key` before being returned, so output row order stays identical to today's serial run. `_analyze_group` itself is unchanged.

**Tech Stack:** Python stdlib `concurrent.futures.ThreadPoolExecutor` (subprocess-bound work releases the GIL, so threads — not processes — are the right tool here; this sidesteps Windows `spawn` overhead and pickling constraints entirely). No new dependencies.

**Spec:** `docs/superpowers/specs/2026-09-10-phylogeny-group-parallelism-design.md`

## Global Constraints

- `--phylo-threads`' existing default formula (`max(1, (os.cpu_count() or 4) - 1)`) stays unchanged — no auto-division by `--phylo-workers`. This is an explicit, accepted-risk decision from the spec (§4.3), not an oversight to "fix" mid-implementation.
- Output row order (`ortho_reassortment_table.csv`) must stay sorted by `group_key`, identical to a serial run, regardless of `--phylo-workers` value or actual thread-completion order.
- Any within-group speedup (MAFFT strategy, extra alignment threads) is explicitly out of scope (spec §3) — do not touch `_analyze_group`'s internals.
- Existing tests must keep passing unmodified in behavior — only the `_args()` fixture helper needs a new field (`phylo_workers=1`) to keep exercising today's exact single-worker path.
- This repo's test convention: stdlib `unittest`, run with `python -m unittest phylogeny -v`. Follow existing patterns in the file (see Task 2's cited line numbers for exact style to match).

---

## Task 1: `--phylo-workers` CLI flag

**Files:**
- Modify: `phylogeny.py:947-982` (`parse_args`)
- Test: `phylogeny.py` (new `ParseArgsPhyloWorkersTest` class, inserted after `GroupCandidatesTest` ends at line 1318, before the `@unittest.skipUnless` decorator at line 1320)

**Interfaces:**
- Consumes: nothing new.
- Produces: `args.phylo_workers` (`int`, default `4`, guaranteed `>= 1` by the time any caller sees it) — Task 2 reads this directly.

- [ ] **Step 1: Write the failing tests**

Insert into `phylogeny.py` right after `GroupCandidatesTest`'s closing (currently ends at line 1318 with a blank line, then `@unittest.skipUnless(...)` at 1320 — insert the new class between them):

```python
class ParseArgsPhyloWorkersTest(unittest.TestCase):
    """--phylo-workers: new flag from docs/superpowers/specs/2026-09-10-phylogeny-
    group-parallelism-design.md. Defaults to 4, rejects non-positive values at parse
    time rather than failing later inside ThreadPoolExecutor with a confusing error.
    """

    def test_defaults_to_four(self):
        args = parse_args([])
        self.assertEqual(args.phylo_workers, 4)

    def test_accepts_explicit_value(self):
        args = parse_args(["--phylo-workers", "8"])
        self.assertEqual(args.phylo_workers, 8)

    def test_rejects_zero(self):
        with self.assertRaises(SystemExit):
            parse_args(["--phylo-workers", "0"])

    def test_rejects_negative(self):
        with self.assertRaises(SystemExit):
            parse_args(["--phylo-workers", "-1"])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m unittest phylogeny.ParseArgsPhyloWorkersTest -v`
Expected: all 4 FAIL with `AttributeError: 'Namespace' object has no attribute 'phylo_workers'` (the flag doesn't exist yet).

- [ ] **Step 3: Add the flag and validation**

In `phylogeny.py`, find this exact block (currently lines 980-982):

```python
    parser.add_argument("--phylo-outgroup", default=None,
                        help="Accession or isolate_key to root trees on, via IQ-TREE's -o.")
    return parser.parse_args(argv)
```

Replace it with:

```python
    parser.add_argument("--phylo-outgroup", default=None,
                        help="Accession or isolate_key to root trees on, via IQ-TREE's -o.")
    parser.add_argument("--phylo-workers", type=int, default=4,
                        help="Number of groups (species/phylogroup trees) to align + build "
                             "concurrently. Each worker still uses --phylo-threads threads for "
                             "its own mafft/iqtree calls, so total thread demand is roughly "
                             "workers x --phylo-threads -- lower --phylo-threads manually if "
                             "you raise this and see CPU oversubscription.")
    parsed = parser.parse_args(argv)
    if parsed.phylo_workers < 1:
        parser.error("--phylo-workers must be >= 1")
    return parsed
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m unittest phylogeny.ParseArgsPhyloWorkersTest -v`
Expected: all 4 PASS.

- [ ] **Step 5: Run the full phylogeny.py suite to confirm no regression**

Run: `python -m unittest phylogeny -v 2>&1 | tail -20`
Expected: same pass/skip count as before this task (27 tests before this task per the branch's existing state; now 31 with the 4 new ones), 0 failures.

- [ ] **Step 6: Commit**

```bash
git add phylogeny.py
git commit -m "phylogeny.py: add --phylo-workers CLI flag (default 4, validated >=1)"
```

---

## Task 2: Parallelize the group loop

**Files:**
- Modify: `phylogeny.py:17-31` (imports)
- Modify: `phylogeny.py:799-832` (`_run_phylogenetics_scope`)
- Modify: `phylogeny.py:1011-1016` (test-section imports)
- Modify: `phylogeny.py:1412-1414` (`PhylogeneticsScopeEndToEndTest._args`)
- Test: `phylogeny.py` (one new method on `PhylogeneticsScopeEndToEndTest`; two new test classes, `GroupSubmissionOrderTest` and `GroupFailureIsolationTest`, inserted in the same location as Task 1's new class — after `GroupCandidatesTest`, before the `@unittest.skipUnless` line. Order within that gap: `ParseArgsPhyloWorkersTest` (Task 1), then `GroupSubmissionOrderTest`, then `GroupFailureIsolationTest`.)

**Interfaces:**
- Consumes: `args.phylo_workers` from Task 1.
- Produces: `_run_phylogenetics_scope`'s public signature and return type are unchanged (`List[Dict[str, Any]]`, same row shape) — no other task depends on new names here.

- [ ] **Step 1: Write the failing tests**

First, the determinism test. Add this method to `PhylogeneticsScopeEndToEndTest` (phylogeny.py, the class starting at line 1321), placed after `test_outgroup_applied_only_where_present_and_logged_when_absent` (its last method):

```python
    def test_worker_count_does_not_change_output_rows(self):
        """Regression: docs/superpowers/specs/2026-09-10-phylogeny-group-parallelism-
        design.md's correctness contract -- parallelizing WHEN groups run must not
        change WHAT rows come out or their order. Runs against separate output
        directories so file writes don't collide; recorded paths are relative to each
        call's own output_dir, so they come out identical either way."""
        serial_out = self.output_dir / "serial"
        serial_out.mkdir()
        parallel_out = self.output_dir / "parallel"
        parallel_out.mkdir()

        serial_rows = _run_phylogenetics_scope(
            "species", self.candidates, serial_out, self._args("species"),
            str(self.mafft_bin), str(self.iqtree_bin))

        parallel_args = self._args("species")
        parallel_args.phylo_workers = 4
        parallel_rows = _run_phylogenetics_scope(
            "species", self.candidates, parallel_out, parallel_args,
            str(self.mafft_bin), str(self.iqtree_bin))

        self.assertEqual(serial_rows, parallel_rows)
```

Next, update `_args` (currently lines 1412-1414) so every existing caller of it also exercises the new field — find:

```python
    def _args(self, scope, outgroup=None):
        import argparse
        return argparse.Namespace(
            phylo_min_taxa=4, bootstrap=100, min_bootstrap_support=95.0, phylo_seed=42,
            phylo_threads=1, phylo_scope=scope, phylo_outgroup=outgroup,
        )
```

Replace with:

```python
    def _args(self, scope, outgroup=None):
        import argparse
        return argparse.Namespace(
            phylo_min_taxa=4, bootstrap=100, min_bootstrap_support=95.0, phylo_seed=42,
            phylo_threads=1, phylo_scope=scope, phylo_outgroup=outgroup, phylo_workers=1,
        )
```

Now the two new standalone test classes. Insert both after `ParseArgsPhyloWorkersTest` (added in Task 1), still before the `@unittest.skipUnless(_HAVE_PHYLO_DEPS, ...)` line that starts `PhylogeneticsScopeEndToEndTest`:

```python
class GroupSubmissionOrderTest(unittest.TestCase):
    """_run_phylogenetics_scope must submit groups largest-taxa-first (the "longest
    processing time first" scheduling heuristic -- see docs/superpowers/specs/
    2026-09-10-phylogeny-group-parallelism-design.md section 4.1), not
    group_candidates' incidental alphabetical order. Uses phylo_workers=1 so execution
    order is deterministically submission order -- a single-worker ThreadPoolExecutor
    runs queued tasks strictly FIFO -- making this a scheduling-order test, not a
    timing-based concurrency test. Mocks _analyze_group entirely, so it needs neither
    Bio.Phylo/dendropy nor real mafft/iqtree.
    """

    def test_groups_analyzed_largest_taxa_first(self):
        candidates = (
            _make_candidates(["Orthohantavirus hantanense"], n_per_species=4, seed=1)
            + _make_candidates(["Orthohantavirus dobravaense"], n_per_species=9, seed=2)
            + _make_candidates(["Orthohantavirus puumalaense"], n_per_species=6, seed=3)
        )
        call_order: List[str] = []

        def fake_analyze_group(group_key, scope, members, *_a, **_kw):
            call_order.append(group_key)
            return []

        args = argparse.Namespace(phylo_min_taxa=4, phylo_workers=1)
        with mock.patch("phylogeny._analyze_group", side_effect=fake_analyze_group):
            _run_phylogenetics_scope("species", candidates, Path("unused"), args, "mafft", "iqtree")

        self.assertEqual(call_order, [
            "Orthohantavirus dobravaense",   # 9 candidates -- largest
            "Orthohantavirus puumalaense",   # 6 candidates
            "Orthohantavirus hantanense",    # 4 candidates -- smallest
        ])


class GroupFailureIsolationTest(unittest.TestCase):
    """One group's failure must not lose other groups' results under the concurrent
    path -- extends _analyze_group's existing per-group error isolation (its own
    try/except, unchanged) with a defensive layer at the future.result() collection
    point, for whatever could still raise before that try block starts (e.g. the
    labels.json write). Mocks _analyze_group entirely.
    """

    def test_one_groups_exception_does_not_lose_other_groups_rows(self):
        candidates = (
            _make_candidates(["Orthohantavirus hantanense"], n_per_species=4, seed=1)
            + _make_candidates(["Orthohantavirus dobravaense"], n_per_species=4, seed=2)
        )

        def fake_analyze_group(group_key, scope, members, *_a, **_kw):
            if group_key == "Orthohantavirus dobravaense":
                raise RuntimeError("simulated failure before _analyze_group's own try block")
            return [{"isolate_key": c.isolate.isolate_key, "phylo_status": "analyzed"}
                    for c in members]

        args = argparse.Namespace(phylo_min_taxa=4, phylo_workers=2)
        with mock.patch("phylogeny._analyze_group", side_effect=fake_analyze_group):
            rows = _run_phylogenetics_scope(
                "species", candidates, Path("unused"), args, "mafft", "iqtree")

        by_key = {row.get("isolate_key"): row for row in rows}
        hantanense_keys = [c.isolate.isolate_key for c in candidates
                           if c.isolate.species_name == "Orthohantavirus hantanense"]
        dobravaense_keys = [c.isolate.isolate_key for c in candidates
                            if c.isolate.species_name == "Orthohantavirus dobravaense"]

        for key in hantanense_keys:
            self.assertEqual(by_key[key]["phylo_status"], "analyzed")
        for key in dobravaense_keys:
            self.assertEqual(by_key[key]["phylo_status"], "failed")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m unittest phylogeny.GroupSubmissionOrderTest phylogeny.GroupFailureIsolationTest -v`
Expected: both FAIL — `GroupSubmissionOrderTest` with `AttributeError: 'Namespace' object has no attribute 'phylo_workers'` (or a `NameError`/`AttributeError` on `mock`, since it isn't imported yet); `GroupFailureIsolationTest` likewise.

Run: `python -m unittest phylogeny.PhylogeneticsScopeEndToEndTest.test_worker_count_does_not_change_output_rows -v` (skipped automatically if biopython/dendropy aren't installed in this environment — that's expected, not a failure; if it does run, expect it to fail with the same `phylo_workers` AttributeError).

- [ ] **Step 3: Add `ThreadPoolExecutor`/`as_completed` and `mock` imports**

In the production imports block (phylogeny.py, currently lines 17-31), find:

```python
from collections import Counter, defaultdict
from dataclasses import dataclass
```

Replace with:

```python
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
```

In the test-section imports block (currently lines 1011-1016), find:

```python
import json
import random
import sys
import tempfile
import unittest
from io import StringIO
```

Replace with:

```python
import json
import random
import sys
import tempfile
import unittest
from io import StringIO
from unittest import mock
```

- [ ] **Step 4: Replace the serial loop with the parallel implementation**

In `phylogeny.py`, find `_run_phylogenetics_scope`'s exact current body (lines 799-832):

```python
def _run_phylogenetics_scope(
    scope: str, candidates: Sequence[PhyloCandidate], output_dir: Path, args: argparse.Namespace,
    mafft_bin: str, iqtree_bin: str,
) -> List[Dict[str, Any]]:
    """One scope's worth of grouping + tree-building + classification. Factored out of
    run_phylogenetics_stage so --phylo-scope=both is just this called twice; the
    scope='species' call alone is exactly what run_phylogenetics_stage did before this
    function existed.
    """
    by_group, group_drop_reasons = group_candidates(candidates, args.phylo_min_taxa, scope)
    for reason in group_drop_reasons:
        logging.info("Phylogenetics: excluded group below --phylo-min-taxa: %s", reason)
    logging.info("Phylogenetics (%s): analyzing %s group(s): %s", scope, len(by_group), sorted(by_group))

    grouped_keys = {c.isolate.isolate_key for members in by_group.values() for c in members}
    candidates_per_group: DefaultDict[str, int] = defaultdict(int)
    for c in candidates:
        candidates_per_group[_candidate_group_key(c, scope)] += 1
    rows: List[Dict[str, Any]] = [
        _phylo_row(c, c.isolate.species_name, candidates_per_group[_candidate_group_key(c, scope)],
                  "species_below_min_taxa", "", scope, _assignment_basis(c, scope))
        for c in candidates if c.isolate.isolate_key not in grouped_keys
    ]

    # Phylogroup output lives under a sibling directory, never inside phylogenetics/, so a
    # species-scope-only run's directory tree (and therefore the alignment/tree path
    # columns it writes) stays exactly what it was before scope existed.
    group_dir_root = output_dir / "phylogenetics" if scope == "species" else output_dir / "phylogenetics_phylogroup"
    for group_key, members in sorted(by_group.items()):
        logging.info("Phylogenetics (%s): %s (%s taxa) -- aligning + building trees",
                     scope, group_key, len(members))
        rows.extend(_analyze_group(
            group_key, scope, members, group_dir_root, output_dir, args, mafft_bin, iqtree_bin))
    return rows
```

Replace it with:

```python
def _run_phylogenetics_scope(
    scope: str, candidates: Sequence[PhyloCandidate], output_dir: Path, args: argparse.Namespace,
    mafft_bin: str, iqtree_bin: str,
) -> List[Dict[str, Any]]:
    """One scope's worth of grouping + tree-building + classification. Factored out of
    run_phylogenetics_stage so --phylo-scope=both is just this called twice; the
    scope='species' call alone is exactly what run_phylogenetics_stage did before this
    function existed.

    Groups are analyzed concurrently, up to args.phylo_workers at once (see
    docs/superpowers/specs/2026-09-10-phylogeny-group-parallelism-design.md) --
    submitted largest-taxa-first to minimize total makespan given this dataset's
    extreme group-size skew, but always returned sorted by group_key so output row
    order stays identical to a serial run regardless of worker count or completion
    order.
    """
    by_group, group_drop_reasons = group_candidates(candidates, args.phylo_min_taxa, scope)
    for reason in group_drop_reasons:
        logging.info("Phylogenetics: excluded group below --phylo-min-taxa: %s", reason)
    logging.info("Phylogenetics (%s): analyzing %s group(s): %s", scope, len(by_group), sorted(by_group))

    grouped_keys = {c.isolate.isolate_key for members in by_group.values() for c in members}
    candidates_per_group: DefaultDict[str, int] = defaultdict(int)
    for c in candidates:
        candidates_per_group[_candidate_group_key(c, scope)] += 1
    rows: List[Dict[str, Any]] = [
        _phylo_row(c, c.isolate.species_name, candidates_per_group[_candidate_group_key(c, scope)],
                  "species_below_min_taxa", "", scope, _assignment_basis(c, scope))
        for c in candidates if c.isolate.isolate_key not in grouped_keys
    ]

    # Phylogroup output lives under a sibling directory, never inside phylogenetics/, so a
    # species-scope-only run's directory tree (and therefore the alignment/tree path
    # columns it writes) stays exactly what it was before scope existed.
    group_dir_root = output_dir / "phylogenetics" if scope == "species" else output_dir / "phylogenetics_phylogroup"

    # Largest-taxa-first submission order (the "longest processing time first"
    # scheduling heuristic): with a fixed worker count and this dataset's extreme
    # group-size skew, starting the biggest groups immediately -- rather than possibly
    # queuing them behind small ones -- minimizes total makespan. Ties broken by
    # group_key for determinism.
    submission_order = sorted(by_group.items(), key=lambda kv: (-len(kv[1]), kv[0]))

    def _run_one_group(group_key: str, members: List[PhyloCandidate]) -> List[Dict[str, Any]]:
        logging.info("Phylogenetics (%s): %s (%s taxa) -- aligning + building trees",
                     scope, group_key, len(members))
        return _analyze_group(
            group_key, scope, members, group_dir_root, output_dir, args, mafft_bin, iqtree_bin)

    group_rows: Dict[str, List[Dict[str, Any]]] = {}
    with ThreadPoolExecutor(max_workers=args.phylo_workers) as executor:
        future_to_group = {
            executor.submit(_run_one_group, group_key, members): (group_key, members)
            for group_key, members in submission_order
        }
        for future in as_completed(future_to_group):
            group_key, members = future_to_group[future]
            try:
                group_rows[group_key] = future.result()
            except Exception as exc:  # noqa: BLE001 - one group's failure must not kill the run
                # _analyze_group already catches everything inside its own try/except and
                # returns 'failed' rows instead of raising (see its docstring) -- this only
                # fires for something outside that guard, e.g. a labels.json write failure
                # before its try block starts. Same isolation guarantee, one layer out.
                logging.exception("Phylogenetics failed for %s (%s scope)", group_key, scope)
                error = f"{type(exc).__name__}:{exc}"
                group_rows[group_key] = [
                    _phylo_row(c, c.isolate.species_name, len(members), "failed", error, scope,
                              _assignment_basis(c, scope))
                    for c in members]

    for group_key in sorted(group_rows):
        rows.extend(group_rows[group_key])
    return rows
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python -m unittest phylogeny.GroupSubmissionOrderTest phylogeny.GroupFailureIsolationTest -v`
Expected: both PASS.

Run: `python -m unittest phylogeny.PhylogeneticsScopeEndToEndTest -v` (auto-skips if biopython/dendropy aren't installed — check the output either says all PASS or all SKIP, not a mix; a mix would mean something else broke).

- [ ] **Step 6: Run the full phylogeny.py suite**

Run: `python -m unittest phylogeny -v 2>&1 | tail -40`
Expected: 0 failures, 0 errors. Pass/skip counts should be exactly Task 1's post-count plus 3 (the determinism test on `PhylogeneticsScopeEndToEndTest`, plus `GroupSubmissionOrderTest` and `GroupFailureIsolationTest`'s one test each) — 34 total (31 after Task 1, some skipped if biopython/dendropy missing).

- [ ] **Step 7: Commit**

```bash
git add phylogeny.py
git commit -m "phylogeny.py: parallelize group loop via ThreadPoolExecutor, largest-first submission"
```

---

## Task 3: Benchmark against the 2026-09-06 serial baseline (requires mafft/iqtree)

This task needs real `mafft`/`iqtree` on PATH and the real dataset — it does not run in
a dev environment without those tools (this repo's Windows dev checkout currently lacks
them; see `.superpowers/sdd/2026-09-09-stage1-reference-snapshot/progress.md` Task 15).
Run it in the Linux/conda environment described in `environment.yml`. This task is
**not** part of the automated TDD loop Tasks 1-2 use — it is a one-time measurement that
decides whether `--phylo-workers=4` ships as the default (spec §5/§7).

**Files:**
- Create (scratch, not committed): a temporary `monitor_rss.py` on the machine running
  the benchmark — not part of this repo's tracked file set, since it is a one-off
  measurement tool, not a pipeline stage.
- Modify (only if the rollback condition fires): `phylogeny.py:973` (the
  `--phylo-workers` default) and this plan's own spec, documenting the measured
  reason.

- [ ] **Step 1: Verify the 2026-09-06 log is still a valid baseline**

`output/ortho_phylogeny.log` is gitignored (`output/*.log` in `.gitignore`) — it is
**not** in git history, it only exists as a local file on whatever machine still has it
from the 2026-09-06 run. Before doing anything that would run Stage 2 again and
overwrite it, copy it aside:

```bash
cp output/ortho_phylogeny.log /tmp/ortho_phylogeny_2026-09-06_baseline.log
```

If this machine does not have that file (e.g. it's a fresh checkout), locate it on
whatever machine ran the 2026-09-06 benchmark before proceeding — without it, §2's
baseline numbers cannot be verified as still valid, only re-derived from scratch (see
the fallback below).

Then run Stage 1 fresh (if not already current — check `output/ortho_summary.json`'s
`generated_at_utc` is recent) and Stage 2 with `--phylo-workers 1` (serial-equivalent,
matching the old log's own behavior) against the same input, and compare group taxa
counts in the new log output against the saved baseline:

```bash
python prepare_dataset.py --input . --output output
python phylogeny.py --output output --phylo-workers 1
grep "aligning + building trees" output/ortho_phylogeny.log | sort > /tmp/new_groups.txt
grep "aligning + building trees" /tmp/ortho_phylogeny_2026-09-06_baseline.log | sort > /tmp/old_groups.txt
diff /tmp/old_groups.txt /tmp/new_groups.txt
```

If `diff` shows any group's taxa count changed, re-derive this plan's §2 baseline numbers
from this fresh run instead of trusting the old log, and update the spec (§2) with the
new figures before proceeding — do not compare the parallel run against stale numbers.
If the 2026-09-06 log could not be located at all, this fresh serial run's own timings
become the baseline instead (same diffing/analysis method as spec §2 used).

- [ ] **Step 2: Write the RSS monitor script**

Save this as `monitor_rss.py` in a scratch location (e.g. `/tmp/monitor_rss.py`) on the
benchmark machine — Linux-only (reads `/proc`), matching this repo's `environment.yml`
deployment target:

```python
#!/usr/bin/env python3
"""Sample total RSS (a root process + all its descendants) at a fixed interval, so
peak RSS during a specific time window can be read off afterward by cross-referencing
sample timestamps against another log's own timestamps. Linux only (reads /proc).

Usage: python3 monitor_rss.py <root_pid> <interval_seconds> <out_csv>
Start phylogeny.py separately, then run this pointed at its PID. Stops on its own once
the root PID exits.
"""
import csv
import os
import sys
import time


def _children_of(pid: int) -> list:
    children = []
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open(f"/proc/{entry}/status") as f:
                for line in f:
                    if line.startswith("PPid:"):
                        if int(line.split()[1]) == pid:
                            children.append(int(entry))
                        break
        except (FileNotFoundError, ProcessLookupError):
            continue  # process exited between listdir() and open()
    return children


def _process_tree(root_pid: int) -> list:
    tree = [root_pid]
    frontier = [root_pid]
    while frontier:
        next_frontier = []
        for pid in frontier:
            next_frontier.extend(_children_of(pid))
        tree.extend(next_frontier)
        frontier = next_frontier
    return tree


def _rss_kb(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])  # kB
    except (FileNotFoundError, ProcessLookupError):
        pass
    return 0


def main():
    root_pid, interval, out_csv = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["timestamp_utc", "total_rss_kb", "process_count"])
        while os.path.exists(f"/proc/{root_pid}"):
            tree = _process_tree(root_pid)
            total_kb = sum(_rss_kb(pid) for pid in tree)
            writer.writerow([time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), total_kb, len(tree)])
            f.flush()
            time.sleep(interval)


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Sanity-check the monitor script**

Before trusting it on a multi-hour real run, verify it against a short dummy process:

```bash
sleep 30 &
python3 /tmp/monitor_rss.py $! 2 /tmp/monitor_test.csv &
wait
cat /tmp/monitor_test.csv
```

Expected: a handful of rows with `process_count=1` and a stable, nonzero `total_rss_kb`
for `sleep`'s duration, then the monitor exits once `sleep` does. If `process_count` or
`total_rss_kb` looks wrong (e.g. always 0), fix the script before using it for real —
don't proceed to Step 4 on a monitor you haven't verified reports sane numbers.

- [ ] **Step 4: Run the parallel benchmark**

`--output` doubles as both the directory phylogeny.py reads Stage 1's reference
snapshot from *and* where it writes its own results (see `phylogeny.py:953`) — so the
benchmark needs a scratch directory seeded with Stage 1's snapshot first, not a bare
empty one, and not the tracked `output/` itself (so this benchmark run doesn't clobber
it):

```bash
mkdir -p output_bench_parallel
cp output/ortho_reference_isolates.json output/ortho_reference_genomic.json \
   output/ortho_reference_cds.json output/ortho_reference_snapshot_manifest.json \
   output_bench_parallel/

python3 phylogeny.py --output output_bench_parallel --phylo-workers 4 &
PHYLO_PID=$!
python3 /tmp/monitor_rss.py $PHYLO_PID 30 /tmp/rss_samples.csv &
wait $PHYLO_PID
```

- [ ] **Step 5: Extract wall-clock and peak-RSS-during-overlap numbers**

Wall-clock: same method as §2's original derivation — diff timestamps between the
"analyzing N group(s)" line and the final "Wrote N rows" line in
`output_bench_parallel/ortho_phylogeny.log`, per scope.

Peak RSS during the giant-overlap windows (spec §5): cross-reference
`/tmp/rss_samples.csv`'s timestamps against the "aligning + building trees" log lines
for puumalaense/hantanense/tulaense (species scope) and arvicolinae_borne/murinae_borne
(phylogroup scope) — find each group's start time from the log, and its end time from
the next group's start time (or the scope's final write line for the last group to
finish). The overlap window is the union of the giants' [start, end) intervals; report
`max(total_rss_kb)` for samples falling inside it.

- [ ] **Step 6: Compare against the spec's predictions and decide**

Compare measured wall-clock against §2's ~98min (species)/~212min (phylogroup)/~310min
(total) predictions, and measured peak RSS during overlap against available RAM on the
benchmark machine.

- If wall-clock is close to predictions and peak RSS during overlap stays well under
  available RAM: no code change needed. Record the measured numbers in this plan's Task
  3 completion note (or the branch's SDD ledger, if one exists for this plan) for future
  reference.
- If the machine OOMs, swaps heavily, or peak RSS during overlap is uncomfortably close
  to available RAM (the rollback condition, spec §5): lower the shipped default in
  `phylogeny.py:973`'s `--phylo-workers` `default=4` to the measured-safe value, and
  update the spec's §5 with the actual numbers that drove the decision. This is a
  data-driven change made after Step 5's measurement — do not pre-guess a lower default
  now.

- [ ] **Step 7: Commit (only if Step 6 changed the default)**

```bash
git add phylogeny.py docs/superpowers/specs/2026-09-10-phylogeny-group-parallelism-design.md
git commit -m "phylogeny.py: lower --phylo-workers default to N per benchmark (measured peak RSS during giant-overlap window)"
```

If Step 6 required no change, there is nothing to commit for this task — the benchmark
result itself (numbers, not code) is the deliverable; record it wherever this branch
already tracks such measurements (e.g. an SDD ledger, or a note added to the spec's §2).
