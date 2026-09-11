# phylogeny.py group-level parallelism — design

**Written 2026-09-10. Branch `06-module-final-simplification`.**

## 1. Purpose

`_run_phylogenetics_scope` (phylogeny.py:799) loops over species/phylogroup groups
strictly serially, calling `_analyze_group` (3x mafft + 3x iqtree subprocess calls per
group, one per S/M/L segment) one group at a time. Groups are fully independent —
nothing about one group's alignment/tree depends on another's. This plan parallelizes
that loop across groups to reduce Stage 2 wall-clock time, without changing anything
about what any individual group computes.

## 2. Background: what the 2026-09-06 log actually shows

`output/ortho_phylogeny.log` (from the run specified in
`docs/superpowers/specs/2026-09-06-overnight-full-pipeline-rerun-design.md`) is the only
full-run log with per-group timestamps present in this checkout. Per-group durations
derived from its timestamps:

| Scope | Groups | Total | Dominant groups | Their share |
|---|---|---|---|---|
| species | 14 | ~201 min | puumalaense (425 taxa, 98 min), hantanense (273 taxa, 37 min), tulaense (316 taxa, 37 min) | 172/201 min (86%) — the other 11 groups total ~28 min combined |
| phylogroup | 4 | ~359 min | arvicolinae_borne (755 taxa, 212 min), murinae_borne (550 taxa, 119 min) | 331/359 min (92%) |

This is the load-bearing fact for the whole design: **runtime is dominated by a handful
of giant groups, not by having many groups to churn through.** Two consequences:

- The makespan floor for a scope is bounded below by its single largest group's own
  runtime — cross-group parallelism cannot make puumalaense's own alignment+tree faster;
  only `--phylo-threads` (already near-maxed by default) touches that.
- Given the group-size distribution above, **expected speedup is ~2.0x for species scope
  (201→98 min), ~1.7x for phylogroup scope (359→212 min), ~1.8x overall (560→310 min)**.
  `--phylo-workers=4` should be judged against this ~1.8x ceiling, not against a naive 4x
  — 4 workers is already enough to saturate the available parallelism in this dataset
  (phylogroup scope has exactly 4 groups; species scope's small groups total only ~28 min
  of work, not enough to usefully occupy a 5th+ worker). **Workers beyond 4 buy ~nothing
  measured on this dataset** — document this as a measured property, not a guess, once
  the benchmark task (§7) confirms it.

## 3. Non-goals

- **Any within-group speedup** (MAFFT strategy changes, additional threads beyond what
  `--phylo-threads` already provides, alternate aligners). The makespan floor is a single
  group's own runtime — shrinking it further would require changing what MAFFT/IQ-TREE
  actually compute, which risks the byte-identity contract (§5) for no gain this plan is
  scoped to capture. Out of scope by design, not an oversight.
- Cross-scope parallelism (species scope and phylogroup scope running concurrently).
  Scopes stay sequential; only the group loop *within* each scope's call parallelizes.
- Auto-dividing `--phylo-threads` by `--phylo-workers`. Considered and explicitly
  rejected (see §4.3) — the user chose to leave `--phylo-threads`' default formula
  unchanged and accept the oversubscription risk, tuning manually if needed.
- Recalibration, code simplification elsewhere in phylogeny.py, or any change to
  `_analyze_group`'s internals. This plan touches exactly the submission/collection
  mechanics around the group loop.

## 4. Design

### 4.1 Architecture

Replace the serial loop in `_run_phylogenetics_scope`:

```python
for group_key, members in sorted(by_group.items()):
    rows.extend(_analyze_group(group_key, scope, members, ...))
```

with a `ThreadPoolExecutor(max_workers=args.phylo_workers)`. `_analyze_group` itself is
unchanged — it already wraps its entire body (mafft/iqtree calls, tree parsing,
RF-distance computation, classification) in one try/except, so submitting it as a task
preserves its existing per-group error isolation for free.

**Submission order: largest-taxa-first**, not the current alphabetical
`sorted(by_group.items())`. This is the standard "longest processing time first"
scheduling heuristic: with a fixed worker count and the extreme size skew measured in
§2, submitting the biggest groups first means they start immediately rather than
possibly queuing behind small ones, which minimizes total makespan. Concretely: build
the submission list sorted by `len(members)` descending (ties broken by `group_key` for
determinism), submit in that order to the executor.

**Output order: still sorted by `group_key`.** Collect results via `as_completed`, then
sort the accumulated rows by `group_key` before returning — output row order (and
therefore `ortho_reassortment_table.csv`'s row order) stays identical to today's serial
run regardless of worker count or actual completion order. Submission order (§ above,
performance-only) and output order (this paragraph, determinism-only) are independent
concerns; keep them that way in the implementation.

Scopes stay sequential: `--phylo-scope=both` runs all species-scope groups to
completion (in parallel), then all phylogroup-scope groups (in parallel).

### 4.2 New CLI flag

`--phylo-workers`, `type=int`, default `4`. Validated `>=1` at parse time — reject
0/negative with a clear argparse error, not a runtime crash. Help text states the
oversubscription tradeoff with `--phylo-threads` explicitly (see §4.3), since this plan
does not auto-mitigate it.

### 4.3 `--phylo-threads` interaction (accepted risk)

`--phylo-threads` already defaults to `max(1, (os.cpu_count() or 4) - 1)` — sized
assuming one group runs at a time. With `--phylo-workers=4`, up to 4 groups' mafft/iqtree
calls can run concurrently, each requesting close to that many threads — real
oversubscription on typical machines. Considered auto-dividing
(`--phylo-threads` default becomes `cpu_count // phylo_workers` when not explicitly
set) and explicitly **rejected**: the user chose to leave `--phylo-threads`' formula
untouched, accepting the risk, and tune `--phylo-threads` down manually when running with
`--phylo-workers > 1` if needed. This is documented in the flag's help text, not silently
handled.

This risk is **not hypothetical on this dataset** (§2): the groups most likely to
overlap under `--phylo-workers=4` are exactly the giants (puumalaense, hantanense,
tulaense in species scope; arvicolinae_borne, murinae_borne in phylogroup scope), since
small groups free their worker within a minute while giants keep running for tens of
minutes to hours — meaning oversubscription would hit precisely the long-running,
expensive alignments, not the quick ones. This is why the benchmark (§7) specifically
instruments the giant-overlap windows rather than treating the whole run as one uniform
measurement.

### 4.4 Correctness contract: same output regardless of scheduling

Per-group mafft/iqtree invocations receive identical arguments whether they run serially
or in parallel — parallelizing changes *when* a group runs, never *what command* runs for
it. On a fixed machine/library-version/seed, this plan claims **byte-identical** output
between a serial (`--phylo-workers 1`) and parallel (`--phylo-workers N`) run on the same
input, following this repo's existing convention (see
`docs/superpowers/specs/2026-08-29-rpi-fitness-function-design.md`: byte-identity is
claimed only within a machine/version, not across OSes/builds) and its established
verification method (`docs/superpowers/specs/2026-09-06-overnight-full-pipeline-rerun-design.md`
§4: checksum diff first, semantic diff — values not just bytes — as fallback on anything
not byte-identical, since IQ-TREE could in principle differ run-to-run even under a fixed
seed for reasons unrelated to this change).

**Unit test coverage note (added after implementation):** `test_worker_count_does_not_change_output_rows`
covers row plumbing only — same rows, correctly sorted by `group_key` — using the fake
mafft/iqtree stubs. It cannot detect numeric violations of this contract (e.g., different
RF-distance values, different alignment columns) because the stubs produce identical
fixed output regardless of worker count. The full byte-identity contract is validated at
Task 3 (§7), which diffs the real `ortho_reassortment_table.csv` and tree files. Note for
Task 3: to isolate concurrent-run variability from other sources, set `OMP_NUM_THREADS=1`
(or `OMP_DYNAMIC=false`) on the benchmark machine before comparing runs.

### 4.5 Error handling

`_analyze_group`'s existing try/except already produces `status='failed'` rows on any
internal failure — unchanged, not retested by this plan (it predates this work). One
defensive addition: if `future.result()` itself raises (something outside
`_analyze_group`'s own guard — a genuinely unexpected bug), that is caught at the
collection point, logged, and recorded as a failed-group row rather than aborting the
whole run — matching the existing "one group's failure never aborts the run" guarantee
under the new concurrent path.

## 5. Safety: OOM rollback condition

Given §4.3's accepted oversubscription risk, and that the highest-risk windows are known
in advance (§2's giant-overlap groups), the benchmark task (§7) must specifically measure
**peak RSS during the giant-overlap windows**:
species scope: the interval where puumalaense, hantanense, and tulaense may run
concurrently; phylogroup scope: the interval where arvicolinae_borne and murinae_borne
may run concurrently. This is where a real OOM condition would fire, not during any
small-group window — measuring average or whole-run peak RSS would dilute the signal.

**Rollback condition:** if the benchmark shows the machine actually OOMs, swaps
heavily, or peak RSS during a giant-overlap window exceeds a safe margin of available
RAM, the plan does not ship `--phylo-workers=4` as the default — the finishing task
lowers the shipped default (to whatever the benchmark shows is safe, e.g. 2, or 1 with a
documented recommendation to raise it manually on larger machines) rather than merging
a default that OOMs the reference environment. This is a measured decision made by the
benchmark task, not a pre-committed number in this spec.

## 6. Testing plan (unit-level, no external tools required)

Existing tests (`PhylogeneticsScopeEndToEndTest` et al.) pin `phylo_threads=1` via a fake
mafft/iqtree stub already; default their fixture's `phylo_workers` to `1` too, so today's
tests keep testing today's exact behavior (the degenerate single-worker case through the
new code path) without new flakiness.

New unit tests:
- `phylo_workers=4` produces output identical (same rows, correctly sorted by
  `group_key`) to `phylo_workers=1` on the same fixture input — determinism check using
  the fake mafft/iqtree stubs (no real tool timing involved).
- Largest-first submission order: given groups of known differing sizes, assert the
  executor receives them in descending-size order (mock/spy on the submission call, not a
  timing-based test — timing-based concurrency tests are flaky by nature).
- One group's induced failure (fake mafft/iqtree stub that raises) doesn't abort other
  groups' results — extends the existing per-group error-isolation guarantee to the
  concurrent path.
- `--phylo-workers 0` and negative values rejected at argparse time with a clear error.

## 7. Benchmark task (real tools, real data — separate from unit tests)

This task runs against actual mafft/iqtree and the real dataset, so it requires an
environment where those are installed (this development environment currently does not
have them — see Task 15's finding in
`.superpowers/sdd/2026-09-09-stage1-reference-snapshot/progress.md`). It is a distinct
task from §6's unit tests and should be scheduled for an environment that has mafft/
iqtree on PATH.

1. **Verify the log is still a valid baseline before trusting it.** The 2026-09-06 log
   predates this session's Stage 1 snapshot refactor (`prepare_dataset.py`'s
   `ortho_reference_isolates/genomic/cds.json` snapshot, `wild_classified_pool`, the I1-I6
   review fixes). Confirm the current pipeline still produces the same group taxa counts
   the log reports (425 puumalaense, 273 hantanense, 316 tulaense, 755 arvicolinae_borne,
   550 murinae_borne, and the rest) before treating the log's timings as the serial
   baseline to compare against. If counts differ, re-derive the baseline numbers in §2
   from a fresh serial run rather than trusting the stale log.
2. **Run the parallel version only** (`--phylo-workers 4`) against the same input the
   2026-09-06 run used, and compare its wall-clock and output against the log/baseline —
   do not re-run the serial path from scratch (the log already is that baseline, once
   step 1 confirms it's still valid).
3. Measure wall-clock for species scope, phylogroup scope, and total; compare against
   §2's ~98/~212/~310-minute floor predictions.
4. Measure peak RSS specifically during the two giant-overlap windows (§5) — this is
   the evidence the rollback condition acts on.
5. Diff the parallel run's `ortho_reassortment_table.csv` and tree files against the
   serial baseline's: checksum diff first, semantic diff (values, not just bytes) as
   fallback per §4.4's contract.
6. If the rollback condition (§5) fires, lower the shipped default and document why with
   the actual measured numbers, rather than merging `--phylo-workers=4` unconditionally.

## 8. Acceptance criteria

- [ ] `--phylo-workers` flag exists, defaults to 4 (or the benchmark-adjusted value per
      §5/§7.6), validated `>=1`.
- [ ] Groups submitted largest-taxa-first; output rows still sorted by `group_key`.
- [ ] `--phylo-threads`' default formula is unchanged; help text documents the
      oversubscription tradeoff.
- [ ] One group's failure doesn't abort the run (existing guarantee, extended to the
      concurrent path with a test).
- [ ] Unit tests (§6) pass; existing test suite (`phylogeny.py`, currently 27 tests)
      stays green with `phylo_workers=1` fixtures.
- [ ] Benchmark task (§7) run in an environment with mafft/iqtree, confirming measured
      speedup and peak-RSS-during-overlap numbers, with the rollback condition evaluated
      against real measurements before the default ships.
