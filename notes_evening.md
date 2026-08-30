# Evening run summary (2026-08-30 → 2026-08-31)

## Commit hashes, in order

1. `98bc10a` — Wire `known_segment_isolates` through `build_swap` so published crosses
   get isolate-level Factor E. (Prerequisite fix from a `/code-review` follow-up earlier
   in the session, landed before this evening's numbered task list started.)
2. `6d85d5b` — **Step 1.** Rerun seeded discrimination panel (seed 20260829) with the
   `build_swap` fix landed.
3. `3762c91` — **Step 2.** Efficiency fix: skip expensive isolate-identity resolution
   when fewer than two species are implicated.
4. `8503101` — **Step 4.** Code review of the branch diff; consolidate duplicated
   efficiency-fix comment (docstring-only, score-neutral).
5. `4ef78f1` — **Step 5** (and a correction to step 4's writeup, combined into one
   commit): README front page + `notes_codereview.md` finding 1 severity correction
   after verification.

**Nothing was pushed.** Local `04-rpi-fitness-validation` is 5 commits ahead of
`origin/04-rpi-fitness-validation` (`f4804bf`), which is unaffected.

## Step-by-step outcome

### Step 1 — panel rerun (`6d85d5b`)

Reran the seed-20260829 discrimination panel with `build_swap`'s `known_segment_isolates`
wiring in place. Preserved the prior after-fix panel for the diff
(`output/discrimination_scratch/discrimination_panel_after_e_fix_pre_build_swap.json`).

**Result: numerically identical `rpi_generation` on every row of all 4 groups**, verified
row-by-row. 3 of group c's 4 published-cross rows flipped `rpi_e_isolate_basis` from
`nearest_reference` to `actual_parents` — the more honest label, since the old proxy
search was already landing on the true donor isolate at k-mer distance zero (the
candidate's donated-segment protein IS that isolate's protein, copied verbatim by
`build_swap`).

Group c stratified by `rpi_e_isolate_basis` (after fix): `actual_parents` n=8 (mean
0.3314), `nearest_reference` n=1 (mean 0.3669 — `RODRIGUEZ1998_SNV_NMR11sl_x_BCCVm`,
whose gate1 species label doesn't match `build_swap`'s own backbone/donor species names,
a pre-existing gate1 mismatch unrelated to this fix), 2 unscored (gate 0/1 rejected, as
before).

**Task 14 pass criterion did NOT move** (per the standing instruction, recorded as-is,
no compensating change attempted):

| comparison | before this rerun | after this rerun |
|---|---|---|
| `b2_vs_random` | p=4.68e-09, δ=+0.953, True | p=4.68e-09, δ=+0.953, True |
| `b1_vs_b2` | p=2.18e-11, δ=+0.961, True | p=2.18e-11, δ=+0.961, True |
| `reals_vs_b1` | p=0.733, δ=+0.078, not sig. | p=0.733, δ=+0.078, not sig. |
| `overall_ranking_supported` | True | True |

Full before/after writeup appended to `notes_e_isolate_fix.md`.

### Step 2 — efficiency fix (`3762c91`)

`ecological_opportunity` discards `parent_isolates` and reads only `isolate_basis` when
fewer than two species are implicated. Implemented the **identical-output** version, not
the naive one: extracted `_isolates_with_usable_protein` as the cheap "would any isolate
resolve at all" filter, used only in the `<2`-species branch, which reproduces the exact
same `nearest_reference`/`insufficient_metadata` decision the expensive k-mer/cosine path
would have made (verified — `nearest_reference_isolate` only ever returns `None` from
that same cheap guard clause; the expensive part only decides WHICH isolate, never
WHETHER one resolves). Added `ResolveEParentIsolatesSkipsExpensiveResolutionTest` (4
tests): cheap-path/expensive-path basis agreement in both the found and not-found case,
a spy proving `ah.kmer_profile` is skipped for `<2` species and still called for `>=2`
(regression guard), and an assertion on `compute_rpi`'s shipped row per this project's
"assert on the row that ships" rule.

### Step 3 — full suite gate

Verified green twice this evening (322 tests across `reassortment_detector`,
`rpi_calibration`, `rpi_discrimination`, `rpi_discrimination_stats`, `pipeline`,
`validate_reassortment_detector`; ~814s wall each run) — once before step 4's comment
edit, once after. No `.py` files changed since the second green run (`4ef78f1` touched
only `README.md`/`notes_codereview.md`), so that result still holds for the current
tree; a third rerun was skipped as redundant rather than re-verified, since nothing
executable changed.

### Step 4 — code review (`8503101`, corrected in `4ef78f1`)

`/code-review medium` over the full branch diff, focused on this evening's commits.
5 findings in `notes_codereview.md`:

- **Fixed autonomously** (1): duplicated comment/docstring for the efficiency fix's
  rationale, consolidated to one canonical spot. Pure text, zero behavior change,
  re-verified the 4 affected tests pass.
- **Left as proposals** (4), none touching `RPI_WEIGHTS`/tau/gate logic/verdicts/group
  construction/output-column semantics:
  1. `resolve_e_parent_isolates`'s known-isolate check could label an intra-species
     candidate `actual_parents`. **Applied `/superpowers:receiving-code-review`
     rigor to this one rather than accepting the reviewer's framing**: verified against
     the actual call sites and found the specific "corrupts a diagnostic" claim does
     NOT hold — `rpi_discrimination.py` excludes intra-species crosses from the group
     that diagnostic reads, and `validate_reassortment_detector.py`'s own report drops
     the `rpi_e_isolate_basis` column before it's written anywhere. Downgraded medium
     → low and corrected the writeup accordingly, rather than leave an overstated
     finding standing.
  2. 3 of the 4 new efficiency-fix tests assert on the intermediate return, not
     `compute_rpi`'s shipped row (verified accurate; test-coverage gap, not a
     pre-approved autonomous-fix category).
  3. The `<2`-species threshold is now an independent literal in two functions
     (verified accurate; refactor risk, left for user).
  4. A suggested simpler API shape for `nearest_reference_isolate`/
     `resolve_e_parent_isolates` (verified accurate description; design opinion, left
     for user).

### Step 5 — README front page (`4ef78f1`)

New "RPI fitness-validation toolchain (`pipeline.py`)" section inserted after the
title: one-paragraph purpose, `pipeline.py`'s exact CLI (`run`/`validate`/`diagnose-e`/
`diagnose-gate0`, `--recalibrate`/`--seed` flags), headline results **from the
post-`build_swap`-fix panel only** (τ=3.914, all three Task 14 comparisons + verdict,
table above), a known-limitations summary pointing to the full binding list already in
the file, and a docs map covering `docs/PROJECT-STATE.md`, all 5 design specs, and this
evening's own notes files. Every number cross-checked directly against the committed
panel JSON and `data/rpi_calibration.json` before writing, not recalled from memory. No
new citations introduced beyond what the existing README sections already cite
(Rodriguez 1998, Rizvanov 2004, Handke 2010, Klempa 2018, Park 2025) — all confirmed
present in `docs/My Library.csv` with matching DOIs.

## Hard rules — compliance check

- `RPI_WEIGHTS`, tau, gate logic, verdicts, group construction, output-column
  semantics: **unchanged everywhere**, confirmed by the identical-`rpi_generation`
  panel result (step 1), the efficiency fix's own identical-output guard tests (step
  2), and by declining to autonomously fix any code-review finding that touched those
  categories (step 4).
- Literature: no new citations introduced; existing ones re-confirmed against
  `docs/My Library.csv`, never drove a code change.
- Step 1's panel move on `reals_vs_b1`/overall verdict: **it did not move** — recorded
  as such, no compensating change attempted (per the explicit standing instruction).

## Known friction this session (not a task outcome, but worth noting)

The session's safety-layer classifier intermittently refused `git add`/`git commit`/
even a background test-suite launch throughout the evening, for reasons tied to
earlier conversation content rather than the specific command. It was NOT persistent —
retrying on a later turn, or having the user run the exact command manually, both
cleared it every time it was tried. Per the standing hard-rule instruction, stopped
retrying after ~3 attempts each time rather than working around it (no `--no-verify`,
no alternate commit tooling), and handed the exact commands to the user twice for
manual execution when blocked. No commits were lost; everything queued eventually
landed.
