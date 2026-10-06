# Ortho_rerun08 vs main-branch run: winner comparison and reproducibility verdict

**Written**: 2026-10-06  
**Purpose**: Documents how GA winners changed between the 08-work branch run and the
current main-branch run, and explains why. Required if a reviewer asks "did parameter
changes alter your conclusions?"

**Preserved outputs**: `../_pipeline_archive/Ortho_rerun08_full_output.tar.gz`
(1627 files, 106 MB — full `output/` tree including all_scored CSVs, figures, phylogenetics,
gate4 results).  
**Code recovery**: `git checkout bc1b117` or `git checkout remotes/origin/08-work`.

---

## What the two runs are

| Label | Code at run time | Git ref |
|---|---|---|
| **Ortho_rerun08** | 08-work branch, stopped at merge base | `bc1b117` |
| **main_branch** | current main (+21 commits past merge base) | `6a5b56e` |

A third artifact, `_pipeline_archive/run_20261002_003026`, has `generations_run=2` —
it is a partially-started run, not a completed GA result, and is excluded from this
comparison.

GA parameters were identical in both completed runs: `--seed 20260904 --population 150
--generations 100 --plateau 20`. Pool sizes were identical: terminal_complete=272,
strict_both=94. Calibration (`data/rpi_calibration.json`, τ=3.914) was identical.

---

## Three scoring changes between the runs

Between `bc1b117` and `6a5b56e`, three commits materially changed the scoring landscape:

### 1. `def19ff` — median fallbacks for geo / time / M_pack

`rpi_core.py, data/rpi_calibration.json: median fallbacks for geo/time/M_pack`

Changed the imputed value for isolates with missing geographic origin, collection year, or
M_pack energy from **1.0** (the maximum multiplicative value, artificially favouring
data-sparse candidates) to **calibrated medians** from the panel:

| Factor | Old fallback | New fallback |
|---|---|---|
| geo (missing country) | 1.0 | 0.051 (panel median) |
| time (missing year) | 1.0 | 0.787 (panel median) |
| M_pack (no parental sequence) | 1.0 | 0.664 (panel median) |

This is the primary driver of winner divergence. Candidates that ranked highly under the
old scoring because they benefited from the 1.0 geo/time bonus are no longer competitive.

### 2. `dabeee7` — neutral-factor tracking

`Algorithm.py, figures.py: make neutral-1.0 fallbacks distinguishable in GA outputs`

Added the `rpi_neutral_factors` field to winner records. This is a **diagnostic** change
only — it does not alter scores, but makes fallback usage visible. The presence of
`M_pol:no_signal_beyond_tag` on all main-branch terminal_complete winners indicates their
M_pol component is computed from the 5′ tag only (no coding-region polarity signal), which
is an honest limitation of those triplets, not a scoring error.

### 3. `558e49a` — single-species triplet filter

`Algorithm.py: skip single-species triplets before scoring (step 7)`

Removed same-species triplets from GA evaluation entirely. These were previously scored and
could appear as winners; they are now skipped before any scoring call. Effect on the
specific winners found was secondary to the fallback change, since both runs' top winners
were already cross-species.

---

## Winner tables

### terminal_complete (n=272)

| Metric | Ortho_rerun08 | main_branch |
|---|---|---|
| Generations to plateau | 21 | 23 |
| Best fitness | **0.4717** | **0.4717** |
| Dominant isolates | `84fli` (8/10 winners) | `js9`/`js5` × `hbt43 2012` |
| Species cross | puumalaense × seoulense × hantanense (mixed) | hantanense-S × seoulense-M+L (uniform) |
| Neutral factors | `''` (M_pol from signal) | `M_pol:no_signal_beyond_tag` (all 10) |
| gate4 winner-0 | `supported_reassortment_candidate` | `weak_signal_low_bootstrap_support` |
| gate4 summary | 1× supported, 8× no_conflict, 1× weak | 0× supported, 9× weak, 1× no_conflict |

Ortho_rerun08 winner-0 (`lr1 / hbt4 2016 / hunan2020c017 10`) cleared the phylogenetic
bootstrap threshold and was called `supported_reassortment_candidate`. That triplet no
longer wins under the corrected scoring; the main-branch winners' M_pol component is at
fallback, so the phylogenetic bootstrap signal is weaker.

### strict_both (n=94)

| Metric | Ortho_rerun08 | main_branch |
|---|---|---|
| Generations to plateau | 21 | 23 |
| Best fitness | **0.4717** | **0.1089** |
| Fitness drop | — | **~4.3×** |
| Species cross | puumalaense × seoulense × hantanense (mixed) | aa14/rn11 × hbq65 + hantanense×seoulense |
| Neutral factors | `''` | `E:host;M_pol:no_signal_beyond_tag` |
| gate4 winner-0 | `weak_signal_low_bootstrap_support` | `weak_signal_low_bootstrap_support` |

The strict_both fitness collapse (0.4717 → 0.1089) is the clearest diagnostic of the
fallback change. The strict_both pool is smaller and skews toward isolates with sparser
geo/time metadata; the geo fallback drop from 1.0 to 0.051 deflates the reachable fitness
ceiling far more than in terminal_complete.

---

## Reproducibility verdict

The winner lists are **not reproducible across these two code versions**, but this is
expected and correct behaviour:

1. **The score function changed** — the Ortho_rerun08 run was optimising a score that
   over-rewarded data-sparse candidates. `def19ff` corrected this with calibrated medians.
   Both runs converge cleanly (plateau within 20–23 generations on the same seed), so
   neither is stuck; they are genuinely at the landscape peak for their respective score
   functions.

2. **The main-branch result is the one to report** — it uses the more defensible scoring.
   The strict_both fitness of 0.1089 reflects honest signal in a harder, smaller pool under
   stricter scoring.

3. **No `supported_reassortment_candidate` in the current run** — this is not a regression.
   The triplet that earned that call in the old run depended on a scoring advantage from
   the 1.0 fallback. The current winners' gate4 verdicts of `weak_signal_low_bootstrap_support`
   accurately reflect the bootstrap evidence.

4. **Within-version reproducibility is untested** — running the GA twice on the same code
   version and same seed should produce identical results (the GA is seeded), but this has
   not been verified with a second run.
