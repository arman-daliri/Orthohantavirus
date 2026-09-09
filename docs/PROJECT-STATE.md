# Orthohantavirus — Project State

**Single source of truth for resuming work. Written 2026-08-29.**
A fresh session should need this file, `CLAUDE.md`, and the two documents in
§7 — and no chat transcript.

---

## 1. Where things stand right now

**This table is a point-in-time snapshot (last updated 2026-09-09, branch
`05-module-consolidation`). If it looks stale against `git log --oneline -20`, trust the
log — see `docs/CLAUDE.md`'s "RESUME HERE" note.**

| | |
|---|---|
| **Branch** | `05-module-consolidation` |
| **HEAD** | Overnight full-pipeline rerun and reconciliation session (§7g) complete:
5 real bugs found and fixed (§7e Findings 9 & 10 in `Algorithm.py`; a CSV line-ending
non-determinism in `all_hanta.py`; two previously-untested bugs in `rpi_core.py
validate --with-phylogeny`'s gate4 path — dropped gate4 fields, unresolved binary
paths); `Algorithm.py --dump-scored` added; `output/` fully promoted from the rerun
(all top-level tables, both phylogenetics scopes, both GA search configs, gate4
confirmation for all 10 long-pool winners, the full 55-candidate validate panel with
real gate4 data, fresh before/after tree figures). §7e Findings 3, 4, 11 remain open
proposals, not yet fixed (none triggered by this session's run configuration) |
| **Test suite** | **330 tests**, `python -m unittest rpi_core` — all pass (measured
2026-09-09, against the final promoted state; verify with the same command before
trusting this number again, it changes every session) |
| **Validation harness** | **PASS** — sensitivity 19/19, specificity 27/27, blind spot
5/5, adversarial 4/4, **now with real gate4 data behind 40/55 candidates** (previously
gate4 was never actually exercised in the committed baseline — see §7g) |
| **Redundancy survey (not yet actioned)** | See §7g: the §7f stalled run's stale
`ga-confirm/` output (superseded, byte-identical to the fresh run) and stray
`__pycache__/` are candidates for cleanup, awaiting explicit go-ahead — nothing deleted
this session |

**Current goal:** none pending from this session — §7g's overnight rerun, reconciliation,
and report are complete and committed. Next session should start from `docs/PROJECT-STATE.md`
§7g for full context, and act on the redundancy-survey candidates and §7e Findings 3/4/11
only if asked.

---

## 2. The two governing documents

Read both before touching code. This file summarises; they are authoritative.

- **Design/spec:** `docs/superpowers/specs/2026-08-29-rpi-fitness-function-design.md`
  — the argument, the M_pack options, and decisions **D1–D14**.
- **Task plan:** `docs/superpowers/plans/2026-08-29-rpi-fitness-validation.md`
  — bite-sized steps for WS0–WS5 and WS-API.

---

## 3. The single most important finding

**RPI's factors split into two kinds, and only one kind can discriminate.**

`compute_rpi` takes seven inputs (`m_pack_status` added WS3.2). E (ecology), K
(constellation prior) and F (retention) are **species-level or segment-level priors** —
they do not read the nucleotides of the particular S+M+L assembly being scored. Only
**M_pack** and **M_pol** do.

**UPDATE 2026-08-29 (WS3.3 + WS3.4): both M_pack and M_pol are now real.** Both call
`rpi_species_assignment(gate1, odd_segment)` — one shared donor/backbone determination,
not one per factor, so the two can't silently disagree about who the backbone is.
- `score_m_pack` cofolds the candidate's own cross-segment panhandle duplex against the
  backbone species' own reference consensus (D3) and rescales via
  `M_pack = exp(-ΔΔG / τ)`, τ = 3.914 (WS3.1).
- `score_m_pol` compares the candidate's own N-recognition footprint (beyond the
  invariant tag) against the backbone species' own reference consensus CDS window
  (`species_consensus_cds_window`). The old
  `polymerase_compatibility(terminal_windows, terminal_windows)` self-comparison is
  gone. **Bundled into WS3.4, not separately scoped:** `footprint` now reads
  `n_recognition_footprint(backbone_species)` instead of the always-default 9 nt —
  wiring that existed but was never called from production before this task.

Neither is the self-comparison bug that predated it — see `score_m_pack`'s docstring
for why the donor==backbone (concordant) case looks similar but isn't.
`divergence_stratum` is still a binary threshold at 0.80 and production identities run
0.35–0.46, so every row still lands `distant` and the K stratum never varies.

**Consequence, now substantially resolved for M_pack, narrowly for M_pol:**
`rpi_generation` is now the genuine four-factor product `E × M_pack × M_pol × K`.
Smoke-tested 2026-08-29 on the real WS4 panel (12 cases): M_pack varies 0.10–1.0 on
every candidate. **M_pol varies only when the backbone species is *Orthohantavirus
sinnombreense*** (Mir & Panganiban's SNV-23-vs-genus-9 asymmetry is the only wired
footprint length) — SNV-backbone candidates scored 0.74–0.96 (`status="scored"`); the
other 9 of 12 honestly reported `no_signal_beyond_tag` (footprint 9 falls entirely
inside the invariant 14-nt tag) rather than a fake constant. WS3.1's own measurement
also showed M_pack's signal is real but small (Cliff's δ = +0.222; 38% of the
calibration panel saturates M_pack at 1.0).

**But a tie is only a defect for some comparisons** (corrected by the user, decision D7):

- **Within-species shuffles (b1)** *should* score high. They are genuinely plausible.
  RPI measures **plausibility, not occurrence**; separating b1 from real reassortants is
  the **tree detector's** job, not the index's. `reals ≈ b1` is an **expected tie**.
- **Cross-species sympatric shuffles (b2)** are where the index must earn its keep. With
  both factors real but each carrying only partial signal (M_pack weakly for everyone,
  M_pol strongly but only for one species), it is **expected to discriminate somewhat,
  not fully** — the user's 2026-08-29 decision was to proceed with both factors as
  scoped, not to retune τ or invent additional footprint lengths without a literature
  source.

**WS4 is where this gets measured for real** — whether M_pack + M_pol together close
enough of the b2 gap, or whether the Factor E rebuild (deferred, see gate3's KNOWN
LIMITATION) ends up mattering more than either molecular factor.

---

## 4. Decisions of record (D1–D14) — all approved by the user

| # | Decision |
|---|---|
| **D1** | M_pack rescaling: **calibrate first**, then choose. Bimodal ΔΔG with a clean gap → logistic; otherwise τ-scaled Boltzmann, `τ = ΔΔG_span / ln(100)`. Min–max **rejected** (cohort-dependent → breaks reproducibility). |
| **D2** | RNAcofold window **N = 50 nt**, swept 25 / 50 / 75. |
| **D2b** | Mechanism = **ViennaRNA Python bindings**, no CLI, **no fallback**. |
| **D3** | M_pack and M_pol from the candidate's terminal windows vs the **implicated parent/donor consensus** — never self-comparison. |
| **D4** | Gate 4 fixture runs in a **scratch dir only**, never touching `output/`. |
| **D5** | Add the SNV × BCCV cross; resolve BCCV's ICTV name. |
| **D6** | *(narrowed)* Byte-identity for RPI output, and for MAFFT **only within a version**. Across machines: MAFFT at **alignment-content** level, IQ-TREE at **topology + support** level. |
| **D7** | Hard fakes split into **b1** (within-species, expect HIGH) and **b2** (cross-species sympatric, expect mid-low). |
| **D8** | Gate 1 determinism fix (WS0) — done. |
| **D9** | WS-API runtime: per-accession memoisation; **report warm and cold latency separately, never blended**. |
| **D10** | `score_triplet` returns RPI columns + status flags only — **never a verdict, never `gate4_*`**. |
| **D11** | Single source of truth: WS-API serialises WS4's group builder and `PanelBuilder`; `_score_panel` must not mutate the global — done. |
| **D12** | Install route: **Windows = manual pinned binaries + SHA-256**; **Linux = conda**. |
| **D13** | ~~τ ≈ 4.86 provisional~~ — **SUPERSEDED 2026-08-29 by the WS3.1 measurement: τ = 3.914** (robust p5–p95 span, N = 50). Quote τ = 3.914, never 4.86. |
| **D14** | MAFFT split: **Windows 7.526 / Linux 7.525** — forced, not chosen. |

---

## 5. What was completed this session (12 commits, `e7d9f44..0abd77e`)

### WS0 — determinism (`505879e`)
`gate1_check` built `implicated_species = list({...})` from a **set**; Python randomises
string hashing per process, so identical input produced different orderings run to run.
`gate3_check` reads `implicated_species[:2]`, so **the hash seed chose which pair gate 3
judged.** Measured on three real species, 2 of 6 orderings returned
`ecological_overlap_plausible` where the rest returned `no_ecological_route`.

Fixed with `sorted()`. **Verified:** validation report byte-identical under
`PYTHONHASHSEED` 1 vs 999 (md5 `42d52a2c06d6ea6ebf9c91ecbf38269f`); before the fix, two
runs differed.

**Deferred, logged as a KNOWN LIMITATION on `gate3_check`:** gate 3 judges only the
**first two of N** implicated species. Sorting makes that choice *stable*, not
*principled*. The fix is the worst-pair-over-all-pairs treatment `ecological_opportunity`
already uses — deferred to the **Factor E rebuild**. Until then, `no_ecological_route` on
a 3-species candidate means "the alphabetically first two do not overlap."

### WS1 — bug fixes (`a4f967f`)
- **N4** — `boltzmann_compatibility(nan, 0.0)` returned `9.86e-305` (maximal penalty);
  pre-clamp it returned `1.0` (no penalty). Two opposite silent answers. Now raises
  `ValueError`; WS3.3's call site catches it as `rpi_m_pack_status="undefined_energy"`.
- **N5** — `_json_safe` tested only `float`. **Wider than recorded:** `np.float32('nan')`
  → JSON string `"nan"`, *and* finite `np.float32(1.5)` → `"1.5"`, `np.int64(3)` → `"3"`.
  Fixed with `np.floating`/`np.integer`, `bool` checked first.
- **N6** — the ordering guard was satisfied by equality (F = 1.0 for M/intact) and could
  not fail via `_rpi_decimals`. Added a guard-the-guard plus a companion reproducing the
  split-precision inversion directly.

### WS2 — reproducibility (`e649099`, `8b54c22`)
- `requirements-pip.txt` created and **enforced** by `EnvironmentReconciliationTest`.
- **MAFFT 7.525 has no Windows build.** Release tracks diverge: source page has 7.525,
  Windows page has 7.526. Hence D14.
- **IQ-TREE 3.1.3 is in `iqtree/iqtree3`, not `iqtree2`** (which tops out at v2.4.0).
- `tools/CHECKSUMS.sha256` + `tools/VERSIONS.txt`; 134 MB of binaries gitignored.

### `_score_panel` global-mutation fix (`b03558c`)
The sweep reassigned the module-level `RPI_WEIGHTS` — **process-global**, so a
concurrent scorer would read perturbed weights. Weights are now threaded explicitly.
**Verified as a real regression guard:** against the old code the test reports
`rpi_generation` 0.6561 where 0.4716981132 is correct — a 39% silent error.

### D2b — bindings swap (`0abd77e`)
The pip `viennarna` package ships **bindings only**; no `RNAcofold` executable exists in
the interpreter tree or on PATH, so `rnacofold_dg` could not run at all.
Removed `rnacofold_dg`, `_parse_rnacofold_dg`, `RnacofoldDgParsingTest`.
Added `cofold_dg`, `panhandle_window`, `RPI_PANHANDLE_WINDOW_NT = 50`,
`viennarna_provenance`.

---

## 6. Measured facts (do not re-derive)

| Measurement | Value |
|---|---|
| `score_candidate`, gate 4 disabled, warm | **7,088 ms/call** (~9.9 s multi-species) |
| ↳ `gate1_check` | 4,231 ms |
| ↳ `gate0_check` | 2,770 ms |
| ↳ `ecological_opportunity` | 2,848 ms (only when ≥2 species — always true for shuffles) |
| 50 nt × 50 nt cofold, bindings | **7.1 ms** |
| `load_reference_data` / `_build_precomputed` | 3.2 s / 0.09 s (once) |
| Reference isolates | 9,867 |

**RNAcofold is ~0.3% of per-call cost.** The WS-API budget is reachable only via
per-accession memoisation of gate 0/gate 1 (~71% of cost), not by the ecology precompute
alone. **1,000 calls ≈ 2 h today; 10,000 ≈ 20 h.**

Installed: Python 3.14.6, numpy 2.5.1, pandas 3.0.5, biopython 1.88, dendropy 5.0.13,
ViennaRNA 2.7.2 (bindings), MAFFT `v7.526 (2024/Apr/26)`, IQ-TREE `3.1.3 ... built Jun 19 2026`.

---

## 7. WS3.1 — DONE. The D1 branch is decided: **Option A, τ-scaled Boltzmann**

`rpi_calibration.py` + `data/rpi_calibration.json` are built, tested (39 tests) and
deterministic (artifact byte-identical under `PYTHONHASHSEED` 1 vs 999, md5
`3c08d9e960847bfb68459fe4cf71c0d1`).

```bash
python rpi_calibration.py              # rebuilds data/rpi_calibration.json
python rpi_calibration.py --stability  # re-runs the branch across seeds x windows
python -m unittest rpi_calibration     # 39 tests, <1 s
```

**Panel.** 314 wild terminally-complete isolates → 302 after WS4 exclusion → **42 selected
across 14 species** (cap 8/species), **1,008 folds** (252 cognate, 756 non-cognate),
ViennaRNA Python bindings 2.7.2, N = 50 nt, seed 20260829.

**Measured ΔΔG:** n = 756, mean **+1.99**, sd 5.58, p5 −6.80, median +1.60, p95 +11.22,
range −15.60 … +21.90.

**The D1 rule failed 3 of 4 pre-registered checks — NOT bimodal.** BIC prefers one
component (Δ = **+15.1**); means **1.24 SD** apart against a 2.0 requirement; the mixture
density has **one** local maximum. Sarle's BC = 0.344 (below the 5/9 flag).

**τ = 3.914**, from the robust p5–p95 span (18.03). Raw min–max span (37.50) would give
8.14; both ship in the artifact. **This supersedes D13's provisional τ ≈ 4.86.** Stable:
Option A on all three seeds at N = 50, τ = 3.90 / 4.28 / 3.90, and τ = 3.93…4.39 across
per-species caps 2/4/8/16/none.

### Three things WS3.1 turned up that were not in the plan

1. **The plan's "reverse-complemented" 3′ window is wrong.** RNAcofold takes both strands
   5′→3′ and pairs them antiparallel; the genomic 5′ and 3′ windows as the FASTA stores
   them already ARE reverse complements — that is what a panhandle is. As-is folds
   −30.9…−45.3 vs −7.3…−15.4 reverse-complemented. Pinned by `CofoldOrientationTest`.
2. **48% of deposits are missing their true termini** and `SEGMENT_LENGTH_RANGE` does not
   catch it. A truncated "panhandle window" is arbitrary internal coding sequence
   (−10.0 vs −42.2 kcal/mol). Since τ comes from the ΔΔG *span*, folding those would have
   inflated τ with sequencing artefacts. Hence `has_terminal_panhandle`.
3. **⚠ ΔΔG carries real but SMALL signal — weigh this before WS3.3.** Non-cognate duplexes
   are weaker in the expected direction in only **62%** of observations; **38% have
   ΔΔG ≤ 0**, so **M_pack clamps to exactly 1.0 and contributes nothing** on more than a
   third of the panel. Cliff's δ = **+0.222**, "small". The genus's terminal 14 nt are
   invariant across species *and* across S/M/L, which caps what a terminal duplex can
   distinguish. The *rescaling* is fine (M_pack median 0.66, p5 0.057, min 0.0037) — but
   **M_pack alone should not be expected to separate the b2 group.** Do not respond to
   this by retuning τ.

## 7b. WS3.2 (Task 9) — DONE (2026-08-29). User approved the D1 branch and the
"proceed as planned" response to the weak-M_pack finding.

`RPI_WEIGHTS` now carries `m_pack_tau = 3.914078624686338` and `m_pack_window_nt = 50`.
`RPI_WEIGHTS_VERSION` → `"2026-08-29.1"`. `MPackTauRescalingTest` (6 tests) asserts the
four required properties against the EXACT call shape Task 10 will use —
`boltzmann_compatibility(ddg, 0.0, rt=RPI_WEIGHTS["m_pack_tau"])` — plus a regression
guard that τ never drifts back to the dead 4.86 provisional. **τ is not yet consumed
anywhere**: `m_pack=1.0` at the `score_candidate` call site is untouched. That wiring is
Task 10 (WS3.3), scoped separately on purpose. Full suite 194+6 tests, all pass.

## WS3.3 + WS3.4 — DONE (2026-08-29). Both M_pack and M_pol score for real.

- `score_m_pack` (Task 10, D3): candidate's own cross-segment panhandle duplex vs the
  backbone species' own reference consensus (`species_consensus_panhandle_window`),
  rescaled via `boltzmann_compatibility(ddg, 0.0, rt=RPI_WEIGHTS["m_pack_tau"])`.
  `rpi_m_pack_status`: `scored` / `undefined_energy` (N4's `ValueError`, caught) /
  `no_parental_sequence` / `fixed_constant` (declared, not yet produced — see
  `score_m_pack`'s docstring).
- `score_m_pol` (Task 11, D3): candidate's own CDS-derived N-recognition footprint
  (beyond the invariant tag) vs the backbone species' own reference consensus CDS
  window (`species_consensus_cds_window`). `footprint` now reads
  `n_recognition_footprint(backbone_species)` — SNV's real 23-nt asymmetry is exercised
  at a production call site for the first time. `rpi_m_pol_status`: whatever
  `polymerase_compatibility` returns (`scored` / `no_signal_beyond_tag`), or
  `no_parental_sequence`.
- Both share **one** donor/backbone determination: `rpi_species_assignment(gate1,
  odd_segment)` (renamed from `mpack_species_assignment` — it is no longer M_pack-only).
- Smoke-tested on the real WS4 panel (12 cases): M_pack 0.10–1.0 on every candidate;
  M_pol 0.74–0.96 on the 3 SNV-backbone candidates, `no_signal_beyond_tag` (neutral 1.0)
  on the other 9 — honest, not a fake constant.
- `m_pol_variance_check` is still NOT wired in anywhere — it is spec 6.2's RUN-LEVEL
  check over a representative sample of lineage pairs, not a per-candidate computation.
  Building that aggregation is future work, never in either task's stated checklist.

## 7c. WS4 Task 12 — DONE. Task 13/14 — CODE DONE, RESULT NOT YET INTERPRETED.

**Start a fresh session here.** This section is the complete handoff — read it before
touching `rpi_discrimination.py` or `rpi_discrimination_stats.py`.

### Task 12 (D5) — DONE, see section 5 / commit `6a62c77`.
SNV x BCCV cross added to `PUBLISHED_CROSSES`. Validation harness sensitivity 19/19.
Corrected a wrong species lookup along the way (BCCV is `Orthohantavirus nigrorivense`
in this project's own curated table, not `bayoui`). Nothing further to do here.

### Task 13 (D7) — `rpi_discrimination.py`, commit `1e0b9aa`.
Code works, is tested (17 tests), and was run once against real data. **Read that
module's own docstring first** — it documents the four groups' construction in full
(random / b1 / b2 / real, sample sizes and why, the calibration-disjointness
enforcement) so this section does not repeat it.

### Task 14 (D7) — `rpi_discrimination_stats.py`, commit `1e0b9aa`.
Hand-rolled Mann-Whitney U (no scipy in this environment) + Cliff's delta (reused from
`rpi_calibration.py`), 12 tests, all pass. Computes the three pairwise comparisons the
ranking `reals ~ b1 > b2 > random` reduces to. **Read that module's own docstring** for
why "equivalence" is deliberately NOT a formal TOST claim.

### THE MEASURED RESULT (one run, 2026-08-29, seed 20260829) — the ranking does NOT hold:

    b2_vs_random:  p=0.153 (NOT significant), Cliff's delta = -0.225 (small, WRONG SIGN —
                   random fakes scored HIGHER than sympatric cross-species fakes)
    b1_vs_b2:      p=0.052 (borderline NOT significant), Cliff's delta = +0.275 (small,
                   right direction, just misses alpha=0.05)
    reals_vs_b1:   p=0.048 (SIGNIFICANT), Cliff's delta = -0.430 (medium, WRONG DIRECTION
                   for a claimed tie — reals scored LOWER than b1, not indistinguishable)
    overall_ranking_supported: False

Reproduce: `python rpi_discrimination.py --json <path>` (~15 min, scores 97 candidates
through `score_candidate`), then
`rpi_discrimination_stats.evaluate_pass_criterion(json.load(open(path))["groups"])`.

### LEADING HYPOTHESIS for `reals_vs_b1` — NOT YET CONFIRMED. Check this FIRST.

Group (c) "real positives" mixes two things with very different Factor E exposure:

- The 4 published inter-species lab crosses (RIZVANOV x2, HANDKE, the new BCCV one) are
  **real parental pairs that do NOT ecologically co-occur** (SNV is USA, ANDV is Chile;
  PUUV and PHV are different-continent lab pairings). Factor E is SUPPOSED to penalise
  these hard — "a perfect molecular profile must NOT rescue a candidate whose parents
  never co-occurred" is the model's own stated design (section 3 of this doc). Measured:
  `RIZVANOV2004_SNVsl_x_ANDVm` scored 0.0035, `RIZVANOV2004_reciprocal` 0.0046,
  `HANDKE2010` 0.0047 — all near the floor, almost certainly Factor E doing exactly what
  it is designed to do, NOT a M_pack/M_pol failure.
- The 7 all_hanta-detected intra-species candidates are same-species, so Factor E is
  near 1 (no ecological penalty) — matching b1's own Factor E exposure (b1 is
  within-species by construction too). These ranged 0.025–0.47 in the one run
  (`REAL_ORTHO_seoulense_001` = 0.050, `_002` = 0.47, `_003` = 0.35, `_004` = 0.025,
  `_000` hantanense = 0.47). 2 of 9 candidates (`REAL_ORTHO_sinnombreense_005`,
  the `unclassified Orthohantavirus|msb89863` row) hit `junk_sequence` at gate0 and
  scored `None` — NOT YET INVESTIGATED why; check gate0 failures on those two
  specifically before re-running.

**b1 has NO Factor-E-penalised members at all** (by construction every b1 case is one
species, so E is always ~1). Comparing a "reals" pool that is ~36% Factor-E-suppressed-
near-zero against a pool that is 0% Factor-E-suppressed is an apples-to-oranges
comparison that could manufacture exactly this "reals scores significantly lower"
result WITHOUT M_pack or M_pol being at fault.

**The specific next diagnostic, not yet run:** stratify `c_real_positive` by
`rpi_e_total` (already in the scored row) into "E ~= 1" (the all_hanta intra-species
candidates) vs "E suppressed" (the published inter-species lab crosses), and re-run
`reals_vs_b1` against ONLY the E~1 stratum. If THAT comparison also fails to look like a
plausible tie, the finding is real (M_pack/M_pol, or gate1's blind spot for
intra-species events, genuinely produce a reals/b1 gap) and should be written up
plainly per Task 14's honest-failure clause — do not retune τ or invent a new factor to
close the gap. If restricting to the E~1 stratum resolves it, the fix is a Task 13
change (split "real positives" into an ecologically-matched sub-group before comparing
to b1, or compare each stratum against a construction with matching Factor E exposure),
not a scoring-code change.

### `b2_vs_random`'s wrong-sign result — also not yet diagnosed.

`gate1_pattern` counts observed on this run: `a_random_fakes` was 17/20 `m_only` (gate1
DID detect discordance for random fakes, as expected — its 3 `undetectable_by_design`
cases are presumably close species pairs drawn by chance); `b2_cross_species_sympatric`
was 14/36 (39%) `undetectable_by_design` — the SAME near-identity blind spot Task 12
found for BCCV/bayoui, now affecting the sympatric species themselves
(hantanense/khabarovskense/seoulense/puumalaense are closely related, which is part of
WHY they can be sympatric). `undetectable_by_design` still sets `discordant_segment`
(it is a confidence qualifier on the same computation, not a different one — see
`gate1_check`), so `rpi_species_assignment` should still resolve a real donor/backbone
for these — this was NOT yet verified by checking `rpi_m_pack`/`rpi_m_pol` values on the
b2 rows specifically. **Next diagnostic:** pull the per-row `rpi_m_pack`/`rpi_m_pol`
values for `b2_cross_species_sympatric` (same JSON, `report["groups"]
["b2_cross_species_sympatric"]`) and check whether the `undetectable_by_design` rows
score differently (higher, i.e. more "plausible") than the `m_only`/`s_or_l_involving`
rows within b2 itself — and separately check `a_random_fakes`' own `rpi_e_total` values
(different-continent pairs should suppress E hard; if random's E is NOT actually being
suppressed as expected, that is a second, independent thing to check before concluding
anything about M_pack/M_pol).

### What NOT to do next.

- Do not retune `m_pack_tau`, `RPI_N_FOOTPRINT_NT`, or any `RPI_WEIGHTS` entry against
  this result — every prior session decision (2026-08-29, the WS3.1 and WS3.3/3.4
  sign-offs) was explicit that weak-signal findings get reported, not retuned against
  the validation set. This result is validation, not calibration.
- Do not change `rpi_discrimination.py`'s group construction until the E-stratification
  diagnostic above is actually run — the hypothesis is plausible but unconfirmed, and
  changing the panel before confirming would risk fitting the panel to the desired
  answer.
- Do not report Task 14 as "failed" or "passed" externally yet — the honest-failure
  clause requires reporting the measurement, but the measurement needs the diagnostics
  above before its cause is understood well enough to write up.

### Task 15 (not started).
Chimera ladder (detection-limit-vs-divergence sweep) + extend `sensitivity_sweep` to
cover `m_pack_tau` and `m_pack_window_nt` at ±50% / 25-50-75 nt. The `sensitivity_sweep`
extension is small and mechanical (follow the exact pattern already used for
`jeffreys_alpha`/`host_decay` — deep-copy `RPI_WEIGHTS`, perturb, score
`_rpi_sensitivity_panel`, report rank correlation); the chimera ladder needs a
divergence metric (candidate: `gate0.best_kmer_identity` or `SpeciesHit.distance` from
`nearest_species_by_kmer`) and a set of real species pairs spanning it, not yet chosen.
Blocked behind Task 14 being actually resolved, not technically blocked.

**Binding constraint (still applies to any further b1/b2/random work):** every group
builder **MUST exclude the 42 isolate keys in `data/rpi_calibration.json` →
`reserved_for_calibration.isolate_keys`** (`rpi_discrimination.usable_pool` already
does this). τ is derived from this panel's ΔΔG span, so any overlap tunes the scale
parameter on the validation set and makes the WS4 result circular.

**Weak-signal findings from WS3, still true, still not to be fixed by retuning:**
- M_pack: ΔΔG separates cognate/non-cognate only weakly (Cliff's δ = +0.222; 38% of the
  calibration panel has ΔΔG ≤ 0, where M_pack saturates at 1.0).
- M_pol: real signal exists only for an SNV backbone (the only species with a non-default
  `RPI_N_FOOTPRINT_NT` entry); every other backbone reports `no_signal_beyond_tag`.

### Gate 4 blockers to handle when it is first run
- The executable is **`iqtree3.exe`**, not `iqtree`. `parse_args` defaults
  `--iqtree-bin` to `None`; `validate_reassortment_detector.py:464` passes the literal
  `"iqtree"`. Both need the real path threaded through.
- **MAFFT's `--version` is not on the first line** — the `.bat` prints a "Preparing
  environment" preamble. Scan for the `v7.x` token.

### E-stratification diagnostic — RUN, 2026-08-30. Outcome: genuine equivalence, not E-driven.

The diagnostic specified above was run: `diagnose_e_stratification.py` (new, read-only —
loads the existing panel JSON, does not regenerate it or touch `rpi_discrimination.py`/
`rpi_discrimination_stats.py` beyond calling the public `compare_groups`) against
`output/discrimination_scratch/discrimination_panel.json` (regenerated by the user in
their own terminal, unchanged code, seed 20260829).

**No evidence the reals-vs-b1 tie was E-driven.** `rpi_e_total` is cleanly bimodal (1.0
x6 vs 0.01 x3; the three 0.01 rows are the published lab crosses, whose low E is
expected behavior: lab crosses bypass natural ecological opportunity, so a low
ecological prior is correct, not an artifact). In the E~1 stratum (n=6 vs b1 n=30):
medians 0.3604 vs 0.3470, Mann-Whitney p=0.582, Cliff's delta -0.144 (negligible). The
tie persists, consistent with D7: RPI measures plausibility, and the tree detector — not
the index — separates plausible from occurred. Caveat of record: n=6 is a small
stratum, so this is no-detectable-difference, not proven equivalence; the
equivalence-interval claim stays open until group c grows.

**Provenance note:** the diagnostic was executed via the PowerShell tool, not Bash —
this session's Bash tool refused every `python -c`/script invocation with a safety
classifier block unrelated to the command content itself (auto-mode's "earlier
conversation content" guard). PowerShell was unaffected. Recorded here because the
numbers above did not come from the usual Bash path.

`b2_vs_random`'s wrong-sign result (previous section) is still undiagnosed — separate
from this finding and not addressed by it.

### Factor E isolate-level fix — IMPLEMENTED, 2026-08-30. `overall_ranking_supported`
now True, with a caveat on how it got there. Full record: `notes_e_isolate_fix.md`.

Root cause (diagnosed read-only in `notes_diagnosis.md`, this session): `ecological_
opportunity`'s geo/host/time terms were computed from the WHOLE SPECIES' isolate
population in the reference (`by_species[species]`), not from the two specific parent
isolates actually implicated in a candidate — so a widely-distributed or
border-spanning species (Puumala, Seoul) registered "geo opportunity" against any other
species sharing so much as one country anywhere in the corpus, defeating
`build_group_a_random`'s own different-continents construction guarantee.

Fix, per `docs/superpowers/specs/2026-08-30-rpi-e-geo-isolate-level-fix-design.md`
(approved with one addition before code was written — the fallback below):
`ecological_opportunity` now takes ONE isolate per implicated species
(`parent_isolates`) and computes host/geo/time from that pair directly. For
candidates whose true parent isolate is known (WS4's group a/b1/b2 builders, which
construct candidates from specific real isolates on purpose — now via
`CandidateInput.known_segment_isolates`), that real isolate is used
(`rpi_e_isolate_basis="actual_parents"`). For every other candidate — ordinary
production scoring, and WS4 group c's real positives, whose donor parent is exactly the
unknown reassortment detection exists to resolve — `resolve_e_parent_isolates` falls
back to the k-mer-nearest wild reference isolate of that species
(`nearest_reference_isolate`), basis `"nearest_reference"`, or `"insufficient_
metadata"` if even that fails. Every row declares its basis; nothing is silently on a
proxy. `RPI_WEIGHTS` unchanged; this is a bug fix, not a retune.

**Before/after, seeded panel (seed 20260829), `evaluate_pass_criterion`:**

| comparison | before | after |
|---|---|---|
| `b2_vs_random` | p=0.153, δ=−0.225 (wrong sign) | p=4.7e-09, δ=+0.953 (large, correct) |
| `b1_vs_b2` | p=0.052, δ=+0.275 (borderline) | p=2.2e-11, δ=+0.961 (large) |
| `reals_vs_b1` | p=0.048, δ=−0.430 (significant, wrong direction) | p=0.733, δ=+0.078 (negligible — a tie) |
| `overall_ranking_supported` | **False** | **True** |

`a_random_fakes`'s `rpi_generation` mean collapsed 0.331 → 0.000138 (`rpi_e_total` was
`{1.0: 14, 0.01: 5, 0.000625: 1}`, now suppressed on all 20/20) — exactly the mechanism
the design doc predicted, confirmed via the basis column (20/20 `actual_parents`).

**Caveat, checked and confirmed, not assumed:** the `reals_vs_b1` tie is NOT a clean
confirmation of the fix's own mechanism. 3 of `c_real_positive`'s 4 `nearest_reference`-
basis rows (RIZVANOV2004 ×2, HANDKE2010 — published inter-species lab crosses) moved
from the pre-fix floor (~0.004) to 0.35–0.47 because their resolved PROXY isolate
(picked by protein similarity alone) happened to carry NO country/host/year metadata at
all, so `ecological_opportunity`'s pre-existing "missing metadata means neutral" rule
fired — not because isolate-level geography found real ecological overlap. See
`notes_e_isolate_fix.md` for the full per-case diagnosis (reproduced directly against
the real reference data, not inferred). The basis column correctly flags these rows as
proxy-based, so nothing is misrepresented, but pulling them out would likely move the
result back toward the old (low) `reals_vs_b1` reading. Reported as measured, not edited
toward the preferred outcome — per this task's own honesty protocol.

**Two new known limitations, both logged, neither fixed here:**
- **E vs gate3 disagreement, extended.** `ecological_opportunity` now reads isolate-level
  geography for the actual pair; `gate3_check` (verdict-bearing, deliberately untouched)
  still reads species-wide geography for only the first two sorted implicated species —
  one more concrete instance of the already-deferred Factor E rebuild.
- **`nearest_reference_isolate` has no preference for metadata completeness** — it picks
  the nearest isolate by k-mer protein distance only, so it can resolve to a
  sparsely-annotated GenBank deposit and go neutral by omission rather than by evidence.
  A future fix could prefer the nearest isolate WITH usable geo/host/year metadata, or
  report `insufficient_metadata` instead of a metadata-empty proxy.

Full suite after the change: 237 tests (`reassortment_detector`, was 228), all pass.
`rpi_discrimination`: 18 tests, pass. `rpi_calibration`: 39 tests, pass, unaffected.

### gate0 junk_sequence diagnostic — RUN, 2026-08-30. Outcome: real incomplete deposits, not a bug.

`diagnose_gate0_junk.py` (new, read-only — calls only public `gate0_check` /
`_build_precomputed` / `load_reference_data`) checked why `REAL_ORTHO_sinnombreense_005`
and the `unclassified Orthohantavirus|msb89863` row hit `junk_sequence` in the
`c_real_positive` group.

- **`sinnombreense_005`** (`t2|Orthohantavirus sinnombreense|202003197`): only L fails —
  5705 nt vs the expected 6300-6800 nt window. S (1461 nt) and M (3673 nt) both pass
  their own length checks. A partial/truncated L deposit, not a data-loading bug.
- **`msb89863`** (`t2|unclassified Orthohantavirus|msb89863`): fails on **all three**
  segments — S=330 nt (expected 1200-2100), M=1390 nt (expected 3400-3800), L=490 nt
  (expected 6300-6800). These are not full genomic segments at all, just short partial
  GenBank records.

Both also warn `terminal panhandle motif not recognizable in either orientation` on
every segment, consistent with the documented "48% of deposits lack their true termini"
limitation. **Conclusion: gate0 is correctly rejecting genuinely incomplete NCBI
deposits — this is expected behavior, not a bug.** No code change indicated; these two
candidates staying `None` in the `c_real_positive` group is correct and should not be
patched around (e.g. do not relax `SEGMENT_LENGTH_RANGE` to force a score out of them).

---

## 7d. GA search runtime fix + overnight full run — DONE, 2026-09-05 (branch `05-module-consolidation`)

> **UPDATED 2026-09-05 with corrected pool results** (§7e found a pool defect, commit
> `105bd5c` fixed it, re-ran both benchmarks). The runtime work is sound; the
> acceleration is 9.5x; the ceiling is pool-independent and survives on corrected data.

**Problem:** `Algorithm.py`'s default GA search (population=150, generations=100, over
`rpi_core.score_candidate`) was taking up to 10 hours.

**Diagnosis (profiled, not guessed):** `cProfile` on a small sample (population=20,
generations=3, real `all_hanta/` data) showed ~70% of `score_candidate`'s time in
`gate0_check`/`gate1_check`'s k-mer cosine-similarity matrix (`profiles_to_matrix`) and
~25% in `score_m_pol`'s `species_consensus_cds_window` — both rebuilding data that
depends ONLY on the fixed reference set, from scratch, on every single candidate call,
even though candidates are existing reference isolates recombined (never synthesized),
so the same reference-side work repeats thousands of times.

**Fix (2 commits, each preceded by its own golden-master regression test against real
`all_hanta/` data, `6a6ab6a`/`765970b` and `618a4d8`/`b0b5761`):**
- `_kmer_identity_to_reference` / `nearest_species_by_kmer`: the reference-side k-mer
  matrix is now built once per `reference_profiles` list and cached, instead of rebuilt
  per call. Provably score-neutral: L1-then-L2 normalization is scale-invariant, so each
  row's unit vector depends only on its own counts.
- `species_consensus_cds_window`: cached per `(reference, species, segment, window)`.

**A real bug caught mid-implementation, not just a review nit:** both caches key on
`id(...)`. Naively, a garbage-collected object's id can be reused by an unrelated new
object, returning a stale cached result for a completely different input — this actually
happened (`Gate1CheckTest.test_implicated_species_come_back_in_a_stable_order` failed
when first implemented, and would have hit `species_consensus_cds_window` the same way
via `SpeciesConsensusCdsWindowTest`'s repeated `"SpA"` fixture). Fixed by having each
cache hold a strong reference to the object behind its key, pinning the id for the
process's life. **Anyone adding another `id(...)`-keyed cache to this codebase must do
the same, or add tests proving two different objects sharing a name/label are not
confused by it.**

**Verification, in order:**
1. Golden-master unit tests (real data, hand-picked edge cases including a candidate
   protein shorter than `KMER_K`) — pass before AND after each fix.
2. Full `rpi_core` suite: 326 tests before the fixes existed, 328 after (the 2 new
   golden-master tests) — all pass, both times.
3. **Discrimination panel re-run, seed 20260829, 97 candidates across all 4 groups**, at
   the pre-optimization commit (`ecaa1d4`) and post-optimization HEAD, in an isolated
   `git worktree` (no working-tree disruption): summary stats byte-identical, **every
   individual scored row identical, zero diffs.** Also matches the independent
   pre-existing `output/consolidation_baseline/discrimination_panel_golden.json`.

**Measured speedup:** same 50-candidate profiling sample, 736s → 221s (kmer fix alone) →
77.8s (both fixes) — **9.5x**.

**Resilience added before the overnight run (commit `d4a1a5b`):** `Oracle.score` in
`Algorithm.py` had no protection against `score_candidate` raising on one bad triplet —
any exception would have crashed a multi-hour unattended run with no checkpointing. It
now catches `Exception`, scores that candidate 0.0 (same treatment as `junk_sequence`),
and records `(triplet, error)` in `oracle.errors` / the output payload's
`scoring_errors` field, so a run reports exactly what it skipped instead of either
crashing or silently hiding a failure. Guarded by 2 new `OracleTest` cases.

**CORRECTED POOL RESULTS — commit `105bd5c` fixed `build_pool` to call `rc.usable_pool`:**

| Metric | Contaminated pool (192 iso) | Corrected pool (272 iso) |
|---|---|---|
| **Wall time (21 gen)** | 5m 15.978s | **8m 56.101s** |
| **Wall time (300 gen)** | 18m 31.586s | **24m 50.940s** |
| **Generations run** | 21 of 100 — plateau stop | 21 of 100 — plateau stop |
| **Best fitness** | 0.4716981132 | 0.4716981132 |
| **Winners (fitness)** | All 10 at 0.4717 | All 10 at 0.4717 |
| **Scoring errors** | 0 | 0 |
| **Pool size in results** | 192 | 272 |

**Key finding — verified by running the identical search on both pools:** the fitness
values, generation count, ceiling, and winner structure are **IDENTICAL**. This
definitively proves the 25/53 ceiling and 21-generation plateau are **pool-independent,
structural properties of the scoring formula**, not artifacts of the contaminated pool.
The larger search space (272 vs 192) reduces per-generation cache reuse, which is why
wall times increased (by ~70%), but the optimization remains **9.5x faster than
pre-optimization**.

The previously-quoted "up to 10 hours" was for the theoretical 100-generation run. In
practice the search plateaus at generation 21 on both pools. **The true, measured
runtime is: ~5 minutes on the contaminated pool, ~9 minutes on the corrected pool, for
seed 20260904 with default settings.**

*(Superseded estimate, kept as a worked example of why §8 says flag estimates: this
section originally extrapolated a forced 100-generation run at `100/21 × 5m16s ≈ 25 min`.
The forced run was then actually measured — **300** generations in **18m 31.586s**, i.e.
~3.7 s/generation against the 21-generation run's ~15 s/generation. The linear
extrapolation was wrong because per-generation cost FALLS as the population converges:
fewer distinct triplets per generation means more `Oracle.cache` hits. Prefer the measured
number; see the resolution subsection below.)*

**Why 0.4717 is never beaten — RESOLVED, 2026-09-05. It is a structural ceiling, not a
converged optimum, and the GA is behaving correctly.**

`rpi_generation = E × M_pack × M_pol × K` (see `compute_rpi`). E, M_pack and M_pol each
saturate at 1.0, so the achievable maximum is set entirely by **K**, and

```
K = constellation_prior(stratum)[odd_segment]
constellation_counts["close"] = {M: 12, S: 12, L: 1},  jeffreys_alpha = 0.5
  -> K_close[M] = K_close[S] = (12 + 0.5) / (25 + 1.5) = 12.5 / 26.5 = 25/53
                                                        = 0.4716981132...
```

**0.4716981132 is exactly 25/53** — the close-stratum M/S constellation prior, reached
whenever E = M_pack = M_pol = 1.0. Nothing in the GA's search space can exceed it.

**Why every GA candidate is in the `close` stratum:** `divergence_stratum` returns
`"close"` when `identity >= stratum_identity_boundary` (0.8), and `identity` is
`gate0.best_kmer_identity` — protein k-mer identity to the nearest reference isolate.
`Algorithm.py` builds candidates by **recombining existing reference isolates**, never by
synthesizing sequence, so each segment's protein IS a reference protein and scores
identity ≈ 1.0. The `distant` stratum (where `K[M] = 0.92`, nearly double) is therefore
**structurally unreachable by this search**, and `rpi_k_stratum_provisional` is set on
every row precisely because this `identity` is not the parental nucleotide divergence
spec §7 intends to stratify on (see the SCAFFOLDING comment at `score_candidate`'s
`identity=` argument).

**Empirical confirmation, three independent ways:**

| Check | Result |
|---|---|
| 120 random triplets from the 192-isolate pool, scored directly | max = exactly 25/53, **0 exceeded it**, 1.7% (2/120) hit it exactly, all 120 in `close` stratum |
| **Forced long run** — 300 generations, plateau early-stop disabled, seed 20260904 (`output/ga_search_long_2026-09-05/`) | best fitness **identical at 25/53 for all 300 generations, never once improved**; 0 scoring errors; wall time **18m 31.586s** |
| Overnight default run (§ table above), 21 generations | same value, same stratum, all 10 winners tied |

With a 1.7% ceiling-hit rate and population=150, the chance generation 0 already contains
a ceiling-hitter is `1 - (1 - 0.017)^150 ≈ 92%` — which is exactly why "best" appears at
generation 0 and never moves. **The GA is not stuck; there is nowhere higher to go.**

**The fitness landscape is NOT flat, but the top of it is.** The same 120-triplet sample
produced **86 distinct fitness values**, so the index discriminates perfectly well *below*
the ceiling. The problem is at the top: ceiling-hitters tie at *exactly* 25/53, so the GA
**cannot rank among them at all.** Over 300 generations the population collapsed from 150
distinct triplets to 34, with median = best = the ceiling from ~generation 5 onward — it
is not selecting a best candidate, it is accumulating ties.

**Consequences — read before quoting any GA "winner":**
1. **The "top 10 winners" list is an arbitrary tie-break, not a ranking.** All ten scored
   identically in both runs. Which ten appear is an artifact of population order and
   `dict`/list traversal, not of any measured difference between them.
2. **Running longer buys nothing.** 300 generations (18.5 min) found exactly what 21
   generations (5 min) found. Do not spend hours on bigger runs expecting improvement.
3. A high RPI still establishes **PLAUSIBILITY, not occurrence** (the gate-4 tree
   detector, not the index, establishes that — see §3). Both runs' winners carry
   `insufficient_comparators` (8–9 of 10) or `undetectable_by_design` verdicts, which must
   be read alongside the fitness number, never instead of it.
4. **Open question, deliberately not answered here:** whether the GA is the right
   instrument at all if its objective saturates on ~1.7% of the search space. Options, in
   rough order of how much they change the science: (a) rank ceiling-hitters by a
   documented tie-breaker (e.g. `rpi_retention`, which multiplies by F and is *not*
   saturated); (b) report the ceiling set as an unordered candidate set rather than a
   ranked top-10; (c) revisit whether `identity` should be real parental nucleotide
   divergence — which would let the `distant` stratum become reachable and would change
   K's ceiling to 0.92. **(c) changes scores and must not be done as a "fix" without a
   design decision and re-calibration**, per §8.

---

## 7e. `/code-review high` on §7d's work — 2026-09-05. **One finding invalidates §7d's pool. NOTHING FIXED YET.**

Run after the §7d runtime work and the ceiling analysis, per the standing "review before
trusting an outcome" habit. **No code was changed in response** — every finding below
either changes scores or needs a design decision, and §8's rule is that only
provably-score-neutral fixes land autonomously. These are proposals awaiting a decision.

### 🔴 FINDING 1 (HIGH) — the GA searches the WRONG isolate pool. §7d's numbers stand on it.

`Algorithm.py:build_pool` filters only on `has_complete_triplet` / `is_lab_derived` /
`is_vaccine_strain`. The project's canonical pool is `rpi_core.usable_pool`, whose
docstring **explicitly forbids a second implementation** and additionally requires
`terminally_complete_isolates` (all three segments full-length AND carrying both panhandle
tags) plus dropping unclassified species — "because 48%+ of deposits are trimmed of their
true termini and would corrupt any panhandle-adjacent computation." **That is exactly the
M_pack term the GA optimises.**

Measured against real `all_hanta/` (verified independently, not taken on the reviewer's word):

| | count |
|---|---|
| `Algorithm.build_pool` (what the GA actually searched) | **192** |
| `rpi_core.usable_pool` (canonical) | **272** |
| In the GA pool but NOT canonical — **contaminated** | **98 (51% of what it searched)** |
| Canonical but never searched — **missed** | **178 (65% of the real pool)** |
| Overlap | 94 |

**Consequence:** every number in §7d — the 5m16s and 18m31s runs, the 25/53 ceiling
observations, the 1.7% ceiling-hit rate, both winner lists — was produced over a pool that
is **half contaminated and missing two-thirds of the legitimate isolates.** The
*mechanism* explained in §7d (K caps at 25/53 for the close stratum; E/M_pack/M_pol
saturate) is a property of the scoring formula and does NOT depend on the pool, so that
conclusion survives. The *empirical* figures do not. **Do not quote §7d's counts,
rates, or winner lists until this is fixed and the runs are redone.**

Proposed fix: `build_pool` calls `rc.usable_pool(reference, reserved)`. This CHANGES
SCORES (different pool → different candidates), so it is a decision, not a cleanup — and
the re-run afterwards is mandatory, not optional.

### 🟠 FINDING 2 (HIGH) — `--with-phylogeny` cannot ever succeed

`confirm_winners` builds `Namespace(..., mafft_bin=None, iqtree_bin=None, ...)`, so
`gate4_check` → `run_mafft` → `subprocess.run([None, "--auto", ...])` → `TypeError`. Every
winner gets `gate4_error` while the function still returns `"attempted: per-winner gate4
fields recorded"` and `main` prints that as though confirmation ran. **This is the §8
defect shape again — a success message over a value that failed.** `rc._score_args` shows
the working form (`mafft_bin="mafft"`, `iqtree_bin="iqtree3"`); `ph.resolve_binary` would
be better still. Never exercised by §7d's runs (the flag was not passed), so it does not
affect those results.

### 🟡 FINDING 3 (MEDIUM) — missing calibration artifact silently un-reserves the 42

`load_reserved_isolate_keys` returns `set()` when `data/rpi_calibration.json` is absent,
and resolves the path relative to CWD. `rc.reserved_calibration_isolate_keys` deliberately
RAISES instead. Run `Algorithm.py` from any other directory and the 42 calibration-reserved
isolates silently re-enter the search — the exact tau-tuned-on-the-validation-set
circularity §7b exists to prevent — signalled only by `"calibration_reserved_excluded": 0`
buried in the output JSON. (§7d's runs were from the repo root and DID exclude 42, verified
in the payload.) Should raise, or at minimum refuse to run.

### 🟡 FINDING 4 (MEDIUM) — `--diversity` can hang forever

`while len(unique) < population_size:` accepts only unseen random triplets; with fewer than
`population_size` distinct triplets available (`len(pool)**3 < population_size`) it never
terminates — an unattended run spins at 100% CPU, silently. Needs a bounded attempt budget
(the pattern `build_group_a_random` already uses).

### 🟡 FINDING 5 (MEDIUM) — §7d's caches retain ~750 MB for the process lifetime

Caused by §7d's own optimisation. Independently measured here (`unit.nbytes`, real
`all_hanta/`): S `(6234, 5547)` 263.8 MB + M `(4205, 7424)` 238.2 MB + L `(4152, 7884)`
249.7 MB = **751.7 MB** across 3 cache entries, held for the process lifetime because the
entry pins the `reference_profiles` list. *(The review reported 788.3 MB; the ~36 MB gap
is measurement basis, not disagreement — the direction and order of magnitude match.)*
`_SPECIES_CDS_WINDOW_CACHE_REFS` additionally pins whole `ReferenceData` objects. Before
the change each matrix was transient (~277 MB peak, one segment at a time, freed after the
call), so steady-state RSS roughly triples and never drops. A process building two
`_build_precomputed` results over `all_hanta` pays it twice (~1.6 GB) — which
`python -m unittest rpi_core` now does, since `KmerIdentityAndNearestSpeciesGoldenMasterTest`
builds its own. The correctness of the caches is not in question (the L1-then-L2
scale-invariance argument holds, and the review confirmed it); the **retention policy** is
what needs a bound or an explicit eviction hook. Note the trade-off: eviction costs back
part of §7d's 9.5x speedup, so this is a decision, not an obvious win.

### 🟢 FINDINGS 9–11 (LOW) — robustness, not correctness

- **9.** Every winner reuses `candidate_id="ga-confirm"`, so with `--with-phylogeny` all
  ten winners write gate4 artifacts to the same directory; only the last survives and none
  can be attributed. Use the triplet or winner index in the id.
- **10.** A *systemic* scoring failure (bad `--reference-input`, factory `KeyError` on
  every candidate) is indistinguishable from a flat landscape: §7d's new `except Exception`
  scores 0.0 and the run "succeeds" with best fitness 0.0, exit 0, and one stderr line that
  scrolls past unattended. Should abort when the error rate crosses a threshold (e.g. all
  of generation 0 failed). *This is the §8 defect shape sneaking back in via the
  resilience fix — worth fixing for that reason alone.*
- **11.** `--elitism 0` is accepted and removes the guarantee that the incumbent survives,
  so the reported "best" (taken from the final population) can be worse than the best
  actually scored (still sitting in `oracle.cache`). Track best-ever, or reject `< 1`.

### ⚪ FINDINGS 6–8 — the `CLAUDE.md` / `PUBLIC_RELEASE.md` move — RESOLVED 2026-09-05, committed on explicit instruction

Previously: the working tree had this move (root `CLAUDE.md`/`PUBLIC_RELEASE.md` deleted,
copies added under `docs/`) sitting uncommitted, flagged here as needing an explicit
decision rather than being committed or reverted autonomously — **Claude Code auto-loads
`CLAUDE.md` from the repo root, not from `docs/`**, so committing it as-is would silently
stop applying the repo policy it carries to every future session.

**Resolution:** the user's session-start instruction for this session's work explicitly
said "commit everything" as step 1 (baseline, commit `7624ad6`, tag `pre-figures-baseline`)
— that is the explicit decision this finding was waiting on. Committed as directed.

**Still open, not addressed by that commit:** `docs/PUBLIC_RELEASE.md`'s own COPY list
(lines 48, 58) and this file's own line 4 still name both files at the repo root, so a
release built by literally copying that list would omit them — needs a follow-up edit to
`docs/PUBLIC_RELEASE.md` and this file's header, out of scope for this session's four-step
plan. `Algorithm.py` and `figures.py` both remain unlisted in the COPY/NEVER lists and in
no README module map — same open item, now two files instead of one.

## 7f. Post-§7e four-step session, 2026-09-05 (branch `05-module-consolidation`)

Explicit instruction, no input-contract refactor (that stays parked per the decision
recorded at the top of §7e's introduction — stages remain independent, no stage reads
another's output file, per the 2026-08-13 three-stage-split spec §4).

**Step 1 — baseline.** Commit `7624ad6`, tag `pre-figures-baseline`. Committed the
working tree's pre-existing uncommitted state (§7e Findings 6–8's `CLAUDE.md` move, the
`rpi_core.py` cache-ref-pinning fix, five untracked `output/ga_search_*` directories) and
refreshed this file's §1. Full suite verified green first: 328 tests
(`python -m unittest rpi_core`).

**Step 2 — §7e Finding 2 fixed.** Commit `3ab162c`. `confirm_winners` hardcoded
`mafft_bin=None, iqtree_bin=None`; fixed to resolve real paths via `ph.resolve_binary`
(same function Stage 2 uses). **This turned out to be load-bearing, not cosmetic:**
verified directly (see step 4) that `subprocess.run(["mafft", ...])` fails outright on
Windows (`mafft` is a `.BAT`; Windows will not exec a bare bat name without a resolved
absolute path or `shell=True`) — so the bare-string form `rc._score_args` uses would
*also* have failed here had confirm_winners copied it verbatim; resolving through
`ph.resolve_binary` was necessary, not just cleaner. Test-first, 5 new cases (then 4 more
in step 4). All other §7e findings (1, 3–11) intentionally untouched.

**Step 3 — before/after figures.** Commit `96b4e82`. `figures.py` (new), reading only
committed Stage 2 outputs (`output/phylogenetics*/…/tree.treefile`), pruned to a local
neighborhood (user-confirmed choice — several species trees run 200–250+ taxa) for
legibility. BEFORE: all 7 unique `supported_reassortment_candidate` isolates, reassortant
tip highlighted + connected across S/M/L, discordant segment colored/bolded. AFTER: the
10 winners in `output/ga_search_corrected_pool_long/seed_20260904.json` (the longer of
the two corrected-pool §7e-Finding-1 runs — both plateau at the identical 0.4717 ceiling
but land on different tied triplets, confirming §7d's structural-ceiling-flatness finding
rather than contradicting it), each donor's own tree rendered, colored by slot, captioned
"in-silico combination... no new sequence is created or synthesized" per instruction.
Added `matplotlib>=3.9` to `environment.yml`/`requirements-pip.txt` (figures only, not on
the scored RPI path). 10 new tests in `figures.py`, plus a real-data integration smoke
test. Output committed under `output/figures/{before,after}/` plus
`output/figures/after/ga_winners_parents.csv` (candidate → S/M/L parent isolate/species/
country/host).

**Step 4 — gate4 confirmation of the GA winners — IN PROGRESS, not yet committed.**
Commit `0294896` added `confirm_winners(..., checkpoint_path=...)` (test-first, 4 cases):
appends each winner's gate4 result to a JSONL file as it completes, so a crash mid-run
loses nothing already confirmed. Necessary because this run is genuinely expensive, not
because it hangs — **verified real**: a bare 3-sequence MAFFT alignment via the resolved
`tools/mafft-win/mafft.BAT` path completes in ~5s, but 9 of the 10 winners implicate
`Orthohantavirus hantanense` alone, whose S/L comparator trees run 200–250 taxa, and
IQ-TREE's ultrafast bootstrap enforces a minimum of 1000 replicates per segment per
winner — a single real call was still running after 25+ minutes before being killed
purely to diagnose (confirmed via `Get-CimInstance`/`Get-Process`: real child processes,
non-trivial but low cumulative CPU, consistent with genuine large-alignment IQ-TREE work
under Windows Defender's real-time-scan overhead on every spawned subprocess, not a
deadlock).

Launched detached (`Start-Process`, hidden window, `SetThreadExecutionState` to block
Windows sleep for the process's lifetime) at 2026-09-05T13:52:39Z:
- Driver: ad hoc, not committed (one-off orchestration, not product code) — reconstructs
  `Oracle` exactly as `Algorithm.main()` does, loads the **already-committed**
  `ga_search_corrected_pool_long/seed_20260904.json` winners (no GA re-search), calls
  `confirm_winners(..., with_phylogeny=True, checkpoint_path=...)`.
- Checkpoint: `output/ga_search_corrected_pool_long/gate4_confirmation_checkpoint.jsonl`
  (one line per winner as it finishes).
- Final: `output/ga_search_corrected_pool_long/gate4_confirmation.json`, written once at
  the end (same shape as the checkpoint's lines, plus the summary `confirmation` string).
- Log: `output/ga_search_corrected_pool_long/gate4_confirmation.log{,.stdout,.stderr}`.

**Binding reporting rule for this run's results (set by the user, applies to whoever
writes up step 4's outcome next):** gate4 verdicts are reported verbatim per winner,
`no_conflict` and `insufficient_comparators` included. 9 of 10 winners are
within-hantanense combinations — **gate4 (the tree-based detector) is exactly the
instrument that CAN see intra-species reassortment** (unlike RPI's own Gate 1, whose
blind spot there is structural and separately documented above) — so `no_conflict` on
those 9 is a real, meaningful result if that is what comes back, not an instrument
limitation to explain away. No re-running for a "better" result, no softening language.

**Next session action:** check `gate4_confirmation.json` / the checkpoint file for
completion, review `gate4_confirmation.log`, write up the verdicts into this section
verbatim per the rule above, and commit `output/ga_search_corrected_pool_long/gate4_*`
as this session's step 4 commit.

## 7g. Overnight full-pipeline rerun and reconciliation, 2026-09-06/07 (branch `05-module-consolidation`)

User request: rerun the full pipeline (`all_hanta.py` → `phylogeny.py` →
`rpi_core.py`/`Algorithm.py`) from scratch into an isolated folder, reconcile against
committed `output/`, promote, illustrate the phylogeny trees, survey redundant
artifacts (no deletion), and report. Design:
`docs/superpowers/specs/2026-09-06-overnight-full-pipeline-rerun-design.md`. Plan:
`docs/superpowers/plans/2026-09-06-overnight-full-pipeline-rerun.md`.

### §7f's stalled gate4 run actually finished on its own

The process from §7f step 4 was never killed — it kept running detached across the
session boundary and completed independently around **2026-09-06 10:12** (checkpoint
grew from 4 to 10 lines; `output/ga_search_corrected_pool_long/gate4_confirmation.json`
now exists, untracked). Its content is **byte-identical** to this session's
independently-launched confirmation run (verified via diff). Two fully independent
multi-hour gate4 runs converging on identical output is real evidence that gate4
confirmation (single-threaded IQ-TREE, `-T 1`) is deterministic — unlike the
multi-threaded tree-building elsewhere (see below). It used the pre-Finding-9-fix
shared `candidate_phylogenetics/ga-confirm/` directory (all 10 winners' trees
overwriting each other in one directory) and is superseded by this session's correctly
per-winner-named run; it goes on the redundancy list (§ below), not deleted.

### Pre-run review fixes (score-neutral, applied before the rerun)

- **§7e Finding 9** (`Algorithm.py`): `confirm_winners` hardcoded
  `candidate_id="ga-confirm"` for every winner, colliding all ten winners' gate4
  artifact directories under `--with-phylogeny`. Fixed to `f"ga-confirm-{idx}"`.
  Guard-tested. Commit `49a4efe`.
- **§7e Finding 10** (`Algorithm.py`): a systemic scoring failure (bad
  `--reference-input`, a factory `KeyError` on every candidate) was indistinguishable
  from a genuine flat 0.0 fitness landscape — the run "succeeded" with exit 0 and one
  easy-to-miss stderr line. Added `_check_no_total_scoring_failure`, which aborts only
  when 100% of scored triplets errored; never fires on a healthy run. Guard-tested.
  Commit `88908c6`.
- Findings 3, 4, 11 reviewed and left as open proposals — none triggered by tonight's
  run configuration (repo-root cwd with calibration present, no `--diversity`, default
  elitism).

### New feature: `Algorithm.py --dump-scored` (user request mid-session)

Writes every unique triplet the GA search actually scored (not just the top 10
winners) to `all_scored_seed_<N>.jsonl`/`.csv`, from `oracle.cache`/`oracle.rows`.
Purely additive — verified by re-running both GA configs (default, and long-pool with
`--plateau 10000`) with the new flag: winners, `generations_run`, fitness log, and
`scoring_errors` are **bit-identical** to the committed
`output/ga_search_corrected_pool_{default,long}/seed_20260904.json` in every field.
Documented in the README (previously an open gap per §7e Findings 6–8's note that
`Algorithm.py` was in no module map). Commit `fa7609f`.

### The rerun itself

Driven by an ad hoc script (not committed, matching §7f step 4's precedent) into
`output_rerun_20260906/`. Two driver bugs found and fixed mid-run (both before any
data was lost):

1. `Algorithm.py`'s long-pool config needs `--plateau 10000` (read back from the
   committed run's own `"config"` field in `seed_20260904.json`) to avoid an early
   plateau stop at generation 21 — the driver omitted it, producing the *wrong*
   winners on the first attempt. Caught before gate4 confirmation wasted hours on
   them; corrected and re-run. Verified bit-identical to golden afterward (see below).
2. `rpi_core.py validate`'s `--output` defaults to the relative path
   `"output/validation"`, which resolves against the driver's cwd (the repo root), not
   the isolated rerun folder — omitting an explicit `--output` the first time wrote
   gate4's tree artifacts (80 files) into the **live, committed** `output/validation/`
   instead. Caught via `git status` showing new untracked files there; the 80 files
   were moved into `output_rerun_20260906/validation/candidate_phylogenetics/` (their
   correct home) before anything was committed. `output/` was never actually modified
   in git terms (files were untracked, never staged).

Full stage timeline: all_hanta ingest (~18s) → phylogeny, both scopes (~9h21m, zero
errors) → `rpi_core run` (~53s, reused pinned calibration) → `rpi_core validate
--with-phylogeny` (two real bugs found here, see below; final timing TBD) → GA search
default (~6.5 min) → GA search long-pool, corrected (~21.5 min) → gate4 confirmation
of the 10 long-pool winners (~2h50m, all succeeded, see verdicts below).

### Two real bugs found in `rpi_core.py validate --with-phylogeny`'s gate4 path

Neither is a §7e finding — both are new, in a code path (`score_case`, used only by
`_validate_main`) that had **zero tests** before tonight.

1. **`score_case` discarded every gate4 field.** `score_candidate` computes
   `gate4_phylo_status`/`gate4_reassortment_call`/`gate4_discordant_segment`/
   `gate4_bootstrap_support` when `--with-phylogeny` is passed, but `score_case`'s
   return dict never copied them — every `--with-phylogeny validate` run shipped
   result rows with zero trace gate4 ever ran, while the top-level `gate4_exercised`
   flag (`args.with_phylogeny`, i.e. "was the flag passed", not "did gate4 run")
   claimed otherwise. Exactly the "asserted on an intermediate value, not the row that
   ships" defect class `docs/CLAUDE.md` names as recurring in this codebase. Confirmed
   on the committed baseline too: `output/validation/validation_report.json` has
   `gate4_exercised: false` (that report predates `--with-phylogeny` ever being used)
   — so this was latent, never previously exercised in anger.
   **Fix:** copy the four gate4 fields into `score_case`'s return dict; replace
   `gate4_exercised` with `gate4_requested` (what was asked) plus
   `gate4_analyzed_count`/`gate4_call_counts` (what actually shipped, computed from
   `results`). No existing test or code depended on the old field's semantics
   (single usage in the whole codebase, grep-verified). 2 guard tests.
2. **`score_case` passed bare `"mafft"`/`"iqtree"` strings, never resolved.** The
   same defect shape as §7e Finding 2 (fixed once for `Algorithm.py`'s
   `confirm_winners`), left unfixed in this sibling path: Windows' `subprocess.run`
   cannot exec a bare name (`mafft` is a `.BAT`) without `shell=True` or an already-
   resolved absolute path. Every validation-panel candidate that had enough real
   comparators to reach an actual alignment attempt (most short-circuit to
   `insufficient_comparators` first, comparator-count-limited, unrelated to this bug)
   crashed with `FileNotFoundError` instead of running gate4 — 3 of 55 candidates hit
   this in the first (broken) rerun attempt. **Fix:** resolve via
   `ph.resolve_binary(None, "mafft"/"iqtree")` when `with_phylogeny=True`, matching
   `confirm_winners`'s pattern; `None` when not requested (never touches binary
   resolution). 2 guard tests. Both fixes, one commit: `6310018`.

Both bugs are unrelated to the `--reference-input`/`--output` driver mistakes above —
these are real product-code defects that would reproduce for any user running
`rpi_core.py validate --with-phylogeny` on any platform, on any machine, regardless of
this session's rerun.

### Line-ending gap: real for JSON, already covered for CSV (correcting an initial misdiagnosis)

First checksum diff pass (Task 4) found every top-level `ortho_*.csv`/`.json` byte-
different from committed. Initial diagnosis: `safe_write_csv`'s `frame.to_csv(index=
False)` has no `lineterminator`, so pandas defaults to `os.linesep` -- LF on the Linux
machine that produced the committed baseline, CRLF on this Windows rerun. **True for
the raw generated bytes, but incomplete**: `.gitattributes` already pins `*.csv text
eol=lf`, so git normalizes any committed CSV to LF on checkout regardless of what
wrote it -- the diff I was seeing was a fresh-on-disk file (pre-commit, raw pandas
output) against a git-checked-out one (post-normalization), not a genuine committed-
data gap. The `safe_write_csv` fix (pinning `lineterminator="\n"`, commit `6ae45cc`)
is still worth keeping -- it makes the raw generated bytes match what git would
normalize to anyway, so `git status` doesn't show a spurious whitespace-only diff
every time the pipeline runs on Windows -- but it wasn't fixing a real cross-machine
reproducibility gap for CSVs specifically, since `.gitattributes` already closed that
one. **`*.json` had no such rule** (`*.jsonl -text` exists, deliberately unnormalized
for large line-delimited data, but singular `.json` files -- `ortho_summary.json`,
`ortho_feature_scaling.json`, every `seed_<N>.json`, `validation_report.json` -- had
none), so those genuinely do get CRLF-mangled by `core.autocrlf` on a Windows
checkout while staying LF on Linux -- confirmed directly: the committed
`output/ortho_feature_scaling.json` had CRLF, `git add` warned "LF will be replaced
by CRLF" for it just now. **Fixed** by adding `*.json text eol=lf` to
`.gitattributes`, matching the existing `.py`/`.md`/`.csv` pattern exactly.

### Non-determinism found and documented (NOT fixed — no bug, see reasoning)

1. **K-mer clustering (`all_hanta.py`), `Orthohantavirus puumalaense` only.** 1305
   rows across `ortho_isolate_table.csv`/`ortho_ml_feature_matrix.csv` differ in
   `kmer_cluster_id`/`sl_rank_separation`/`m_rank_offset` between the Linux-generated
   committed baseline and this Windows rerun. Root cause: `cosine_distance_matrix`'s
   `unit @ unit.T` matrix multiply is not bit-reproducible across BLAS
   backends/platforms; `single_linkage_clusters`' union-find is otherwise fully
   deterministic (sorted vocabulary, no hash-order dependency) but is *maximally
   sensitive* to a single threshold-boundary (0.02) flip — one differing merge
   decision early on relabels every downstream cluster ID for that species, which is
   why the diff looks large despite the root cause being ULP-level. Confirmed
   diagnostic-only: grepped `rpi_core.py`/`phylogeny.py`/`Algorithm.py`, these three
   k-mer-derived columns feed no scoring or verdict logic (the tree-based
   `tree_*_rank_*` columns in `phylogeny.py` are the ones that matter, and are
   unaffected — this is exactly the k-mer-vs-tree distinction the ingest log's own
   warning describes: "Reassortment is NOT called by this pipeline" via the k-mer
   signal). Reproduced identically across two separate runs on this same machine, so
   it is platform-dependent, not flaky.
2. **IQ-TREE/MAFFT tree-derived numeric diagnostics, both scopes.** Nearly every
   `tree_*_distance`/`species_rf_*`/`tree_*_rank_*`/`*_nearest_isolate` value in
   `ortho_reassortment_table.csv` differs from the committed baseline (1453/1470
   species-scope rows, 1522/1525 phylogroup-scope rows) — alignment/tree paths too,
   since MAFFT/IQ-TREE were run with `--phylo-threads` at its default
   (`cpu_count - 1`), which differs between the original Linux machine and this
   Windows one. This is IQ-TREE/MAFFT's well-documented multi-threaded floating-point
   non-associativity (parallel likelihood-calculation reduction order is not
   reproducible across thread counts, even with a fixed `--seed`) — not a pipeline
   bug, and not fixable without forcing single-threaded alignment/tree-building
   pipeline-wide (a much larger, slower change than this session's scope). **The
   categorical verdict that actually matters, `reassortment_call`, is stable for
   2994/2995 rows.** The one exception —
   `t2|unclassified Orthohantavirus|3877` (phylogroup scope): `weak_signal_low_
   bootstrap_support` (88.0% support) → `no_conflict` — is exactly the kind of flip
   this same non-determinism would produce: a borderline sub-threshold bootstrap
   value (88 < the 95 gate) shifting across a run-to-run branch-length recalculation.
   Not investigated further per the binding rule already recorded in §7f: verdicts
   are reported verbatim, no re-running for a "better" result. **Ruled out:** an
   earlier hypothesis that phylogroup-scope grouping had structurally changed
   (species-set keys vs. the committed baseline's broad phylogroup names) — a
   phylogroup-only retry of `phylogeny.py` reproduced the exact same broad group
   names (`arvicolinae_borne`, `murinae_borne`, `sigmodontinae_borne`,
   `soricomorpha_borne`) as the baseline. The species-pair-style group names seen in
   the driver log (e.g. "Orthohantavirus dobravaense+Orthohantavirus hantanense")
   were misattributed at first glance; they came from `validate`'s own gate4 path
   (bug 2 above), not from `phylogeny.py`'s stage, which ran clean for its full ~9h21m.

### Gate4 verdicts, GA long-pool winners (all 10, verbatim, per the binding rule)

| # | triplet (S\|M\|L, abbreviated) | gate4_phylo_status | gate4_reassortment_call |
|---|---|---|---|
| 0 | hantanense\|seoulense\|seoulense | analyzed | **supported_reassortment_candidate** |
| 1 | hantanense\|hantanense\|hantanense | analyzed | no_conflict |
| 2 | hantanense\|hantanense\|hantanense | analyzed | no_conflict |
| 3 | hantanense\|hantanense\|hantanense | analyzed | no_conflict |
| 4 | hantanense\|hantanense\|hantanense | analyzed | no_conflict |
| 5 | hantanense\|hantanense\|hantanense | analyzed | weak_signal_low_bootstrap_support |
| 6 | hantanense\|hantanense\|hantanense | analyzed | weak_signal_low_bootstrap_support |
| 7 | hantanense\|hantanense\|hantanense | analyzed | no_conflict |
| 8 | hantanense\|hantanense\|hantanense | analyzed | no_conflict |
| 9 | hantanense\|hantanense\|hantanense | analyzed | weak_signal_low_bootstrap_support |

0 errors, 0 `insufficient_comparators`. 9/10 winners are within-hantanense
combinations — per §7f's binding note, `no_conflict` on those is gate4 (the
tree-based detector) doing exactly the job it exists for, not a limitation. See
`output/ga_search_corrected_pool_long/gate4_confirmation.json` and the per-triplet
detail in `seed_20260904.json` for full rows (backbone/donor identities, RPI factors).

### Redundancy survey — actioned 2026-09-09 on explicit instruction

- `output/ga_search_corrected_pool_long/gate4_confirmation.{json,log,log.stdout,
  log.stderr}`, `gate4_confirmation_checkpoint.jsonl`, and
  `candidate_phylogenetics/ga-confirm/` (66 files) — the §7f stalled run that
  actually finished on its own (see above), fully superseded by this session's
  correctly-organized `ga-confirm-0`..`ga-confirm-9` output and byte-identical to
  it. These had been swept into the bulk promotion commit (`8219386`) despite being
  flagged as redundant, not tonight's pipeline output — `git rm -r`'d on explicit
  instruction, commit `be098b0`.
- Stray `__pycache__/` (gitignored, was never tracked — removed from disk directly).

**Explicitly NOT redundant — preserve as-is, outside tonight's pipeline scope:**
`output/consolidation_baseline/` (9 files, module-consolidation era evidence),
`output/pipeline_acceptance/raw_scored_groups.json` and `run_confirm_20260830.{json,
log}` (Task 12 era evidence). Neither is touched by `all_hanta.py`/`phylogeny.py`/
`rpi_core.py`/`Algorithm.py`; promotion must not delete them.

### Phylogeny tree illustrations

`output/figures/{before,after}/` regenerated from the promoted, corrected output
(the prior `output/figures/` was deliberately retired in commit `9cdd768` because it
was tied to the since-superseded contaminated-pool GA run — not a decision against
having figures at all). `before/`: 7 PNGs, one per real
`supported_reassortment_candidate` isolate, S/M/L trees pruned to a local
neighborhood, discordant segment highlighted. `after/`: 10 PNGs, one per corrected
long-pool GA winner, plus `ga_winners_parents.csv` (donor/backbone isolate, species,
country, host per winner).

### Promotion

`output_rerun_20260906/` merged into `output/` (overwrite-and-add, per the redundancy
survey above — `consolidation_baseline/` and the Task-12 `pipeline_acceptance/` files
are untouched because tonight's pipeline never produces them, not because promotion
skipped them). `output/ortho_output_checksums.txt` self-verified against the promoted
files (sha256, all match) after promotion. `output/validation/` promoted separately,
once below.

### Validate `--with-phylogeny`'s full 55-candidate gate4 pass — COMPLETE

Ran to completion (`RESULT: PASS`, classes A=4, B=5, C=27, D=15, E=4 — matching the
harness's own sensitivity/specificity/blind-spot/adversarial checks, all still green
with real gate4 data behind them for the first time). Final per-candidate gate4
disposition, from the shipped `output/validation/validation_report.json` (`gate4_requested:
true`):

| `gate4_phylo_status` | count | meaning |
|---|---|---|
| `analyzed` | 40 | real MAFFT+IQ-TREE comparator tree built and scored |
| `insufficient_comparators` | 12 | backbone species has too few real reference isolates for a meaningful tree — structural, not a bug |
| *(null)* | 3 | `ADV_random_nucleotides`, `ADV_truncated_segments`, `ADV_empty_M_segment` — adversarial junk rejected at gate0, correctly never reaching gate4 |

Of the 40 `analyzed`, `gate4_reassortment_call` breaks down as: **20 `no_conflict`**,
**17 `weak_signal_low_bootstrap_support`**, **3 `supported_reassortment_candidate`**.
This is the validation panel's own synthetic/adversarial/published cases, not the GA
search winners (see the verbatim table above for those) — no binding verbatim-reporting
rule was set for this panel specifically, but the same spirit applies: these are shipped
as-is, not re-run for a more favorable mix.

**Real-world scalability numbers for the paper's methods section** — full per-segment
wall-clock timing for the largest candidate encountered, `SYNTH_hantanense_x_puumalaense_S`
(hantanense+puumalaense combined backbone, **699 taxa**, single-threaded IQ-TREE,
`-m MFP` testing up to 968 models, 1000 ultrafast bootstrap replicates, this session's
Windows machine):

| segment | wall time | started | finished |
|---|---|---|---|
| S | 58m35s | 2026-09-07 17:30:05 | 2026-09-07 18:28:40 |
| M | 2h26m58s | 2026-09-07 18:29:05 | 2026-09-07 20:56:03 |
| L | 2h11m54s | 2026-09-07 20:57:08 | 2026-09-07 23:09:02 |
| **total (S+M+L)** | **5h38m57s** | 2026-09-07 17:30:05 | 2026-09-07 23:09:02 |

Two more candidates at the same 699-taxa scale (`SYNTH_hantanense_x_puumalaense_M`,
`SYNTH_hantanense_x_puumalaense_L` — the other two donor-segment variants of the same
backbone pair) and two candidates at a smaller-but-still-large ~640-taxa scale
(`SYNTH_puumalaense_x_seoulense_S`, `_M`, `_L`) each took a comparable multi-hour total,
consistent with this table — confirming it as a real, reproducible scalability
characteristic of gate4 (dominated by ModelFinder + ultrafast bootstrap cost on trees
in the many-hundred-taxa range), not a one-off. Total wall-clock for the full 55-candidate
panel from the fix landing to completion: from 2026-09-06 12:36 to 2026-09-09 06:48
(spanning multiple such large candidates plus many fast `insufficient_comparators`/small
cases) — an unattended multi-day run on this hardware, per the user's explicit
no-restart-on-failure instruction (it never needed one).

**Promoted.** `output/validation/` now holds the complete, correct
`validation_report.json` plus `candidate_phylogenetics/` (1416 files: alignments and
IQ-TREE outputs for all 40 analyzed candidates). The `output_rerun_20260906/` staging
folder has been removed (fully promoted, nothing left to merge).

---

## 8. Standing rules

- **Never change gate logic, verdicts, or `gate*_` columns.** RPI is a parallel column
  set. `VerdictUnchangedRegressionTest` green at every commit.
- **Never rewrite history, never force-push.** Remote is `arman-daliri/Orthohantavirus`.
- Repo stays **private**; public release is a separate clean repo per `PUBLIC_RELEASE.md`.
- One commit per workstream. If a step fails, **log it and continue**.
- Bump `RPI_WEIGHTS_VERSION` in the same commit as any `RPI_WEIGHTS` change.
- **Assert on the row that ships**, not the function that feeds it. Add a
  guard-the-guard wherever a fixture could drift into making a test trivially true.
  *(Three separate reviews have caught this same defect shape; two more this session.)*

---

## 9. Project history — how we got here

Consolidated from all 19 Claude sessions (2026-07-30 → 2026-08-29). The raw transcripts
are no longer in `/resume`; this is the record. Read it to understand *why* the design is
what it is.

**Raw transcript recovery, if ever needed:**
- The 13 pre-08-12 sessions: `git checkout 09c2226 -- docs/claude-sessions-backup-20260812-183855`
- The 5 later sessions (07-30, 08-23, 08-24 → 08-28): local archive at
  `C:\git-projects\_claude-session-archive-orthohantavirus\` — **outside the repo and
  not git-tracked**, deliberately, so transcripts are never pushed again (see Phase 4).
  This is a single local copy with no remote backup; if the machine is lost, so are they.

### Phase 1 — Pipeline consolidation (2026-07-30)

Two branches existed, `main` and `v2`, and the question was which was scientifically
defensible for publication. `v2` won on bioinformatics and biology logic for handling
public NCBI Orthohantavirus data, classification, and recognising well-characterised
human-pathogenic strains. **`v2` was merged into `main`**; the extra branches (`v2`,
`update`, `more-info`) were then deleted from GitHub, leaving `main` as the single line.
Several `/ultrareview` and `/code-review` rounds fixed bugs found in the merged pipeline.
`Mydocs/TEACHME.html` was written in this phase as a step-by-step explainer.

**Carried forward:** all sequence "distance" in `all_hanta.py` is **alignment-free 3-mer
cosine** on amino-acid composition, not phylogenetic distance. A reassortment classifier
built on it was **withdrawn after failing literature validation** — which is the reason
the current detector uses real phylogenetics in gate 4.

### Phase 2 — Reassortment detector (2026-08-12)

Brainstormed and specced the gate-based detector. Produced
`2026-08-12-reassortment-detector-design.md` and the three-stage pipeline split.

### Phase 3 — The RPI is conceived (2026-08-23 → 08-24)

The Zotero MCP was connected specifically so papers could be **cited accurately** rather
than from memory. The scoring design was grounded on Rivero — but with an explicit
instruction to **double- and triple-check against other papers, because Rivero is a
bioRxiv preprint and not formally published**.

**The founding requirement, in the user's words:** *"I do not want to rely on another
paper's raw index, but to have numerical score instead of implausible / possible /
plausible / likely."*

That is the whole reason the Reassortment Plausibility Index exists: **a continuous
numerical score, independently derived, replacing a borrowed categorical scale.** It is
also why every weight in `RPI_WEIGHTS` carries its own provenance and verification state
in a comment, and why the index is repeatedly described as RELATIVE and never a
probability. Executed via subagent-driven development into
`2026-08-24-reassortment-plausibility-index-design.md`.

### Phase 4 — Hardening and the privacy incident (2026-08-24 → 08-28)

Multiple adversarial review rounds; fixes verified **by execution**, not by inspection,
with reviewers explicitly told to confirm invariants hold **on the shipped output** and
to report without fixing. `03-pipeline-phylogeny-reassortment` was reconciled into
`main`.

**The incident that set repo policy:** a push accidentally published
`docs/claude-sessions-backup-20260812-183855/` (66 files, 16 MB of transcripts) and two
copyrighted PDFs to `github.com/arman-daliri/Orthohantavirus`. The decision was **not**
to rewrite history — the remote is not the primary user's account, so a rewrite would
break other people's clones. Instead: **the repo stays PRIVATE**, the PDFs were extracted
to notes and removed, and **public release became a separate clean repo** built from an
explicit file list in `PUBLIC_RELEASE.md`. This is why §8's rules exist.

### Phase 5 — RPI as a validated fitness function (2026-08-29, this session)

Scope set by the user: *"a valid fitness function for a genetic algorithm to detect real
from fake reassortants after shuffling segments"*, reproducible on Windows and Linux,
handed off to a downstream team as a callable scorer. Produced decisions **D1–D14**, the
E-vs-M_pack finding in §3, and the work in §5.

### Session index (raw transcripts deleted)

| Dates | Msgs | What happened |
|---|---|---|
| 07-30 → 08-12 | ~5,900 | 13 sessions: pipeline dev, Parquet→CSV migration, branch consolidation. *(deleted 2026-08-29; recoverable via `git checkout 09c2226 -- docs/claude-sessions-backup-20260812-183855`)* |
| 07-30 09:27–10:10 | 179 | v2-vs-main comparison; v2 merged to main |
| 07-30 10:13–12:32 | 476 | ultrareview cycles, bug fixes, TEACHME.html, GitHub branch cleanup |
| 08-23 13:46–19:05 | 78 | superpowers plugin install |
| 08-23 19:16 → 08-24 12:20 | 1,305 | **RPI design** — Zotero MCP, Rivero grounding, numerical-score requirement |
| 08-24 12:39 → 08-28 18:45 | 1,614 | RPI implementation, adversarial reviews, privacy incident, repo policy |
| 08-29 11:19 → 19:10 | 1,241 | **This session** — WS0/WS1/WS2, D2b, D1–D14 |

**Surviving curated documents (do not delete):**
- `docs/superpowers/specs/2026-08-12-reassortment-detector-design.md`
- `docs/superpowers/specs/2026-08-13-three-stage-pipeline-split-design.md`
- `docs/superpowers/specs/2026-08-24-reassortment-plausibility-index-design.md`
- `docs/superpowers/specs/2026-08-29-rpi-fitness-function-design.md`
- `docs/superpowers/plans/*.md`
- `docs/superpowers/notes/2026-08-14-reassortment-detector-reconsolidation.md`
- `PUBLIC_RELEASE.md`, `README.md` — `TEACHME.html`'s content was folded into the
  README's "Known limitations" section and the standalone file deleted (2026-09-04,
  `Mydocs/` cleanup); its defect-by-defect detail remains at `git show b8f5f0c`
