# Orthohantavirus — Project State

**Single source of truth for resuming work. Written 2026-08-29.**
A fresh session should need this file, `CLAUDE.md`, and the two documents in
§7 — and no chat transcript.

---

## 1. Where things stand right now

| | |
|---|---|
| **Branch** | `04-rpi-fitness-validation`, pushed, tracking `origin/` |
| **HEAD** | `0abd77e` |
| **`main`** | `e7d9f44`, pushed (was 4 commits behind at session start) |
| **Test suite** | **194 tests, all pass**, ~13.5 min (`python -m unittest reassortment_detector`) |
| **Validation harness** | **PASS** — sensitivity 18/18, specificity 27/27, blind spot 5/5, adversarial 4/4 |
| **Working tree** | Clean except `Mydocs/My Library.csv`, deleted by an external process; user said leave it |
| **Unpushed** | Nothing |

**Current goal:** make the Reassortment Plausibility Index (RPI) work end-to-end and
demonstrate it separates real reassortants from shuffled fakes, reproducibly on Windows
and Linux, then hand a scoring function to a downstream team.

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

`compute_rpi` takes six inputs. E (ecology), K (constellation prior) and F (retention)
are **species-level or segment-level priors** — they do not read the nucleotides of the
particular S+M+L assembly being scored. Only **M_pack** and **M_pol** do.

Today `m_pack` is hardcoded `1.0`, and `m_pol` is computed by passing
`terminal_windows` as **both** arguments — a self-comparison, so it is exactly 1.0 on
every candidate. `divergence_stratum` is a binary threshold at 0.80 and production
identities run 0.35–0.46, so every row lands `distant` and the stratum never varies.

**Consequence:** for a fixed species set and odd segment, `rpi_generation` has exactly
one possible value. Two assemblies matching on species/place/year get identical scores,
necessarily.

**But this is only a defect for some comparisons** (corrected by the user, decision D7):

- **Within-species shuffles (b1)** *should* score high. They are genuinely plausible.
  RPI measures **plausibility, not occurrence**; separating b1 from real reassortants is
  the **tree detector's** job, not the index's. `reals ≈ b1` is an **expected tie**.
- **Cross-species sympatric shuffles (b2)** are where the index must earn its keep, and
  where it provably cannot today.

**This is why WS3.3 and WS3.4 are mandatory, not optional.**

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
| **D13** | **τ ≈ 4.86 is PROVISIONAL** pending WS3.1. **Do not quote as settled anywhere.** |
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

## 7. NEXT ACTION — WS3.1 calibration panel

**Approved and ready to start. No open questions.**

Build `rpi_calibration.py` → `data/rpi_calibration.json`.

**Cognate definition (user-approved this session):**
- **Cognate = within-isolate** — S/M/L from the *same* isolate. The one reading nobody
  disputes: segments that demonstrably travelled together in a real virion.
- **Non-cognate = across species**, at matched sampling depth.
- **Both restricted to isolates not used in any WS4 group**, so the disjointness
  guarantee is real and τ is not tuned on the validation set.

**Then:**
1. Fold via `cofold_dg` (bindings, D2b) over the 50-nt panhandle windows.
2. Record **mechanism + version** via `viennarna_provenance()` — the whole τ derivation
   depends on which implementation produced the energies.
3. Test bimodality (dip statistic or 2-component GMM BIC) and report the gap.
4. **Branch (D1):** bimodal with a clean gap → logistic; otherwise τ-scaled Boltzmann.
5. **Report the branch outcome to the user before implementing WS3.2.**

**Then WS3.2/3.3/3.4 → WS4 → WS-API → WS5.** WS2 Task 7 (cross-OS fixture) also remains.

### Gate 4 blockers to handle when it is first run
- The executable is **`iqtree3.exe`**, not `iqtree`. `parse_args` defaults
  `--iqtree-bin` to `None`; `validate_reassortment_detector.py:464` passes the literal
  `"iqtree"`. Both need the real path threaded through.
- **MAFFT's `--version` is not on the first line** — the `.bat` prints a "Preparing
  environment" preamble. Scan for the `v7.x` token.

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

## 9. Archived chat sessions (deleted 2026-08-29, recoverable from git)

`docs/claude-sessions-backup-20260812-183855/` held **13 raw JSONL transcripts, 66
files, 16 MB**, spanning **2026-07-30 → 2026-08-12**. They were raw logs, not curated
notes, and their outcomes are already captured in the design docs listed below. Deleted
to keep the working tree clean; recover with
`git checkout <commit-before-deletion> -- docs/claude-sessions-backup-20260812-183855`.

| Dates | Msgs | Topic | Superseded by |
|---|---|---|---|
| 07-30 | 10 | Protein-count discrepancy (15504 vs NCBI 15568) | — |
| 07-30 → 07-31 | 1797 | Main pipeline development | `all_hanta.py`, README |
| 07-30 | 25 | Parquet file troubles | Parquet→CSV migration |
| 07-31 | 147 | Tooling/model questions | — |
| 07-31 | 27 | Claude-desktop setup | — |
| 07-31 | 24+15+32 | ML framework / pipeline / curation design briefs | `2026-08-12-reassortment-detector-design.md` |
| 08-01 | 460 | Pipeline Q&A, output differences | README |
| 08-01 | 72 | Parquet output issues | Parquet→CSV migration |
| 08-01 → 08-12 | 1824 | Parquet→CSV conversion, pipeline-complete branch | `all_hanta.py` |
| 08-02 → 08-12 | 1102 | Pipeline explanation and extension | `2026-08-13-three-stage-pipeline-split-design.md` |
| 08-12 | 337 | Brainstorming → reassortment detector | `2026-08-12-reassortment-detector-design.md` |

**Surviving curated documents (do not delete):**
- `docs/superpowers/specs/2026-08-12-reassortment-detector-design.md`
- `docs/superpowers/specs/2026-08-13-three-stage-pipeline-split-design.md`
- `docs/superpowers/specs/2026-08-24-reassortment-plausibility-index-design.md`
- `docs/superpowers/specs/2026-08-29-rpi-fitness-function-design.md`
- `docs/superpowers/plans/*.md`
- `docs/superpowers/notes/2026-08-14-reassortment-detector-reconsolidation.md`
