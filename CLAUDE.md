# Orthohantavirus pipeline — working notes for Claude sessions

## RESUME HERE (last session: 2026-08-29)

**READ `docs/PROJECT-STATE.md` FIRST.** It is the single source of truth for resuming:
current state, all decisions D1-D14, measured performance numbers, what was completed,
and the exact next action. This section is only a pointer.

**State: green, WS3.1–WS3.4 committed locally, not yet pushed.** Branch
`04-rpi-fitness-validation`. Suite **228 tests, all pass** (~14 min) + **39 tests, <1 s**
(`rpi_calibration`). Validation harness **PASS**. `main` pushed, unaffected by this branch.

**WS3.1–WS3.2 ARE DONE.** τ = 3.914 (WS3.1's D1 branch: Option A, τ-scaled Boltzmann —
ΔΔG is NOT bimodal). `RPI_WEIGHTS["m_pack_tau"]`/`["m_pack_window_nt"]` set, version
bumped to `"2026-08-29.1"`. **User approved this branch 2026-08-29.**

**WS3.3 + WS3.4 ARE DONE — M_pack and M_pol both score for real now.** Both call
`rpi_species_assignment(gate1, odd_segment)` (renamed from `mpack_species_assignment` —
shared, not M_pack-only) for ONE donor/backbone determination per candidate.
- `score_m_pack`: candidate's own cross-segment panhandle duplex vs the backbone
  species' own reference consensus, `M_pack = exp(-ΔΔG/τ)`.
- `score_m_pol`: candidate's own N-recognition footprint (beyond the invariant tag) vs
  the backbone species' own reference consensus CDS window. Self-comparison
  (`polymerase_compatibility(terminal_windows, terminal_windows)`) is gone. Bundled
  in: `footprint` now reads `n_recognition_footprint(backbone_species)` — SNV's real
  23-nt asymmetry wired to production for the first time.

**NEXT ACTION: WS4 (Task 12+) — validation, the core deliverable.** WS3 is done; WS2
Task 7 (cross-OS fixture) also remains, independently. See `docs/PROJECT-STATE.md` §7b.

**Two weak-signal findings, user-acknowledged 2026-08-29 — proceed as planned, do not
retune.**
- M_pack: ΔΔG separates cognate/non-cognate only weakly (Cliff's δ = +0.222; 38% of the
  calibration panel has ΔΔG ≤ 0, where M_pack saturates at 1.0).
- M_pol: real signal exists ONLY for an SNV backbone (the sole species with a
  non-default `RPI_N_FOOTPRINT_NT` entry). Smoke-tested on the WS4 panel: 3 of 12
  candidates (SNV-backbone) scored 0.74–0.96; the other 9 honestly reported
  `no_signal_beyond_tag` rather than a fake constant.

Neither is a bug to fix by retuning τ or inventing footprint lengths without a
literature source — **WS4 decides whether M_pack + M_pol together are enough** to
separate the b2 cross-species sympatric group, or whether the deferred Factor E rebuild
(gate3's KNOWN LIMITATION) matters more.

**Binding on WS4 Task 13:** the b1/b2/random builders MUST exclude the 42 isolate keys in
`data/rpi_calibration.json` → `reserved_for_calibration.isolate_keys`, or τ is tuned on
the validation set.

**`m_pol_variance_check` is still NOT wired in anywhere.** It is spec 6.2's RUN-LEVEL
check over a representative sample of lineage pairs, not a per-candidate computation —
building that aggregation was never in WS3.4's stated checklist and remains future work.

**The finding that shaped WS3:** E, K and F are species/segment-level priors and cannot
see which particular S+M+L assembly is being scored. Only M_pack and M_pol can, and
both are now real (above). Within-species shuffles (b1) scoring high is CORRECT, not a
false positive — RPI measures plausibility; occurrence is the tree detector's layer.

**Do not redo:** N4/N5/N6 fixed; gate1 set-ordering nondeterminism fixed (`sorted()` at
`gate1_check`); `_score_panel` no longer mutates the `RPI_WEIGHTS` global; RNAcofold
subprocess path replaced by the ViennaRNA Python bindings (no CLI fallback); MAFFT 7.526
and IQ-TREE 3.1.3 installed and checksummed in `tools/`.

**Known limitation, deliberately deferred:** `gate3_check` judges only the FIRST TWO of N
implicated species. Sorting made that choice stable, not principled. Deferred to the
Factor E rebuild — see the KNOWN LIMITATION comment in `gate3_check`.

**τ ≈ 4.86 is DEAD — quote τ = 3.914.** WS3.1 measured it under N = 50 (D13 superseded).
The old figure came from a ΔΔG span of 22.4 taken before the window was chosen; the
measured robust p5–p95 span is 18.03.

**Two corrections WS3.1 forced, do not reintroduce:**
- **The 3′ panhandle window is NOT reverse-complemented** (the plan said it was).
  RNAcofold pairs both 5′→3′ strands antiparallel, and the genomic 5′/3′ windows already
  are reverse complements — that is what a panhandle is. As-is folds −30.9…−45.3 vs
  −7.3…−15.4 reverse-complemented. Pinned by `CofoldOrientationTest`.
- **48% of deposits lack their true termini** and `SEGMENT_LENGTH_RANGE` does not catch
  it. Any panhandle work must filter on `has_terminal_panhandle` or it folds arbitrary
  internal coding sequence (−10.0 vs −42.2 kcal/mol).

## Repository policy

**This repo stays PRIVATE.** It is the single working repo: code, data, docs, notes and
session backups all live here together. Do not propose splitting them out.

**Never rewrite history. Never force-push.** `git filter-repo`, `filter-branch`,
`rebase` onto pushed commits, and `push --force` are all off the table. The remote is
`arman-daliri/Orthohantavirus` — not the primary user's account — so a rewrite would
break other people's clones. If something sensitive lands in history, the answer is
"the repo is private," not a rewrite.

**Public release is a SEPARATE, CLEAN repository**, built by copying an explicit file
list into a fresh repo. Never by flipping this repo's visibility. The list lives in
`PUBLIC_RELEASE.md` — read it before any release work, and keep it current as the
project grows.

**Saving session backups and working notes under `docs/` is fine** while the repo is
private. They are excluded from public release by `PUBLIC_RELEASE.md`, not by deletion.

**If unsure whether a file belongs in the public release, ask.** Do not guess.

## Environment

Full detail in `tools/VERSIONS.txt` and `requirements-pip.txt`. Summary:

- Active interpreter is **system Python 3.14.6**, not the `hanta_phylo` conda env.
  `conda`/`mamba`/`micromamba` are not on PATH here.
- **RECONCILED 2026-08-29.** `requirements-pip.txt` documents the system-Python path
  (numpy 2.5.1, pandas 3.0.5, biopython 1.88, dendropy 5.0.13, viennarna 2.7.2), and
  `EnvironmentReconciliationTest` asserts it does not contradict `environment.yml`.
  `environment.yml` now describes the **Linux/conda machine only**.
- **`mafft` and `iqtree` ARE NOW INSTALLED** under `tools/` (gitignored binaries;
  checksums in `tools/CHECKSUMS.sha256`). Gate 4 has still never been exercised against
  them — two things block it: the executable is **`iqtree3.exe`**, not `iqtree`, and
  `validate_reassortment_detector.py:464` passes the literal string `"iqtree"`.
  Tests still substitute fake binaries; on Windows those need a `.cmd` shim
  (`_write_fake_binary`) because `CreateProcess` cannot exec a `.py` directly.
- **MAFFT version split (D14):** Windows **7.526**, Linux **7.525**. Forced, not chosen —
  MAFFT publishes no 7.525 Windows build and no 7.526 source release.
- **ViennaRNA 2.7.2 IS installed**, as **Python bindings only** — the pip package ships
  no `RNAcofold` executable. The production path uses `RNA.fold_compound`/`mfe_dimer`
  via `cofold_dg` (D2b); the subprocess/CLI path was removed, with no fallback.
  `M_pack` is still hardcoded to 1.0 at the call site until WS3.3 lands.

## Running things

```bash
python -m unittest rpi_calibration                # 39 tests, <1 s
python rpi_calibration.py                         # rebuilds data/rpi_calibration.json (~12 s)
python rpi_calibration.py --stability             # branch across seeds x windows (~2 min)
python -m unittest reassortment_detector          # 194 tests, ~13.5 min (gate 4 is slow)
python validate_reassortment_detector.py --json output/validation/validation_report.json
python reassortment_detector.py --rpi-report ...  # writes rpi_validation_report.json
```

## Known limitations to keep in mind

- **Gate 1 cannot see intra-species reassortment.** It assigns each segment to its
  nearest reference *species*; same-species parents all assign to that species, so it
  reports `concordant`. This is the class Rodriguez 1998, Klempa 2018 and Park 2025 all
  describe. Gate 4 is the layer that would catch it, and gate 4 has never been validated
  against real binaries. Do not read "0 false positives" as "detects all reassortment."
- **`rpi_generation` is currently `E × K`**, not the four-factor product the schema
  advertises — `M_pack` and `M_pol` are both scaffolded at 1.0. **This is why RPI cannot
  yet discriminate cross-species sympatric shuffles: E, K and F are species/segment-level
  priors, so only M_pack and M_pol read the actual assembly being scored.** WS3.3/WS3.4.
- **`rpi_k_stratum` is derived from a protein k-mer proxy**, not parental nucleotide
  divergence. Every row carries `rpi_k_stratum_provisional = true`.
- See the README's "Reassortment Plausibility Index — limitations" section (items 1–10)
  before quoting any number.

## A recurring defect class in this codebase

Three separate reviews have caught the same shape: **an invariant asserted on an
intermediate value while the value that actually ships quietly violates it.**

- `rpi_m_pol_status` reported `"scored"` for a constant 1.0
- the sensitivity sweep reported rank correlations that were arithmetic identities
- `rpi_e_total` rounded to a literal `0.0` while the test asserted on the raw dict

**Assert on the row that ships, not on the function that feeds it.** When adding a test
for a guarantee, check it against the CSV/JSON artifact, not the intermediate return.
Also add a guard-the-guard assertion where a fixture could drift into making the test
trivially true.
