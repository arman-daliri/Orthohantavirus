# Orthohantavirus Reassortment Detection Pipeline — Consolidated Project Reference

**Single source of truth. Last consolidated 2026-09-18.** This file replaces the former
per-session `docs/PROJECT-STATE.md` narrative, the design specs and task plans under
`docs/superpowers/{specs,plans,notes}/`, and the local `.superpowers/sdd/` task ledgers —
all folded into one document and then removed from the working tree. It compresses
~19,000 lines of design specs, task plans, and session ledgers into one narrative.
Granular per-number detail (exact test counts from three sessions ago, byte-for-byte
diffs, individual commit hashes) is deliberately dropped here — it remains fully
recoverable via `git log`/`git show` on the commit that removed those files. Read this
file first; go to git history only when you need a number this file omits.

**Progress log extended 2026-09-18 (evening)** to cover 2026-09-10 → 2026-09-18 (§4 Phases
11–13, §5). That pass folded in the four remaining Claude Code chat sessions for this
project (2026-09-09 → 2026-09-18) and then deleted them; this file, `README.md` and
`git log` are now the only record. Numbers in Phases 11–13 are sourced as follows: run
logs and JSON under `output/` (read directly), `README.md` at `4df4135`, `git log`, and the
sessions' own compaction summaries (marked *(session summary)* where nothing on disk
re-confirms them).

---

## 1. Project aim

This is a **retrospective, in-silico surveillance-analysis pipeline** over public NCBI
*Orthohantavirus* genome data (genus *Orthohantavirus*, family *Hantaviridae*). Hantaviruses
have a tripartite genome (S, M, L segments) and can **reassort** — swap whole segments when
two strains co-infect the same cell. Reassortment can produce novel, potentially more
pathogenic strains, so detecting which observed or hypothetical S+M+L combinations are
biologically plausible is a real surveillance question.

The pipeline does two things, deliberately kept separate:

1. **Detect** real reassortment signal in the actual NCBI dataset, using real phylogenetic
   trees (never by calling reassortment from raw sequence similarity — an earlier k-mer-based
   classifier was built, validated against the literature, found to produce the *inverse* of
   the known species ranking, and withdrawn; see §4 Phase 1).
2. **Score plausibility** of any S+M+L combination — real, hypothetical, or GA-searched — with
   the **Reassortment Plausibility Index (RPI)**, a continuous, literature-anchored numerical
   score replacing a vaguer categorical scale. The founding requirement, in the project
   owner's words: *"I do not want to rely on another paper's raw index, but to have a numerical
   score instead of implausible / possible / plausible / likely."*

**Publication goal.** Every output, limitation, and validation result is written to the
standard a virology/bioinformatics peer reviewer would accept — hence the insistence on
primary-source citations (tiered by how directly they're relied on), explicit uncertainty
where evidence is thin, and honest reporting of negative/weak-signal findings rather than
tuning them away.

**Two standing epistemic rules that shape every design decision:**
- **RPI measures plausibility, not occurrence.** A within-species shuffle scoring as high as
  a real reassortant is *correct*, not a bug — plausibility-vs-occurrence is exactly the
  question the phylogenetic tree detector (gate 4), not the index, is built to answer.
- **All alignment-free distances (3-mer cosine k-mer composition) are a proxy, not
  phylogenetic distance.** They're retained as informational/diagnostic signal but are never
  used to call reassortment, and every RPI row that depends on this proxy (`rpi_k_stratum`)
  is marked `provisional`.

---

## 2. Architecture — the three-stage pipeline

```
all_hanta/ (raw NCBI package, gitignored)
        │
        ▼
Stage 1  prepare_dataset.py   (~20-26s)   ── the ONLY file that reads all_hanta/ncbi_dataset/
        │  writes: reference snapshot (ortho_reference_isolates/genomic/cds.json,
        │           gitignored, ~75MB derived, regenerable),
        │           ga_pool_*.csv, ortho_*.csv/json tables
        ▼
   ┌────┴─────────────────────────┐
   ▼                              ▼
Stage 2  phylogeny.py       Stage 3a  rpi_core.py         Stage 3b  Algorithm.py
(~5-9h, MAFFT+IQ-TREE,      (RPI scoring engine;          (GA search over triplets;
 real per-segment trees,     ViennaRNA bindings for        imports rpi_core as a frozen
 group-parallel since        M_pack; minutes)               scorer library; ~14 min,
 2026-09-10)                                                 gate4 confirm +hours)
        │                          │                              │
        ▼                          ▼                              ▼
ortho_reassortment_table.csv  RPI scores, calibration      ga_search/ winners,
phylogenetics/<species>/      artifact, discrimination      fitness logs, gate4
  <segment>/                  panel, validation report       confirmations
```

Stage 1 must run first; Stages 2 and 3 are independent of each other and of Stage 1's
*downstream* consumers (no stage reads another stage's output file except through the Stage 1
snapshot — this "no false dependency" design was an explicit, evaluated decision).

### Module map (current, after the 2026-08-31 single-file consolidation)

| File | Role |
|---|---|
| `prepare_dataset.py` (was `all_hanta.py`) | Stage 1 ingest: taxonomy resolution, isolate grouping, segment/completeness resolution, host/pathogenicity classification, k-mer clustering (diagnostic only), pool construction. Writes the reference snapshot every downstream stage reads. |
| `phylogeny.py` | Stage 2: real MAFFT alignment + IQ-TREE ML trees per species/phylogroup/segment, Robinson-Foulds topology comparison, bootstrap-gated reassortment calls. Group-level thread-pool parallelism added 2026-09-10. |
| `rpi_core.py` | Stage 3a: the **single consolidated module** (was 6 files: `reassortment_detector.py`, `rpi_calibration.py`, `rpi_discrimination.py`, `rpi_discrimination_stats.py`, `validate_reassortment_detector.py`, `pipeline.py`) — the 4-gate detector (gate0 sanity, gate1 k-mer pre-screen, gate2 structural annotation, gate3 ecological veto) plus gate4 (real phylogenetics, called from the candidate-scoring path) plus the entire RPI toolchain: calibration, discrimination-panel scoring, statistics, validation harness, CLI. One CLI entry point: `python rpi_core.py {run,calibrate,validate,score,diagnose-e,diagnose-gate0}`. |
| `Algorithm.py` | Stage 3b: genetic-algorithm search over segment-donor triplets, treating `rpi_core.py` as a frozen scoring oracle. Never creates/edits sequence — only recombines existing reference isolates. Runs **two pools sequentially** (`terminal_complete`, `strict_both`; §4 Phase 11), writes `ga_search/<pool>/seed_<N>.json` plus `ga_pool_comparison.json`; `--dump-scored` writes every scored triplet; `--with-phylogeny` gate4-confirms the 10 winners per pool (resolves `mafft`/`iqtree3` from **PATH** — there is no `--mafft-bin`/`--iqtree-bin` flag on this script). |
| `figures.py` | Renders before/after phylogenetic tree figures from Stage 2's committed outputs and the GA winners. |

### The gate framework (inside `rpi_core.py`, unchanged by RPI work)

| Gate | Purpose | Can it veto/fail a candidate? |
|---|---|---|
| **Gate 0** | Sanity/junk screen: segment length, ORF plausibility, terminal panhandle motif (informational only — checked in both orientations since many real deposits are reverse-complement or trimmed), k-mer identity floor | Yes — `junk_sequence` terminates scoring |
| **Gate 1** | Fast k-mer nearest-species pre-screen per segment; top-3 hits reported, not just winner | No — only decides Gate 4's comparator set; cannot itself produce `supported_reassortment_candidate` |
| **Gate 2** | Structural/functional annotation (WAASA cleavage motif, TMIS fusion residues, pH-switch residues, aligned not fixed-position; lineage-specific NP/GP markers; S/L-mismatch and receptor-class flags) | No — annotation only, never a veto, per the source report's own stated limitation (no real structural modeling) |
| **Gate 3** | Ecological opportunity veto: host/geography/time overlap of implicated parent species | Yes — `no_ecological_route` on positive evidence of non-overlap; missing metadata is `ecology_unknown`, never a fail |
| **Gate 4** | Real phylogenetic evidence: candidate injected into a live comparator tree (MAFFT + IQ-TREE, bootstrap-gated), the methodologically load-bearing gate | Yes — the only gate that can produce `supported_reassortment_candidate` |

**RPI is a parallel column set, never merged into the verdict.** `combine_verdict`'s closed
verdict list is frozen; `rpi_generation`/`rpi_retention` are computed for every candidate,
including vetoed ones, specifically so the categorical gate-3 veto and the continuous Factor E
can be cross-checked against each other (disagreement is expected and informative, not a bug —
gate3 only ever examines the first two of N implicated species, a documented, deferred
limitation).

### RPI formula

```
rpi_generation = E × M_pack × M_pol × K
rpi_retention  = rpi_generation × F
```

Multiplicative by design (necessary-condition composition, McDonald et al.'s cognate-
competition framing, and to prevent molecular elegance from rescuing ecological
impossibility) — a deliberate departure from the pipeline's other, additive scoring
(`surveillance_priority`).

| Factor | What it measures | Status as shipped |
|---|---|---|
| **E** (ecological opportunity) | host/geo/time overlap of implicated parent lineages, computed from the **specific parent isolates** (fixed 2026-08-30; previously computed from whole-species populations, which defeated the "different continents" guarantee for widely-distributed species) | Real, isolate-level where the true parent is known; falls back to k-mer-nearest reference isolate (`rpi_e_isolate_basis` flags which) |
| **M_pack** (RNA–RNA co-packaging) | Boltzmann-style rescaling of ViennaRNA-cofolded ΔΔG between the candidate's own cross-segment panhandle duplex and the backbone species' reference consensus | Real but small signal (Cliff's δ = +0.222; 38% of the calibration panel saturates at 1.0) |
| **M_pol** (protein–RNA polymerase/N recognition) | N-recognition footprint beyond the genus-invariant 14-nt terminal tag, candidate vs. backbone-species consensus | Real signal only for an *O. sinnombreense* (SNV) backbone (only species with a non-default footprint length); every other backbone reports `no_signal_beyond_tag` |
| **K** (constellation prior) | Empirical prior over which segment is "odd one out," from ~35 total published in-vitro observations, close/distant divergence-stratified, Jeffreys-smoothed | Real, but `divergence_stratum` is fed by a k-mer identity proxy, not real parental divergence — every row marked `rpi_k_stratum_provisional=true` |
| **F** (retention) | Generation→retention shift (M:S:L 12:12:1 in vitro vs. 14:2:1 in nature) — immediate post-reassortment fitness only; compensatory adaptation (McDonald's mode c) is out of scope | Weakly determined (small natural-data denominators) |

---

## 3. Key design decisions still in force

### From the RPI fitness-function decisions (D1–D14, 2026-08-29, all approved)

| # | Decision |
|---|---|
| D1 | M_pack rescaling: **calibrate first, then choose.** Calibration (WS3.1) found the ΔΔG distribution is *not* bimodal (failed 3 of 4 pre-registered bimodality checks) → **Option A, τ-scaled Boltzmann** (`M_pack = exp(-ΔΔG/τ)`) wins over the logistic runner-up. Min–max rejected outright (cohort-dependent, breaks reproducibility). |
| D2 / D2b | RNAcofold window = terminal **50 nt** (swept 25/50/75). Mechanism = **ViennaRNA Python bindings** (`RNA.fold_compound`), no CLI — the pip package ships bindings only, no executable exists. |
| D3 | M_pack/M_pol always compare the candidate's own sequence against the **implicated parent/donor's own reference consensus** — never self-comparison (an earlier bug scored the candidate against itself, trivially returning 1.0). |
| D4 | Gate 4 fixture runs in a scratch dir only, never touching `output/`. |
| D5 | SNV × BCCV cross added to published validation crosses. |
| D6 | Byte-identity claimed for RPI output and for MAFFT **only within a version**; across the two machines' different MAFFT builds (D14), only alignment-content equivalence is claimed; IQ-TREE compared at topology+support level only (never byte-identity — floating-point non-associativity across threads/platforms is real and documented, not a bug). |
| D7 | Hard fakes split into **b1** (within-species shuffle — expected to score as high as real positives; a plausibility-ceiling control) and **b2** (cross-species sympatric — expected mid-low, the group that actually tests the index). |
| D8 | Gate 1 determinism fix (`sorted()` on a set-built list that was silently hash-seed-dependent). |
| D9 | WS-API runtime: per-accession memoization; warm/cold latency always reported separately. |
| D10 | `score_triplet`/RPI returns columns and status flags only — **never a verdict, never `gate4_*` fields** in that path. |
| D11 | Single source of truth for weights: `RPI_WEIGHTS` passed explicitly, never mutated as a process global (a real bug: a sweep once left perturbed weights live for a concurrent scorer, producing a 39%-wrong score). |
| D12 | Install route: Windows = manual pinned binaries + SHA-256 checksums (bioconda has no native Windows build); Linux = conda. |
| D13 | (superseded) τ ≈ 4.86 provisional → **τ = 3.914**, the actual WS3.1 measurement. Never quote 4.86. |
| D14 | MAFFT version forced apart by platform: Windows 7.526, Linux 7.525 (no single version ships on both). |

### Other load-bearing decisions

- **Public release is a separate, clean repository** (§7 below) — never a visibility flip on
  this repo.
- **`all_hanta.py` → `prepare_dataset.py` rename (2026-09-09)** forced a repo-wide dependency
  inversion: Stage 1 became the *only* file allowed to touch raw `all_hanta/ncbi_dataset/`
  data; Stages 2/3 now read a JSON reference snapshot Stage 1 writes, instead of each
  independently re-parsing ~15 MB of raw package data. The snapshot is gitignored
  (~75 MB derived blob whose ungitignored source doesn't exist in the repo either) and
  regenerated on demand.
- **`Algorithm.py`'s GA pools must derive from `rpi_core.usable_pool`**, not a second bespoke
  filter — a code-review finding showed an earlier version searched a pool that was 51%
  contaminated and missing 65% of the legitimate isolates (fixed `105bd5c`; the empirical
  numbers from before that fix were discarded and redone). Since 2026-09-11 there are **two**
  pools by design: `terminal_complete` (= `usable_pool`, 314 rows) and `strict_both`
  (terminally-complete ∩ coding-complete, 108 rows). Stage 1 writes both **unfiltered**; the
  42 calibration-reserved isolates are removed inside `Algorithm.py` at load time, giving 272
  and 94 searchable isolates.
- **Output files are named for their consumer** (2026-09-11): `phylo_rpi_ga_*` = read by
  phylogeny + rpi_core + Algorithm (the reference snapshot), `ga_*` = Algorithm only, `rpi_*`
  = rpi_core, `ortho_*` = the Stage 2 result tables. Stage 1 writes nothing no downstream
  stage reads (6b21bb4 removed the rest).
- **The prose uses "Step N", not "Stage N"** (README convention set 2026-09-13; Step 3a =
  `rpi_core.py`, Step 3b = `Algorithm.py`). This file still says "Stage" in older sections.
- **Never change gate logic, verdicts, or `gate*_` columns.** RPI stays a parallel column set,
  enforced by a standing regression test.
- **One commit per workstream; if a step fails, log it and continue** rather than aborting an
  unattended run.
- **Assert on the row that ships, not the function that feeds it.** A recurring defect class
  (found independently at least 4 separate times across sessions): an invariant holds on an
  intermediate return value while the value that actually reaches the CSV/JSON quietly
  violates it (e.g. `rpi_e_total` rounding a real nonzero product to a literal `0.0`; a
  sensitivity sweep reporting a tautological rank correlation of 1.0 instead of a real one;
  `m_pol_variance_check` never actually being called in production despite a status field
  claiming otherwise).

---

## 4. Progress timeline

**Phase 1 — Pipeline consolidation (2026-07-30).** Two branches (`main`, `v2`) existed; `v2`
won on bioinformatics/biology-handling grounds and was merged into `main`. Several
code-review rounds followed. Carried forward as a standing limitation: all sequence
"distance" is alignment-free 3-mer cosine on amino-acid composition, not phylogenetic
distance — the reassortment classifier built on it was later withdrawn after failing
literature validation (its species-level reassortment-frequency ranking was the *inverse*
of the published literature — e.g. flagging ANDV, where reassortment is rare, at 72%, and
PUUV, where it's frequent, at only 7.7%). This withdrawal is the direct reason the current
detector's load-bearing gate (gate 4) uses real phylogenetics instead.

**Phase 2 — Reassortment detector design (2026-08-12).** Specced the 4(+1)-gate framework
described in §2, turning a human-applied heuristic into a runnable tool. Included
primary-source-verified structural anchors (WAASA/TMIS/pH-switch residues, cross-checked
against two independently solved structures that turned out to share numbering) for Gate 2.

**Phase 2b — Three-stage pipeline split (2026-08-13) + reconsolidation (2026-08-14).** Split
the monolithic `all_hanta.py` into `all_hanta.py` / `phylogeny.py` / `reassortment_detector.py`
so each stage has an explicit, documented, independently-runnable CLI. A day later, the
detector's own further split into 4 files was explicitly reversed back to one file on direct
request — recorded as a legitimate opposite judgment call, not a correctness fix.

**Phase 3 — RPI conceived (2026-08-23→24).** Zotero MCP connected specifically to verify
citations rather than rely on memory; scoring grounded in Rodriguez 1998, Rizvanov 2004,
Klempa 2018, Lowen 2018, McDonald 2016 (Tier 1, read in full), with the Rivero 2026 bioRxiv
preprint used only for structural corroboration, never for coefficients (explicitly flagged:
unreviewed, and self-reports poor cross-species transfer). Design work produced the
`E × M_pack × M_pol × K` architecture and the entire limitations list reproduced in README
(§6 below). Implementation (9 subagent-driven tasks) surfaced and fixed several defects along
the way, most notably at final whole-branch review: `m_pol_variance_check` was declared but
never invoked in production (status field lied); the leave-one-pair-out validation test held
nothing out (the prior was read from a fixed constant regardless); `modal_agrees` could not
structurally return `False` on the benchmark panel. All three fixed in one dedicated fix wave.

**Phase 4 — Hardening and the privacy incident (2026-08-24→28).** Multiple adversarial review
rounds, reviewers explicitly told to verify invariants **on the shipped output**, not by
inspection. **The incident that set standing repo policy:** a push accidentally published
66 files (16 MB) of session transcripts plus two copyrighted PDFs to the shared GitHub remote.
Decision: **never rewrite history** (the remote isn't the primary user's own account; a
rewrite would break other clones) — instead the repo stays **private permanently**, the PDFs
were extracted to notes and removed, and **public release became a separate, clean,
copy-built repository**. This is the origin of §7 below.

**Phase 5 — RPI as a validated fitness function (2026-08-29).** The single most important
finding of the whole project: as originally wired, `m_pack` was a hardcoded `1.0` and `m_pol`
compared the candidate to *itself* (always 1.0) — meaning `rpi_generation` was mathematically
incapable of distinguishing a real reassortant from a same-species hard fake; only
species-level priors (E, K) varied at all. Produced decisions D1–D14 (§3) and drove the work
that made both molecular factors score for real:

- **WS3.1 calibration** (42 isolates, 14 species, 756 non-cognate folds): measured τ = 3.914,
  confirmed the ΔΔG distribution is *not* bimodal, and found the terminal panhandle window
  needs its own genuine-terminus check (~48% of "complete" NCBI deposits lack true termini —
  a structural fact that recurs throughout the project as `has_terminal_panhandle`/pool
  filtering).
- **WS3.2–3.4**: τ wired in; M_pack and M_pol scored against real backbone-species consensus
  sequences via one shared donor/backbone determination.
- **WS4 discrimination panel** (groups a=random, b1=within-species shuffle, b2=cross-species
  sympatric shuffle, c=real positives): first run showed the pass criterion
  (`reals ~ b1 > b2 > random`) **failing** — `b2_vs_random` had the wrong sign, `reals_vs_b1`
  showed reals scoring significantly *lower* than b1 (not the expected tie). Root-caused
  (2026-08-30) to Factor E computing host/geo/time overlap from **whole-species populations**
  instead of the two specific parent isolates — defeating the random group's
  different-continent construction guarantee for widely-distributed species (Seoul, Puumala).
  Fixed at the isolate level, with an explicit proxy fallback (`rpi_e_isolate_basis`) for
  candidates whose true parent is unknown. **After the fix: `overall_ranking_supported`
  flipped False → True**, with an honestly-reported caveat (3 of the "real positive" rows'
  post-fix scores moved because their resolved proxy isolate happened to carry no
  host/geo/year metadata, triggering the "missing metadata = neutral" rule rather than
  genuine isolate-level ecological evidence — flagged, not hidden).

**Phase 6 — Single-file consolidation (2026-08-31).** Collapsed 6 RPI-related modules into
one `rpi_core.py` with one CLI (`run`/`validate`/`calibrate`/`score`/`diagnose-e`/
`diagnose-gate0`), after an earlier intermediate state left 12 files at repo root (6 shims +
2 generator scripts + 4 outside scripts) that the project owner found confusing. Discovered
mid-plan: `reassortment_detector.py`'s own CLI (`score`, Stage 3's single-candidate scoring
path) had never actually been wired into the shim dispatcher at all — fixed as part of this
same pass.

**Phase 7 — GA combinatorial search (2026-09-04→06).** `Algorithm.py` added: a genetic
algorithm over segment-donor triplets (existing reference isolates only, never synthesized
sequence), treating `rpi_core.py` as a frozen oracle. Profiling found ~70% of per-candidate
cost was rebuilding the same reference-side k-mer matrix from scratch on every call (since
candidates recombine existing isolates, that data never changes within a run) — caching it
gave a **9.5x speedup** (measured, not estimated). Along the way discovered and fixed a real
id()-reuse cache-collision hazard (garbage-collected objects can have their id recycled by an
unrelated new object, serving a stale cached result). A subsequent code-review pass found the
GA had actually been searching a **contaminated, incomplete pool** (51% wrong isolates, 65%
of the legitimate pool never searched) — fixed by routing through `rpi_core.usable_pool`.
Also resolved *why* the GA's best fitness plateaus at exactly 25/53 (0.4716981132...): it is
the exact `close`-stratum M/S constellation-prior ceiling reached whenever E=M_pack=M_pol=1.0,
and every GA candidate lands in the `close` stratum by construction (candidates are
recombined *existing* isolates, so their k-mer identity to the reference is always ≈1.0) — a
**structural ceiling of the scoring formula**, not a stuck search; confirmed three independent
ways (random sampling, a forced 300-generation run, the standard run), and explicitly flagged
that the "top-10 winners" from a run at the ceiling are an arbitrary tie-break, not a ranking.

**Phase 8 — Overnight full-pipeline rerun and reconciliation (2026-09-06→09).** Full
`all_hanta.py → phylogeny.py → rpi_core.py/Algorithm.py` rerun from scratch into an isolated
folder, reconciled against committed `output/`, and promoted. Found and fixed two real,
previously-untested bugs in `rpi_core.py validate --with-phylogeny`'s gate4 path (it silently
discarded every gate4 field it computed, and it passed unresolved binary names that crash on
Windows) — both instances of the same "assert on the row that ships" defect class. Documented,
not fixed, two genuine sources of cross-machine/cross-run non-determinism: k-mer clustering
(BLAS floating-point non-associativity, affects only diagnostic columns, not scoring/verdicts)
and multi-threaded IQ-TREE/MAFFT numeric diagnostics (well-documented parallel
non-associativity; the *categorical* `reassortment_call` stayed stable for 2994/2995 rows).
Completed the first-ever full validation run with **real gate4 data behind 40 of the 55
validation candidates** (12 returned `insufficient_comparators` — backbone species with too
few reference isolates for a tree, structural not a bug; 3 adversarial-junk cases were
correctly rejected at gate0 and never reached gate4; of the 40 analyzed: 20 `no_conflict`,
17 `weak_signal_low_bootstrap_support`, 3 `supported_reassortment_candidate`) and real gate4
confirmation on all 10 GA long-pool winners — 9 of 10 implicated
a single species (*O. hantanense*) alone, so `no_conflict` on those is gate4 correctly doing
its intra-species job (gate 1's structural blind spot), reported verbatim per a standing
binding rule (no re-running for a "nicer" result).

**Phase 9 — Stage 1 reference snapshot (2026-09-09).** Fixed a broken-import cascade caused
by an earlier rename (`all_hanta.py` → a non-importable `01-prepare_dataset.py`). Used the
occasion to fix the underlying design smell it exposed: three separate files were each
independently re-parsing the raw ~15 MB NCBI package. `prepare_dataset.py` now writes a JSON
reference snapshot that every downstream stage loads instead. A 15-task, review-gated
implementation surfaced and fixed several real bugs along the way (a genomic-index dual-keying
mismatch between the writer and reader; the snapshot's genomic-accession count double-counting
because of that same dual-keying; a completeness-overlap calculation mixing filtered and
unfiltered isolate universes). A final whole-branch review flagged one Critical (the tracked
`output/` had no reference snapshot, so every documented default command would
`FileNotFoundError` post-merge) — resolved by generating the snapshot locally to prove it
works, then gitignoring it (regenerable derived data, its raw source is gitignored too) and
documenting the required Stage-1-first step in the README, rather than committing a ~75 MB
blob.

**Phase 10 — Phylogeny group-level parallelism (2026-09-10).** `phylogeny.py`'s per-group
tree-building loop (previously strictly serial) parallelized with a thread pool
(`--phylo-workers`, default 4), submitting largest-taxa-first to minimize makespan, while
keeping output row order independent of completion order (sorted by `group_key`) for
determinism. Measured/predicted from the 2026-09-06 log: ~1.8x overall speedup
(560→310 min), since runtime is dominated by a handful of giant groups (puumalaense 425 taxa,
arvicolinae_borne 755 taxa) rather than by having many groups to churn through — workers
beyond 4 buy essentially nothing on this dataset. The default=4 shipped without its
real-tool OOM-rollback benchmark actually being run (the development environment lacks
mafft/iqtree) — an explicit, adjudicated risk acceptance, not an oversight, with the
benchmark task left open for whenever a Linux/conda environment is available.

**Phase 11 — Windows full reruns, dual GA pools, branch merge (2026-09-09 → 09-11).**
- *Pools.* The GA pool was first narrowed to a strict-merged set (`bf65f78`: terminally-
  complete ∩ coding-complete, 94 searchable isolates, down from 272), and the three
  completeness definitions (coding-complete / terminally-complete / strict merged) were made
  to compare over one shared universe (`783359a`). Best fitness stayed 0.4717 — the optimum
  was already inside the strict set. On 09-11 the single pool became **two** (`c90ccfc`):
  `Algorithm.py` runs `terminal_complete` and `strict_both` back to back and writes a
  side-by-side `ga_pool_comparison.json`.
- *First full rerun on this Windows machine (2026-09-10).* Step 1 ≈ 26 s. Step 2 with
  `--phylo-workers 4`: log start 12:01, species scope done 15:46, phylogroup scope done
  22:57 — **≈ 11 h** *(session summary of that run's log, which no longer exists on disk;
  it quoted 687 min)*, i.e. **not** faster than the
  2026-09-06 serial run on the same machine (≈ 9 h 21 m, pre-consolidation §7g). See §5
  open item on the unreproduced speed-up. Algorithm.py default run: 669 s, 1,734 unique
  triplets scored, plateau at generation 21, best 0.4717.
- *Why "the numbers differ between runs".* Answered 2026-09-10: the pool changed (272 → 94)
  and the earlier `--with-phylogeny` run had no `candidate_phylogenetics/` because the flag
  was not passed; that earlier gate4 run took 48 h+ and had stalled at 4 of 10 winners
  *(session summary)*.
- *Step 1 trimmed (`6b21bb4`)*: output writes no downstream step reads were removed.
- *Naming + scrub (`c90ccfc`, `fb5c0ae`)*: consumer-prefix output names (§3), references to
  removed files scrubbed. *Code-review fixes (`4d1cb8e`)* on `Algorithm.py`.
- *Merge (`d071055`, 2026-09-11)*: `06-module-final-simplification` merged into `main`; the
  branch was then deleted. `main` is the only line.

**Phase 12 — Clean-slate rerun and output/ reset (2026-09-11 → 09-13).** Requirement set by
the project owner: behave as if no previous output ever existed — no test may silently skip
because a folder is missing, nothing may consume an old artifact. The owner renamed `output/`
to `output_old/` (kept as a read-only historical record, never used as input) and
everything was regenerated:
- *Step 1* fresh; a second Step 1 run produced byte-identical pool CSVs and identical SHA-256
  digests for the three snapshot files (only `generated_at_utc` differs) *(session summary)*.
  Pools: 314 / 108 rows raw; 272 / 94 after removing the 42 reserved-calibration isolates.
  (The old README's "227 / 206" was a misreading: 272 is post-exclusion terminal-complete,
  206 is the old `ortho_ml_sequence_table.csv` row count. Nothing changed in the data.)
- *Step 2*: the `--phylo-workers 4` attempt was **OOM-killed**; the rerun with
  `--phylo-workers 1 --phylo-threads 2` finished **2026-09-11 20:03:51 → 09-12 06:15:57
  (10 h 12 min)**, per `output/ortho_phylogeny.log`. 9,867 isolates, 223 complete triplets,
  2,995 rows × 38 cols: 2,975 analyzed, 20 `species_below_min_taxa`, 0 failed; 8
  `supported_reassortment_candidate`, 6 `weak_signal…`, 2,961 `no_conflict`.
- *Step 3, dual GA, no gate4*: 390.4 s (`terminal_complete`, pool 272) and 406.8 s
  (`strict_both`, pool 94); both plateau-stop at generation 21 at best fitness
  0.4716981132 (= 25/53, the structural ceiling of §6 #13); 0 scoring errors; 1,526 / 1,734
  unique triplets scored (`output/ga_search/ga_pool_comparison.json`, `all_scored_*.jsonl`).
- *Gate4 confirmation of the winners* (`Algorithm.py --with-phylogeny`): attempt 1 failed
  because `Algorithm.py` has no tool-path flags (fix: put `tools/` dirs on PATH); attempt 2
  was **OOM-killed** on the large cross-species trees (≈ 486 and ≈ 698 taxa) until memory was
  freed; the from-scratch rerun completed on 2026-09-13 — ≈ 23 h wall-clock for both pools
  *(request → completion timestamps in the session, so an upper bound)*. 60 gate4 treefiles
  (10 winners × 3 segments × 2 pools), 20 winners analyzed, 0 errors *(session summary)*. **Verdicts, verbatim** (10 winners per pool):
  `terminal_complete` — 1 `supported_reassortment_candidate` (winner 0, hantanense +
  seoulense, gate1 `s_or_l_involving`), 4 `weak_signal_low_bootstrap_support`, 5
  `no_conflict`; `strict_both` — 3 `weak_signal_low_bootstrap_support`, 7 `no_conflict`.
  (README briefly carried 1 / 9 for `strict_both`; corrected in `4df4135` against the
  seed JSON.) As in §4 Phase 8, `no_conflict` on same-species winners is gate4 doing its job.
- *Tests, all five suites, 0 skipped* *(session summary, 2026-09-11/12)*: `rpi_core` 334,
  `Algorithm` 27, `phylogeny` 35, `prepare_dataset` 18, `figures` 10 (12 after `d6f9243`).
  They are `unittest` classes inside the modules — `pytest` collects 0; use
  `python -m unittest <module>`.
- *Old-vs-new comparison, then `output_old/` deleted (2026-09-13, explicit instruction).*
  Phylogeny results were **not** byte-identical (IQ-TREE thread-count non-determinism,
  as already documented in Phase 8); the GA output layout/schema differs by design; the
  `strict_both` scored-triplet JSONL **is** byte-identical, and pool row counts match.
  Caveat of record: "isolates are classified identically" was **inferred** from that scoring
  identity plus the two-run Step 1 digests, not from a direct byte diff of the old and new
  snapshot files. `output_old/` was untracked by git, so its deletion is not recoverable.
- *Docs and repo hygiene.* README restructured into ordered Steps 0–3b with Linux and
  Windows run instructions, a reproducibility section, and ASD-STE100-style tables
  (`9b43985`, `287e8da`, `974e7e7`, `fbe036a`); pipeline artefacts removed from git tracking
  (`7bf6906`) and `.gitignore` extended (`304bdbd`, `6d90786`) after `/ultrareview` refused a
  2,379-file / 4.8 M-line diff; tracked `output/` is now `ortho_reassortment_table.csv`, the
  `phylogenetics*/` trees and the two figure folders. Before/after tree figures are rendered
  per pool into `output/figures_terminal_complete/` and `output/figures_terminal_cds_complete/`
  (`b7c5a66`, `9b19954`).

**Phase 13 — Consolidation and review sweep (2026-09-18).**
- Planning docs, task ledgers and (in the evening pass) the chat sessions folded into this
  file (`7d4fbc0`, `8786071`).
- *Code-review coverage.* `rpi_core.py`, `phylogeny.py`, `Algorithm.py` have formal review
  history. `prepare_dataset.py` had none until today: reviewed, **no findings**. `figures.py`:
  two findings — `local_neighborhood` silently returned the whole unpruned tree when a tip
  was missing, and `render_before_figures` let case/punctuation-colliding isolate keys
  overwrite each other's PNG. Neither was triggered by the committed output (all 7 "before"
  keys stay distinct; "after" figures are index-named); both fixed in `d6f9243` with two
  guard tests, re-run reproduces the same 7 filenames. Two `/code-review /ultrareview` runs
  on the clean tree returned no diff and no findings.
- *Methods-draft cross-check.* A local draft ("method - sep 17", untracked, not in the working
  tree now) was checked against the code: the MAFFT / IQ-TREE / phylogeny / RPI-formula / GA
  sections (1.2–1.7) matched the implementation; the host-classification table in §1.1 was
  scrambled (category labels and counts do not match README's host table; a 338-vs-33
  inconsistency) and `rpi_retention`'s F factor was undocumented. Not edited — still open.

---

## 5. Current state (as of 2026-09-18 evening; `git log` outranks this section)

- **Branch:** `main` only, in sync with `origin/main` at `4df4135` (checked with
  `git status -sb`). Branches `05-module-consolidation` and `06-module-final-simplification`
  are merged and gone.
- **Pipeline last run end-to-end:** the 2026-09-11 → 09-13 clean-slate run of §4 Phase 12
  (Step 1 → Step 2 → dual GA → gate4 confirmation), all outputs regenerated from scratch.
  Committed results are the gate4 verdicts in Phase 12 and the tracked `output/` items named
  there; everything else under `output/` is gitignored and regenerable.
- **Test suite:** per-module counts in Phase 12 (334 / 27 / 35 / 18 / 10→12), measured
  2026-09-11/12 and **not re-run in the 09-18 pass** — re-run before quoting. **Known
  environment quirk:** a whole-module `python -m unittest rpi_core` can OOM-kill in
  constrained environments; the workaround is one test class per process.
- **Validation harness:** last full run **2026-09-09**, `PASS` — sensitivity 19/19,
  specificity 27/27, blind spot 5/5, adversarial 4/4, with real gate4 data behind 40 of 55
  candidates (§4 Phase 8). `output/validation/validation_report.json` was untracked in
  `7bf6906` and is **not on disk**; it was not regenerated in the clean-slate run.
  Regenerating it (`python rpi_core.py validate --with-phylogeny`) is a multi-day job on this
  hardware (largest candidate ≈ 5 h 39 min alone).
- **Reference snapshot:** gitignored; a fresh checkout needs
  `python prepare_dataset.py --input . --output output` before any downstream command.
- **Repo root:** `prepare_dataset.py`, `phylogeny.py`, `rpi_core.py`, `Algorithm.py`,
  `figures.py` — no shim or generator scripts.
- **Redundancy survey:** actioned 2026-09-09 (`be098b0`); nothing pending.

### Open items (none blocks the others; all need an explicit decision or a run)

1. **`Algorithm.py` §7e findings still open**: (3) a missing `data/rpi_calibration.json`
   silently returns an empty reserved set, un-reserving the 42 calibration isolates when run
   from another cwd; (4) `--diversity` can loop forever when fewer than `population` distinct
   triplets exist; (11) `--elitism 0` drops the incumbent-survives guarantee; (5) the k-mer
   and CDS caches pin ≈ 750 MB for the process lifetime (a trade-off against the 9.5× speed-up,
   not an obvious fix). A fix for 3 and 4 was drafted on 2026-09-18 and **discarded, not
   committed**, because it was made while reading a stale copy of this file. Per the standing
   rule only provably score-neutral fixes with a guard test land autonomously.
2. **README gate4 timing is internally inconsistent**: four places say gate4 adds "~11 hr per
   pool", while the Step 3b runtime paragraph says ~328 min (~5.5 h) per pool. The 09-12→13
   rerun took ≈ 23 h for **both** pools together, which fits ~11 h per pool. Reconcile before
   quoting either.
3. **Step 2 parallel speed-up is not reproduced on Windows.** README's table projects
   ~310 min for `--phylo-workers 4` (from the 09-06 log); measured on this machine: ≈ 561 min
   serial (09-06), ≈ 11 h with 4 workers (09-10), 10 h 12 min with 1 worker × 2 threads
   (09-11, after the 4-worker attempt OOM-killed). Cause not investigated (memory pressure and
   4 × (cores − 1) thread oversubscription are candidates — unverified). The real-tool
   OOM-rollback benchmark from Phase 10 is still deferred.
4. **Methods draft §1.1 host table** needs correcting against README's host table (Phase 13).
5. **`--with-phylogeny` needs `mafft`/`iqtree3` on PATH** — `Algorithm.py` has no flag for it
   (asymmetric with `phylogeny.py`'s `--mafft-bin`/`--iqtree-bin`).
6. **Three GA-winner donors (`hantanense|lr1` in `terminal_complete` winners 00/06/07) have
   no tree path in `ortho_reassortment_table.csv`**; the "after" figure prints
   "not in a committed tree" for them by design.
7. **Public release** has not been built; `docs/PUBLIC_RELEASE.md` is the list to follow.
8. **Open question carried from §6 #13:** whether a GA whose objective saturates on ~1.7% of
   the search space is the right instrument — options are a documented tie-breaker (e.g.
   `rpi_retention`), reporting the ceiling set as unordered, or making `divergence_stratum`
   use real parental nucleotide divergence (which changes scores and needs re-calibration).
9. **Working-tree hazard — read before committing anything.** On 2026-09-18 `README.md`,
   `docs/PROJECT-STATE.md`, `docs/PUBLIC_RELEASE.md` and `figures.py` repeatedly reappeared
   as "modified" holding **pre-consolidation content** (older PROJECT-STATE, a `figures.py`
   without the `d6f9243` fix); the four files carried identical-to-the-millisecond mtimes,
   i.e. one batch write. Committing that state would revert the consolidation and the
   `figures.py` fixes. Always `git diff` first; restore with `git restore <file>`. The writer
   was not identified (an editor holding stale buffers is the first thing to check).

---

## 6. Known limitations / open items (binding — reproduce before any Methods paragraph)

Reproduced in full in `README.md`'s own limitations section — this is a compressed index:

1. **The index is relative, not an absolute probability.** Ranking is meaningful; the
   absolute value is not.
2. **Every in-vitro anchor is Vero E6** (monkey kidney cell line), not reservoir rodent cells.
3. **Total observed reassortants across all direct measurements ≈ 35** — the constellation
   prior (K) rests on this thin base; no smoothing manufactures real statistical power.
4. **MOI is confounded across the two foundational in-vitro studies** (Rizvanov MOI 3 vs.
   Rodriguez MOI 1) — cross-study reassortment *rates* are never compared; only Rodriguez's
   internal contrast is controlled.
5. **Compatibility is directional in reality; the metrics are largely symmetric** (documented,
   not solved — an open question).
6. **No gold-standard label set exists** — 8 `supported_reassortment_candidate` calls out of
   ~2,995 rows in one snapshot; nothing can be fitted to 8 positives. This is why the index is
   literature-anchored rather than learned.
7. **`rpi_retention` models immediate post-reassortment fitness only** — compensatory
   adaptation over time is out of scope by design, not an oversight.
8. **M_pack and M_pol are real but weak signals** — M_pack: Cliff's δ = +0.222 ("small"), 38%
   of the calibration panel saturates at 1.0 (contributes nothing). M_pol: real signal exists
   only when the backbone species is *O. sinnombreense* (SNV); every other backbone reports
   `no_signal_beyond_tag`.
9. **The K stratum boundary is fed by a protein k-mer identity proxy, not real parental
   nucleotide divergence** — every scored row is marked `rpi_k_stratum_provisional=true` for
   this reason. This is the same alignment-free-distance limitation named in §1, recurring
   here at the K-factor level.
10. **`modal_agrees` saturates and must never be quoted alone** — every benchmark
    observation's top predicted set contains M, so it cannot return `False` on this panel;
    `composition_distance` (quantitative, can fail) is the primary agreement measure.
11. **Alignment-free 3-mer cosine k-mer distance is not phylogenetic/PED distance anywhere in
    this pipeline.** Real per-segment ML trees exist (Stage 2), but no DEmARC-style
    concatenated pairwise-evolutionary-distance value is computed — Stage 2 targets the
    reassortment question, not ICTV species re-demarcation.
12. **Gate 1 has a structural, validated blind spot for intra-species reassortment.** When
    both parents are the same species, all three segments assign to that species and gate 1
    reports `concordant` — correctly, by construction. Gate 4's per-segment phylogenetics is
    the only layer that can see intra-species reassortment; confirmed against 5 published
    intra-species crosses, all silent at gate 1 as predicted, and confirmed as *not* a test
    failure.
13. **GA search winners at the fitness ceiling are an unordered tie set, not a ranking.**
    Nothing in the search space exceeds the close-stratum constellation-prior ceiling
    (25/53 ≈ 0.4717) when E=M_pack=M_pol=1.0, which every GA-recombined candidate reaches by
    construction (recombined isolates always score k-mer identity ≈1.0 to the reference).
14. **~48% of "complete" NCBI segment deposits lack their true 5'/3' genomic termini** —
    trimmed at submission. This corrupts any panhandle-adjacent computation if not filtered
    out; `has_terminal_panhandle`/the terminally-complete pool exists specifically to exclude
    these.
15. **Cross-platform reproducibility is scoped, not absolute.** RPI arithmetic (pure
    Python/numpy) is byte-reproducible everywhere. MAFFT is byte-identical only within a
    pinned version (Windows and Linux are forced onto *different* MAFFT point releases,
    7.526 vs 7.525, because neither ships for both platforms) — cross-machine, only
    alignment-content equivalence is claimed. IQ-TREE is never claimed byte-identical across
    OS/thread-count, only topology+support equivalence (documented, real floating-point
    non-associativity in parallel likelihood calculation, not a bug to chase).
16. **Several citations remain at a lower verification tier and must be resolved before
    publication** — most notably the ViennaRNA papers (Lorenz 2011, Bernhart 2006, gate the RT
    constant) and Mir & Panganiban 2006 (gates the entire M_pol footprint-length rationale).

---

## 7. Repo policy (binding, unchanged since the 2026-08-24→28 privacy incident)

- **This repository stays PRIVATE, permanently.** It is the single working repo: code, data,
  docs, session notes all live here together. Never proposed to be split out.
- **Never rewrite history. Never force-push.** `git filter-repo`, `filter-branch`, rebase
  onto pushed commits, and `push --force` are all off the table — the remote is not the
  primary user's own account, so a rewrite would break other people's clones. If something
  sensitive is ever found in history, the answer is "the repo is private," never a rewrite.
- **Public release is a SEPARATE, CLEAN repository**, built by *copying* an explicit file
  list into a fresh repo — never by flipping this repo's visibility. The authoritative,
  current list lives in `docs/PUBLIC_RELEASE.md` — read it in full before any release work; it
  enumerates exactly what to copy (the three stage scripts, `data/rpi_calibration.json`,
  README, CLAUDE.md, environment files, small example data only) and what to never copy
  (session backups, the full `all_hanta/` package, `output/`, `.superpowers/`, absolute local
  paths, any non-CC-BY PDF).
- **This file is itself an internal working document**, not part of the public release list
  unless explicitly added.
- **If unsure whether a file belongs in the public release: ask. Do not guess.**

---

## Provenance of this consolidation

On 2026-09-18, at the project owner's request, this file absorbed and replaced:
`docs/superpowers/plans/*.md` (8 files), `docs/superpowers/specs/*.md` (10 files),
`docs/superpowers/notes/*.md` (1 file), and the local (gitignored, never tracked)
`.superpowers/sdd/*/` task-execution ledgers from four prior sessions. All of the above
remain recoverable from git history (the tracked files) up to the commit immediately
before this one; the gitignored SDD ledgers were local-only working material and are not
recoverable once deleted. `README.md` and `docs/PUBLIC_RELEASE.md` were updated in the same
commit to stop pointing at the removed files.

**Chat sessions, 2026-09-18 evening.** At the owner's request the remaining Claude Code
sessions for this project were read (user prompts, assistant outcomes and each session's
compaction summaries), their pipeline progress was folded into §4 Phases 11–13 and §5, and
they were then deleted, leaving only the session that performed the consolidation. Deleted:
`7a348444` (09-09, an earlier consolidation + Stage 1 snapshot session), `d9309959`
(09-09 → 09-13, phylogeny parallelism, reruns, README restructure — the large one), and two
empty Remote-Control bridge stubs (`1dd9bd43`, `68c5589a`). Raw transcripts are not
recoverable; anything not written here or in git is gone. The manual copy at
`~/.claude/projects/.../backup_sessions_2026-09-05/` (eight older transcripts) was **not**
touched.
