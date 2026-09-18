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
| `Algorithm.py` | Stage 3b: genetic-algorithm search over segment-donor triplets, treating `rpi_core.py` as a frozen scoring oracle. Never creates/edits sequence — only recombines existing reference isolates. |
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
- **`Algorithm.py`'s GA search pool must be `rpi_core.usable_pool`**, not a second bespoke
  filter — a code-review finding showed an earlier version searched a pool that was 51%
  contaminated and missing 65% of the legitimate isolates.
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
Completed the first-ever full validation run with **real gate4 data behind all 55 validation
candidates** and real gate4 confirmation on all 10 GA long-pool winners — 9 of 10 implicated
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

---

## 5. Current state (as of the latest recorded session, 2026-09-10 branch work; see
`git log` for anything after)

- **Branch:** work has progressed through `05-module-consolidation` (merged) into
  `06-module-final-simplification` for the phylogeny-parallelism work. Always trust
  `git log --oneline -20` over any stale narrative, including this one.
- **Test suite:** `python -m unittest rpi_core` has run as high as **330+ tests**, all
  passing, at various checkpoints (test count grows session to session — re-run before
  quoting). **Known environment quirk:** the whole-module invocation reliably OOM-kills in
  constrained environments; the established workaround is running each of `rpi_core.py`'s
  ~80 test classes as its own process.
- **Validation harness:** last full run (with real gate4 data behind all 55 candidates) —
  **PASS**: sensitivity 19/19, specificity 27/27, blind spot 5/5, adversarial 4/4.
- **Reference snapshot:** gitignored, regenerated on demand by running Stage 1 first; a
  fresh checkout needs `python prepare_dataset.py --input . --output output` before any
  downstream command.
- **Repo root is down to the three stage files** (`prepare_dataset.py`, `phylogeny.py`,
  `rpi_core.py`) plus `Algorithm.py` and `figures.py` — no more shim/generator scripts.
- **Open, not yet actioned:** a redundancy survey of superseded output artifacts (flagged,
  awaiting explicit go-ahead, nothing deleted); the phylogeny-parallelism real-tool OOM
  benchmark (blocked on an environment with mafft/iqtree); several LOW-severity
  code-review findings deliberately parked as not provably score-neutral (`Algorithm.py`'s
  `--diversity` infinite-loop risk under a tiny pool; `--elitism 0` losing the
  incumbent-survives guarantee; a systemic scoring failure being indistinguishable from a
  genuinely flat 0.0 landscape without an error-rate abort threshold).

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
