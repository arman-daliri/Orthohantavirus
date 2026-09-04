# Combinatorial search over segment-donor assignments — design spec

Date: 2026-09-04. Branch: new branch off `main` (post-merge of `05-module-consolidation`). Status: approved design, pending implementation.

## 1. Purpose and scope

This repo is a retrospective, in-silico surveillance-analysis project over public NCBI
orthohantavirus sequences. `rpi_core.py` is the validated plausibility scorer (RPI): it
ranks how plausible a combination of three existing genome segments is.

This spec adds a classic population-based combinatorial search (genetic algorithm, per
Holland 1992, summarized in `GA.md`) that searches the space of **segment-donor
assignments** and ranks them with the existing scorer.

Hard scope rule: an individual is a triplet of **existing reference isolates** assigned
to segment slots. Nothing here creates, synthesizes, edits, or modifies any sequence.
The output is a ranked list of plausible constellations for further study. Motivation:
the full triplet space (pool³) is far too large to score exhaustively; guided search
samples thousands of candidates instead of billions.

## 2. Encoding

- Individual: `(s_key, m_key, l_key)` — three isolate keys drawn from the usable pool.
- Pool: same definition as the discrimination panel's usable pool — coding-complete
  triplet, not lab-derived, not vaccine strain — and **excluding the calibration-reserved
  isolate keys** recorded in `data/rpi_calibration.json` (the repo's anti-circularity
  rule applies to this search too).
- Fixed-length, 3-slot, integer-label encoding (Holland's bit strings become per-slot
  donor labels; a one-slot change produces a meaningful score change, satisfying
  Holland's representation requirement).

## 3. Operators

- **Initialization:** uniform random triplets from the pool, seeded RNG.
- **Selection:** tournament (size 3). Not fitness-proportionate: RPI scores cluster near
  floors and ceilings, where proportionate selection degenerates.
- **Crossover:** single-point over the 3 slots (S|ML or SM|L swap), two children per
  pair. This mirrors the segment-swap the scorer was built to evaluate.
- **Mutation:** per-slot probability (default 0.02) replaces one donor with a uniform
  random pool isolate. Never touches sequence content.
- **Diversity maintenance (optional, default OFF):** deduplicate triplets per generation
  and niche by (backbone species, donor-species pair). OFF when the goal is the single
  most plausible constellation; ON when the goal is a diverse watchlist of distinct
  high-scoring constellations.

## 4. Oracle wiring (fitness evaluation)

- `import rpi_core as rc`; `rc.load_reference_data(Path("allhanta"))` once;
  `rc._build_precomputed(reference)` once; `args = rc._score_args()` — the existing
  fast mode with the phylogenetic gate disabled, the same convention the validation
  harness uses.
- Per individual: `rc.score_candidate(rc.CandidateInput(candidate_id=..., s_genomic=...,
  m_genomic=..., l_genomic=...), reference, precomputed, args)`; fitness =
  `row["rpi_generation"]`; a `junk_sequence` verdict maps to fitness 0.0.
- Memoize on `(s_key, m_key, l_key)` — the scorer is deterministic (frozen
  `RPI_WEIGHTS`, pinned τ), so a triplet is scored at most once per run.
- Frozen configuration: no weight, τ, gate, seed, or panel artifact is modified anywhere
  in the loop. Validation artifacts (frozen panel, published crosses) are never read
  inside the loop.

## 5. Loop and budget

- Defaults: population 150, generations 100, elitism 2 (best two carried unchanged).
- Stop: generation budget reached, or best score unchanged for 20 generations.
- Log per generation: best / mean / median fitness and count of distinct triplets.
- CLI: `python Algorithm.py --seed 20260904 --population 150 --generations 100
  --mutation-rate 0.02 [--diversity]`. Every run embeds its full config and seed in its
  output.
- Output: `output/Algorithm/seed_<N>.json` (per-generation log + final winners with
  their complete score rows) plus a short human-readable summary.

## 6. Expected outcome and evaluation (stated up front)

- **Expectation:** convergence toward ecologically coherent triplets resembling panel
  group b1 (within-species-like). That is the RPI landscape's designed ceiling, not a
  failure of the search.
- **Evaluation protocol (outside the loop):** compare final winners against the frozen
  seed-20260829 discrimination panel; report whether winners show the segment-exchange
  pattern the literature reports most (M most often donated, S+L co-inherited), and any
  overlap with the panel's real-positive group.
- **Re-confirmation:** final winners go once through the full gate path including
  phylogenetics. If the Windows MAFFT/IQ-TREE blockers in PROJECT-STATE.md are still
  open, record "confirmation deferred" — never substitute stub binaries.

## 7. Implementation constraints (repo rules)

- New file `Algorithm.py` importing `rpi_core` as a library. **`rpi_core.py` is frozen
  and validated — do not modify it.**
- Tests in-file, in this repo's style: seeded determinism (same seed → identical run),
  memoization actually hits, junk verdict → 0, crossover swaps whole slots only,
  mutation replaces donors only, pool excludes calibration-reserved keys.
- Full test suite green (324 existing + new) before any commit. The user runs all git
  commands.

## 8. Out of scope

- Base-level sequence generation or editing of any kind; GA.md's inversion and
  classifier-system sections; re-tuning any RPI weight, τ, or threshold; modifying the
  frozen panel artifacts.

## 9. Open items

- Pool choice: the panel's usable pool (recommended, keeps results comparable with the
  frozen panel) vs. the smaller coding-complete ML pool.
- Where the final gate-4 confirmation runs, given the open Windows binary blockers.
