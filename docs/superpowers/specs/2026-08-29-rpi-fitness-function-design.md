# RPI as a Fitness Function — Design Doc

**Status: APPROVED 2026-08-29.** D1, D2, WS3.4 scope, MAFFT/IQ-TREE install, SNV x BCCV and
the cross-OS scoping are all settled — see "Decisions of record" below. Implementation
proceeds; the task plan is at `docs/superpowers/plans/2026-08-29-rpi-fitness-validation.md`.

**Branch:** `04-rpi-fitness-validation`
**Date:** 2026-08-29
**Supersedes nothing.** Extends `2026-08-24-reassortment-plausibility-index-design.md`.

**Goal:** Make the Reassortment Plausibility Index work end-to-end and demonstrate that it
separates real reassortants from randomly shuffled fakes — with a reproducible pipeline on
both Windows and Linux.

---

## Global Constraints

- **Never change existing gate logic, verdicts, or output columns.** RPI stays a parallel
  column set. `combine_verdict` and every `gate*_` column are frozen.
- **Commits per workstream.**
- **If a step fails, log it and continue** — no aborting the run on one bad case.
- **No code until the plan AND the M_pack rescaling choice are approved.**
- Rescaled M_pack must be: **bounded in (0, 1]**, **monotonic** in ddG, **no exact zeros**,
  **sweepable at ±50%** in the sensitivity analysis.
- Reproducible: pinned MAFFT / IQ-TREE / ViennaRNA, seeded everything, hash check on a small
  fixture, identical across OSes.

---

## Decisions of record (2026-08-29)

| # | Decision | Resolution |
|---|---|---|
| **D1** | M_pack rescaling | **Calibrate first, then choose.** Run the WS3.1 calibration panel; if the ddG distribution is bimodal with a clean gap -> logistic (Option C); otherwise -> tau-scaled Boltzmann (Option A) with `tau = ddG_span / ln(100)` ~= 4.86. Min-max (Option B) rejected: cohort-dependence breaks reproducibility. |
| **D2** | RNAcofold window | **N = 50 nt**, swept 25 / 50 / 75. Whole-segment folding is both computationally infeasible and biologically wrong — the compatibility signal is terminal (panhandle), not distributed. |
| **D3** | WS3.4 (M_pol scored for real) | **In scope and required.** Compute M_pack and M_pol from the candidate's actual terminal windows against the **implicated parent/donor consensus** — never self-comparison. This is what makes the score sequence-dependent. |
| **D4** | MAFFT + IQ-TREE install | **Go.** Pinned versions. Gate 4 fixture runs in a **scratch dir only**, never touching existing outputs. |
| **D5** | SNV x BCCV | Add the cross spec; resolve the ICTV species name against the reference index. |
| **D6** | Cross-OS reproducibility | Byte-identity for MAFFT and RPI outputs; **topology + support-level agreement** for IQ-TREE. |
| **D7** | Hard-fake taxonomy | Split into **b1** (within-species shuffles) and **b2** (cross-species sympatric). See §4/WS4 — this changed the pass criterion. |

---

## 0. The finding that reorders this whole plan

**As the code stands today, RPI is mathematically incapable of separating a real reassortant
from a "hard fake" (same species, same place, same year). This is not a tuning problem. It is
an identifiability problem, and it makes M_pack the load-bearing element of the deliverable
rather than a nice-to-have.**

Here is the proof, read off the production call site.

`compute_rpi` (`reassortment_detector.py:1484`) computes:

```
rpi_generation = E_total * M_pack * M_pol * K
rpi_retention  = rpi_generation * F
```

Its inputs, and what each actually depends on:

| Input | Source at the call site (`reassortment_detector.py:2351`, `score_candidate`) | Varies with the *specific sequences*? |
|---|---|---|
| `implicated_species`, `isolates` -> **E** | `gate1.implicated_species` + the whole reference set | **No** — species-level only |
| `odd_segment` -> **K**, **F** | `gate1.discordant_segment or "M"` | **No** — one of S/M/L |
| `identity` -> stratum -> **K** | `gate0.best_kmer_identity` | Only via a **binary** threshold |
| `m_pack` | hardcoded `1.0` (`reassortment_detector.py:2456`) | **No** — constant |
| `m_pol` | `polymerase_compatibility(w, w)` — **same dict passed twice** | **No** — always exactly 1.0 |
| `structural_intact` | gate2a booleans | Boolean only |

`divergence_stratum` (`reassortment_detector.py:1124`) is a binary threshold at 0.80, and the
production `identity` runs ~0.35–0.46 in this dataset — so in practice **every row lands in
`distant` and the stratum does not vary at all.**

`polymerase_compatibility` is handed `terminal_windows` as both arguments, so its identity is
1.0 by construction on every candidate.

**Consequence.** For a fixed implicated-species set and a fixed odd segment, `rpi_generation`
takes exactly **one** value, and `rpi_retention` takes at most **two** (structural intact or
not). A hard fake is *defined* by matching the real reassortant on species, place and year —
which is precisely what fixes E. So the hard fake and the real reassortant receive **identical
RPI scores**, necessarily, today.

**Only sequence-level factors read the actual combination being scored.** E, K and F are all
species-level or segment-level priors. M_pack and M_pol are the only two terms in the model
that look at the nucleotides of *this* S+M+L assembly.

### What this does and does not mean (corrected 2026-08-29, per D7)

An earlier draft of this section treated "RPI gives a real reassortant and a same-species
shuffle the same score" as a defect to be engineered away. **That was a misframing, and it is
worth stating plainly because it would have set the wrong success criterion for the whole
validation.**

RPI measures **plausibility**, not **occurrence**. A within-species shuffle *is* plausible:
its parents co-circulate, co-infect, and have compatible termini, so a high RPI is the
*correct* answer for it. Asking RPI to score it below a real reassortant would be asking the
index to report something it does not measure — and the layer that separates "could happen"
from "did happen" is the phylogenetic detector, not the index. **Plausibility is not evidence
of occurrence.**

So the finding splits in two:

- **For within-species shuffles (b1):** identical-to-real scores are the *intended* behaviour.
  These are a **plausibility-ceiling control** — they tell us the top of the scale is where it
  should be. Gate 4's tree calls, not RPI, must separate them from real positives.
- **For cross-species sympatric shuffles (b2):** E is high by construction (same host order,
  country, sampling years), so E cannot discriminate — and with `m_pack = 1.0` and a
  self-compared `m_pol = 1.0`, **nothing else can either.** b2 is where the index must earn
  its keep, and today it provably cannot.

> **The operative conclusion is unchanged, only its target:** WS3.3 and WS3.4 are mandatory,
> because b2 discrimination is impossible without them. Rescaling alone fixes the annihilation
> problem but leaves both factors constant at 1.0.

This is why the plan puts "score M_pack against real parental sequences" and "score M_pol
against real parental termini" inside the M_pack workstream rather than treating the rescaling
as a standalone arithmetic change. **Turning M_pack and M_pol on is the deliverable; the
rescaling is one decision inside it.**

---

## 1. Decision D1 — M_pack rescaling (**needs your approval**)

### The problem, with real numbers

`boltzmann_compatibility` (`reassortment_detector.py:1265`) returns `exp(-ddG / RT)` with
`RT = 0.0019872 * 310.15 = 0.6163` kcal/mol.

Measured energies from the previous session: cognate panhandle dG ~= **-22.4**, non-cognate
dG ~= **-3.6 to 0.0**, so ddG ~= **18.8 to 22.4** kcal/mol.

| ddG (kcal/mol) | raw `exp(-ddG/RT)` |
|---|---|
| 0.0 | 1.0 |
| 3.6 | 2.9e-03 |
| 18.8 | **5.7e-14** |
| 22.4 | **1.6e-16** |

Because `rpi_generation` is a product, switching this on as written multiplies every
non-cognate row by ~1e-14 and annihilates the index. The other factors span far narrower
ranges — `host_decay` spans 100x (1.0 -> 0.01), `F` spans ~13x (1.0 -> 0.075), `K` spans ~8x.
**The Boltzmann weight is not wrong physics; it is simply not commensurate with the other
factors.** A 1e-14 equilibrium occupancy ratio and a 0.5 host-proximity weight are not
quantities that should be multiplied as peers.

### Option A — scale parameter tau replacing RT (**recommended**)

`M_pack = exp(-ddG / tau)`, with tau a declared free scale parameter, **not** a physical
temperature, calibrated so the observed ddG span maps onto a target dynamic range `D`:

```
tau = ddG_span / ln(D)
```

With `D = 100` (matching `host_decay`, the widest existing factor) and ddG_span = 22.4:
**tau ~= 4.864**.

| ddG | M_pack |
|---|---|
| 0.0 | 1.0000 |
| 3.6 | 0.4771 |
| 10.0 | 0.1280 |
| 18.8 | 0.0210 |
| 22.4 | **0.0100** — lands exactly on `eco_floor` |

**Why I recommend this:**

- One parameter. Keeps `boltzmann_compatibility` — an already-implemented, already-unit-tested
  function with a validated sign convention — intact; tau enters through its existing `rt`
  argument, so the diff is small and the existing tests stay meaningful.
- Bounded (0, 1], strictly monotonic decreasing, no exact zeros (the existing -700 exponent
  clamp already guarantees this), trivially sweepable at ±50% (tau in {2.43, 4.86, 7.30}).
- **Cohort-independent and therefore reproducible**: a candidate's score does not change when
  another candidate is added to the run. This is the property that kills Option B.
- tau is derived from a *pinned, versioned* calibration panel shipped in the repo, not picked
  by eye — so the number is auditable.
- Honest Methods sentence: *"M_pack is the relative Boltzmann weight of the non-cognate versus
  cognate panhandle duplex, rescaled by a free parameter tau so that its dynamic range is
  commensurate with the other index factors. tau is not a temperature and the result is not an
  occupancy probability; the index is explicitly relative."*

**The cost, stated plainly:** it forfeits the claim that M_pack is a *physical* Boltzmann
weight. That claim is arguably already unavailable — `compute_rpi`'s docstring says the index
"is RELATIVE. It is not a probability, and must never be reported as one" — so we would be
giving up something we do not actually have. But this is the one point a reviewer will press
on, and `RATIONALE.md` must meet it head-on rather than bury it.

### Option B — min-max normalisation over the run's ddG distribution

`M_pack = floor + (1 - floor) * (ddG_max - ddG) / (ddG_max - ddG_min)`.

**Rejected.** It is cohort-dependent: adding one candidate rescales every other candidate's
score. That directly violates the reproducibility requirement and would make two runs over
different candidate sets non-comparable — fatal for a fitness function meant to rank across
groups (a), (b) and (c). Listed here so the rejection is on the record.

### Option C — logistic on ddG (runner-up)

`M_pack = 1 / (1 + exp((ddG - ddG_50) / s))`, e.g. ddG_50 = 11.2 (midpoint of the observed
regimes), s = 3.0:

| ddG | M_pack |
|---|---|
| 0.0 | 0.9766 |
| 3.6 | 0.9264 |
| 10.0 | 0.5987 |
| 18.8 | 0.0736 |
| 22.4 | 0.0234 |

**Pros:** does not pretend to be a Boltzmann weight at all, so it dodges the physics objection
entirely; the sigmoid's saturation at both ends is arguably the right shape (a duplex much
worse than threshold is not meaningfully worse than one merely worse). Two interpretable
parameters. Bounded, monotonic, never exactly 0 or 1.

**Cons:** two parameters instead of one, both requiring calibration and both entering the ±50%
sweep (four extra sweep cells). It replaces rather than reuses the tested
`boltzmann_compatibility`. And the saturation cuts both ways — it compresses exactly the
cognate end of the range where we most need resolution to separate real from hard-fake.

**Option C is worth taking only if the calibration panel shows the ddG distribution is
genuinely bimodal with a clean gap** — in which case the sigmoid describes a real threshold
rather than imposing one. The calibration task (WS3.1) produces exactly that evidence.
**If you want the most defensible answer rather than the fastest one, that is the
recommendation: approve WS3.1 now, decide A vs C on the measured distribution.**

### D1 — what I need from you

Pick one:

1. **Option A, tau from D = 100** — proceed as specified. *(fastest defensible path)*
2. **Run WS3.1 calibration first, then choose A or C on the measured ddG distribution.**
   *(most defensible; costs one extra review cycle)*
3. Option C with stated parameters.
4. Something else.

---

## 2. Decision D2 — the RNAcofold window (**needs your approval, lower stakes**)

`rnacofold_dg`'s current docstring and `boltzmann_compatibility`'s clamp comment both refer to
**"whole-segment duplexes whose MFEs run to hundreds of kcal/mol."** That is not runnable.
RNAcofold is O(n^3) in time and O(n^2) in memory; hantavirus S ~= 1.8 kb, M ~= 3.6 kb,
L ~= 6.5 kb, so an L x M duplex is a ~10 kb fold — on the order of 1e12 operations and tens of
GB. It will not complete.

The measured energies quoted above (~= -22.4 kcal/mol) are consistent with **panhandle termini
only**, which is also the biologically correct window: hantavirus segment packaging is mediated
by the complementary 5'/3' termini, not the whole genome.

**Proposal:** fix the window explicitly and pin it as a versioned constant — the terminal
**N = 50 nt** of the 5' end of one segment against the terminal 50 nt of the 3' end of the
other, reverse-complemented as the panhandle requires. 50 nt comfortably contains the conserved
14-nt tag plus the lineage-variable region beyond it (`HANTAVIRUS_TERMINAL_TAG_NT = 14`,
`RPI_N_FOOTPRINT_NT` up to 23), folds in milliseconds, and gives dG magnitudes in the observed
range. **N enters the ±50% sweep** (25 / 50 / 75 nt).

This also requires fixing the stale docstring on `rnacofold_dg` and the clamp comment in
`boltzmann_compatibility`, both of which currently justify the -700 clamp by reference to
whole-segment folding. The clamp should stay (still correct, still cheap) but its stated
rationale needs to match reality.

**D2:** approve N = 50, or name a different window.

---

## 3. File structure

| File | Change | Responsibility |
|---|---|---|
| `reassortment_detector.py` | Modify | N4/N5/N6 fixes; tau/window constants into `RPI_WEIGHTS`; real M_pack + M_pol call sites in `score_candidate`; bump `RPI_WEIGHTS_VERSION` |
| `rpi_calibration.py` | **Create** | Builds the pinned ddG calibration panel; derives tau; writes `data/rpi_calibration.json`. Standalone so calibration is reproducible and auditable, and so `reassortment_detector.py` reads a pinned artifact rather than recomputing. |
| `rpi_discrimination.py` | **Create** | The group (a)/(b)/(c) discrimination harness, planted-chimera ladder, and the ±50% sweep driver. **Separate from `validate_reassortment_detector.py`**, which validates *verdicts*; this validates *scores*. Different question, different artifact. |
| `validate_reassortment_detector.py` | Modify | Add the missing SNV x BCCV cross to `PUBLISHED_CROSSES`; expose `PanelBuilder` for reuse |
| `requirements-pip.txt` | **Create** | The system-Python path, reconciled against `environment.yml` |
| `environment.yml` | Modify | Confirm/adjust pins; drop the "could not be reconciled" caveat once verified |
| `RATIONALE.md` | **Create** | Why each design decision was made — the defence document |
| `docs/methods-paragraph.md` | **Create** | Corrected Methods text |

Two new modules rather than growing `reassortment_detector.py`: it is already 4,444 lines with
the test suite inline. Calibration and discrimination are distinct concerns with distinct
artifacts, and keeping them out means the frozen scoring path stays small enough to review.

---

## 4. Workstreams

Ordered as specified. Each ends in its own commit.

### WS1 — Bug fixes (N4, N5, N6)

No dependencies. Can start immediately, before D1/D2 are decided.

- **N4 — nan ddG sign flip.** `boltzmann_compatibility(float('nan'), 0.0)` returns
  `9.86e-305`: `max(-700.0, nan)` yields `-700.0`, so a nan ddG silently becomes *maximal
  penalty*. Neither `1.0` (the pre-clamp behaviour) nor `9.86e-305` is a defensible answer for
  "the energy is undefined." **Fix: raise `ValueError` on non-finite ddG.** A nan here means
  RNAcofold returned garbage, and per the module's own no-silent-sentinel rule that must
  surface, not propagate. Per the "log it and continue" constraint, the *caller* in
  `score_candidate` catches it, records `rpi_m_pack_status = "undefined_energy"`, and falls
  back to the declared-constant path — the same honesty pattern `rpi_m_pol_status` already uses.
- **N5 — `_json_safe` np.float32 hole.** `np.float32('nan')` is not `isinstance(x, float)`, so
  it survives sanitisation and `default=str` serialises it as the JSON *string* `"nan"` before
  `allow_nan=False` can fire. **Fix: test `numbers.Real` / `np.floating` rather than `float`.**
  Not reachable today — but WS3 introduces numpy-typed energies into the report, which is
  exactly what makes it reachable.
- **N6 — thin ordering test.** `test_rounding_never_reports_retention_above_generation`
  (`reassortment_detector.py:3947`) passes `floor**5` to both columns at the call site, so it
  cannot fail via any change to `_rpi_decimals`. **Fix: parametrise the two precisions
  independently, and add the guard-the-guard assertion its N1 sibling has** — assert that at
  deliberately mismatched precisions the test *does* fail, so the guard is proven live.

**Regression tests for all three, asserted on the shipped CSV/JSON artifact, not on the
intermediate return** — per the repo's recurring-defect note in `CLAUDE.md`.

### WS2 — Reproducibility

Depends on nothing. Runs in parallel with WS1.

- Pin and **verify** MAFFT, IQ-TREE, ViennaRNA. ViennaRNA 2.7.2 is already installed from PyPI
  and matches the pin. MAFFT and IQ-TREE are not on PyPI — this needs conda, or manual Windows
  binaries. **This is the one workstream that cannot complete without your go-ahead on
  installing them**; everything else proceeds regardless.
- `requirements-pip.txt` reconciling the system-Python path (`dendropy`, `viennarna`) against
  `environment.yml`, with a header stating which interpreter each file describes.
- Seed audit: `random.Random(seed)` in `validate_reassortment_detector.py:396` is seeded;
  `--phylo-seed` defaults to 42. Add explicit seeding to any new RNG in WS3/WS4 and record the
  seed in every artifact.
- **Cross-OS hash check — with an honest scope limit.** MAFFT alignments (fixed options) and
  all RPI arithmetic (pure Python/numpy) are byte-reproducible. **IQ-TREE likelihood values are
  not guaranteed bit-identical across OS, compiler, or thread count**, even with a fixed seed.
  So the fixture check hashes the deterministic artifacts byte-for-byte, and compares IQ-TREE
  output on **topology + support values**, not bytes. Claiming byte-identity for IQ-TREE across
  OSes would be a claim we cannot keep.

### WS3 — M_pack (**gated on D1 and D2**)

- **WS3.1 — Calibration panel.** Build a pinned set of cognate and non-cognate panhandle pairs
  from the reference index, run RNAcofold over all of them, write `data/rpi_calibration.json`
  with the ddG distribution and derived tau. **This artifact is the evidence for the A-vs-C
  decision** if you choose D1 option 2.
- **WS3.2 — Rescaling.** Implement the approved form. tau (or ddG_50 / s) and the window N go
  into `RPI_WEIGHTS`; bump `RPI_WEIGHTS_VERSION`.
- **WS3.3 — Score M_pack for real.** Replace `m_pack=1.0` at `reassortment_detector.py:2456`
  with a real cognate-vs-non-cognate comparison against parental sequences sourced from
  `reference.genomic_index` per implicated species. Add `rpi_m_pack_status` mirroring
  `rpi_m_pol_status` (`scored` / `fixed_constant` / `undefined_energy` / `no_parental_sequence`).
- **WS3.4 — Score M_pol for real.** Same fix, same source: stop passing `terminal_windows`
  twice, source real parental termini, and delete the forced `m_pol_status="fixed_constant"`
  override so `m_pol_variance_check` can do its job. **Without this, M_pol stays a constant and
  contributes nothing to discrimination** — see §0.
- Fix the stale whole-segment docstrings on `rnacofold_dg` and `boltzmann_compatibility`, and
  remove the now-obsolete "ASSUMED FORMAT, NOT VALIDATED" caveat on `_parse_rnacofold_dg`
  (the regex was validated last session).

### WS4 — Validation (the core deliverable)

Depends on WS3.

**Four groups** (per D7), scored on `rpi_generation` and `rpi_retention`:

- **(a) Random fakes** — segments drawn from different species *and* different continents.
  Expect near-floor. E collapses (`host_decay` different_order = 0.01, geo floor = 0.01) *and*
  M_pack collapses. This group separates on E alone, so **it is the weak test** and must be
  reported as such rather than as a success.
- **(b1) Within-species shuffles** — segments from different isolates of the *same* species.
  **Expected to score HIGH, at or near the real positives.** This is a
  **plausibility-ceiling control**, not a negative: these assemblies genuinely are plausible,
  and RPI scoring them highly is correct behaviour. **RPI is not required to separate b1 from
  the real positives** — the tree detector is the layer that does that. What b1 tests is that
  the top of the RPI scale sits where it should.
- **(b2) Cross-species sympatric shuffles** — different species, but same host order, same
  country, overlapping sampling years. **E is high by construction, so M_pack and M_pol must
  do all of the discriminating.** Expected mid-to-low. **This is the group that actually tests
  the index**, and per §0 it is impossible without WS3.3 + WS3.4.
- **(c) Real positives** — the 8 supported candidates plus the published lab pairs.
  **Gap found:** `PUBLISHED_CROSSES` currently has SNV x ANDV (both orientations), PUUV x PHV,
  and the Rodriguez SNV NMR11 x CC107 intra-species pair — but **SNV x BCCV is not there**. It
  exists only as counts in `RPI_INVITRO_BENCHMARK`. Per D5, WS4 adds it as a cross spec, which
  requires resolving BCCV's ICTV species name against the reference index and confirming a
  complete triplet exists. If it does not, that is logged as a skip, not silently dropped.

### Pass criterion (two parts, both required)

1. **RPI ranking:** median `rpi_retention` ordered **reals ~= b1 > b2 > random fakes**.
   `reals ~= b1` is an *expected tie*, asserted as an equivalence (overlapping distributions),
   **not** as a separation. Reporting a real-over-b1 gap as a success would be claiming the
   index measures occurrence, which it does not.
2. **Detector separation:** the detector's **tree calls (gate 4) separate reals from b1.**
   This is the half that establishes occurrence, and it is the first real exercise of gate 4
   against live binaries (D4).

Gaps that must separate (b1 vs b2, b2 vs random) are reported with a Mann-Whitney U and an
effect size. Not a t-test — the scores are bounded products and nowhere near normal. The
expected tie (reals vs b1) is reported as an equivalence interval, since a non-significant
U-test is absence of evidence, not evidence of the tie we are claiming.

- **Planted chimeras at increasing parental divergence** to find the detection limit: sweep
  donor/backbone pairs across the observed identity range and report the divergence at which
  separation vanishes.
- **±50% sensitivity sweep on all weights**, extending the existing `sensitivity_sweep`
  (`reassortment_detector.py:1932`) to cover tau and N.

**Honest-failure clause:** if group (b2) does not separate from b1, the deliverable is the
measurement plus a stated limitation — *not* a retuned tau chosen to make it separate. Tuning
the scale parameter against the validation set would make the validation circular, and it is
the single easiest thing for a reviewer to catch.

### WS5 — Outputs

- Results tables (CSV) and plots (matplotlib; distribution per group, chimera ladder, sweep
  tornado).
- Corrected Methods paragraph in `docs/methods-paragraph.md`.
- `RATIONALE.md` — every design decision with its reasoning, explicitly including: why tau is
  not a temperature, why min-max was rejected, why the window is 50 nt, why IQ-TREE is compared
  on topology rather than bytes, and what group (b) actually measures.

---

## 5. Risks

| Risk | Mitigation |
|---|---|
| **Group (b2) does not separate from b1** — the predicted-from-source outcome if WS3.3/3.4 are skipped | WS3.3 and WS3.4 are both mandatory, not optional. If it still does not separate with both scored, report it. |
| **b1 read as a false positive** | b1 scoring high is *correct*. The pass criterion asserts `reals ~= b1` as an equivalence, and `RATIONALE.md` states why. Guarding against a future reader "fixing" this. |
| Gate 4 fixture run clobbers existing outputs | Per D4, the fixture runs in a scratch dir only. Asserted in the WS2 test, not just intended. |
| tau tuned against the validation set (circularity) | tau is derived in WS3.1 from a calibration panel **disjoint from the WS4 validation groups**, and pinned before WS4 runs. |
| RNAcofold energies are parameter-set dependent | ViennaRNA 2.7.2 pinned and now installed; the calibration artifact records the version. |
| Cross-OS bit-identity claimed too broadly | Scoped in WS2: bytes for MAFFT + RPI, topology/support for IQ-TREE. |

---

## 6. Status

All open questions resolved — see "Decisions of record". The bite-sized task plan is at
`docs/superpowers/plans/2026-08-29-rpi-fitness-validation.md`. Implementation starts at WS1.

One decision remains *deferred by design*, not unanswered: **D1's A-vs-C choice is made from
the WS3.1 calibration artifact**, on the bimodality test stated in the table. That is a
data-driven branch inside the plan, not a pending approval.
