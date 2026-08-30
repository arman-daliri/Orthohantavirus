# Design note: fix `rpi_e_geo` to compare the two parent isolates, not species-level unions

**Written 2026-08-30. No code changed yet — this is the plan, for review before
implementation.** Prerequisite reading: `docs/PROJECT-STATE.md` §7c and
`notes_diagnosis.md` (this session), which found and confirmed the root cause below.

## 1. The exact change

`ecological_opportunity` (reassortment_detector.py:1258) currently computes the
geography term like this, per implicated-species pair `(species_a, species_b)`:

```python
members_a, members_b = by_species.get(species_a), by_species.get(species_b)
...
countries_a = {i.country for i in members_a if i.country}
countries_b = {i.country for i in members_b if i.country}
...
elif countries_a & countries_b:
    geo = 1.0
```

`by_species[species]` is every non-lab, non-vaccine isolate of that **species** present
anywhere in the reference dataset — not the two isolates that actually contributed the
candidate's segments. `geo` is 1.0 whenever the two species' full country sets intersect
anywhere in the corpus, regardless of where the specific parent isolates being scored
were collected.

**The fix:** `geo` (and, by the same reasoning, `host` and `time`) must be computed from
the **specific parent isolates implicated in this candidate** — the actual S/M/L source
isolates identified for the two sides of the discordant split — not from the species-wide
`by_species` union. Concretely, `ecological_opportunity` needs the caller to pass (or
needs to derive from `case`/`result` context already available at the `score_candidate`
call site) the **two actual isolate records** whose segments are in tension, and compare
`isolate_a.country` vs `isolate_b.country` (and `host_order`, and year) directly, instead
of first fanning out to "every isolate this species has ever had sequenced." The
worst-pair-over-all-implicated-species-pairs structure (`worst = min(...)` across pairs)
stays — only the *inputs* to each pairwise host/geo/time computation change, from
"every isolate of species A vs every isolate of species B" to "isolate A vs isolate B,
the ones actually in this candidate."

This is a plausibility-not-occurrence-respecting fix, not a stricter gate: E still asks
"could these two lineages have met," it just asks it about the two lineages that are
actually present in the candidate being scored, which is what Factor E's own docstring
already claims to measure ("relative opportunity for the implicated lineages to
co-infect one cell") — `implicated_species` is the species names, but nothing in the
current signature carries which *isolate* stood for that species in this particular
candidate, and that is the gap.

## 2. Does `gate3_check` have the same pattern? Is it covered by this fix?

**Yes, `gate3_check` (reassortment_detector.py:970) has the exact same species-level
union pattern** — `by_species` built the same way, `countries[a] & countries[b]` on the
full per-species country sets, not the two parent isolates.

**No, this fix does not touch it, and should not.** Two independent reasons, both
already on record:

- `gate3_check` is verdict-bearing (`docs/PROJECT-STATE.md` §8: "Never change gate
  logic, verdicts, or `gate*_` columns... RPI is a parallel column set"). This fix is
  scoped to `ecological_opportunity`, which feeds `rpi_e_*` columns only.
- `gate3_check`'s own docstring already flags this exact class of problem (there, framed
  as "only the first two of N implicated species," a related but distinct symptom of the
  same underlying species-level-not-isolate-level design) and explicitly defers the fix:
  *"The right answer is the worst-pair-over-all-pairs treatment that
  `ecological_opportunity` already implements for Factor E, and moving gate 3 onto it
  changes a verdict-bearing gate — out of scope for the RPI work. Deferred to the Factor
  E rebuild."*

So after this fix, `ecological_opportunity` and `gate3_check` will disagree in a new,
specific way: `ecological_opportunity` will read isolate-level geography for the pair
actually in the candidate, while `gate3_check` will still read species-level geography
for (only) the first two sorted species. Per `ecological_opportunity`'s own docstring
("Disagreement between the two is expected and informative"), this is consistent with
the model's existing design, not a new inconsistency introduced by this change — it is
one more concrete instance of the already-deferred Factor E rebuild, to be picked up
together with gate 3's "first two of N" limitation, not fixed piecemeal here.

## 3. Expected consequence

- **`a_random_fakes`'s `rpi_e_total` distribution collapses toward the floor.** Measured
  today: 14/20 rows at `E=1.0` (species-level country overlap defeats the
  construction-time "≥2 continents" guarantee for widely-distributed or
  border-spanning species like Seoul virus or the China/Russia-overlapping
  Hantaan/Puumala/Seoul trio), 6/20 already correctly suppressed. After the fix, all 20
  should suppress toward `eco_floor` (0.01) on the geo term alone, since by
  construction every `a_random_fakes` candidate's two contributing isolates are on
  different continents (`build_group_a_random` already enforces this at the isolate
  level — the fix makes `ecological_opportunity` finally *see* that same guarantee).
- **`b2_vs_random`'s sign should correct.** `b2_cross_species_sympatric` is
  ecologically sympatric by construction (`E=1.0`, 36/36, unaffected by this fix — its
  isolates genuinely do share country/host/time). With `a_random_fakes`'s `E` correctly
  suppressed instead of sitting at the ceiling for 70% of the group, `random`'s
  `rpi_generation` should drop well below `b2`'s, which is the expected direction
  (`b2 > random`) per Task 13's D7 design.
- **Every downstream panel number changes, not just `a_random_fakes`.** `E` also feeds
  `b1_within_species_shuffle` (single-species, `len(implicated_species) < 2` in the
  common case, so likely unaffected — worth confirming, not assuming) and
  `c_real_positive` (the published inter-species lab crosses currently score
  `E≈0.01`, already at the floor via the host/time terms per PROJECT-STATE §7c's
  E-stratification finding — these may or may not move depending on whether geo was
  already the binding term for them or host/time already floored them independently).
  **The `reals_vs_b1` and `b1_vs_b2` statistics must be re-run and re-read after this
  fix**, not assumed unaffected — Task 14's whole pass criterion is computed over
  `rpi_generation`, which is `E × M_pack × M_pol × K`, and `E` is changing.
- `tau`, `M_pack`, `M_pol`, and `K` are untouched by this fix — nothing about the
  calibration panel (`data/rpi_calibration.json`) or `RPI_WEIGHTS` changes, so no
  recalibration step is implied.

## 4. Honesty protocol

- **The before-panel output is preserved, not overwritten, for a mandatory
  before/after diff.** `output/discrimination_scratch/discrimination_panel.json` and
  `output/pipeline_acceptance/raw_scored_groups.json` (both seed 20260829, current
  code) stay in place untouched when this fix lands. The fix's own commit (or the very
  next one) must include a fresh seeded run's output alongside them, and the commit
  message and/or a diagnostic note must show the actual before/after numbers for
  `a_random_fakes`'s `rpi_e_total` distribution and the three `evaluate_pass_criterion`
  comparisons (`b2_vs_random`, `b1_vs_b2`, `reals_vs_b1`) — not just an assertion that
  they changed, the actual old-vs-new values side by side.
- **This is framed as a bug fix, never as tuning.** `ecological_opportunity`'s own
  docstring already states the intended semantics ("relative opportunity for the
  implicated LINEAGES to co-infect one cell") and `build_group_a_random`'s docstring
  already states its intended construction ("different species AND continents,"
  enforced at the isolate level). The current code fails to implement either
  specification correctly — it substitutes "has this species ever been seen near that
  species anywhere in the whole reference corpus" for "could these two specific
  lineages have met." Fixing that substitution makes the implementation match the
  already-documented design; it does not change what E is supposed to measure, adjust
  any weight or threshold, or respond to Task 13/14's measured result by moving a
  target. Per the standing instruction in PROJECT-STATE.md ("do not retune `m_pack_tau`,
  `RPI_N_FOOTPRINT_NT`, or any `RPI_WEIGHTS` entry against this result"), this fix
  changes zero entries in `RPI_WEIGHTS` and zero calibration parameters — it corrects
  which isolates a pre-existing formula is applied to.
- **The commit message for the eventual implementation must say both things
  explicitly**: (a) which specific lines changed and why, framed as "implementation
  did not match documented behavior," and (b) that the panel numbers moved as a
  *consequence* of that correction, with the before/after values inline or linked to
  the diagnostic note that carries them — never phrased as "improved the score" or "made
  the ranking pass," since neither framing would be honest about what changed and why.
- **If the fix does not make `overall_ranking_supported` become `True`, that is not a
  failure of the fix and must not trigger a second, compensating change.** The fix's
  job is to make `ecological_opportunity` measure what it already claims to measure;
  whether that happens to move Task 14's pass criterion over the line is a separate,
  honestly-reported fact, not a success condition for this change.
