# Factor E isolate-level fix -- implementation record (2026-08-30)

Implements `docs/superpowers/specs/2026-08-30-rpi-e-geo-isolate-level-fix-design.md`
(approved, with the `rpi_e_isolate_basis` addition folded into the note before this code
was written). Root cause and mechanism were already confirmed read-only in
`notes_diagnosis.md` -- this note is the fix + before/after record, not a re-diagnosis.

## What changed

- `ecological_opportunity` (reassortment_detector.py) now takes `parent_isolates: Dict[
  species, isolate]` and `isolate_basis: str` (both required, no default -- same
  "state your own status" rule `m_pack_status`/`m_pol_status` already follow) and
  computes host/geo/time from the ONE isolate that stands for each implicated species,
  never from that species' full country/host/year union across the reference.
- `nearest_reference_isolate` + `resolve_e_parent_isolates` (new): the fallback for
  candidates whose true parent isolate is unknown -- the k-mer-nearest wild reference
  isolate of each implicated species, on the segment gate1 top-matched it against. Basis
  is `actual_parents` / `nearest_reference` / `insufficient_metadata`, echoed onto every
  row as `rpi_e_isolate_basis` -- a row on the proxy is marked, never silent.
- `CandidateInput.known_segment_isolates` (new, optional, defaults to `None`): a caller
  that already knows which real isolate contributed each segment can say so.
  `rpi_discrimination.py`'s `_candidate_from_isolates` now sets it for every group it
  builds (a/b1/b2/c) -- group c's real candidates still correctly fall back to the proxy
  for the (unknown) donor side, because `resolve_e_parent_isolates` only credits
  `actual_parents` when EVERY implicated species is covered by the known map, and a
  genuinely discordant real candidate's donor species is exactly the thing that map
  does not contain.
- `score_candidate` now calls `resolve_e_parent_isolates` before `compute_rpi` and wires
  its output through. `compute_rpi` gained the two new required params and passes them
  straight to `ecological_opportunity`.
- `_score_panel` (sensitivity sweep) derives `parent_isolates` from its own synthetic
  panel's `isolates` list (already one real isolate per species) -- basis
  `actual_parents` throughout, unaffected in kind by this fix.

No `RPI_WEIGHTS` entry changed; `RPI_WEIGHTS_VERSION` unchanged. This is a bug fix, not
a retune -- see the design doc's honesty protocol section.

## Regression tests added (reassortment_detector.py)

- `EcologicalOpportunityTest`: rewritten to pass `parent_isolates`/`isolate_basis`
  explicitly; added `test_same_species_country_pool_but_different_parent_isolates_
  still_separates`, which puts two species in overlapping country POOLS but scores a
  specific parent pair that never co-occurs, and a same-pool pair that does -- the two
  return opposite `rpi_e_geo` answers from the same wider isolate list.
- `NearestReferenceIsolateTest`: correctness of the k-mer-nearest fallback (picks the
  closer protein; `None` on empty candidate protein, empty pool, or no isolate carrying
  a protein for the requested segment).
- `ResolveEParentIsolatesTest`: `actual_parents` when every species is covered by
  `known_isolates`; falls back to the proxy for the WHOLE row (never a silent per-species
  mix) when one species is missing from the known map, matching group c's real
  candidates; `insufficient_metadata` when even the fallback can't resolve a species; and
  `test_proxy_basis_is_flagged_on_the_shipped_row_not_silent`, which asserts
  `rpi_e_isolate_basis` on `compute_rpi`'s RETURNED row (the row that ships, not an
  intermediate) per the project's own "assert on the row that ships" rule.
- `rpi_discrimination.py`: `CandidateFromIsolatesTest` -- `_candidate_from_isolates` sets
  `known_segment_isolates` to exactly the sources it was given.

Full suite after the change: 237 tests (`reassortment_detector`, was 228), all pass, 815s
wall (gate 4 dominates, unaffected by this fix). `rpi_discrimination`: 18 tests, all
pass. `rpi_calibration`: 39 tests, all pass, unaffected (Factor E is not part of the
calibration panel).

## Before/after: the seeded discrimination panel (seed 20260829)

Before = `output/discrimination_scratch/discrimination_panel.json` (preserved, untouched,
same code that produced the numbers in PROJECT-STATE.md section 7c). After =
`output/discrimination_scratch/discrimination_panel_after_fix.json` (this fix, same
seed, same reference data).

### `rpi_generation` summary

| group | n | mean before | mean after | min before | min after | max before | max after |
|---|---|---|---|---|---|---|---|
| a_random_fakes | 20 | 0.3314 | 0.000138 | 0.000295 | 1.31e-06 | 0.4717 | 0.001179 |
| b1_within_species_shuffle | 30 | 0.3238 | 0.2885 | 0.0908 | 0.0146 | 0.4717 | 0.4717 |
| b2_cross_species_sympatric | 36 | 0.2519 | 0.02145 | 0.01176 | 0.000118 | 0.4717 | 0.1179 |
| c_real_positive | 9 scored / 11 | 0.1947 | 0.3353 | 0.00349 | 0.02514 | 0.4717 | 0.4717 |

`a_random_fakes` collapses almost three orders of magnitude, as the design doc
predicted: `rpi_e_total` was `{1.0: 14, 0.01: 5, 0.000625: 1}` before (14/20 sitting at
the ungated ceiling because two widely-distributed/border-spanning species' full country
sets overlapped somewhere in the corpus) and is now suppressed on every single row (15
distinct small values, all well under 0.01) because the fix reads the two ACTUAL parent
isolates' countries, which `build_group_a_random` already guarantees are on different
continents.

`b2_cross_species_sympatric` also drops (mean 0.2519 -> 0.02145). This was NOT predicted
to be geo-suppressed -- b2 is genuinely sympatric by construction -- but M_pack/M_pol's
own (already-documented, weak-but-real) discrimination is now visible without random's
inflated ceiling drowning it out in the comparison; b2's own median (0.00398) sits well
above random's (3.07e-05), which is the point.

`b1_within_species_shuffle` moved more than "likely unaffected, worth confirming, not
assuming" (design doc section 3) turned out to mean: mean 0.3238 -> 0.2885, not
identical. Basis breakdown explains it: 27/30 b1 rows are `actual_parents` (gate1
correctly reads a single implicated species, E stays neutral exactly as before), but
3/30 resolved via `nearest_reference` -- rows where gate1's own per-segment k-mer match
implicated a SECOND, sibling species even though every segment genuinely came from the
same true species (a pre-existing, documented gate1 blind spot, not something this fix
changes). For those 3 rows only, E is no longer trivially neutral, and the fallback
isolate for the phantom second species can be ecologically separated from the real one.
This is a real, isolate-level-correct consequence, not a bug in this fix.

`c_real_positive` rises (mean 0.1947 -> 0.3353). Basis breakdown: 5/9 scored rows are
`actual_parents` (gate1 saw a single, concordant species -- the isolate's own -- so E was
always neutral and unaffected), 4/9 are `nearest_reference` (the genuinely discordant
real candidates, whose donor species is unknown by construction).

**Correction to an assumption checked and found wrong before writing this note down:**
the 4 `nearest_reference`-basis rows do NOT all stay near the floor. Diagnosed directly
(not assumed) by reproducing `resolve_e_parent_isolates` against the real reference data
for all four published crosses:

| case | before | after | resolved isolates' geo/host/year |
|---|---|---|---|
| `RIZVANOV2004_SNVsl_x_ANDVm` | 0.0035 | 0.3494 | ANDV proxy isolate `t2\|...andesense\|chile 9717869` has EMPTY country/host_order/host_family/host_name/year -- every field blank |
| `RIZVANOV2004_reciprocal_ANDVsl_x_SNVm` | 0.0046 | 0.4569 | same ANDV proxy (empty), plus the SNV proxy this direction resolved to also has empty country/year (only host fields populated, from a human clinical record) |
| `HANDKE2010_PUUVsl_x_PHVm` | 0.0047 | 0.4717 | PHV proxy isolate `acc\|NC_038940.1` has EMPTY country/host_order/host_family/host_name/year -- every field blank |
| `RODRIGUEZ1998_SNV_NMR11sl_x_BCCVm` | 0.3669 | 0.3669 | unchanged -- both proxies resolve to `country="USA"` (a real, non-missing value); already `nearest_reference`-basis and near-neutral before this specific rerun's baseline too |

Three of the four rows went from floor to mid/high not because the isolate-level fix
found real ecological overlap, but because `ecological_opportunity`'s pre-existing
"missing metadata means neutral, not separated" rule (deliberate, and correct on its
own terms -- see `test_species_absent_from_reference_entirely_stays_neutral` and
`test_unknown_host_metadata_stays_neutral_not_maximal_separation`) fired on a proxy
isolate that `nearest_reference_isolate` picked by PROTEIN SIMILARITY ALONE, with no
preference for metadata completeness. `t2|...andesense|chile 9717869` and
`acc|NC_038940.1` both happen to be sparse GenBank deposits with no country/host/year
recorded at all -- the nearest isolate by k-mer distance, not the best-annotated one.

This is a REAL, newly-surfaced limitation of the `nearest_reference` fallback, distinct
from the fix's own correctness: `rpi_e_isolate_basis="nearest_reference"` on these rows
IS an honest, accurate flag (the row is not on real known parents), so a reader checking
that column before trusting the row is not misled -- the honesty mechanism the design
doc asked for is working exactly as intended. But the PROXY QUALITY itself has a gap:
nothing in `nearest_reference_isolate` prefers isolates with usable ecological metadata
over ones without it. Logged as a known limitation below and in PROJECT-STATE.md;
per the design doc's own honesty protocol and the standing "do not chase a compensating
change" instruction, this is NOT fixed here.

The fourth row (RODRIGUEZ1998 SNV x BCCV) is different in kind: its resolved proxies
both carry a REAL, non-missing `country="USA"` -- not a metadata gap. At country-level
granularity (this codebase's existing geo resolution, unchanged by this fix) SNV
(continent-wide US) and BCCV (Louisiana, US) are simply not distinguishable, which is a
pre-existing coarseness of the geo term, not something this fix introduced or could
have corrected.

### Task 14 pass criterion (`rpi_discrimination_stats.evaluate_pass_criterion`)

| comparison | before | after |
|---|---|---|
| `b2_vs_random` | p=0.153, delta=-0.225 (small, WRONG SIGN), met_expectation=False | p=4.68e-09, delta=+0.953 (large, correct sign), met_expectation=**True** |
| `b1_vs_b2` | p=0.052, delta=+0.275 (small, borderline), met_expectation=False | p=2.18e-11, delta=+0.961 (large), met_expectation=**True** |
| `reals_vs_b1` | p=0.048, delta=-0.430 (medium, WRONG DIRECTION for a claimed tie), significant | p=0.733, delta=+0.078 (negligible), **not significant** -- a plausible tie, as D7 predicts |
| `overall_ranking_supported` | **False** | **True** |

This is reported as the honest consequence of a bug fix, not a tuning result -- per the
design doc's honesty protocol, nothing in `RPI_WEIGHTS` or the group-construction code
changed between these two runs. `b2_vs_random`'s wrong sign (notes_diagnosis.md section
1) is now resolved by the mechanism the fix targeted: `a_random_fakes`'s 20/20 rows are
all `actual_parents` basis, and the isolate-level geo term now correctly reads the real,
different-continent parent pair each candidate was built from. That part of the result
is exactly what the design doc predicted, confirmed by the basis column, not assumed.

`reals_vs_b1`'s prior significant-wrong-direction result is now a negligible-effect tie,
but **this part of the result is NOT a clean confirmation of the fix's own mechanism.**
Three of `c_real_positive`'s four `nearest_reference`-basis rows (the RIZVANOV x2 and
HANDKE published crosses) moved from the floor to mid/high because the fallback's proxy
isolate happened to have no geo/host/year metadata at all (see the table above), not
because isolate-level geography found real ecological overlap -- a genuine, separately-
logged limitation of `nearest_reference_isolate`, not evidence that these three
candidates are ecologically plausible. Pulling those three rows out of `c_real_positive`
would very likely restore something closer to the OLD result (reals scoring low) for the
remaining rows, since the mechanism that raised them is a metadata gap specific to those
three isolates, not a property of "reals" as a class. This was NOT re-run with those
three rows excluded -- doing so would be constructing a third, cherry-picked group
after seeing which one produces the preferred pass/fail outcome, which the design doc's
honesty protocol and the standing "do not chase a compensating result" instruction both
rule out. The number reported above (`overall_ranking_supported: True`) is the actual,
unedited output of the actual, unedited fix and group construction -- report it as
measured, with this caveat attached, not as a clean validation of D7's plausibility-tie
claim.

## New known limitations (logged in PROJECT-STATE.md)

1. **E vs gate3 disagreement.** `ecological_opportunity` (Factor E) now reads
   isolate-level geography for the actual parent pair; `gate3_check` (verdict-bearing,
   deliberately untouched -- see design doc section 2) still reads species-wide
   geography for only the first two sorted implicated species. The two can now disagree
   in a new, specific way beyond the pre-existing "first two of N" limitation. Both are
   already-documented, deferred to the same Factor E rebuild -- this is one more
   concrete instance, not a new problem.

2. **`nearest_reference_isolate` picks by protein similarity only, with no preference
   for metadata completeness.** Measured directly (see the table above): 3 of 4
   `c_real_positive` rows using the `nearest_reference` fallback resolved to a proxy
   isolate with EMPTY country/host_order/host_family/host_name/year -- a sparse GenBank
   deposit that happened to be the closest match by k-mer distance -- so
   `ecological_opportunity`'s pre-existing (deliberate, correct-on-its-own-terms)
   "missing metadata means neutral" rule fired, and E went to 1.0 by omission rather
   than by evidence of real ecological overlap. The `rpi_e_isolate_basis` column
   correctly flags these rows as proxy-based, so no row is silently misrepresented --
   but the proxy's quality itself has a gap a future change could close (e.g. prefer the
   nearest isolate WITH usable metadata over the nearest isolate period, or fall back to
   `insufficient_metadata` rather than a metadata-empty proxy). Not fixed here --
   flagged for whoever picks up the Factor E rebuild next, per the "log it and continue"
   standing rule (PROJECT-STATE.md section 8) and the instruction not to chase a
   compensating change against this task's own measured result.
