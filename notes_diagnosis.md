# WS4 Task 13/14 — read-only diagnosis (2026-08-30)

No code edits, no parameter changes. All numbers pulled from the existing
`output/discrimination_scratch/discrimination_panel.json` (seed 20260829, unchanged
code) plus one read-only recomputation of `ecological_opportunity` against the loaded
reference data, run interactively (not committed as a script).

## 1. `b2_vs_random` wrong sign — root cause identified

**Which E sub-factor fails to suppress `a_random_fakes`: `rpi_e_geo` (geography).**

`build_group_a_random` (rpi_discrimination.py) requires the three chosen S/M/L
**isolates** to span >=2 continents via `COUNTRY_CONTINENT` on each isolate's own
`country`. That guarantee is real at construction time. But the score itself does not
re-check geography at the isolate level — `ecological_opportunity` (reassortment_detector.py:1258)
computes `geo` from the **union of countries across every isolate of the whole species**
in the reference dataset (`by_species[species]`), not from the two specific parent
isolates chosen for this fake. `rpi_e_total = rpi_e_host * rpi_e_geo * rpi_e_time`, and
`rpi_e_geo = 1.0` whenever `countries_a & countries_b` is non-empty at the species level.

Measured distribution, `a_random_fakes` (n=20):

    rpi_e_total: 1.0 x14, 0.01 x5, 0.000625 x1

The 14 rows at `E=1.0` are every `['Asia','Europe']`-labeled case built from the trio
`{hantanense, puumalaense, seoulense}` (or similar closely related Palearctic species).
Recomputing `ecological_opportunity(['hantanense','puumalaense','seoulense'], all_isolates)`
directly returns `{host: 1.0, geo: 1.0, time: 1.0, total: 1.0}` — confirmed by pulling each
species' full country set from the reference:

    hantanense:  China, India, Russia, South Korea, Viet Nam
    puumalaense: Austria, Belarus, ..., China, Russia, South Korea, ... (26 countries)
    seoulense:   Argentina, Belgium, ..., China, Russia, South Korea, ..., USA (22 countries)

`puumalaense` (nominally European) and `seoulense` (globally invasive via rat trade) both
have real reference isolates from China/Russia/South Korea, so every pairwise
species-level country-set intersection is non-empty even though the two ACTUAL parent
isolates used to build a given fake candidate were sampled continents apart. `host` and
`time` are also 1.0 for this trio (same host order, broad multi-decade year spans), so
`geo` isn't the only factor at 1.0, but it is the one whose failure mode is specific and
diagnosable: it is checking "does this species ever appear near that species anywhere in
the corpus," not "did these two specific parent isolates ever co-occur," which silently
defeats the "different continents" intent the group-a builder encodes at construction
time. The 6 rows that DO get suppressed (`E=0.01` / `0.000625`) are pairs whose species
have genuinely disjoint reference country sets (e.g. Asia/North America, Asia/South
America, Africa/Asia) — no such widely-distributed or border-spanning species is
involved.

This is a species-level-prior limitation of the same shape already documented in
`docs/PROJECT-STATE.md` §3 ("E, K and F...do not read the nucleotides of the particular
assembly being scored") — extended here to say E also does not read the *isolate-level
geography* of the particular assembly being scored, only species presence/absence.

**Separately checked and ruled out:** `undetectable_by_design` vs `m_only`/
`s_or_l_involving` rows *within* `b2_cross_species_sympatric` do NOT differ meaningfully.

    b2, undetectable_by_design (n=14):        rpi_generation median 0.2270, m_pack median 0.884, m_pol median 1.0
    b2, m_only/s_or_l_involving (n=22):        rpi_generation median 0.2549, m_pack median 0.833, m_pol median 1.0

Gate1's confidence qualifier does not depress or inflate scores within b2 — this rules
out "the undetectable_by_design blind spot is itself dragging b2's median down" as a
contributor to the wrong-sign result. The wrong sign is fully explained by (1) above:
`a_random_fakes`'s construction-time continent guarantee is defeated by scoring-time
species-level geography, so ~70% of random fakes score at the ungated ceiling
(`rpi_generation = 0.4716981132`, i.e. `E=1.0` times whatever K/M_pack/M_pol give), while
`b2_cross_species_sympatric` is ALL `E=1.0` by construction (36/36) but has real
M_pack/M_pol competition depressing its median below random's inflated one.

## 2. `msb89863` PARTIAL-record check

Confirmed via the isolate's own metadata in the loaded reference data (not re-derived
from length heuristics this time):

    s_completeness = partial       m_completeness = partial       l_completeness = partial
    triplet_status = all_segments_present_not_all_complete
    has_complete_triplet = False
    segments_present = ['S', 'M', 'L']
    coding_complete_segments = []
    warnings = ['species_unclassified']

NCBI's own per-segment completeness annotation says `partial` on all three segments —
this is authoritative, not inferred from `SEGMENT_LENGTH_RANGE` alone (S=330 nt vs the
gate0 diagnostic's already-noted 1200-2100 nt expectation, M=1390 vs 3400-3800, L=490 vs
6300-6800, all previously logged in PROJECT-STATE.md's gate0 diagnostic). This corroborates
rather than contradicts the earlier finding: `msb89863` is a genuinely partial GenBank
deposit by NCBI's own record, and gate0's `junk_sequence` rejection of it in
`c_real_positive` is correct, not a bug. No code change indicated.

## Bottom line

Neither finding calls for retuning `RPI_WEIGHTS` or changing group construction code.
The `b2_vs_random` wrong sign has a confirmed, specific mechanism (species-level geo
lookup in `ecological_opportunity` vs isolate-level continent guarantee in
`build_group_a_random`) but per the standing instruction this diagnosis does not
recommend or make a fix — that decision belongs to whoever picks up Task 13/14 next.
