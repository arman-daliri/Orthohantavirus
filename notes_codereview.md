# Code review — branch `04-rpi-fitness-validation` vs `main` (2026-08-30, evening)

`/code-review medium` over the full branch diff, with focused scrutiny on the three
commits landed this session: `98bc10a` (wire `known_segment_isolates` through
`build_swap`), `6d85d5b` (seeded panel rerun), `3762c91` (efficiency fix in
`resolve_e_parent_isolates`). Independently confirmed by the reviewing agent: the panel
rerun is byte-for-byte identical in `rpi_generation` across all 97 rows pre/post
`build_swap` fix, matching `notes_e_isolate_fix.md`'s own claim.

Per the evening task list's binding rule, autonomous fixes are limited to **provably
score-neutral** changes (dead code, docstrings, pure efficiency with an identical-output
guard). Anything that could move a number, or touches output-column semantics, gate
logic, verdicts, or group construction, is written up here and left for the user.

## Findings

### 1. `resolve_e_parent_isolates`'s known-isolate check may mislabel intra-species published crosses (`reassortment_detector.py:1437`) — PROPOSAL, downgraded after verification

Since `98bc10a`, `build_swap` always populates `known_segment_isolates`. For the
`B_published_intra_species` entries in `PUBLISHED_CROSSES`
(`validate_reassortment_detector.py`) — backbone and donor are the SAME species,
different isolates — gate1 cannot see intra-species reassortment (a documented,
pre-existing blind spot), so `implicated_species` collapses to that one species.
`resolve_e_parent_isolates`'s `all(species in by_species for species in
implicated_species)` check is then true over a single-element list, so these rows
WOULD get `rpi_e_isolate_basis="actual_parents"` if scored.

**Verified against the codebase (this is not currently reachable, contrary to the
reviewer's original "why it matters" claim):** the reviewer's finding said this
"corrupt[s] exactly the 'stratify real positives by Factor E basis' diagnostic." Checked
both places `score_candidate` is actually called on `PUBLISHED_CROSSES` rows:
- `rpi_discrimination.py`'s `load_published_inter_species_crosses` — the function that
  feeds Task 14's group-c stats and the basis-stratification diagnostic — filters to
  `spec["klass"] != "A_published_inter_species": continue` (line 440-441).
  `B_published_intra_species` rows never reach `resolve_e_parent_isolates` through this
  path at all; the module's own docstring on this function states this exclusion is
  deliberate ("Class B (intra-species) is EXCLUDED here").
- `validate_reassortment_detector.py`'s `score_case` (line 565-601), which DOES call
  `score_candidate` on every `PUBLISHED_CROSSES` entry including class B, builds its
  output row from an explicit fixed key list that does not include
  `rpi_e_isolate_basis` at all (only `rpi_generation`/`rpi_retention`/`rpi_e_total`/
  `rpi_k_stratum` survive into the report).

So the vacuous-truth branch is real code, and would fire if a caller ever built a
single-species candidate with `known_segment_isolates` set and asked for its basis
label directly — but no current caller does that in a way that reaches any diagnostic
or report. `rpi_e_total` is unaffected either way (stays neutral 1.0, per
`ecological_opportunity`'s own early return). Keeping this as a proposal rather than
autonomously touching it (output-column semantics is still not mine to decide), but
correcting the severity: this is a latent code-shape observation, not a current data-
quality bug in anything the pipeline reports today.

**Severity:** low (was reported as medium; downgraded — no reachable path to a
misleading output currently exists, verified directly rather than assumed).

### 2. Three of four new efficiency-fix tests assert on the intermediate, not the shipped row (`reassortment_detector.py:4563,4575,4586`) — PROPOSAL, not fixed

`test_single_species_that_resolves_gets_nearest_reference_basis_without_kmer_work`,
`test_single_species_with_no_usable_protein_gets_insufficient_metadata`, and
`test_two_species_still_pays_for_full_resolution` all call `resolve_e_parent_isolates`
directly and assert on `resolved`/`basis`, never on `compute_rpi`'s returned row — the
exact defect class this project's own `CLAUDE.md` names ("assert on the row that ships,
not on the function that feeds it"). Only the fourth test in that class
(`test_ecological_opportunity_row_identical_whether_or_not_identity_was_resolved`)
checks a shipped row. If a future change makes `compute_rpi`/`ecological_opportunity`
read `resolved` differently than these three tests assume, none of them would catch it.

Not autonomously fixed: strengthening test assertions is a genuine improvement but
isn't one of the three pre-approved score-neutral categories (dead code / docstrings /
pure efficiency with a guard) — it's a test-coverage change, left for the user's call on
priority.

**Severity:** low (test-quality gap; the fourth test in the same class already covers
the row-that-ships guarantee for the one case that matters for the skip's own safety).

### 3. `skip_identity_resolution`'s threshold duplicates `ecological_opportunity`'s own, as an independent literal (`reassortment_detector.py:1451` vs `1296`) — PROPOSAL, not fixed

`len(implicated_species) < 2` is now written in two separate functions with nothing but
a cross-referencing comment keeping them in sync. If `ecological_opportunity`'s neutral-
return condition ever changes without a matching update here, the skip would silently
diverge — skipping isolate resolution in cases where Factor E now actually needs
`parent_isolates`, with no compiler/type check to catch the drift.

Not autonomously fixed: any restructuring to share the threshold (e.g. a shared
constant or a query function) is a refactor of production logic, which carries
functional risk that can't be "proven" neutral without redoing the same verification
rigor as the efficiency fix itself — left for the user.

**Severity:** low (latent maintenance risk, not a present bug).

### 4. Comment/docstring for the efficiency fix duplicated near-verbatim across two spots — FIXED

The rationale ("`ecological_opportunity` discards `parent_isolates` when fewer than two
species implicated, so resolving WHICH isolate is wasted") was written out in full both
in `resolve_e_parent_isolates`'s inline comment and in
`ResolveEParentIsolatesSkipsExpensiveResolutionTest`'s class docstring. Consolidated:
the inline comment stays canonical; the test class docstring now references it instead
of repeating it. Pure comment text, zero behavior change — `ResolveEParentIsolatesSkipsExpensiveResolutionTest`'s 4 tests re-verified green after the edit.

**Severity:** trivial (readability only). **Disposition: fixed autonomously** (docstring
category, provably score-neutral).

### 5. `_isolates_with_usable_protein` + a three-way branch is a heavier shape than a single-flag API on `nearest_reference_isolate` — PROPOSAL, not fixed

Suggested alternative: give `nearest_reference_isolate` an `identity_needed: bool` flag
that internally short-circuits before the k-mer/cosine work, rather than exposing a
second helper function and branching in the caller. Would reduce the number of
contracts a reader has to hold in mind and avoid a future third caller having to choose
between calling the private helper directly or re-duplicating the caller's ternary.

Not autonomously fixed: a genuine simplification opinion, but restructuring the
function signature is a production-code refactor with the same "can't prove neutral
without re-verifying" reasoning as finding 3 — left for the user to weigh against
`nearest_reference_isolate`'s existing callers and its own test suite.

**Severity:** low (design/maintainability preference, not a defect).

## Summary

| # | finding | severity | disposition |
|---|---|---|---|
| 1 | intra-species crosses may get `actual_parents` basis with no comparison behind it | low (downgraded from medium — verified not currently reachable via any diagnostic or report) | proposal — left for user (output-column semantics) |
| 2 | 3 new tests assert on intermediate, not shipped row | low | proposal — left for user (test coverage) |
| 3 | `<2`-species threshold duplicated as independent literal | low | proposal — left for user (refactor risk) |
| 4 | efficiency-fix rationale duplicated across comment + docstring | trivial | **fixed** (docstring consolidation) |
| 5 | helper-function shape heavier than a flag-based API | low | proposal — left for user (refactor risk) |

No `RPI_WEIGHTS`, `tau`, gate logic, verdict, group-construction, or output-column
semantics were changed by this review's autonomous fix. Finding 4's edit is comment-only
(no code lines changed); `ResolveEParentIsolatesSkipsExpensiveResolutionTest`'s 4 tests
were re-run directly and pass. A full-suite rerun follows as part of closing out this
step, before proceeding to the README.
