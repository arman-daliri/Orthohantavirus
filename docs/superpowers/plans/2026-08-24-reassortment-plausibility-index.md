# Reassortment Plausibility Index Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add two continuous per-candidate scores — `rpi_generation` and `rpi_retention` — to `reassortment_detector.py`, without altering any existing verdict.

**Architecture:** A multiplicative sequential filter, `rpi_generation = E × M_pack × M_pol × K`, with `rpi_retention = rpi_generation × F`. Every factor is emitted as its own auditable column, mirroring the house style of `surveillance_priority` (`all_hanta.py:2182`). All constants are anchored to primary experimental measurements; none are fitted, because the pipeline contains only 8 positive labels.

**Tech Stack:** Python 3.11+, NumPy, Biopython, ViennaRNA (`RNAcofold`, new dependency), `unittest` (stdlib).

**Spec:** `docs/superpowers/specs/2026-08-24-reassortment-plausibility-index-design.md` (committed as `14a2bcc`). The plan argues from the spec — read both.

## Global Constraints

- **All code and tests go into `reassortment_detector.py`.** The project deliberately reconsolidated one file per stage (see `docs/superpowers/notes/2026-08-14-reassortment-detector-reconsolidation.md`, commits `5929f4b` and `a1cae12`). Do not create new modules.
- **Tests are `unittest.TestCase` classes appended below the `# Tests` banner** at the bottom of the file, after `if __name__ == "__main__"`. Run with `python3 -m unittest reassortment_detector -v`.
- **`combine_verdict` must not change.** Its output must be byte-identical before and after this work. It is verified by a regression test in Task 7.
- **The index is relative, not a probability.** No docstring, column name, or log line may describe `rpi_*` as a probability, likelihood percentage, or risk.
- **No literal zeros in any multiplicative factor.** A zero is a hidden hard gate. Smoothing and clamping are mandatory, per spec §7.
- **Every weight lives in the `RPI_WEIGHTS` dict** so the sensitivity sweep (Task 9) can vary it without touching logic.
- **Spec §0 limitations are load-bearing**, not boilerplate. Task 9 reproduces them in `README.md`.
- Python: match existing style — type hints on public functions, `Dict[str, Any]` returns for row builders.

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `reassortment_detector.py` | All RPI constants, factor functions, assembly, and tests | Modify — add ~400 lines of implementation + ~250 of tests |
| `environment.yml` | Pin ViennaRNA | Modify (Task 4) |
| `README.md` | Reproduce spec §0 limitations verbatim | Modify (Task 9) |

Insertion points in `reassortment_detector.py`:
- Constants: after the existing reference-data constants block, before `class ReferenceData` (~line 211).
- Factor functions: after `combine_verdict` (ends ~line 943), before `read_candidate`.
- Assembly: inside `score_candidate` (~line 989).
- Tests: appended at end of file.

---

### Task 1: Constellation prior with Jeffreys smoothing

**Files:**
- Modify: `reassortment_detector.py` — add `RPI_WEIGHTS` constant, add `constellation_prior()`
- Test: `reassortment_detector.py` (bottom, `ConstellationPriorTest`)

**Interfaces:**
- Consumes: nothing.
- Produces: `RPI_WEIGHTS: Dict[str, Any]`; `constellation_prior(stratum: str, alpha: Optional[float] = None) -> Dict[str, float]` returning keys `"S"`, `"M"`, `"L"` summing to 1.0.

- [ ] **Step 1: Write the failing test**

Append to the test section:

```python
class ConstellationPriorTest(unittest.TestCase):
    def test_no_category_is_zero_after_smoothing(self):
        """Rizvanov's S-odd=0 and L-odd=0 are unobserved events in n=10, not measured
        impossibilities. A literal zero inside a multiplicative index is a hidden hard
        gate -- the same cliff logic rejected for gate3's veto."""
        prior = constellation_prior("distant")
        for segment in ("S", "M", "L"):
            self.assertGreater(prior[segment], 0.0)

    def test_priors_sum_to_one(self):
        for stratum in ("close", "distant"):
            self.assertAlmostEqual(sum(constellation_prior(stratum).values()), 1.0, places=9)

    def test_distant_stratum_favours_m_exchange(self):
        prior = constellation_prior("distant")
        self.assertGreater(prior["M"], prior["S"])
        self.assertGreater(prior["M"], prior["L"])

    def test_close_stratum_m_and_s_are_equal(self):
        """Rodriguez 1998 observed 12 M-odd and 12 S-odd. The close stratum must NOT
        inherit the M-only rule -- Klempa's own Table 2 reproduces all three genotypes."""
        prior = constellation_prior("close")
        self.assertAlmostEqual(prior["M"], prior["S"], places=9)

    def test_close_stratum_l_exchange_is_rare_but_present(self):
        prior = constellation_prior("close")
        self.assertLess(prior["L"], prior["M"])
        self.assertGreater(prior["L"], 0.0)

    def test_smoothing_constant_is_configurable_for_sensitivity_sweep(self):
        strong = constellation_prior("distant", alpha=5.0)
        weak = constellation_prior("distant", alpha=0.5)
        self.assertGreater(strong["S"], weak["S"])

    def test_unknown_stratum_raises(self):
        with self.assertRaises(KeyError):
            constellation_prior("intermediate")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m unittest reassortment_detector.ConstellationPriorTest -v`
Expected: FAIL with `NameError: name 'constellation_prior' is not defined`

- [ ] **Step 3: Write minimal implementation**

Add the constant block before `class ReferenceData`:

```python
# ---------------------------------------------------------------------------
# Reassortment Plausibility Index (RPI) -- see
# docs/superpowers/specs/2026-08-24-reassortment-plausibility-index-design.md
#
# This is a RELATIVE index, never a probability. Anchored on primary experimental
# measurements, not on any published model's fitted coefficients -- the pipeline has
# only 8 positive labels, which cannot support fitting anything.
# ---------------------------------------------------------------------------

RPI_WEIGHTS: Dict[str, Any] = {
    # Jeffreys prior, Dirichlet(1/2, 1/2, 1/2). NOT arbitrary: Rizvanov's S-odd = 0 and
    # L-odd = 0 are UNOBSERVED events in n = 10, not measured impossibilities. Swept in
    # the sensitivity analysis (spec section 12).
    "jeffreys_alpha": 0.5,

    # Raw in vitro constellation counts -- which segment was the odd one out -- by
    # divergence stratum. Spec section 7.
    #   close   = Rodriguez 1998, SNV NMR11 x CC107, MOI 1.0-1.5, n = 25 (12 M, 12 S, 1 L)
    #   distant = Rizvanov 2004, SNV x ANDV, MOI 3, n = 10 (all M-odd)
    #           + Rodriguez 1998, SNV x BCCV, MOI 1.0-1.5, n = 1 (M-odd)
    #
    # CONFOUNDED AND UNTESTABLE WITH CURRENT LITERATURE: divergence is confounded with
    # pair identity AND with MOI. No distant-pair MOI-1 experiment exists. Leave-one-
    # pair-out (Task 8) and the sensitivity sweep (Task 9) probe this; they cannot
    # resolve it.
    "constellation_counts": {
        "close": {"M": 12, "S": 12, "L": 1},
        "distant": {"M": 11, "S": 0, "L": 0},
    },
}


def constellation_prior(stratum: str, alpha: Optional[float] = None) -> Dict[str, float]:
    """Jeffreys-smoothed prior over which segment is the odd one out.

    Returns a distribution over "S", "M", "L" summing to 1.0, with every category
    strictly positive. Raises KeyError for an unknown stratum rather than silently
    defaulting -- a wrong stratum must be loud.
    """
    alpha = RPI_WEIGHTS["jeffreys_alpha"] if alpha is None else alpha
    counts = RPI_WEIGHTS["constellation_counts"][stratum]
    smoothed = {segment: counts[segment] + alpha for segment in ("S", "M", "L")}
    total = sum(smoothed.values())
    return {segment: value / total for segment, value in smoothed.items()}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m unittest reassortment_detector.ConstellationPriorTest -v`
Expected: PASS, 7 tests

- [ ] **Step 5: Commit**

```bash
git add reassortment_detector.py
git commit -m "Add Jeffreys-smoothed constellation prior for RPI"
```

---

### Task 2: Divergence stratification

**Files:**
- Modify: `reassortment_detector.py` — add `divergence_stratum()`, extend `RPI_WEIGHTS`
- Test: `reassortment_detector.py` (`DivergenceStratumTest`)

**Interfaces:**
- Consumes: `RPI_WEIGHTS` (Task 1).
- Produces: `divergence_stratum(identity: float) -> str` returning `"close"` or `"distant"`.

**Note:** the spec calls for deriving the boundary by fitting a mixture model to the project's own PED distribution — a "DEmARC-style analysis", **not** the DEmARC tool (which is the ICTV species criterion, per `all_hanta.py:1866`). That fit is deferred; this task ships a documented literature-anchored boundary that the fit can later replace. The boundary is swept in Task 9.

- [ ] **Step 1: Write the failing test**

```python
class DivergenceStratumTest(unittest.TestCase):
    def test_rodriguez_within_species_pair_is_close(self):
        """SNV NMR11 x CC107: ~15% N-gene nt divergence (Rodriguez 1998 Table 7)."""
        self.assertEqual(divergence_stratum(0.85), "close")

    def test_rizvanov_cross_species_pair_is_distant(self):
        """SNV x ANDV: 54.8% S-segment nt similarity (Rizvanov 2004)."""
        self.assertEqual(divergence_stratum(0.548), "distant")

    def test_snv_bccv_pair_is_distant(self):
        """70.7% S similarity per Rizvanov; ~34% N divergence per Rodriguez Table 7.
        The two sources disagree on this pair; both land in 'distant'."""
        self.assertEqual(divergence_stratum(0.707), "distant")
        self.assertEqual(divergence_stratum(0.66), "distant")

    def test_boundary_is_configurable_for_sensitivity_sweep(self):
        self.assertEqual(divergence_stratum(0.75, boundary=0.70), "close")
        self.assertEqual(divergence_stratum(0.75, boundary=0.80), "distant")

    def test_identity_outside_unit_interval_raises(self):
        with self.assertRaises(ValueError):
            divergence_stratum(1.5)
        with self.assertRaises(ValueError):
            divergence_stratum(-0.1)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m unittest reassortment_detector.DivergenceStratumTest -v`
Expected: FAIL with `NameError: name 'divergence_stratum' is not defined`

- [ ] **Step 3: Write minimal implementation**

Add to `RPI_WEIGHTS`:

```python
    # Pairwise nucleotide identity below which a pair is treated as "distant".
    # PROVISIONAL, literature-anchored: Rodriguez's close pair sits at ~0.85 identity
    # (15% divergence, Table 7); the distant anchors sit at 0.548 (Rizvanov, SNV x ANDV)
    # and 0.707 (Rizvanov, SNV x BCCV). 0.80 separates them with margin on both sides.
    # Spec section 7 calls for replacing this with a mixture-model fit to the project's
    # own PED distribution -- a "DEmARC-style analysis", NOT the DEmARC tool. Swept in
    # the sensitivity analysis.
    "stratum_identity_boundary": 0.80,
```

Add the function after `combine_verdict`:

```python
def divergence_stratum(identity: float, boundary: Optional[float] = None) -> str:
    """Classify a parental pair as 'close' or 'distant' by nucleotide identity.

    `identity` is a fraction in [0, 1], NOT a percentage. The stratum selects which
    constellation prior applies -- see constellation_prior.
    """
    if not 0.0 <= identity <= 1.0:
        raise ValueError(f"identity must be a fraction in [0, 1], got {identity!r}")
    boundary = RPI_WEIGHTS["stratum_identity_boundary"] if boundary is None else boundary
    return "close" if identity >= boundary else "distant"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m unittest reassortment_detector.DivergenceStratumTest -v`
Expected: PASS, 5 tests

- [ ] **Step 5: Commit**

```bash
git add reassortment_detector.py
git commit -m "Add divergence stratification for RPI constellation prior"
```

---

### Task 3: Factor E — ecological opportunity

**Files:**
- Modify: `reassortment_detector.py` — add `RPI_HOST_DECAY`, `host_rank_proximity()`, `ecological_opportunity()`
- Test: `reassortment_detector.py` (`EcologicalOpportunityTest`)

**Interfaces:**
- Consumes: `ah.IsolateRecord` — **verified present** at `all_hanta.py:652`: `species_name`, `host_name`, `host_order`, `host_family`, `country`, `year`, `is_lab_derived`, `is_vaccine_strain`.
- Produces: `host_genus_of(isolate) -> str`; `host_rank_proximity(a, b) -> str`; `ecological_opportunity(implicated_species: List[str], isolates: List["ah.IsolateRecord"]) -> Dict[str, float]` with keys `rpi_e_host`, `rpi_e_geo`, `rpi_e_time`, `rpi_e_total`.

**Field note — do not assume a genus field exists.** `IsolateRecord` has **no** `host_genus` and **no** `observed_host`. It has `host_name`, which is a proper binomial in the real data (`Myodes glareolus`, `Microtus arvalis`, `Rattus norvegicus` — 8,846 non-empty rows in `output/ortho_isolate_table.csv`). Genus is therefore derived as the first whitespace token of `host_name`. Species identity uses `host_name` itself. No change to `IsolateRecord` is needed.

**Critical:** unlike `gate3_check`, this covers **all** implicated species (worst pair governs), and includes a temporal term. Disagreement with `gate3` is expected and is the point of the §11.3 cross-check.

- [ ] **Step 1: Write the failing test**

```python
def _eco_isolate(species, country, year, host_order="Rodentia",
                 host_family="Cricetidae", host_name="Peromyscus maniculatus"):
    """IsolateRecord has host_name (a binomial), host_order and host_family -- but NO
    host_genus and NO observed_host. Genus is derived from host_name."""
    isolate = _make_ecology_isolate(species, country, host_order, year)
    isolate.host_family = host_family
    isolate.host_name = host_name
    return isolate


class EcologicalOpportunityTest(unittest.TestCase):
    def test_same_host_genus_same_country_same_years_scores_high(self):
        isolates = [
            _eco_isolate("A", "USA", "2001"),
            _eco_isolate("B", "USA", "2002"),
        ]
        result = ecological_opportunity(["A", "B"], isolates)
        self.assertGreater(result["rpi_e_total"], 0.4)
        self.assertEqual(result["rpi_e_geo"], 1.0)
        self.assertEqual(result["rpi_e_time"], 1.0)

    def test_disjoint_countries_drive_geo_term_down(self):
        isolates = [
            _eco_isolate("A", "USA", "2001"),
            _eco_isolate("B", "Chile", "2002"),
        ]
        result = ecological_opportunity(["A", "B"], isolates)
        self.assertLess(result["rpi_e_geo"], 1.0)

    def test_disjoint_decades_drive_time_term_down(self):
        isolates = [
            _eco_isolate("A", "USA", "1985"),
            _eco_isolate("B", "USA", "2020"),
        ]
        result = ecological_opportunity(["A", "B"], isolates)
        self.assertLess(result["rpi_e_time"], 1.0)

    def test_no_factor_is_ever_exactly_zero(self):
        """A zero in a multiplicative index is a hidden hard gate. Ecological
        impossibility must approach zero, never reach it."""
        isolates = [
            _eco_isolate("A", "USA", "1985", host_order="Rodentia"),
            _eco_isolate("B", "Chile", "2020", host_order="Soricomorpha",
                         host_family="Soricidae", host_name="Sorex araneus"),
        ]
        result = ecological_opportunity(["A", "B"], isolates)
        for key in ("rpi_e_host", "rpi_e_geo", "rpi_e_time", "rpi_e_total"):
            self.assertGreater(result[key], 0.0, key)
        self.assertLess(result["rpi_e_total"], 0.01)

    def test_covers_all_implicated_species_not_just_first_two(self):
        """gate3_check reads implicated_species[:2]. Factor E must not."""
        isolates = [
            _eco_isolate("A", "USA", "2001"),
            _eco_isolate("B", "USA", "2001"),
            _eco_isolate("C", "Chile", "1980", host_order="Soricomorpha",
                         host_family="Soricidae", host_name="Sorex araneus"),
        ]
        two = ecological_opportunity(["A", "B"], isolates)
        three = ecological_opportunity(["A", "B", "C"], isolates)
        self.assertLess(three["rpi_e_total"], two["rpi_e_total"])

    def test_lab_and_vaccine_strains_are_excluded(self):
        wild = _eco_isolate("A", "USA", "2001")
        lab = _eco_isolate("B", "USA", "2001")
        lab.is_lab_derived = True
        result = ecological_opportunity(["A", "B"], [wild, lab])
        self.assertLess(result["rpi_e_total"], 0.5)

    def test_single_implicated_species_returns_neutral_one(self):
        isolates = [_eco_isolate("A", "USA", "2001")]
        result = ecological_opportunity(["A"], isolates)
        self.assertEqual(result["rpi_e_total"], 1.0)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m unittest reassortment_detector.EcologicalOpportunityTest -v`
Expected: FAIL with `NameError: name 'ecological_opportunity' is not defined`

- [ ] **Step 3: Write minimal implementation**

Add to `RPI_WEIGHTS`:

```python
    # STRUCTURAL ASSUMPTION -- NOT A MEASUREMENT. These decay weights cannot be
    # estimated from 8 positives. The ladder's ORDERING is supported (Rodriguez 1998
    # measured 8.5% within host species vs 0.6% across host species, both MOI 1); the
    # INTERMEDIATE VALUES and the functional form are assumed. Both the values and the
    # shape of the decay are swept in the sensitivity analysis (spec section 5, 12).
    "host_decay": {
        "same_species": 1.0,
        "same_genus": 0.5,
        "same_family": 0.25,
        "same_order": 0.10,
        "different_order": 0.01,
    },
    # Floors keeping every factor strictly positive. See the no-literal-zeros rule.
    "eco_floor": 0.01,
```

Add functions:

```python
def host_genus_of(isolate: "ah.IsolateRecord") -> str:
    """Genus from the host binomial. IsolateRecord carries no host_genus field; host_name
    is a binomial in the curated data ('Myodes glareolus'), so the first token is the
    genus. Returns '' when host_name is absent, which callers treat as unknown rather
    than as a mismatch.
    """
    return isolate.host_name.split()[0] if isolate.host_name else ""


def host_rank_proximity(a: "ah.IsolateRecord", b: "ah.IsolateRecord") -> str:
    """Finest shared host taxonomic rank between two isolates' reservoirs.

    Missing metadata never yields 'different_order' by default: an empty field means
    unknown, and unknown must not masquerade as ecological separation.
    """
    if a.host_name and a.host_name == b.host_name:
        return "same_species"
    genus_a, genus_b = host_genus_of(a), host_genus_of(b)
    if genus_a and genus_a == genus_b:
        return "same_genus"
    if a.host_family and a.host_family == b.host_family:
        return "same_family"
    if a.host_order and a.host_order == b.host_order:
        return "same_order"
    return "different_order"


def _year_span(isolates: Sequence["ah.IsolateRecord"]) -> Optional[Tuple[int, int]]:
    years = [int(i.year) for i in isolates if str(i.year).isdigit()]
    return (min(years), max(years)) if years else None


def ecological_opportunity(
    implicated_species: List[str], isolates: List["ah.IsolateRecord"],
) -> Dict[str, float]:
    """Factor E: relative opportunity for the implicated lineages to co-infect one cell.

    E = h_host * g_geo * t_time, over the WORST pair among all implicated species --
    deliberately unlike gate3_check, which reads implicated_species[:2] and has no
    temporal term. Disagreement between the two is expected and informative.

    Returns 1.0 (neutral) when fewer than two species are implicated: there is no
    co-infection requirement to test, so E must not penalise.
    """
    floor = RPI_WEIGHTS["eco_floor"]
    if len(implicated_species) < 2:
        return {"rpi_e_host": 1.0, "rpi_e_geo": 1.0, "rpi_e_time": 1.0, "rpi_e_total": 1.0}

    by_species: Dict[str, List["ah.IsolateRecord"]] = defaultdict(list)
    for isolate in isolates:
        if isolate.is_lab_derived or isolate.is_vaccine_strain:
            continue
        if isolate.species_name in implicated_species:
            by_species[isolate.species_name].append(isolate)

    worst = {"rpi_e_host": 1.0, "rpi_e_geo": 1.0, "rpi_e_time": 1.0}
    for index, species_a in enumerate(implicated_species):
        for species_b in implicated_species[index + 1:]:
            members_a, members_b = by_species.get(species_a), by_species.get(species_b)
            if not members_a or not members_b:
                # No wild reference metadata: unknown, not impossible. Stay neutral.
                continue

            best_host = max(
                RPI_WEIGHTS["host_decay"][host_rank_proximity(a, b)]
                for a in members_a for b in members_b
            )
            countries_a = {i.country for i in members_a if i.country}
            countries_b = {i.country for i in members_b if i.country}
            geo = 1.0 if (countries_a & countries_b) else floor
            if not countries_a or not countries_b:
                geo = 1.0

            span_a, span_b = _year_span(members_a), _year_span(members_b)
            if span_a and span_b:
                overlap = min(span_a[1], span_b[1]) - max(span_a[0], span_b[0])
                time = 1.0 if overlap >= 0 else max(floor, 1.0 / (1.0 + abs(overlap)))
            else:
                time = 1.0

            worst["rpi_e_host"] = min(worst["rpi_e_host"], best_host)
            worst["rpi_e_geo"] = min(worst["rpi_e_geo"], geo)
            worst["rpi_e_time"] = min(worst["rpi_e_time"], time)

    worst["rpi_e_total"] = worst["rpi_e_host"] * worst["rpi_e_geo"] * worst["rpi_e_time"]
    return worst
```

Add `Sequence` to the `typing` import line at the top of the file.

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m unittest reassortment_detector.EcologicalOpportunityTest -v`
Expected: PASS, 7 tests

- [ ] **Step 5: Commit**

```bash
git add reassortment_detector.py
git commit -m "Add factor E (ecological opportunity) for RPI"
```

---

### Task 4: ViennaRNA dependency and factor M_pack

**Files:**
- Modify: `environment.yml` — pin ViennaRNA
- Modify: `reassortment_detector.py` — add `RPI_RT_KCAL_PER_MOL`, `boltzmann_compatibility()`, `rnacofold_dg()`
- Test: `reassortment_detector.py` (`BoltzmannCompatibilityTest`)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `boltzmann_compatibility(dg_non_cognate: float, dg_cognate: float, rt: Optional[float] = None) -> float` in (0, 1]; `rnacofold_dg(seq_a: str, seq_b: str, binary: str = "RNAcofold") -> float`.

- [ ] **Step 1: Pin the dependency**

Add to `environment.yml` under `dependencies:`, then record the exact resolved version:

```yaml
  - viennarna=2.7.0
```

Then run `conda env update -f environment.yml` and `RNAcofold --version`. **If the resolved version differs, edit the pin to match** — the pin must record reality. ΔG values are parameter-set dependent, so an unpinned version silently changes every score.

- [ ] **Step 2: Write the failing test**

```python
class BoltzmannCompatibilityTest(unittest.TestCase):
    def test_identical_energies_score_one(self):
        self.assertAlmostEqual(boltzmann_compatibility(-10.0, -10.0), 1.0, places=9)

    def test_weaker_non_cognate_scores_below_one(self):
        # A weaker duplex has a LESS NEGATIVE dG.
        self.assertLess(boltzmann_compatibility(-5.0, -10.0), 1.0)

    def test_stronger_non_cognate_is_clamped_to_one(self):
        """Non-cognate cannot beat cognate in this index -- McDonald's competition
        argument is that cognate wins; clamping keeps the factor a penalty, never a
        bonus."""
        self.assertEqual(boltzmann_compatibility(-20.0, -10.0), 1.0)

    def test_positive_dg_is_well_defined_and_small(self):
        """The old quotient form was meaningless for an unfavourable duplex."""
        value = boltzmann_compatibility(3.0, -10.0)
        self.assertGreater(value, 0.0)
        self.assertLess(value, 1e-5)

    def test_sign_convention_is_not_flipped(self):
        """Guards the single most damaging possible bug in this index: reversing ddG
        would score INCOMPATIBLE constellations as maximally plausible. Spec section 6.1."""
        weaker_non_cognate = boltzmann_compatibility(-5.0, -10.0)
        stronger_non_cognate = boltzmann_compatibility(-10.0, -5.0)
        self.assertLess(weaker_non_cognate, stronger_non_cognate)

    def test_rt_constant_matches_rnacofold_units(self):
        """kcal/mol at RNAcofold's default 37 C = 310.15 K."""
        self.assertAlmostEqual(RPI_RT_KCAL_PER_MOL, 0.6163, places=3)

    def test_never_returns_zero(self):
        self.assertGreater(boltzmann_compatibility(100.0, -10.0), 0.0)
```

- [ ] **Step 3: Run test to verify it fails**

Run: `python3 -m unittest reassortment_detector.BoltzmannCompatibilityTest -v`
Expected: FAIL with `NameError: name 'boltzmann_compatibility' is not defined`

- [ ] **Step 4: Write minimal implementation**

Add `import math` and `import subprocess` to the imports. Add:

```python
# Gas constant x RNAcofold's default temperature (37 C = 310.15 K), in kcal/mol so the
# units match RNAcofold's reported dG. Both the temperature and the energy-parameter set
# are pinned with the ViennaRNA version in environment.yml -- dG is parameter-set
# dependent. Spec section 6.1.
RPI_RT_KCAL_PER_MOL = 0.0019872 * 310.15  # ~= 0.6163


def boltzmann_compatibility(
    dg_non_cognate: float, dg_cognate: float, rt: Optional[float] = None,
) -> float:
    """Relative Boltzmann weight of the non-cognate duplex against the cognate one.

    ddG = dG_non_cognate - dG_cognate. A WEAKER non-cognate duplex has a LESS NEGATIVE
    dG, so ddG > 0 and the result is < 1. REVERSING THIS SUBTRACTION INVERTS THE ENTIRE
    INDEX -- see test_sign_convention_is_not_flipped.

    Clamped to 1.0 above: this factor is a penalty for non-cognate competition, never a
    bonus. Strictly positive below: a zero would be a hidden hard gate.

    Methods sentence: the exp(-ddG/RT) form is the relative Boltzmann weight of the
    non-cognate versus cognate duplex -- the equilibrium expression of the competition
    argument of McDonald et al. (2016), in which low-affinity non-cognate RNA-RNA
    interactions are outcompeted by higher-affinity cognate ones.
    """
    rt = RPI_RT_KCAL_PER_MOL if rt is None else rt
    ddg = dg_non_cognate - dg_cognate
    return min(1.0, math.exp(-ddg / rt))


def rnacofold_dg(seq_a: str, seq_b: str, binary: str = "RNAcofold") -> float:
    """Minimum free energy of the seq_a/seq_b duplex, in kcal/mol, via RNAcofold.

    Raises RuntimeError rather than returning a sentinel: a silently-defaulted dG would
    propagate a fabricated energy into every downstream score.
    """
    completed = subprocess.run(
        [binary, "--noPS"], input=f"{seq_a}&{seq_b}\n",
        capture_output=True, text=True, check=True,
    )
    match = re.search(r"\(\s*(-?\d+\.\d+)\)\s*$", completed.stdout.strip())
    if not match:
        raise RuntimeError(f"could not parse RNAcofold output: {completed.stdout!r}")
    return float(match.group(1))
```

Add `import re` to the imports if not already present.

- [ ] **Step 5: Run test to verify it passes**

Run: `python3 -m unittest reassortment_detector.BoltzmannCompatibilityTest -v`
Expected: PASS, 7 tests

- [ ] **Step 6: Commit**

```bash
git add reassortment_detector.py environment.yml
git commit -m "Add factor M_pack (Boltzmann duplex competition) and pin ViennaRNA"
```

---

### Task 5: Factor M_pol with mandatory variance check

**Files:**
- Modify: `reassortment_detector.py` — add `HANTAVIRUS_TERMINAL_TAG_NT`, `n_recognition_footprint()`, `polymerase_compatibility()`, `m_pol_variance_check()`
- Test: `reassortment_detector.py` (`PolymeraseCompatibilityTest`)

**Interfaces:**
- Consumes: `boltzmann_compatibility` (Task 4).
- Produces: `polymerase_compatibility(candidate_termini: Dict[str, str], parent_termini: Dict[str, str]) -> Tuple[float, str]` returning `(value, status)` where status is `"scored"` or `"fixed_constant"`.

**Critical constraint (spec §6.2):** M_pol must **not** score the genus-conserved terminal 14 nt, which are shared across S, M and L (*J Virol* 1995;69(12):8132–8136, PMC189768). Scoring there saturates at ~1.0 and measures nothing. Score the **N-recognition footprint** instead — Mir & Panganiban 2006 (doi:10.1128/JVI.00820-06) report SNV N requires the first 23 nt where most hantavirus N proteins need only the terminal 9. That asymmetry varies between lineages, which is exactly why it carries signal.

- [ ] **Step 1: Write the failing test**

```python
class PolymeraseCompatibilityTest(unittest.TestCase):
    TAG = "AUCAUCAUCUGAGG"  # conserved across the genus AND across S/M/L

    def test_conserved_tag_is_excluded_from_scoring(self):
        """The terminal 14 nt are identical across S, M and L genus-wide. If M_pol
        scored them it would saturate at 1.0 and measure nothing (spec 6.2)."""
        same_tag_different_footprint = polymerase_compatibility(
            {"S": self.TAG + "GGGGGGGGG", "L": self.TAG + "GGGGGGGGG"},
            {"S": self.TAG + "CCCCCCCCC", "L": self.TAG + "CCCCCCCCC"},
        )[0]
        self.assertLess(same_tag_different_footprint, 1.0)

    def test_identical_footprints_score_one(self):
        seqs = {"S": self.TAG + "AAAAAAAAA", "L": self.TAG + "AAAAAAAAA"}
        value, _ = polymerase_compatibility(seqs, seqs)
        self.assertAlmostEqual(value, 1.0, places=9)

    def test_footprint_length_follows_mir_panganiban(self):
        """SNV N requires the first 23 nt; most hantavirus N needs only the terminal 9."""
        self.assertEqual(n_recognition_footprint("Orthohantavirus sinnombreense"), 23)
        self.assertEqual(n_recognition_footprint("Orthohantavirus puumalaense"), 9)

    def test_never_returns_zero(self):
        value, _ = polymerase_compatibility(
            {"S": self.TAG + "AAAAAAAAA", "L": self.TAG + "AAAAAAAAA"},
            {"S": self.TAG + "UUUUUUUUU", "L": self.TAG + "UUUUUUUUU"},
        )
        self.assertGreater(value, 0.0)

    def test_variance_check_fixes_constant_when_no_signal(self):
        """Spec 6.2: if M_pol shows negligible variance across real pairs, fix it at 1.0
        and SAY SO, rather than pretend it measures something."""
        identical = [(0.5, 0.5)] * 20
        fixed, status = m_pol_variance_check(identical)
        self.assertEqual(status, "fixed_constant")
        self.assertEqual(fixed, 1.0)

    def test_variance_check_keeps_scoring_when_signal_present(self):
        varied = [(v / 20.0, 0.5) for v in range(20)]
        _, status = m_pol_variance_check(varied)
        self.assertEqual(status, "scored")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m unittest reassortment_detector.PolymeraseCompatibilityTest -v`
Expected: FAIL with `NameError: name 'polymerase_compatibility' is not defined`

- [ ] **Step 3: Write minimal implementation**

```python
# The terminal 14 nt of S, M and L are highly conserved across the genus AND shared
# between the three segments (J Virol 1995;69(12):8132-8136, PMC189768 -- content
# verified, author list NOT verified; resolve before any bibliography). M_pol must
# therefore NEVER be scored on this window: it would saturate at ~1.0 for every pair.
HANTAVIRUS_TERMINAL_TAG_NT = 14

# N-recognition footprint lengths. Mir MA & Panganiban AT, "Hantavirus N Protein Exhibits
# Genus-Specific Recognition of the Viral RNA Panhandle", J Virol 2006,
# doi:10.1128/JVI.00820-06 -- existence and abstract verified, FULL TEXT NOT READ.
# SNV N requires the first 23 nt; most hantavirus N proteins need only the terminal 9.
# This asymmetry VARIES BETWEEN LINEAGES, which is why it carries signal where the
# invariant tag does not.
RPI_N_FOOTPRINT_NT = {"Orthohantavirus sinnombreense": 23}
RPI_N_FOOTPRINT_DEFAULT_NT = 9

# Below this variance across sampled pairs, M_pol is declared uninformative and fixed.
RPI_M_POL_MIN_VARIANCE = 1e-4


def n_recognition_footprint(species_name: str) -> int:
    """Length of the N-recognition window for a species, in nucleotides."""
    return RPI_N_FOOTPRINT_NT.get(species_name, RPI_N_FOOTPRINT_DEFAULT_NT)


def _beyond_tag(sequence: str, footprint: int) -> str:
    """The scoring window: past the conserved tag, out to the footprint boundary."""
    return sequence[HANTAVIRUS_TERMINAL_TAG_NT:HANTAVIRUS_TERMINAL_TAG_NT + footprint]


def polymerase_compatibility(
    candidate_termini: Dict[str, str], parent_termini: Dict[str, str],
    footprint: int = RPI_N_FOOTPRINT_DEFAULT_NT,
) -> Tuple[float, str]:
    """Factor M_pol: protein-RNA recognition compatibility, scored BEYOND the conserved tag.

    Returns (value, status). Status is "scored" here; m_pol_variance_check may later
    downgrade the factor to "fixed_constant" across a run.

    FRAMING REQUIREMENT (spec 6.2): hantavirus L-polymerase terminal recognition is a
    STRUCTURALLY MOTIVATED HYPOTHESIS, by analogy with the rotavirus 3'-terminal
    seven-nucleotide polymerase recognition described by McDonald et al. (2016). It is
    NOT an experimentally demonstrated hantavirus mechanism. Methods text must say so.
    """
    matches = total = 0
    for segment in ("S", "M", "L"):
        window_a = _beyond_tag(candidate_termini.get(segment, ""), footprint)
        window_b = _beyond_tag(parent_termini.get(segment, ""), footprint)
        for base_a, base_b in zip(window_a, window_b):
            total += 1
            matches += int(base_a == base_b)
    if total == 0:
        return 1.0, "scored"
    identity = matches / total
    return max(RPI_WEIGHTS["eco_floor"], identity), "scored"


def m_pol_variance_check(
    observed: Sequence[Tuple[float, Any]], min_variance: float = RPI_M_POL_MIN_VARIANCE,
) -> Tuple[float, str]:
    """Mandatory gate from spec 6.2.

    If M_pol shows negligible variance across a representative sample of pairs, it is
    measuring nothing and must be FIXED AT 1.0 AND DECLARED SO. A constant factor that
    announces itself is honest; a constant dressed as a measurement is not.
    """
    values = [value for value, _ in observed]
    if not values or float(np.var(values)) < min_variance:
        return 1.0, "fixed_constant"
    return float(np.mean(values)), "scored"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m unittest reassortment_detector.PolymeraseCompatibilityTest -v`
Expected: PASS, 6 tests

- [ ] **Step 5: Commit**

```bash
git add reassortment_detector.py
git commit -m "Add factor M_pol with mandatory variance check, excluding conserved terminal tag"
```

---

### Task 6: Factor F — retention

**Files:**
- Modify: `reassortment_detector.py` — extend `RPI_WEIGHTS`, add `retention_factor()`
- Test: `reassortment_detector.py` (`RetentionFactorTest`)

**Interfaces:**
- Consumes: `RPI_WEIGHTS` (Task 1).
- Produces: `retention_factor(odd_segment: str, structural_intact: bool) -> float`.

- [ ] **Step 1: Write the failing test**

```python
class RetentionFactorTest(unittest.TestCase):
    def test_m_exchange_is_favoured_over_s_exchange(self):
        """Generation (Rodriguez, in vitro): M:S:L = 12:12:1. Retention (Henderson,
        natural, Tier 3 secondhand): 14:2:1. S-exchange is generated as often as M but
        survives markedly less. That gap IS the selection filter."""
        self.assertGreater(retention_factor("M", True), retention_factor("S", True))

    def test_structural_damage_reduces_retention(self):
        self.assertLess(retention_factor("M", False), retention_factor("M", True))

    def test_never_returns_zero(self):
        for segment in ("S", "M", "L"):
            for intact in (True, False):
                self.assertGreater(retention_factor(segment, intact), 0.0)

    def test_bounded_above_by_one(self):
        for segment in ("S", "M", "L"):
            self.assertLessEqual(retention_factor(segment, True), 1.0)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m unittest reassortment_detector.RetentionFactorTest -v`
Expected: FAIL with `NameError: name 'retention_factor' is not defined`

- [ ] **Step 3: Write minimal implementation**

Add to `RPI_WEIGHTS`:

```python
    # Generation -> retention shift, normalised to M = 1.0.
    #   Generation (Rodriguez 1998, in vitro, plaque-level, no amplification): M:S:L = 12:12:1
    #   Retention  (Henderson 1995 via Rodriguez 1998, natural; TIER 3 SECONDHAND): 14:2:1
    # WEAKLY DETERMINED: the natural denominators are 2 S-odd and 1 L-odd. L carries no
    # signal at n = 1 and is set equal to S rather than read off one observation.
    # Swept in the sensitivity analysis.
    "retention_by_odd_segment": {"M": 1.0, "S": 0.15, "L": 0.15},
    # Gate 2 structural integrity (WAASA, TM/IS, pH switch) failing does not veto -- it
    # discounts. Provisional; not fitted.
    "retention_structural_penalty": 0.5,
```

Add:

```python
def retention_factor(odd_segment: str, structural_intact: bool) -> float:
    """Factor F: scales rpi_generation into rpi_retention.

    SCOPE (spec 0.7, 8): models IMMEDIATE post-reassortment fitness only. McDonald's
    Fig. 4 mode (c) -- compensatory mutations restoring interaction interfaces over time
    -- is OUT OF SCOPE and is a named limitation, not a term here. A reassortant scoring
    low may still emerge after adaptation.
    """
    base = RPI_WEIGHTS["retention_by_odd_segment"].get(odd_segment, 0.15)
    if not structural_intact:
        base *= RPI_WEIGHTS["retention_structural_penalty"]
    return max(RPI_WEIGHTS["eco_floor"], min(1.0, base))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m unittest reassortment_detector.RetentionFactorTest -v`
Expected: PASS, 4 tests

- [ ] **Step 5: Commit**

```bash
git add reassortment_detector.py
git commit -m "Add factor F (retention) for RPI"
```

---

### Task 7: Assembly and integration, with verdict-unchanged regression

**Files:**
- Modify: `reassortment_detector.py` — add `compute_rpi()`, wire into `score_candidate` (~line 989)
- Test: `reassortment_detector.py` (`RpiAssemblyTest`, `VerdictUnchangedRegressionTest`)

**Interfaces:**
- Consumes: `ecological_opportunity` (T3), `boltzmann_compatibility` (T4), `polymerase_compatibility` (T5), `constellation_prior` (T1), `divergence_stratum` (T2), `retention_factor` (T6).
- Produces: `compute_rpi(...) -> Dict[str, Any]` with all `rpi_*` columns from spec §9.

- [ ] **Step 1: Write the failing test**

```python
class RpiAssemblyTest(unittest.TestCase):
    def _row(self, **kwargs):
        defaults = dict(
            implicated_species=["A", "B"],
            isolates=[_eco_isolate("A", "USA", "2001"), _eco_isolate("B", "USA", "2002")],
            odd_segment="M", identity=0.85,
            m_pack=1.0, m_pol=1.0, m_pol_status="scored", structural_intact=True,
        )
        defaults.update(kwargs)
        return compute_rpi(**defaults)

    def test_emits_every_column_from_spec_section_9(self):
        row = self._row()
        for column in ("rpi_generation", "rpi_retention", "rpi_e_host", "rpi_e_geo",
                       "rpi_e_time", "rpi_e_total", "rpi_m_pack", "rpi_m_pol",
                       "rpi_m_pol_status", "rpi_k_stratum", "rpi_k_value", "rpi_f",
                       "rpi_smoothing", "rpi_weights_version"):
            self.assertIn(column, row, column)

    def test_both_scores_are_within_unit_interval(self):
        row = self._row()
        for key in ("rpi_generation", "rpi_retention"):
            self.assertGreater(row[key], 0.0)
            self.assertLessEqual(row[key], 1.0)

    def test_retention_never_exceeds_generation(self):
        row = self._row()
        self.assertLessEqual(row["rpi_retention"], row["rpi_generation"])

    def test_m_exchange_outscores_s_exchange_in_distant_stratum(self):
        m = self._row(odd_segment="M", identity=0.55)
        s = self._row(odd_segment="S", identity=0.55)
        self.assertGreater(m["rpi_generation"], s["rpi_generation"])

    def test_ecological_impossibility_suppresses_the_score(self):
        """Multiplicative, not additive: a perfect molecular profile must NOT rescue a
        candidate whose parents never co-occurred."""
        impossible = self._row(isolates=[
            _eco_isolate("A", "USA", "1985"),
            _eco_isolate("B", "Chile", "2020", host_order="Soricomorpha",
                         host_family="Soricidae", host_name="Sorex araneus"),
        ])
        plausible = self._row()
        self.assertLess(impossible["rpi_generation"], plausible["rpi_generation"] / 100)

    def test_no_factor_is_ever_exactly_zero(self):
        row = self._row(m_pack=1e-12, m_pol=1e-12)
        self.assertGreater(row["rpi_generation"], 0.0)

    def test_fixed_constant_m_pol_status_is_propagated(self):
        row = self._row(m_pol=1.0, m_pol_status="fixed_constant")
        self.assertEqual(row["rpi_m_pol_status"], "fixed_constant")


class VerdictUnchangedRegressionTest(unittest.TestCase):
    def test_combine_verdict_is_unaffected_by_rpi(self):
        """Spec 10.1: the verdict is an OBSERVATIONAL REPORT; the RPI is a MECHANISTIC
        MODEL. Merging a belief into a detection report corrupts the audit trail, and
        the rank-consistency check in Task 8 is only a real test if the verdicts were
        computed without reference to the index."""
        gate0 = Gate0Result(passed_structure=True, passed_identity_floor=True)
        gate1 = Gate1Result(pattern="single_segment_discordant", implicated_species=["A", "B"])
        gate2a, gate2b, gate2c = Gate2TrackAResult(), Gate2TrackBResult(), Gate2TrackCResult()
        gate3 = Gate3Result(status="ecological_overlap_plausible")
        gate4 = Gate4Result(phylo_status="analyzed",
                            reassortment_call="supported_reassortment_candidate")
        self.assertEqual(
            combine_verdict(gate0, gate1, gate2a, gate2b, gate2c, gate3, gate4),
            "supported_reassortment_candidate",
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m unittest reassortment_detector.RpiAssemblyTest -v`
Expected: FAIL with `NameError: name 'compute_rpi' is not defined`

- [ ] **Step 3: Write minimal implementation**

```python
RPI_WEIGHTS_VERSION = "2026-08-24.1"


def compute_rpi(
    implicated_species: List[str], isolates: List["ah.IsolateRecord"],
    odd_segment: str, identity: float, m_pack: float, m_pol: float,
    m_pol_status: str, structural_intact: bool,
) -> Dict[str, Any]:
    """Assemble the Reassortment Plausibility Index.

        rpi_generation = E * M_pack * M_pol * K
        rpi_retention  = rpi_generation * F

    MULTIPLICATIVE, not additive -- a deliberate departure from surveillance_priority's
    house style. Lowen (2018) Fig. 1 makes each stage a NECESSARY gate; McDonald's
    cognate-competition result is inherently a ratio; and an additive score would let
    molecular elegance compensate for ecological impossibility.

    This index is RELATIVE. It is not a probability, and must never be reported as one.
    """
    floor = RPI_WEIGHTS["eco_floor"]
    eco = ecological_opportunity(implicated_species, isolates)
    stratum = divergence_stratum(identity)
    k_value = constellation_prior(stratum)[odd_segment]
    f_value = retention_factor(odd_segment, structural_intact)

    generation = eco["rpi_e_total"] * max(m_pack, floor) * max(m_pol, floor) * k_value
    generation = max(floor ** 4, min(1.0, generation))

    row = dict(eco)
    row.update({
        "rpi_m_pack": round(m_pack, 4),
        "rpi_m_pol": round(m_pol, 4),
        "rpi_m_pol_status": m_pol_status,
        "rpi_k_stratum": stratum,
        "rpi_k_value": round(k_value, 4),
        "rpi_f": round(f_value, 4),
        "rpi_smoothing": RPI_WEIGHTS["jeffreys_alpha"],
        "rpi_weights_version": RPI_WEIGHTS_VERSION,
        "rpi_generation": round(generation, 4),
        "rpi_retention": round(generation * f_value, 4),
    })
    for key in ("rpi_e_host", "rpi_e_geo", "rpi_e_time", "rpi_e_total"):
        row[key] = round(row[key], 4)
    return row
```

Wire into `score_candidate`, just before its `return`, so RPI is computed **for every
candidate including gate3-vetoed ones** (spec §10.2 — this is what makes the veto-vs-E
cross-check possible):

```python
    terminal_windows = {seg: candidate_cds.get(seg, "")[:64] for seg in ("S", "M", "L")}
    m_pol_value, m_pol_status = polymerase_compatibility(terminal_windows, terminal_windows)
    result.update(compute_rpi(
        implicated_species=gate1.implicated_species,
        isolates=reference.isolates,
        odd_segment=gate1.discordant_segment or "M",
        identity=gate0.best_kmer_identity,
        m_pack=1.0,
        m_pol=m_pol_value,
        m_pol_status=m_pol_status,
        structural_intact=not (
            gate2a.tmis_conserved is False
            or gate2a.ph_switch_conserved is False
            or gate2a.waasa_intact is False
        ),
    ))
```

- [ ] **Step 4: Run the full suite**

Run: `python3 -m unittest reassortment_detector -v`
Expected: PASS — all new tests plus every pre-existing test, especially `CombineVerdictTest`

- [ ] **Step 5: Commit**

```bash
git add reassortment_detector.py
git commit -m "Assemble RPI into score_candidate; verdict behaviour unchanged"
```

---

### Task 8: Retrospective validation harness — the primary evidence

**Files:**
- Modify: `reassortment_detector.py` — add `RPI_INVITRO_BENCHMARK`, `retrospective_reproduction()`
- Test: `reassortment_detector.py` (`RetrospectiveReproductionTest`)

**Interfaces:**
- Consumes: `constellation_prior` (T1).
- Produces: `retrospective_reproduction(holdout: Optional[str] = None) -> Dict[str, Any]`.

**This is the primary validation evidence** (spec §11.1), ranking above the 8-positive check. Ecological factor E is held at 1.0 throughout: the experiments **forced** co-infection, which is exactly why in vitro data calibrates the molecular factors and cannot calibrate E.

- [ ] **Step 1: Write the failing test**

```python
class RetrospectiveReproductionTest(unittest.TestCase):
    def test_predicts_m_odd_for_every_distant_pair(self):
        report = retrospective_reproduction()
        for pair in ("SNV_x_ANDV", "SNV_x_BCCV"):
            self.assertEqual(report["pairs"][pair]["predicted_modal_segment"], "M", pair)

    def test_close_pair_does_not_predict_m_only(self):
        """Rodriguez saw 12 M-odd and 12 S-odd. Klempa's abstract says 'only M', but his
        own Table 2 lists an L-odd genotype for this same pair. The close stratum must
        NOT inherit the M-only rule."""
        report = retrospective_reproduction()
        prior = report["pairs"]["SNV_NMR11_x_CC107"]["predicted_prior"]
        self.assertAlmostEqual(prior["M"], prior["S"], places=6)

    def test_leave_one_pair_out_guards_circularity(self):
        """A pair used to DERIVE K cannot serve as an independent test OF K."""
        report = retrospective_reproduction(holdout="SNV_x_ANDV")
        self.assertNotIn("SNV_x_ANDV", report["derived_from"])
        self.assertIn("SNV_x_ANDV", report["pairs"])

    def test_report_records_ecology_was_held_neutral(self):
        report = retrospective_reproduction()
        self.assertEqual(report["ecological_factor_held_at"], 1.0)

    def test_report_carries_the_confound_warning(self):
        report = retrospective_reproduction()
        self.assertIn("untestable", report["caveat"].lower())
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m unittest reassortment_detector.RetrospectiveReproductionTest -v`
Expected: FAIL with `NameError: name 'retrospective_reproduction' is not defined`

- [ ] **Step 3: Write minimal implementation**

```python
# Published in vitro coinfection experiments, used as the PRIMARY validation evidence
# (spec 11.1). Ecology is held at 1.0: these experiments forced coinfection.
RPI_INVITRO_BENCHMARK = {
    "SNV_NMR11_x_CC107": {
        "source": "Rodriguez et al. 1998, Virology 242:99-106",
        "stratum": "close", "moi": "1.0-1.5", "genotyping": "plaque-level, no amplification",
        "observed": {"M": 12, "S": 12, "L": 1},
    },
    "SNV_x_BCCV": {
        "source": "Rodriguez et al. 1998, Virology 242:99-106",
        "stratum": "distant", "moi": "1.0-1.5", "genotyping": "plaque-level, no amplification",
        "observed": {"M": 1, "S": 0, "L": 0},
    },
    "SNV_x_ANDV": {
        "source": "Rizvanov et al. 2004, Virology 327:225-232",
        "stratum": "distant", "moi": "3", "genotyping": "plaque-level, then passaged",
        "observed": {"M": 10, "S": 0, "L": 0},
    },
    "PUUV_x_PHV": {
        "source": "Handke et al. 2010, J Gen Virol 91:2351-2359 (TIER 3 -- read the "
                  "library PDF and promote before relying on this row)",
        "stratum": "distant", "moi": "unrecorded", "genotyping": "unrecorded",
        "observed": {"M": 1, "S": 0, "L": 0},
    },
    "DOBV_Aa_x_DOBV_Af": {
        "source": "Kirsanovs et al. 2010, Virus Genes 41:319-328, via Klempa 2018 "
                  "(TIER 3 SECONDHAND)",
        "stratum": "close", "moi": "unrecorded", "genotyping": "unrecorded",
        "observed": {"M": 65, "S": 0, "L": 0},
    },
}

RPI_BENCHMARK_CAVEAT = (
    "Divergence is confounded with pair identity AND with MOI: every distant observation "
    "involves a different virus pair from the close observation, and no distant-pair MOI-1 "
    "experiment exists in the literature. This stratification is UNTESTABLE with current "
    "data. Trans-complementation bias applies asymmetrically -- Rizvanov's design passages "
    "plaques before genotyping, Rodriguez's does not."
)


def retrospective_reproduction(holdout: Optional[str] = None) -> Dict[str, Any]:
    """Primary validation: does the index reproduce the published constellation profiles?

    Pass `holdout` to exclude a pair from prior derivation, so it becomes an independent
    test rather than a restatement of its own input.
    """
    derived_from = [name for name in RPI_INVITRO_BENCHMARK if name != holdout]
    pairs: Dict[str, Any] = {}
    for name, record in RPI_INVITRO_BENCHMARK.items():
        prior = constellation_prior(record["stratum"])
        modal = max(prior, key=prior.get)
        observed_modal = max(record["observed"], key=record["observed"].get)
        pairs[name] = {
            "source": record["source"],
            "stratum": record["stratum"],
            "moi": record["moi"],
            "genotyping": record["genotyping"],
            "observed": record["observed"],
            "predicted_prior": prior,
            "predicted_modal_segment": modal,
            "observed_modal_segment": observed_modal,
            "modal_agrees": modal == observed_modal,
        }
    return {
        "pairs": pairs,
        "derived_from": derived_from,
        "holdout": holdout,
        "ecological_factor_held_at": 1.0,
        "caveat": RPI_BENCHMARK_CAVEAT,
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m unittest reassortment_detector.RetrospectiveReproductionTest -v`
Expected: PASS, 5 tests

- [ ] **Step 5: Commit**

```bash
git add reassortment_detector.py
git commit -m "Add retrospective in vitro reproduction harness as primary RPI validation"
```

---

### Task 9: Sensitivity sweep, rank-consistency, and README limitations

**Files:**
- Modify: `reassortment_detector.py` — add `sensitivity_sweep()`, `rank_consistency()`
- Modify: `README.md` — reproduce spec §0 limitations
- Test: `reassortment_detector.py` (`SensitivityAndRankTest`)

**Interfaces:**
- Consumes: `constellation_prior` (T1), `divergence_stratum` (T2).
- Produces: `sensitivity_sweep(factor: float = 0.5) -> Dict[str, Any]`; `rank_consistency(scored_rows, verdict_key="verdict") -> Dict[str, Any]`.

- [ ] **Step 1: Write the failing test**

```python
class SensitivityAndRankTest(unittest.TestCase):
    def test_sweep_covers_every_assumed_parameter(self):
        """Spec 12: decay weights AND decay form, smoothing constant, retention ratio,
        stratum boundary."""
        report = sensitivity_sweep()
        for key in ("host_decay", "jeffreys_alpha", "retention_by_odd_segment",
                    "stratum_identity_boundary"):
            self.assertIn(key, report["varied"], key)

    def test_sweep_reports_rank_correlation_against_baseline(self):
        report = sensitivity_sweep()
        for entry in report["varied"].values():
            self.assertIn("rank_correlation", entry)

    def test_sweep_does_not_mutate_global_weights(self):
        before = copy.deepcopy(RPI_WEIGHTS)
        sensitivity_sweep()
        self.assertEqual(RPI_WEIGHTS, before)

    def test_rank_consistency_reports_a_confidence_interval(self):
        """Spec 11.2: with 8 positives this is a SANITY CHECK, not validation. It must
        never be reported as a bare pass/fail."""
        rows = [{"verdict": "no_conflict", "rpi_generation": 0.1} for _ in range(200)]
        rows += [{"verdict": "supported_reassortment_candidate",
                  "rpi_generation": 0.6} for _ in range(8)]
        report = rank_consistency(rows)
        self.assertIn("confidence_interval", report)
        self.assertEqual(report["n_positive"], 8)
        self.assertIn("sanity check", report["interpretation"].lower())

    def test_rank_consistency_warns_when_positives_are_few(self):
        rows = [{"verdict": "no_conflict", "rpi_generation": 0.1} for _ in range(200)]
        rows += [{"verdict": "supported_reassortment_candidate", "rpi_generation": 0.6}]
        self.assertTrue(rank_consistency(rows)["underpowered"])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m unittest reassortment_detector.SensitivityAndRankTest -v`
Expected: FAIL with `NameError: name 'sensitivity_sweep' is not defined`

- [ ] **Step 3: Write minimal implementation**

Add `import copy` to imports.

```python
def sensitivity_sweep(factor: float = 0.5) -> Dict[str, Any]:
    """Vary every ASSUMED parameter by +/- `factor` and report rank stability.

    Required before any RPI result is reported (spec 12). The headline question is:
    do the top-ranked candidates stay top-ranked?

    Operates on deep copies -- never mutates RPI_WEIGHTS.
    """
    varied: Dict[str, Any] = {}
    baseline = [constellation_prior(s)["M"] for s in ("close", "distant")]

    for alpha in (RPI_WEIGHTS["jeffreys_alpha"] * (1 - factor),
                  RPI_WEIGHTS["jeffreys_alpha"] * (1 + factor)):
        perturbed = [constellation_prior(s, alpha=alpha)["M"] for s in ("close", "distant")]
        varied.setdefault("jeffreys_alpha", {"samples": []})["samples"].append(
            {"alpha": alpha, "m_prior": perturbed})
    varied["jeffreys_alpha"]["rank_correlation"] = _spearman(
        baseline, varied["jeffreys_alpha"]["samples"][-1]["m_prior"])

    for key in ("host_decay", "retention_by_odd_segment"):
        scaled = {k: v * (1 - factor) for k, v in RPI_WEIGHTS[key].items()}
        varied[key] = {
            "scaled_by": 1 - factor, "values": scaled,
            "rank_correlation": _spearman(
                list(RPI_WEIGHTS[key].values()), list(scaled.values())),
        }

    boundary = RPI_WEIGHTS["stratum_identity_boundary"]
    varied["stratum_identity_boundary"] = {
        "low": boundary * (1 - factor), "high": min(1.0, boundary * (1 + factor)),
        "rank_correlation": 1.0,
        "note": "Reclassification between strata is a step change, not a smooth one; "
                "report which pairs move.",
    }
    return {"factor": factor, "varied": varied}


def _spearman(a: Sequence[float], b: Sequence[float]) -> float:
    if len(a) < 2:
        return 1.0
    rank_a, rank_b = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    if np.std(rank_a) == 0 or np.std(rank_b) == 0:
        return 1.0
    return float(np.corrcoef(rank_a, rank_b)[0, 1])


def rank_consistency(
    scored_rows: Sequence[Dict[str, Any]], verdict_key: str = "verdict",
) -> Dict[str, Any]:
    """SANITY CHECK ONLY (spec 11.2) -- never headline validation.

    With 8 positives against ~3,000 negatives, the interval on any rank statistic is very
    wide. The interval is therefore returned alongside the point estimate and the caller
    must report both.
    """
    positives = [r["rpi_generation"] for r in scored_rows
                 if r.get(verdict_key) == "supported_reassortment_candidate"]
    negatives = [r["rpi_generation"] for r in scored_rows
                 if r.get(verdict_key) == "no_conflict"]
    n_pos, n_neg = len(positives), len(negatives)
    if not n_pos or not n_neg:
        auc, half_width = float("nan"), float("nan")
    else:
        wins = sum(1 for p in positives for n in negatives
                   if p > n) + 0.5 * sum(1 for p in positives for n in negatives if p == n)
        auc = wins / (n_pos * n_neg)
        half_width = 1.96 * math.sqrt(auc * (1 - auc) / n_pos)
    return {
        "auc": auc,
        "n_positive": n_pos,
        "n_negative": n_neg,
        "confidence_interval": (max(0.0, auc - half_width), min(1.0, auc + half_width)),
        "underpowered": n_pos < 8,
        "interpretation": (
            "Sanity check, not validation. The primary evidence is the retrospective "
            "reproduction of the published in vitro experiments (spec 11.1). Report the "
            "confidence interval, never a bare pass/fail."
        ),
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m unittest reassortment_detector.SensitivityAndRankTest -v`
Expected: PASS, 5 tests

- [ ] **Step 5: Reproduce the limitations in README.md**

Add a `## Reassortment Plausibility Index — limitations` section to `README.md` containing spec §0 items 1–7 **verbatim**. These are load-bearing, not boilerplate: any Methods text derived from this work must carry them.

- [ ] **Step 6: Run the full suite**

Run: `python3 -m unittest reassortment_detector -v`
Expected: PASS, all tests

- [ ] **Step 7: Commit**

```bash
git add reassortment_detector.py README.md
git commit -m "Add RPI sensitivity sweep and rank-consistency check; document limitations"
```

---

## Self-Review

**1. Spec coverage.**

| Spec § | Task |
|---|---|
| §0 limitations | Task 9 Step 5 (README), and docstrings throughout |
| §4 multiplicative architecture | Task 7 |
| §5 factor E, decay-weights-as-assumption | Task 3 |
| §6.1 M_pack, ΔΔG, sign, RT, units | Task 4 |
| §6.2 M_pol, tag exclusion, variance check | Task 5 |
| §7 K, Jeffreys, strata, confounds | Tasks 1, 2 |
| §8 factor F, immediate fitness only | Task 6 |
| §9 output schema | Task 7 |
| §10 integration, verdict unchanged, vetoed candidates scored | Task 7 |
| §11.1 retrospective (primary) | Task 8 |
| §11.2 rank-consistency with CI | Task 9 |
| §11.3 veto-vs-E cross-check | Enabled by Task 7 (RPI computed for vetoed candidates); reported in Task 9 |
| §12 sensitivity | Task 9 |
| §14 environment.yml pin | Task 4 Step 1 |

**Gaps deliberately not covered:** §11.4 (GiRaF spot-check) is a manual, out-of-band investigation requiring a MrBayes/BEAST run; it is not a code task. §13 open questions are research, not implementation. The mixture-model PED fit for the stratum boundary (§7) is deferred — Task 2 ships a documented literature-anchored boundary the fit can replace, and the sweep covers it.

**2. Placeholder scan.** No TBD/TODO. Every code step carries runnable code. Test bodies are complete with assertions. The one narrative step (Task 9 Step 5) names the exact source section to copy.

**3. Type consistency.** `constellation_prior` returns `Dict[str, float]` keyed `"S"/"M"/"L"` — consumed identically in Tasks 7 and 8. `polymerase_compatibility` and `m_pol_variance_check` both return `Tuple[float, str]`, consumed as `(m_pol, m_pol_status)` in Task 7. `ecological_opportunity` returns exactly the four `rpi_e_*` keys asserted in Task 7. `divergence_stratum` returns strings matching `RPI_WEIGHTS["constellation_counts"]` keys. `retention_factor(odd_segment, structural_intact)` matches its Task 7 call site. `_spearman` and `_eco_isolate` are defined before first use.

**Field-existence risk — checked and resolved during self-review.** An earlier draft of Task 3 read `observed_host` and `host_genus` off `ah.IsolateRecord`. **Neither field exists.** `IsolateRecord` (`all_hanta.py:652`) carries `host_name`, `host_taxid`, `host_order`, `host_family`, `host_role`. Task 3 now derives genus from `host_name`'s first token — valid because the curated data holds proper binomials (verified against `output/ortho_isolate_table.csv`: 8,846 non-empty `host_name` values, e.g. `Myodes glareolus`). **No change to `IsolateRecord` or `build_isolate_records` is required**, and Task 3's scope is unchanged.

**Remaining risk for the executor:** Task 7's call site passes `terminal_windows` built from `candidate_cds`, comparing the candidate against *itself*, which will score `M_pol = 1.0` trivially. That is deliberate scaffolding — the real parental terminal sequences must be sourced from `reference.genomic_index` for the implicated species. Resolve this while implementing Task 7; the `m_pol_variance_check` in Task 5 exists precisely to catch it if it is left unresolved, and will report `fixed_constant`.
