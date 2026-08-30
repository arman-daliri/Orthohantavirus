# rpi_discrimination.py
"""WS4 Task 13 (D7) -- the four-group RPI discrimination panel.

Builds four groups of candidates, scores each through `reassortment_detector.score_candidate`,
and reports the resulting `rpi_generation`/`rpi_retention` distributions so Task 14's pass
criterion (`reals ~ b1 > b2 > random`) can be evaluated against real numbers, not intuition.

THE FOUR GROUPS
----------------
  (a) RANDOM FAKES -- a full S+M+L triplet assembled from three different isolates whose
      species AND continents all differ from one another. The weakest possible test:
      Factor E (ecological_opportunity) alone should already reject these on
      geography/timing, before M_pack or M_pol ever get a vote. Included as a floor, not
      a discriminating test -- reported as such, not oversold.

  (b1) WITHIN-SPECIES SHUFFLES -- a full S+M+L triplet assembled from three DIFFERENT
      isolates of the SAME species. This is RPI's plausibility-ceiling control: these
      genuinely are plausible assemblies (real segments that could have circulated
      together), and RPI measures plausibility, not occurrence -- scoring HIGH here is
      the correct, expected answer, not a false positive. Separating these from real
      reassortants is gate 4's job (Task 14), not this index's.

  (b2) CROSS-SPECIES SYMPATRIC -- a full S+M+L triplet built across two DIFFERENT species
      that nonetheless share host order, country, AND overlapping collection years in
      this dataset. THE REAL TEST of M_pack + M_pol: Factor E cannot reject these on
      ecological grounds (they really do co-occur), so any separation from group (c)'s
      real positives has to come from the molecular factors. Measured (2026-08-29): only
      6 such (host_order, country, species_A, species_B) tuples exist in the usable pool,
      all Rodentia hosts in China or South Korea (hantanense / khabarovskense / seoulense
      / puumalaense) -- true ecological sympatry between different hantavirus species is
      a narrow, specific configuration, which is exactly what makes it a real test rather
      than an easy one.

  (c) REAL POSITIVES -- two independently-sourced sets, neither constructed by this
      module: the 8 `supported_reassortment_candidate` isolates from all_hanta.py's own
      genome-wide phylogenetic surveillance (`output/ortho_reassortment_table.csv`,
      scored here on their own real, unmodified S/M/L), plus the published inter-species
      lab crosses already reconstructed in `validate_reassortment_detector.py`'s
      `PUBLISHED_CROSSES` (class A only -- class B's intra-species crosses are gate1's
      known blind spot, not something RPI's E factor can distinguish from group b1
      either, so they are excluded from "real positives" here to keep the comparison to
      inter-species events, matching what b2 and random both are).

DISJOINTNESS FROM CALIBRATION (binding constraint, PROJECT-STATE.md section 7b). Every
group builder here excludes the 42 isolates `data/rpi_calibration.json` ->
`reserved_for_calibration.isolate_keys` reserved for deriving tau. Overlap would let this
panel's b2 result be read as validating a tau that was tuned partly on these same
isolates -- the exact circularity `rpi_calibration.py`'s own docstring warns against, in
the other direction.

WHAT THIS MODULE DOES NOT DO. It does not compute the Task 14 pass criterion (Mann-Whitney
U, effect sizes, the reals~b1 equivalence interval) -- that is a separate, smaller function
over this module's OUTPUT, kept separate so the panel construction (subject to real
data availability) and the statistical test (subject to methodology) can be reviewed and
changed independently. It also does not touch gate 4 (Task 14's "detector separation"
half) -- gate4_check needs mafft/iqtree, and their Windows blockers (iqtree3.exe naming,
MAFFT --version parsing) are still open per PROJECT-STATE.md.

Run:
    python rpi_discrimination.py                      # builds panel, scores it, prints summary
    python rpi_discrimination.py --json report.json   # also writes the full scored panel
    python -m unittest rpi_discrimination              # unit tests
"""

import argparse
import csv
import json
import random
import sys
import unittest
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import all_hanta as ah
import reassortment_detector as rd
import rpi_calibration as cal
import validate_reassortment_detector as vd


# Every RNG in this module takes this explicitly, mirroring rpi_calibration.py's own rule
# -- nothing here reads global random state, so a rerun on the same reference data
# reproduces the same panel.
RPI_DISCRIMINATION_SEED = 20260829

# Sample sizes. None of these is exhaustive except (b2), which uses every sympatric
# tuple this dataset actually contains (see module docstring) -- (a) and (b1) are capped
# because both COULD be enumerated far more densely (hantanense alone has 170 usable
# isolates) and a denser sample would not answer a different question, only take longer.
RANDOM_FAKES_TARGET = 20
WITHIN_SPECIES_SHUFFLES_PER_SPECIES = 10
WITHIN_SPECIES_MIN_ISOLATES = 3  # need 3 distinct isolates to fill S, M, L independently

# Minimal, explicit map -- NOT a general geo library. Covers exactly the countries
# present in the terminally-complete pool (verified 2026-08-29); a KeyError on a new
# country is the intended failure mode (silently guessing a continent would defeat the
# point of the "different continents" criterion), not a gap to paper over with a
# permissive default.
#
# Russia is mapped Europe: every Russian isolate in this dataset with usable metadata is
# a Puumala-lineage sample from European Russia, not the Far East -- checked against the
# panel this module actually builds from, not assumed.
COUNTRY_CONTINENT: Dict[str, str] = {
    "Belgium": "Europe", "China": "Asia", "Finland": "Europe", "Germany": "Europe",
    "Greece": "Europe", "Guinea": "Africa", "Lithuania": "Europe", "Russia": "Europe",
    "South Korea": "Asia", "Sweden": "Europe", "USA": "North America",
    "United Kingdom": "Europe", "Venezuela": "South America", "Viet Nam": "Asia",
}

ORTHO_REASSORTMENT_TABLE_PATH = Path("output") / "ortho_reassortment_table.csv"


# ---------------------------------------------------------------------------
# Shared pool: terminally complete, wild, disjoint from calibration
# ---------------------------------------------------------------------------


def reserved_calibration_isolate_keys(
        artifact_path: Path = cal.DEFAULT_ARTIFACT_PATH) -> Set[str]:
    """The 42 isolate keys `data/rpi_calibration.json` reserved for deriving tau.

    Loaded from the pinned artifact, not re-derived -- re-deriving here would risk this
    module and rpi_calibration.py silently disagreeing about which 42 isolates are
    reserved if either one's selection logic ever changes without the other noticing.
    """
    with open(artifact_path, encoding="utf-8") as f:
        artifact = json.load(f)
    return set(artifact["reserved_for_calibration"]["isolate_keys"])


def usable_pool(
        reference: "rd.ReferenceData",
        reserved_keys: Optional[Set[str]] = None) -> List["ah.IsolateRecord"]:
    """Terminally complete, wild isolates, excluding calibration's reserved 42 (binding
    constraint, PROJECT-STATE.md section 7b). Reuses rpi_calibration.py's own filter
    (`terminally_complete_isolates`) rather than a separate implementation -- both this
    module's groups and the calibration panel need the SAME notion of "a real, complete,
    genuinely-terminal genome" for the same reason (D2/WS3.1's finding that 48%+ of
    deposits are trimmed of their true termini and would corrupt any panhandle-adjacent
    computation).
    """
    if reserved_keys is None:
        reserved_keys = reserved_calibration_isolate_keys()
    pool = cal.terminally_complete_isolates(reference)
    return [i for i in pool if i.isolate_key not in reserved_keys]


def _candidate_from_isolates(
        candidate_id: str, reference: "rd.ReferenceData",
        segment_sources: Dict[str, "ah.IsolateRecord"]) -> "rd.CandidateInput":
    """A CandidateInput whose S/M/L each come from `segment_sources[seg]` -- one isolate
    per segment, which may be the same isolate for all three (a real, unmodified genome)
    or three different ones (a constructed shuffle).

    `known_segment_isolates=segment_sources` (design doc 2026-08-30-rpi-e-geo-isolate-
    level-fix): every case this module builds already knows exactly which real isolate
    contributed each segment, so Factor E's geo/host/time terms should read those real
    isolates directly (rpi_e_isolate_basis="actual_parents") rather than fall back to a
    nearest-reference proxy. For group c's supported_reassortment_candidate rows, the
    SAME isolate stands for all three segments -- resolve_e_parent_isolates only
    recognises this as covering an implicated species when gate1 actually names that
    isolate's own species, so the (unknown) OTHER/donor species in a genuinely
    discordant real candidate still correctly falls through to the proxy fallback.
    """
    genomic = {
        seg: cal.genomic_sequence(reference, cal.segment_accession(isolate, seg))
        for seg, isolate in segment_sources.items()
    }
    return rd.CandidateInput(
        candidate_id=candidate_id, s_genomic=genomic["S"], m_genomic=genomic["M"],
        l_genomic=genomic["L"], known_segment_isolates=dict(segment_sources))


def _label(isolate: "ah.IsolateRecord") -> str:
    return vd.PanelBuilder.label(isolate)


# ---------------------------------------------------------------------------
# (a) Random fakes: different species AND different continents
# ---------------------------------------------------------------------------


def build_group_a_random(
        reference: "rd.ReferenceData", pool: Sequence["ah.IsolateRecord"],
        seed: int = RPI_DISCRIMINATION_SEED, target: int = RANDOM_FAKES_TARGET,
) -> List[Dict[str, Any]]:
    """Full S+M+L triplets, each segment from a different isolate, with all three
    isolates' species AND continents pairwise distinct (D7 group a).

    Raises KeyError up front if any NON-EMPTY country in `pool` is missing from
    COUNTRY_CONTINENT -- deliberate, see that table's own comment. Isolates with no
    country at all are silently excluded (missing metadata, not a stale table -- the
    same distinction gate3_check draws between "ecology_unknown" and "no_ecological_
    route"), but a country the table simply never learned about must fail loudly rather
    than be quietly dropped from the eligible pool.
    """
    unmapped = sorted({i.country for i in pool if i.country
                       and i.country not in COUNTRY_CONTINENT})
    if unmapped:
        raise KeyError(f"COUNTRY_CONTINENT is missing: {unmapped}")
    eligible = [i for i in pool if i.country]
    rng = random.Random(seed)
    cases: List[Dict[str, Any]] = []
    attempts = 0
    seen_combos: Set[Tuple[str, str, str]] = set()
    # Bounded attempt budget: with few species/continents represented, not every draw
    # satisfies the pairwise-distinct constraint, and an unbounded loop over a pool that
    # cannot reach `target` combinations would spin forever.
    max_attempts = target * 50
    while len(cases) < target and attempts < max_attempts:
        attempts += 1
        s_iso, m_iso, l_iso = rng.sample(eligible, 3)
        species = (s_iso.species_name, m_iso.species_name, l_iso.species_name)
        continents = {COUNTRY_CONTINENT[i.country] for i in (s_iso, m_iso, l_iso)}
        if len(set(species)) != 3 or len(continents) < 2:
            continue
        combo = tuple(sorted(i.isolate_key for i in (s_iso, m_iso, l_iso)))
        if combo in seen_combos:
            continue
        seen_combos.add(combo)
        case_id = f"RANDOM_{len(cases):03d}"
        cases.append({
            "case_id": case_id, "klass": "a_random_fakes",
            "candidate": _candidate_from_isolates(
                case_id, reference, {"S": s_iso, "M": m_iso, "L": l_iso}),
            "backbone_species": s_iso.species_name,  # descriptive only -- no single
            "donor_species": f"{m_iso.species_name} (M) / {l_iso.species_name} (L)",
            "donated_segment": "M+L",
            "segment_sources": {
                "S": (s_iso.species_name, _label(s_iso)),
                "M": (m_iso.species_name, _label(m_iso)),
                "L": (l_iso.species_name, _label(l_iso)),
            },
            "continents": sorted(continents),
        })
    return cases


# ---------------------------------------------------------------------------
# (b1) Within-species shuffles: plausibility-ceiling control
# ---------------------------------------------------------------------------


def build_group_b1_within_species(
        reference: "rd.ReferenceData", pool: Sequence["ah.IsolateRecord"],
        seed: int = RPI_DISCRIMINATION_SEED,
        per_species: int = WITHIN_SPECIES_SHUFFLES_PER_SPECIES,
) -> List[Dict[str, Any]]:
    """Full S+M+L triplets, each segment from a DIFFERENT isolate of the SAME species
    (D7 group b1). Expected to score HIGH -- see module docstring.

    Draws are independent across shuffles (may reuse an isolate across different
    shuffles of the same species) but never within one shuffle (S, M, L are always three
    distinct isolates) -- with as few as WITHIN_SPECIES_MIN_ISOLATES qualifying isolates
    for the smallest species, requiring disjoint draws ACROSS shuffles too would leave
    almost no shuffles for that species.
    """
    by_species: Dict[str, List["ah.IsolateRecord"]] = defaultdict(list)
    for isolate in pool:
        by_species[isolate.species_name].append(isolate)

    rng = random.Random(seed)
    cases: List[Dict[str, Any]] = []
    for species in sorted(by_species):
        members = sorted(by_species[species], key=lambda i: i.isolate_key)
        if len(members) < WITHIN_SPECIES_MIN_ISOLATES:
            continue
        seen_combos: Set[Tuple[str, str, str]] = set()
        for _ in range(per_species * 5):  # bounded attempt budget, as in group (a)
            if len(seen_combos) >= per_species:
                break
            s_iso, m_iso, l_iso = rng.sample(members, 3)
            combo = (s_iso.isolate_key, m_iso.isolate_key, l_iso.isolate_key)
            if combo in seen_combos:
                continue
            seen_combos.add(combo)
            case_id = f"B1_{species.split()[-1]}_{len(cases):03d}"
            cases.append({
                "case_id": case_id, "klass": "b1_within_species_shuffle",
                "candidate": _candidate_from_isolates(
                    case_id, reference, {"S": s_iso, "M": m_iso, "L": l_iso}),
                "backbone_species": species, "donor_species": species,
                "donated_segment": "M+L",
                "segment_sources": {
                    "S": (species, _label(s_iso)), "M": (species, _label(m_iso)),
                    "L": (species, _label(l_iso)),
                },
            })
    return cases


# ---------------------------------------------------------------------------
# (b2) Cross-species sympatric: the real test
# ---------------------------------------------------------------------------


def _year_int(year: str) -> Optional[int]:
    try:
        return int(str(year)[:4])
    except (ValueError, TypeError):
        return None


def find_sympatric_species_pairs(
        pool: Sequence["ah.IsolateRecord"],
) -> List[Tuple[str, str, str, str]]:
    """(host_order, country, species_a, species_b) tuples where two DIFFERENT species
    genuinely co-occur in this dataset: same host order, same country, and overlapping
    collection years. D7 group b2's precondition -- see module docstring for why this
    is deliberately NOT a large set.
    """
    by_key: Dict[Tuple[str, str], List["ah.IsolateRecord"]] = defaultdict(list)
    for isolate in pool:
        if isolate.host_order and isolate.country:
            by_key[(isolate.host_order, isolate.country)].append(isolate)

    pairs: Set[Tuple[str, str, str, str]] = set()
    for (host_order, country), isolates in by_key.items():
        years_by_species: Dict[str, List[int]] = defaultdict(list)
        for isolate in isolates:
            year = _year_int(isolate.year)
            if year is not None:
                years_by_species[isolate.species_name].append(year)
        species_list = sorted(years_by_species)
        for a_idx in range(len(species_list)):
            for b_idx in range(a_idx + 1, len(species_list)):
                sp_a, sp_b = species_list[a_idx], species_list[b_idx]
                years_a, years_b = years_by_species[sp_a], years_by_species[sp_b]
                if min(years_a) <= max(years_b) and min(years_b) <= max(years_a):
                    pairs.add((host_order, country, sp_a, sp_b))
    return sorted(pairs)


def build_group_b2_sympatric(
        reference: "rd.ReferenceData", pool: Sequence["ah.IsolateRecord"],
) -> List[Dict[str, Any]]:
    """Full S+M+L triplets built across genuinely sympatric species pairs (D7 group b2).

    EXHAUSTIVE over `find_sympatric_species_pairs`'s real tuples, both backbone/donor
    directions, and all three donated-segment positions -- not sampled, because the
    precondition (real ecological sympatry between different species) is already the
    limiting factor; there is no larger population to subsample from.
    """
    tuples = find_sympatric_species_pairs(pool)
    by_species: Dict[str, List["ah.IsolateRecord"]] = defaultdict(list)
    for isolate in pool:
        by_species[isolate.species_name].append(isolate)

    cases: List[Dict[str, Any]] = []
    for host_order, country, sp_a, sp_b in tuples:
        for backbone_species, donor_species in ((sp_a, sp_b), (sp_b, sp_a)):
            backbone_candidates = [
                i for i in by_species[backbone_species]
                if i.host_order == host_order and i.country == country]
            donor_candidates = [
                i for i in by_species[donor_species]
                if i.host_order == host_order and i.country == country]
            if not backbone_candidates or not donor_candidates:
                continue
            backbone = backbone_candidates[0]
            donor = donor_candidates[0]
            for donated_segment in ("S", "M", "L"):
                sources = {seg: (donor if seg == donated_segment else backbone)
                          for seg in ("S", "M", "L")}
                case_id = (f"B2_{backbone_species.split()[-1]}_x_"
                          f"{donor_species.split()[-1]}_{donated_segment}")
                cases.append({
                    "case_id": case_id, "klass": "b2_cross_species_sympatric",
                    "candidate": _candidate_from_isolates(case_id, reference, sources),
                    "backbone_species": backbone_species, "donor_species": donor_species,
                    "donated_segment": donated_segment,
                    "sympatric_context": {"host_order": host_order, "country": country},
                    "segment_sources": {
                        seg: (iso.species_name, _label(iso))
                        for seg, iso in sources.items()
                    },
                })
    return cases


# ---------------------------------------------------------------------------
# (c) Real positives
# ---------------------------------------------------------------------------


def load_supported_reassortment_candidates(
        reference: "rd.ReferenceData",
        table_path: Path = ORTHO_REASSORTMENT_TABLE_PATH,
) -> List[Dict[str, Any]]:
    """The `supported_reassortment_candidate` rows from all_hanta.py's own genome-wide
    phylogenetic surveillance (`output/ortho_reassortment_table.csv`) -- an
    INDEPENDENTLY-sourced set of real positives, not constructed by this module. Scored
    on the isolate's own real, unmodified S/M/L.

    De-duplicated by isolate_key: the source table can list the same isolate more than
    once across different phylogroup analysis scopes.
    """
    by_key = {i.isolate_key: i for i in reference.isolates}
    seen: Set[str] = set()
    cases: List[Dict[str, Any]] = []
    with open(table_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["reassortment_call"] != "supported_reassortment_candidate":
                continue
            key = row["isolate_key"]
            if key in seen:
                continue
            seen.add(key)
            isolate = by_key.get(key)
            if isolate is None:
                continue
            case_id = f"REAL_ORTHO_{isolate.species_name.split()[-1]}_{len(cases):03d}"
            cases.append({
                "case_id": case_id, "klass": "c_real_positive",
                "candidate": _candidate_from_isolates(
                    case_id, reference,
                    {"S": isolate, "M": isolate, "L": isolate}),
                "backbone_species": isolate.species_name,
                "donor_species": isolate.species_name,
                "donated_segment": row.get("discordant_segment", ""),
                "source": "all_hanta.py genome-wide phylogenetic surveillance "
                          "(output/ortho_reassortment_table.csv)",
                "source_isolate_key": key,
                "source_bootstrap_support": row.get("discordant_clade_bootstrap_support"),
            })
    return cases


def load_published_inter_species_crosses(
        reference: "rd.ReferenceData", builder: "vd.PanelBuilder",
) -> List[Dict[str, Any]]:
    """The class-A (inter-species) entries of `validate_reassortment_detector.py`'s
    `PUBLISHED_CROSSES` -- real lab-published reassortants. Class B (intra-species) is
    EXCLUDED here: gate1 cannot see those (the known blind spot), so they measure the
    same thing group b1 already measures for RPI's purposes, not an independent "real
    positive" signal.
    """
    cases: List[Dict[str, Any]] = []
    for spec in vd.PUBLISHED_CROSSES:
        if spec["klass"] != "A_published_inter_species":
            continue
        built = vd.build_swap(builder, spec)
        if built is None:
            continue
        cases.append({
            "case_id": built["case_id"], "klass": "c_real_positive",
            "candidate": built["candidate"],
            "backbone_species": spec["backbone_species"],
            "donor_species": spec["donor_species"],
            "donated_segment": spec["donated_segment"],
            "source": f"published lab cross: {spec['paper']}",
        })
    return cases


def build_group_c_real_positives(
        reference: "rd.ReferenceData", builder: "vd.PanelBuilder",
) -> List[Dict[str, Any]]:
    return (load_supported_reassortment_candidates(reference)
            + load_published_inter_species_crosses(reference, builder))


# ---------------------------------------------------------------------------
# Scoring and orchestration
# ---------------------------------------------------------------------------


def _score_args(output_dir: str = "output/discrimination_scratch") -> argparse.Namespace:
    """gate4 (phylogenetics) is disabled the same way validate_reassortment_detector.py
    disables it by default: a comparator-count threshold no real species reaches, so
    gate4_check returns "insufficient_comparators" without shelling out to mafft/iqtree.
    This module measures RPI, not gate4's tree-based call -- Task 14's "detector
    separation" half is a separate piece of work, gated on the still-open MAFFT/IQ-TREE
    Windows blockers (PROJECT-STATE.md)."""
    return argparse.Namespace(
        output=output_dir, mafft_bin="mafft", iqtree_bin="iqtree3",
        phylo_min_taxa=vd._PHYLO_DISABLED_MIN_TAXA, bootstrap=0,
        min_bootstrap_support=70, phylo_seed=42, phylo_threads=1)


def score_cases(
        reference: "rd.ReferenceData", precomputed: Dict[str, Any],
        cases: Sequence[Dict[str, Any]], args: Optional[argparse.Namespace] = None,
) -> List[Dict[str, Any]]:
    """Runs score_candidate on every case and folds the RPI columns into it."""
    if args is None:
        args = _score_args()
    scored = []
    for case in cases:
        result = rd.score_candidate(case["candidate"], reference, precomputed, args)
        row = dict(case)
        row.pop("candidate")  # not JSON-serialisable and not needed past this point
        for key in ("rpi_generation", "rpi_retention", "rpi_m_pack", "rpi_m_pack_status",
                    "rpi_m_pol", "rpi_m_pol_status", "rpi_e_total", "rpi_e_isolate_basis",
                    "rpi_k_stratum", "gate1_pattern", "gate1_discordant_segment",
                    "gate1_implicated_species", "verdict"):
            row[key] = result.get(key)
        scored.append(row)
    return scored


def build_discrimination_panel(
        reference: "rd.ReferenceData", seed: int = RPI_DISCRIMINATION_SEED,
) -> Dict[str, List[Dict[str, Any]]]:
    """Builds all four (unscored) groups. Separate from scoring so a caller can inspect
    or filter the panel before paying score_candidate's cost."""
    reserved = reserved_calibration_isolate_keys()
    pool = usable_pool(reference, reserved)
    builder = vd.PanelBuilder(reference)
    return {
        "a_random_fakes": build_group_a_random(reference, pool, seed=seed),
        "b1_within_species_shuffle": build_group_b1_within_species(reference, pool, seed=seed),
        "b2_cross_species_sympatric": build_group_b2_sympatric(reference, pool),
        "c_real_positive": build_group_c_real_positives(reference, builder),
    }


def summarise_groups(scored_groups: Dict[str, List[Dict[str, Any]]]) -> Dict[str, Any]:
    """Descriptive summary only -- Task 14's Mann-Whitney U / effect-size pass criterion
    is a separate function over this same scored output, not computed here."""
    summary = {}
    for name, rows in scored_groups.items():
        values = [r["rpi_generation"] for r in rows if r.get("rpi_generation") is not None]
        summary[name] = {
            "n": len(rows), "n_scored": len(values),
            "rpi_generation_mean": (sum(values) / len(values)) if values else None,
            "rpi_generation_min": min(values) if values else None,
            "rpi_generation_max": max(values) if values else None,
        }
    return summary


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-input", default=".",
                        help="Package root containing all_hanta/ncbi_dataset/data")
    parser.add_argument("--json", help="Write the full scored panel to this path")
    parser.add_argument("--seed", type=int, default=RPI_DISCRIMINATION_SEED)
    args = parser.parse_args(argv)

    reference = rd.load_reference_data(Path(args.reference_input) / "all_hanta")
    precomputed = rd._build_precomputed(reference)
    panel = build_discrimination_panel(reference, seed=args.seed)

    scored = {name: score_cases(reference, precomputed, cases)
             for name, cases in panel.items()}
    summary = summarise_groups(scored)

    print(json.dumps(summary, indent=2, default=str))

    if args.json:
        payload = {"seed": args.seed, "summary": summary, "groups": scored}
        Path(args.json).write_text(
            json.dumps(rd._json_safe(payload), indent=2, default=str, allow_nan=False),
            encoding="utf-8")
    return 0


# ---------------------------------------------------------------------------
# Tests
#
# Cover THIS module's own logic (group construction, sympatric-pair-finding,
# disjointness enforcement, de-duplication). rd.score_candidate and its RPI factors
# already have their own extensive test suites in reassortment_detector.py -- these
# tests do not re-verify that scoring is correct, only that this module hands it the
# right candidates and reports its output correctly.
#
# Run: python -m unittest rpi_discrimination
# ---------------------------------------------------------------------------


def _fixture_isolate(
        species: str, key_suffix: str, host_order: str = "Rodentia",
        country: str = "USA", year: str = "2001", s_seq: str = "ATG" + "A" * 30 + "TAA",
        m_seq: str = "ATG" + "A" * 30 + "TAA", l_seq: str = "ATG" + "A" * 30 + "TAA",
        lab: bool = False, vaccine: bool = False) -> "ah.IsolateRecord":
    """Minimal real-shaped IsolateRecord for this module's own tests -- deliberately
    NOT terminally-complete (no genomic panhandle tags), since these tests exercise
    group-construction LOGIC against a synthetic pool passed in directly, never
    `terminally_complete_isolates` itself (that filter has its own tests in
    rpi_calibration.py)."""
    key = f"{species}:{key_suffix}"
    return ah.IsolateRecord(
        isolate_key=key, species_name=species, species_taxid="1",
        species_status="classified", biosample="", isolate=key_suffix, strain="",
        host_name="Test host", host_taxid="1", host_order=host_order,
        host_family="Cricetidae", host_role="reservoir_rodent",
        collection_date=year, year=year, location="", country=country, region="",
        grouping_tier="tier_2_species_isolate", grouping_confidence=1.0,
        s_accession=f"{key}_S", m_accession=f"{key}_M", l_accession=f"{key}_L",
        s_completeness="complete", m_completeness="complete", l_completeness="complete",
        s_nt_length=len(s_seq), m_nt_length=len(m_seq), l_nt_length=len(l_seq),
        s_protein_aa="", m_protein_aa="", l_protein_aa="", concatenated_protein_aa="",
        segments_present=["S", "M", "L"], coding_complete_segments=["S", "M", "L"],
        triplet_status="coding_complete_triplet", has_complete_triplet=True,
        is_lab_derived=lab, is_vaccine_strain=vaccine, has_nss_orf=False,
        metadata_conflict=False, taxonomic_incongruence=False, cross_species_partners=[],
        kmer_cluster_id="", s_protein_accessions=[], m_protein_accessions=[],
        l_protein_accessions=[], all_accessions=[], warnings=[],
    ), {f"{key}_S": s_seq, f"{key}_M": m_seq, f"{key}_L": l_seq}


def _fixture_reference(*isolate_and_seqs) -> "rd.ReferenceData":
    isolates = [pair[0] for pair in isolate_and_seqs]
    sequences: Dict[str, str] = {}
    for _isolate, seqs in isolate_and_seqs:
        sequences.update(seqs)
    genomic_index = {
        acc: ah.FastaRecord(accession=acc, accession_base=acc, header="", sequence=seq,
                            length=len(seq))
        for acc, seq in sequences.items()
    }
    return rd.ReferenceData(isolates=isolates, cds_index={}, genomic_index=genomic_index)


class CountryContinentTest(unittest.TestCase):
    def test_every_country_in_the_usable_pool_is_mapped(self):
        """Guards against the table going silently stale for the pool build_group_a_random
        ACTUALLY sees (usable_pool, i.e. after the calibration-reserved exclusion) -- not
        the broader terminally_complete_isolates() pool, which currently contains 5
        countries (Czech Republic, Ethiopia, Panama, Peru, Switzerland) that happen to
        belong entirely to isolates reserved for calibration and so never reach this
        module's group builders. A new country in a future all_hanta/ update must be
        added to COUNTRY_CONTINENT explicitly (see its own comment on why a KeyError,
        not a permissive default, is the intended failure when one does)."""
        reference = rd.load_reference_data(Path("all_hanta"))
        pool = usable_pool(reference)
        unmapped = sorted({i.country for i in pool if i.country} - set(COUNTRY_CONTINENT))
        self.assertEqual(unmapped, [], f"unmapped countries: {unmapped}")


class ReservedCalibrationIsolateKeysTest(unittest.TestCase):
    def test_loads_the_pinned_42(self):
        keys = reserved_calibration_isolate_keys()
        self.assertEqual(len(keys), 42)

    def test_usable_pool_excludes_every_reserved_key(self):
        reference = rd.load_reference_data(Path("all_hanta"))
        reserved = reserved_calibration_isolate_keys()
        pool = usable_pool(reference, reserved)
        self.assertTrue(reserved.isdisjoint({i.isolate_key for i in pool}))


class FindSympatricSpeciesPairsTest(unittest.TestCase):
    def test_same_host_order_country_overlapping_years_is_sympatric(self):
        a, seq_a = _fixture_isolate("SpA", "i1", country="China", year="2010")
        b, seq_b = _fixture_isolate("SpB", "i1", country="China", year="2010")
        pairs = find_sympatric_species_pairs([a, b])
        self.assertIn(("Rodentia", "China", "SpA", "SpB"), pairs)

    def test_non_overlapping_years_is_not_sympatric(self):
        a, _ = _fixture_isolate("SpA", "i1", country="China", year="2000")
        b, _ = _fixture_isolate("SpB", "i1", country="China", year="2020")
        self.assertEqual(find_sympatric_species_pairs([a, b]), [])

    def test_different_country_is_not_sympatric(self):
        a, _ = _fixture_isolate("SpA", "i1", country="China", year="2010")
        b, _ = _fixture_isolate("SpB", "i1", country="USA", year="2010")
        self.assertEqual(find_sympatric_species_pairs([a, b]), [])

    def test_different_host_order_is_not_sympatric(self):
        a, _ = _fixture_isolate("SpA", "i1", country="China", year="2010",
                                host_order="Rodentia")
        b, _ = _fixture_isolate("SpB", "i1", country="China", year="2010",
                                host_order="Eulipotyphla")
        self.assertEqual(find_sympatric_species_pairs([a, b]), [])

    def test_same_species_is_not_a_pair(self):
        """Sympatry is a CROSS-species question -- two isolates of the SAME species
        overlapping in place and time is just... a species living where it lives."""
        a, _ = _fixture_isolate("SpA", "i1", country="China", year="2010")
        b, _ = _fixture_isolate("SpA", "i2", country="China", year="2010")
        self.assertEqual(find_sympatric_species_pairs([a, b]), [])


class CandidateFromIsolatesTest(unittest.TestCase):
    def test_known_segment_isolates_is_wired_through_for_factor_e(self):
        """Design doc 2026-08-30-rpi-e-geo-isolate-level-fix: this module's own group
        builders already know the real parent isolate per segment, so the CandidateInput
        they hand to score_candidate must carry it -- otherwise Factor E would fall back
        to the nearest-reference proxy even for candidates whose true parents ARE
        known, which is exactly the silent-basis failure the fix's own honesty
        protocol forbids."""
        s_pair = _fixture_isolate("SpA", "s")
        m_pair = _fixture_isolate("SpB", "m")
        l_pair = _fixture_isolate("SpB", "l")
        reference = _fixture_reference(s_pair, m_pair, l_pair)
        sources = {"S": s_pair[0], "M": m_pair[0], "L": l_pair[0]}
        candidate = _candidate_from_isolates("TEST", reference, sources)
        self.assertEqual(candidate.known_segment_isolates, sources)


class BuildGroupARandomTest(unittest.TestCase):
    def test_every_case_has_pairwise_distinct_species_and_at_least_two_continents(self):
        pool_pairs = [
            _fixture_isolate(f"Sp{i}", "i1", country=c)
            for i, c in enumerate(["USA", "China", "Finland", "Guinea", "Venezuela"])
        ]
        isolates = [p[0] for p in pool_pairs]
        reference = _fixture_reference(*pool_pairs)
        cases = build_group_a_random(reference, isolates, target=5)
        self.assertGreater(len(cases), 0)
        for case in cases:
            species = [v[0] for v in case["segment_sources"].values()]
            self.assertEqual(len(set(species)), 3, case["case_id"])
            self.assertGreaterEqual(len(case["continents"]), 2, case["case_id"])

    def test_raises_on_an_unmapped_country(self):
        """Deliberate -- see COUNTRY_CONTINENT's own comment."""
        pool_pairs = [
            _fixture_isolate(f"Sp{i}", "i1", country=c)
            for i, c in enumerate(["USA", "China", "Nowhereland"])
        ]
        isolates = [p[0] for p in pool_pairs]
        reference = _fixture_reference(*pool_pairs)
        with self.assertRaises(KeyError):
            build_group_a_random(reference, isolates, target=5)


class BuildGroupB1WithinSpeciesTest(unittest.TestCase):
    def test_every_case_is_three_distinct_isolates_of_one_species(self):
        pool_pairs = [_fixture_isolate("SpA", f"i{n}") for n in range(5)]
        isolates = [p[0] for p in pool_pairs]
        reference = _fixture_reference(*pool_pairs)
        cases = build_group_b1_within_species(reference, isolates, per_species=3)
        self.assertGreater(len(cases), 0)
        for case in cases:
            species = {v[0] for v in case["segment_sources"].values()}
            self.assertEqual(species, {"SpA"})
            labels = [v[1] for v in case["segment_sources"].values()]
            self.assertEqual(len(set(labels)), 3, case["case_id"])

    def test_species_below_the_minimum_isolate_count_is_skipped(self):
        pool_pairs = [_fixture_isolate("SpTiny", "i1"), _fixture_isolate("SpTiny", "i2")]
        isolates = [p[0] for p in pool_pairs]
        reference = _fixture_reference(*pool_pairs)
        self.assertEqual(build_group_b1_within_species(reference, isolates), [])


class BuildGroupB2SympatricTest(unittest.TestCase):
    def test_builds_both_directions_and_all_three_segments(self):
        pool_pairs = [
            _fixture_isolate("SpA", "i1", country="China", year="2010"),
            _fixture_isolate("SpB", "i1", country="China", year="2010"),
        ]
        isolates = [p[0] for p in pool_pairs]
        reference = _fixture_reference(*pool_pairs)
        cases = build_group_b2_sympatric(reference, isolates)
        directions = {(c["backbone_species"], c["donor_species"]) for c in cases}
        self.assertEqual(directions, {("SpA", "SpB"), ("SpB", "SpA")})
        segments = {c["donated_segment"] for c in cases}
        self.assertEqual(segments, {"S", "M", "L"})
        self.assertEqual(len(cases), 6)  # 2 directions x 3 segments

    def test_no_sympatric_pairs_gives_an_empty_group(self):
        pool_pairs = [
            _fixture_isolate("SpA", "i1", country="China", year="2010"),
            _fixture_isolate("SpB", "i1", country="USA", year="2010"),
        ]
        isolates = [p[0] for p in pool_pairs]
        reference = _fixture_reference(*pool_pairs)
        self.assertEqual(build_group_b2_sympatric(reference, isolates), [])


class LoadSupportedReassortmentCandidatesTest(unittest.TestCase):
    def test_deduplicates_by_isolate_key(self):
        pool_pairs = [_fixture_isolate("SpA", "i1")]
        isolate = pool_pairs[0][0]
        reference = _fixture_reference(*pool_pairs)
        table_path = Path("output") / "_test_scratch_reassortment_table.csv"
        table_path.parent.mkdir(exist_ok=True)
        with open(table_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f, fieldnames=["isolate_key", "reassortment_call", "discordant_segment",
                              "discordant_clade_bootstrap_support"])
            writer.writeheader()
            for _ in range(2):  # the same isolate_key listed twice
                writer.writerow({
                    "isolate_key": isolate.isolate_key,
                    "reassortment_call": "supported_reassortment_candidate",
                    "discordant_segment": "L",
                    "discordant_clade_bootstrap_support": "100.0",
                })
        try:
            cases = load_supported_reassortment_candidates(reference, table_path)
        finally:
            table_path.unlink()
        self.assertEqual(len(cases), 1)

    def test_rows_not_supported_are_excluded(self):
        pool_pairs = [_fixture_isolate("SpA", "i1")]
        isolate = pool_pairs[0][0]
        reference = _fixture_reference(*pool_pairs)
        table_path = Path("output") / "_test_scratch_reassortment_table2.csv"
        table_path.parent.mkdir(exist_ok=True)
        with open(table_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f, fieldnames=["isolate_key", "reassortment_call", "discordant_segment",
                              "discordant_clade_bootstrap_support"])
            writer.writeheader()
            writer.writerow({
                "isolate_key": isolate.isolate_key, "reassortment_call": "concordant",
                "discordant_segment": "", "discordant_clade_bootstrap_support": "",
            })
        try:
            cases = load_supported_reassortment_candidates(reference, table_path)
        finally:
            table_path.unlink()
        self.assertEqual(cases, [])


class SummariseGroupsTest(unittest.TestCase):
    def test_reports_n_and_rpi_generation_stats(self):
        groups = {
            "a": [{"rpi_generation": 0.1}, {"rpi_generation": 0.3}],
            "b": [{"rpi_generation": None}],
        }
        summary = summarise_groups(groups)
        self.assertEqual(summary["a"]["n"], 2)
        self.assertEqual(summary["a"]["n_scored"], 2)
        self.assertAlmostEqual(summary["a"]["rpi_generation_mean"], 0.2)
        self.assertEqual(summary["b"]["n_scored"], 0)
        self.assertIsNone(summary["b"]["rpi_generation_mean"])


if __name__ == "__main__":
    sys.exit(main())
