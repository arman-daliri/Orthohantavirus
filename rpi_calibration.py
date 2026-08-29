# rpi_calibration.py
"""WS3.1 -- the pinned ddG calibration panel behind the M_pack rescaling (decision D1).

Builds a cognate / non-cognate panhandle-duplex panel from real sequences in
`all_hanta/`, folds every duplex through the ViennaRNA Python bindings (D2b), measures
the ddG distribution under the N = 50 nt window (D2), tests it for bimodality, and
writes `data/rpi_calibration.json`.

WHY THIS IS A SEPARATE MODULE. `reassortment_detector.py` must read a *pinned artifact*
rather than recompute tau at scoring time. Recomputing would make every candidate's score
depend on whatever cohort happened to be loaded -- the exact cohort-dependence that got
Option B rejected in D1. Calibration is therefore a distinct, auditable, seeded step whose
output is committed to the repo.

THE CIRCULARITY THIS MODULE EXISTS TO AVOID. tau is derived from the observed ddG span, so
the panel's composition directly sets the scale parameter every downstream M_pack depends
on. If the panel overlapped the WS4 validation groups, tau would be tuned on the validation
set and the WS4 result would be circular. `ws4_reserved_accessions` enforces the
disjointness against the WS4 groups that exist today, and `calibrate` records the isolates
it consumed so the WS4 groups that do not exist yet can exclude them in turn. See
`reserved_for_calibration` in the artifact.

Run:
    python rpi_calibration.py                       # writes data/rpi_calibration.json
    python -m unittest rpi_calibration              # unit tests
"""

import argparse
import json
import math
import random
import sys
import unittest
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

import all_hanta as ah
import reassortment_detector as rd
import validate_reassortment_detector as vd


# ---------------------------------------------------------------------------
# Panel composition
# ---------------------------------------------------------------------------

# Every RNG in this module takes this explicitly; nothing reads global random state.
RPI_CALIBRATION_SEED = 20260829

DEFAULT_ARTIFACT_PATH = Path("data") / "rpi_calibration.json"

# The genus-conserved genomic 5' terminus. Hantavirus genomic RNA begins 5'-UAGUAGUAU...
# and ends with its reverse complement, and it is that complementarity that forms the
# panhandle. A deposit missing these tags is missing its actual termini, so its "terminal
# 50 nt" is arbitrary internal coding sequence -- see has_terminal_panhandle.
HANTAVIRUS_TERMINAL_5P_TAG = "TAGTAGTA"
HANTAVIRUS_TERMINAL_3P_TAG = "TACTACTA"  # == _reverse_complement of the 5' tag

# Matched sampling depth for the non-cognate side (D1 / plan Task 8 Step 2). The
# terminally-complete pool is dominated by two species -- 178 hantanense and 87 seoulense
# against 10 species with a single isolate each. Uncapped, ~75% of every "across species"
# pair would have a hantanense or seoulense partner and tau would be a measurement of
# those two species' termini rather than of the genus. The cap is applied to the focal
# side and the foreign side alike.
CALIBRATION_MAX_ISOLATES_PER_SPECIES = 8

# Foreign partners drawn per focal isolate, each from a DIFFERENT species, so no single
# species can dominate one focal isolate's non-cognate observations either.
CALIBRATION_FOREIGN_PARTNERS_PER_FOCAL = 3

# Species label used for deposits with no accepted ICTV species. Excluded from the panel:
# "non-cognate = across species" is not a statement that can be made about a partner whose
# species is unknown, and 34 such isolates would otherwise enter the foreign pool.
UNCLASSIFIED_SPECIES_LABEL = "unclassified Orthohantavirus"

SEGMENTS = ("S", "M", "L")

# The six ORDERED inter-segment duplexes. M_pack asks whether a DONATED segment's terminus
# pairs with a BACKBONE segment's terminus as well as the cognate one would, and that
# question is asymmetric: S's 5' end against M's 3' end is a different duplex from M's 5'
# end against S's 3' end. Intra-segment duplexes (X == Y, the classical intramolecular
# panhandle) are deliberately excluded -- they are not the quantity M_pack scores, and a
# reassortant never swaps a segment against itself.
SEGMENT_DUPLEXES: Tuple[Tuple[str, str], ...] = tuple(
    (x, y) for x in SEGMENTS for y in SEGMENTS if x != y
)


# ---------------------------------------------------------------------------
# Pre-registered decision rule for D1
# ---------------------------------------------------------------------------
#
# These thresholds are declared HERE, above the measurement, and are asserted by
# BimodalityDecisionRuleTest against synthetic distributions with known answers. The point
# is that "bimodal with a clean gap" is settled before the data is seen. A rule chosen
# after looking at the histogram is not a decision rule, it is a description of the
# histogram, and it would make the A-vs-C branch unfalsifiable.

# A two-component fit must beat the one-component fit on BIC, which already penalises the
# three extra parameters.
BIMODALITY_REQUIRES_BIC_IMPROVEMENT = True

# Neither component may be a negligible outlier cluster. 0.10 is the conventional floor
# below which a mixture component is treated as describing outliers rather than a mode.
BIMODALITY_MIN_COMPONENT_WEIGHT = 0.10

# Component means must be at least this many pooled SDs apart. Two Gaussians of equal
# weight and equal variance produce a density with a genuine valley only above a
# separation of exactly 2 SD; below it the mixture is unimodal however well it fits. This
# is the "clean gap" clause, and it is a mathematical threshold rather than a taste one.
BIMODALITY_MIN_MEAN_SEPARATION_SD = 2.0

# The mixture density itself must have two local maxima on a fine grid. Redundant with the
# clause above for equal-variance components and NOT redundant for unequal-variance ones,
# where a wide component can swallow a narrow one and leave a single peak. Checking the
# density directly is the assertion that actually matters; the SD clause is the
# interpretable summary of it.
BIMODALITY_DENSITY_GRID_POINTS = 2001

# Boltzmann branch: the dynamic range tau is calibrated to span. D = 100 matches
# host_decay, the widest existing RPI factor (spec section "Option A").
RPI_MPACK_TARGET_DYNAMIC_RANGE = 100.0

# The span percentiles for the SHIPPED tau. See derive_tau for why the raw min-max span is
# reported but not used.
TAU_SPAN_LOW_PERCENTILE = 5.0
TAU_SPAN_HIGH_PERCENTILE = 95.0


# ---------------------------------------------------------------------------
# Sequence access and the terminal-completeness filter
# ---------------------------------------------------------------------------


def segment_accession(isolate: "ah.IsolateRecord", segment: str) -> str:
    if segment not in SEGMENTS:
        raise ValueError(f"segment must be one of {SEGMENTS}, got {segment!r}")
    return {"S": isolate.s_accession, "M": isolate.m_accession,
            "L": isolate.l_accession}[segment]


def genomic_sequence(reference: "rd.ReferenceData", accession: str) -> str:
    """The genomic sequence for `accession`, or "" if the index does not hold it.

    Mirrors PanelBuilder.genomic's versioned/unversioned lookup rather than reimplementing
    it: the isolate records carry versioned accessions and the FASTA index is keyed both
    ways depending on the deposit.
    """
    if not accession:
        return ""
    record = (reference.genomic_index.get(accession)
              or reference.genomic_index.get(accession.split(".")[0]))
    return record.sequence if record else ""


def has_terminal_panhandle(sequence: str) -> bool:
    """Does this deposit actually carry both genuine termini?

    THIS FILTER IS THE DIFFERENCE BETWEEN MEASURING A PANHANDLE AND MEASURING NOISE, and
    it is not in the plan -- it comes from the data. Of the 2,622 full-length segments in
    the wild complete-triplet pool, only 1,190 (45%) both start with TAGTAGTA and end with
    TACTACTA; 1,261 (48%) carry neither tag. GenBank deposits are routinely trimmed of
    their extreme termini, and `SEGMENT_LENGTH_RANGE` does not catch it because a segment
    missing 30 nt at each end is still comfortably inside the ICTV band.

    Measured on one such truncated deposit, the 5' x 3' duplex folded to -10.0 kcal/mol
    against -42.2 on a terminally complete one: the truncated "panhandle window" is
    arbitrary internal coding sequence that happens to be 50 nt from an end. Folding those
    would put a large, purely technical variance component into both the cognate and the
    non-cognate distribution -- and since tau is derived from the ddG SPAN, inflating the
    span with sequencing artefacts would directly and silently inflate tau.
    """
    sequence = sequence.strip().upper()
    return (sequence.startswith(HANTAVIRUS_TERMINAL_5P_TAG)
            and sequence.endswith(HANTAVIRUS_TERMINAL_3P_TAG))


def _is_wild(isolate: "ah.IsolateRecord") -> bool:
    return not isolate.is_lab_derived and not isolate.is_vaccine_strain


def _segment_is_full_length(reference: "rd.ReferenceData",
                            isolate: "ah.IsolateRecord", segment: str) -> bool:
    sequence = genomic_sequence(reference, segment_accession(isolate, segment))
    low, high = ah.SEGMENT_LENGTH_RANGE[segment]
    return bool(sequence) and low <= len(sequence) <= high


def terminally_complete_isolates(
        reference: "rd.ReferenceData") -> List["ah.IsolateRecord"]:
    """Wild, classified isolates whose S, M and L are all full-length AND terminal.

    Returned in `isolate_key` order so the pool is deterministic regardless of the order
    `load_reference_data` happens to yield -- the same class of set-ordering
    nondeterminism that WS0 fixed in gate1_check.
    """
    pool = []
    for isolate in reference.isolates:
        if not _is_wild(isolate):
            continue
        if isolate.species_name == UNCLASSIFIED_SPECIES_LABEL or not isolate.species_name:
            continue
        if not all(_segment_is_full_length(reference, isolate, s) for s in SEGMENTS):
            continue
        if not all(has_terminal_panhandle(
                genomic_sequence(reference, segment_accession(isolate, s)))
                for s in SEGMENTS):
            continue
        pool.append(isolate)
    return sorted(pool, key=lambda i: i.isolate_key)


# ---------------------------------------------------------------------------
# Disjointness from the WS4 validation groups
# ---------------------------------------------------------------------------


def ws4_reserved_accessions(reference: "rd.ReferenceData") -> Set[str]:
    """Every accession consumed by the WS4 validation panel as it stands today.

    Sourced by BUILDING THE ACTUAL PANEL rather than by re-deriving the specs, per D11's
    single-source-of-truth rule: `build_panel` picks backbones and donors through
    `find_backbone`/`find_donor`, whose tie-breaks (`hinted[0]`, `ok[0]`) are order
    dependent, and a reimplementation here could silently pick a different isolate and
    hand back a disjointness guarantee that was not true.

    Returns accessions, not isolate keys, because the panel records accessions and an
    accession identifies the sequence a duplex would actually be folded from.
    """
    builder = vd.PanelBuilder(reference)
    panel, _skipped = vd.build_panel(builder)
    reserved: Set[str] = set()
    for case in panel:
        for accession in (case.get("backbone_accessions") or {}).values():
            if accession:
                reserved.add(accession)
        if case.get("donor_accession"):
            reserved.add(case["donor_accession"])
    return reserved


def _isolate_touches(reference: "rd.ReferenceData", isolate: "ah.IsolateRecord",
                     accessions: Set[str]) -> bool:
    """Is any of this isolate's three segments one of `accessions`?

    Checked segment by segment rather than by isolate identity because an isolate can
    contribute a single donated segment to a WS4 case without being its backbone. Sharing
    even one segment is enough to break disjointness: that segment's terminal window would
    appear in both the calibration panel and the validation set.
    """
    for segment in SEGMENTS:
        accession = segment_accession(isolate, segment)
        if not accession:
            continue
        if accession in accessions or accession.split(".")[0] in accessions:
            return True
    return False


def select_calibration_isolates(
        reference: "rd.ReferenceData",
        max_per_species: int = CALIBRATION_MAX_ISOLATES_PER_SPECIES,
        seed: int = RPI_CALIBRATION_SEED,
) -> Tuple[List["ah.IsolateRecord"], Dict[str, Any]]:
    """The calibration pool: terminally complete, disjoint from WS4, depth-matched."""
    pool = terminally_complete_isolates(reference)
    reserved = ws4_reserved_accessions(reference)
    eligible = [i for i in pool if not _isolate_touches(reference, i, reserved)]

    by_species: Dict[str, List["ah.IsolateRecord"]] = defaultdict(list)
    for isolate in eligible:
        by_species[isolate.species_name].append(isolate)

    rng = random.Random(seed)
    selected: List["ah.IsolateRecord"] = []
    for species in sorted(by_species):
        members = sorted(by_species[species], key=lambda i: i.isolate_key)
        if len(members) > max_per_species:
            members = rng.sample(members, max_per_species)
        selected.extend(members)
    selected.sort(key=lambda i: i.isolate_key)

    provenance = {
        "wild_full_length_terminally_complete": len(pool),
        "ws4_reserved_accessions": len(reserved),
        "excluded_as_ws4_reserved": len(pool) - len(eligible),
        "eligible_after_ws4_exclusion": len(eligible),
        "max_isolates_per_species": max_per_species,
        "selected": len(selected),
        "species_selected": len({i.species_name for i in selected}),
        "isolates_per_species": dict(
            sorted(Counter(i.species_name for i in selected).items())),
    }
    return selected, provenance


# ---------------------------------------------------------------------------
# The duplex panel
# ---------------------------------------------------------------------------


def duplex_dg(reference: "rd.ReferenceData",
              five_prime_isolate: "ah.IsolateRecord", five_prime_segment: str,
              three_prime_isolate: "ah.IsolateRecord", three_prime_segment: str,
              window_nt: int) -> float:
    """dG of the (5' window of one segment) x (3' window of another) duplex.

    NEITHER WINDOW IS REVERSE-COMPLEMENTED, and that is a correction to the plan's Task 8
    Step 1, which specified the 3' window "reverse-complemented". Verified two ways:

      1. Convention control. RNAcofold takes both strands 5'->3' and pairs them
         ANTIPARALLEL, so a duplex forms between sequences that are reverse complements of
         each other as written: cofold_dg(a, revcomp(a)) = -20.2 where cofold_dg(a, a) =
         -16.6 and the parallel arrangement gives -7.9.
      2. The genomic 5' and 3' terminal windows, both written 5'->3' as the FASTA stores
         them, ALREADY are reverse complements of each other -- that is what a panhandle
         is. Reverse-complementing the 3' window therefore destroys the duplex rather than
         forming it. Measured on a terminally complete isolate across all nine segment
         pairings: as-is -30.9 to -45.3, reverse-complemented -7.3 to -15.4.

    Reverse-complementing would have scored noise and, because tau is derived from the
    resulting span, would have propagated silently into every downstream M_pack.
    """
    five = rd.panhandle_window(
        genomic_sequence(reference, segment_accession(five_prime_isolate, five_prime_segment)),
        "5", window_nt)
    three = rd.panhandle_window(
        genomic_sequence(reference, segment_accession(three_prime_isolate, three_prime_segment)),
        "3", window_nt)
    return rd.cofold_dg(five, three)


def build_panel(
        reference: "rd.ReferenceData",
        window_nt: int = rd.RPI_PANHANDLE_WINDOW_NT,
        seed: int = RPI_CALIBRATION_SEED,
        max_per_species: int = CALIBRATION_MAX_ISOLATES_PER_SPECIES,
        foreign_partners: int = CALIBRATION_FOREIGN_PARTNERS_PER_FOCAL,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Matched cognate / non-cognate observations, one ddG each.

    THE PAIRING IS MATCHED, not two independent samples. Each observation holds the 5'
    window of ONE focal segment fixed and swaps only the 3' partner: cognate is that focal
    isolate's own partner segment, non-cognate is a different SPECIES' same partner
    segment. ddG = dG_non_cognate - dG_cognate is therefore a like-for-like difference on a
    shared anchor, which is exactly the quantity `boltzmann_compatibility` consumes and
    exactly what M_pack will ask at scoring time (D3).

    Two independent pools -- all cognate dGs against all non-cognate dGs -- would have
    confounded the contrast with per-segment GC content, since S, M and L termini differ
    systematically and an unmatched difference would partly measure that.

    COGNATE = WITHIN-ISOLATE (approved 2026-08-29): S/M/L from the same isolate, segments
    that demonstrably travelled together in one virion. Deliberately NOT within-species,
    which is what the WS4 b1 group tests.
    """
    selected, pool_provenance = select_calibration_isolates(
        reference, max_per_species=max_per_species, seed=seed)
    if len(selected) < 2:
        raise RuntimeError(
            f"calibration pool has {len(selected)} isolate(s) -- cannot build "
            "cognate/non-cognate pairs")

    by_species: Dict[str, List["ah.IsolateRecord"]] = defaultdict(list)
    for isolate in selected:
        by_species[isolate.species_name].append(isolate)
    all_species = sorted(by_species)

    rng = random.Random(seed + 1)
    observations: List[Dict[str, Any]] = []
    cognate_cache: Dict[Tuple[str, str, str], float] = {}

    for focal in selected:
        # One foreign partner per foreign SPECIES, so a species with eight selected
        # isolates gets no more weight here than one with a single isolate.
        other_species = [s for s in all_species if s != focal.species_name]
        if not other_species:
            continue
        chosen_species = (rng.sample(other_species, foreign_partners)
                          if len(other_species) > foreign_partners else list(other_species))
        partners = [rng.choice(sorted(by_species[s], key=lambda i: i.isolate_key))
                    for s in sorted(chosen_species)]

        for five_segment, three_segment in SEGMENT_DUPLEXES:
            key = (focal.isolate_key, five_segment, three_segment)
            if key not in cognate_cache:
                cognate_cache[key] = duplex_dg(
                    reference, focal, five_segment, focal, three_segment, window_nt)
            dg_cognate = cognate_cache[key]

            for partner in partners:
                dg_non_cognate = duplex_dg(
                    reference, focal, five_segment, partner, three_segment, window_nt)
                observations.append({
                    "focal_isolate_key": focal.isolate_key,
                    "focal_species": focal.species_name,
                    "focal_accession": segment_accession(focal, five_segment),
                    "partner_isolate_key": partner.isolate_key,
                    "partner_species": partner.species_name,
                    "partner_accession": segment_accession(partner, three_segment),
                    "five_prime_segment": five_segment,
                    "three_prime_segment": three_segment,
                    "dg_cognate": dg_cognate,
                    "dg_non_cognate": dg_non_cognate,
                    "ddg": dg_non_cognate - dg_cognate,
                })

    panel_provenance = dict(pool_provenance)
    panel_provenance.update({
        "window_nt": window_nt,
        "foreign_partners_per_focal": foreign_partners,
        "segment_duplexes": ["%s5'x%s3'" % pair for pair in SEGMENT_DUPLEXES],
        "cognate_folds": len(cognate_cache),
        "non_cognate_folds": len(observations),
        "total_folds": len(cognate_cache) + len(observations),
        "observations": len(observations),
    })
    return observations, panel_provenance


# ---------------------------------------------------------------------------
# Distribution statistics
# ---------------------------------------------------------------------------


def _summarise(values: Sequence[float]) -> Dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "n": int(array.size),
        "mean": float(array.mean()),
        "sd": float(array.std(ddof=1)) if array.size > 1 else 0.0,
        "min": float(array.min()),
        "p5": float(np.percentile(array, 5)),
        "p25": float(np.percentile(array, 25)),
        "median": float(np.median(array)),
        "p75": float(np.percentile(array, 75)),
        "p95": float(np.percentile(array, 95)),
        "max": float(array.max()),
    }


def bimodality_coefficient(values: Sequence[float]) -> float:
    """Sarle's bimodality coefficient, BC = (skew^2 + 1) / (kurt + 3(n-1)^2/((n-2)(n-3))).

    Above 5/9 ~= 0.5556 -- the value of a uniform distribution -- is the conventional
    flag for bimodality. Reported as a cheap, assumption-light second opinion alongside
    the mixture fit; it is NOT part of the decision rule, because BC also rises on skewed
    unimodal distributions and would let a long tail masquerade as a second mode.
    """
    array = np.asarray(values, dtype=float)
    n = array.size
    if n < 4:
        raise ValueError(f"bimodality coefficient needs n >= 4, got {n}")
    mean, sd = array.mean(), array.std(ddof=1)
    if sd == 0:
        raise ValueError("bimodality coefficient is undefined for zero variance")
    z = (array - mean) / sd
    skew = float(np.mean(z ** 3)) * math.sqrt(n * (n - 1)) / (n - 2)
    kurt = ((n - 1) * ((n + 1) * float(np.mean(z ** 4)) - 3 * (n - 1))
            / ((n - 2) * (n - 3)))
    return (skew ** 2 + 1.0) / (kurt + 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3)))


def _normal_logpdf(x: np.ndarray, mean: float, sd: float) -> np.ndarray:
    return -0.5 * np.log(2 * np.pi * sd ** 2) - 0.5 * ((x - mean) / sd) ** 2


def fit_two_component_gaussian(values: Sequence[float], max_iter: int = 1000,
                               tol: float = 1e-10) -> Dict[str, float]:
    """1-D two-component Gaussian mixture by EM.

    Hand-rolled rather than sklearn's GaussianMixture because scipy and sklearn are BOTH
    absent from this environment (numpy 2.5.1, pandas 3.0.5, biopython, dendropy,
    viennarna -- see requirements-pip.txt). Adding a dependency to run one 1-D EM would
    have to clear EnvironmentReconciliationTest and the D12 install route on two
    platforms, for about thirty lines of arithmetic.

    Initialised deterministically at the 25th and 75th percentiles -- no RNG, so the fit
    is reproducible without threading yet another seed through.
    """
    x = np.asarray(values, dtype=float)
    if x.size < 4:
        raise ValueError(f"mixture fit needs n >= 4, got {x.size}")
    # Variance floor: a component that collapses onto a single point drives its sd to 0
    # and the log-likelihood to +inf, which would win on BIC while describing nothing.
    variance_floor = max(float(np.var(x)) * 1e-6, 1e-12)

    mu = np.array([np.percentile(x, 25), np.percentile(x, 75)], dtype=float)
    if mu[0] == mu[1]:
        mu = np.array([x.min(), x.max()], dtype=float)
    sd = np.array([max(float(np.std(x)), math.sqrt(variance_floor))] * 2, dtype=float)
    weight = np.array([0.5, 0.5], dtype=float)

    previous = -np.inf
    loglik = -np.inf
    for _ in range(max_iter):
        # E step, in log space so a far-out point cannot underflow both components to 0
        # and produce a 0/0 responsibility.
        log_component = np.vstack([
            np.log(weight[k]) + _normal_logpdf(x, mu[k], sd[k]) for k in (0, 1)])
        peak = log_component.max(axis=0)
        total = peak + np.log(np.exp(log_component - peak).sum(axis=0))
        responsibility = np.exp(log_component - total)
        loglik = float(total.sum())

        # M step
        mass = responsibility.sum(axis=1)
        weight = mass / x.size
        mu = (responsibility @ x) / mass
        for k in (0, 1):
            variance = float((responsibility[k] * (x - mu[k]) ** 2).sum() / mass[k])
            sd[k] = math.sqrt(max(variance, variance_floor))

        if abs(loglik - previous) < tol:
            break
        previous = loglik

    order = np.argsort(mu)
    mu, sd, weight = mu[order], sd[order], weight[order]
    # 5 free parameters: two means, two SDs, one weight (the second is 1 - the first).
    bic = -2.0 * loglik + 5.0 * math.log(x.size)
    return {
        "mean_low": float(mu[0]), "mean_high": float(mu[1]),
        "sd_low": float(sd[0]), "sd_high": float(sd[1]),
        "weight_low": float(weight[0]), "weight_high": float(weight[1]),
        "loglik": loglik, "bic": float(bic),
    }


def fit_one_component_gaussian(values: Sequence[float]) -> Dict[str, float]:
    x = np.asarray(values, dtype=float)
    mean, sd = float(x.mean()), float(x.std())
    if sd <= 0:
        raise ValueError("one-component fit is undefined for zero variance")
    loglik = float(_normal_logpdf(x, mean, sd).sum())
    return {"mean": mean, "sd": sd, "loglik": loglik,
            "bic": float(-2.0 * loglik + 2.0 * math.log(x.size))}


def _mixture_density_modes(fit: Dict[str, float], low: float, high: float,
                           points: int = BIMODALITY_DENSITY_GRID_POINTS
                           ) -> Tuple[int, Optional[float]]:
    """Interior local maxima of the fitted mixture density, and the valley between them.

    Checks the DENSITY rather than inferring bimodality from the parameters, because two
    components of very unequal width produce a single-peaked density however far apart
    their means are -- the wide one simply swallows the narrow one.
    """
    grid = np.linspace(low, high, points)
    density = (fit["weight_low"] * np.exp(_normal_logpdf(grid, fit["mean_low"], fit["sd_low"]))
               + fit["weight_high"] * np.exp(
                   _normal_logpdf(grid, fit["mean_high"], fit["sd_high"])))
    peaks = [i for i in range(1, points - 1)
             if density[i] > density[i - 1] and density[i] >= density[i + 1]]
    if len(peaks) < 2:
        return len(peaks), None
    first, last = peaks[0], peaks[-1]
    valley = first + int(np.argmin(density[first:last + 1]))
    return len(peaks), float(grid[valley])


def assess_bimodality(values: Sequence[float]) -> Dict[str, Any]:
    """Apply the pre-registered D1 rule. Returns the verdict AND every input to it."""
    x = np.asarray(values, dtype=float)
    one = fit_one_component_gaussian(x)
    two = fit_two_component_gaussian(x)

    pooled_sd = math.sqrt(two["weight_low"] * two["sd_low"] ** 2
                          + two["weight_high"] * two["sd_high"] ** 2)
    separation_sd = (abs(two["mean_high"] - two["mean_low"]) / pooled_sd
                     if pooled_sd > 0 else float("inf"))
    min_weight = min(two["weight_low"], two["weight_high"])
    modes, valley = _mixture_density_modes(two, float(x.min()), float(x.max()))

    checks = {
        "bic_improves": bool(two["bic"] < one["bic"]),
        "both_components_substantial": bool(min_weight >= BIMODALITY_MIN_COMPONENT_WEIGHT),
        "means_separated": bool(separation_sd >= BIMODALITY_MIN_MEAN_SEPARATION_SD),
        "density_has_two_modes": bool(modes >= 2),
    }
    return {
        "rule": {
            "requires_bic_improvement": BIMODALITY_REQUIRES_BIC_IMPROVEMENT,
            "min_component_weight": BIMODALITY_MIN_COMPONENT_WEIGHT,
            "min_mean_separation_sd": BIMODALITY_MIN_MEAN_SEPARATION_SD,
            "pre_registered": "declared in rpi_calibration.py above the measurement; "
                              "asserted by BimodalityDecisionRuleTest",
        },
        "one_component": one,
        "two_component": two,
        "pooled_sd": pooled_sd,
        "mean_separation_sd": separation_sd,
        "min_component_weight": min_weight,
        "bic_delta_two_minus_one": two["bic"] - one["bic"],
        "density_local_maxima": modes,
        "density_valley_ddg": valley,
        "bimodality_coefficient": bimodality_coefficient(x),
        "bimodality_coefficient_threshold": 5.0 / 9.0,
        "checks": checks,
        "bimodal_with_clean_gap": bool(all(checks.values())),
    }


# ---------------------------------------------------------------------------
# The D1 branch
# ---------------------------------------------------------------------------


def derive_tau(values: Sequence[float],
               dynamic_range: float = RPI_MPACK_TARGET_DYNAMIC_RANGE) -> Dict[str, float]:
    """tau = ddG_span / ln(D), with span taken robustly.

    DEVIATION FROM THE LITERAL D1 FORMULA, FLAGGED RATHER THAN TAKEN SILENTLY. D1 writes
    `tau = ddG_span / ln(D)` without defining the span. The raw min-max span is a function
    of the two most extreme folds in the panel, so a single outlying duplex moves tau and
    therefore moves every M_pack in the project -- which is the same cohort-sensitivity
    that got Option B rejected. The p5-p95 span is reported as `tau` and the raw min-max
    span is reported alongside it as `tau_raw_span`, so the choice is visible in the
    artifact and reversible by changing one constant.
    """
    x = np.asarray(values, dtype=float)
    raw_span = float(x.max() - x.min())
    robust_span = float(np.percentile(x, TAU_SPAN_HIGH_PERCENTILE)
                        - np.percentile(x, TAU_SPAN_LOW_PERCENTILE))
    log_range = math.log(dynamic_range)
    return {
        "tau": robust_span / log_range,
        "tau_raw_span": raw_span / log_range,
        "ddg_span_robust_p5_p95": robust_span,
        "ddg_span_raw_min_max": raw_span,
        "target_dynamic_range": dynamic_range,
        "span_percentiles": [TAU_SPAN_LOW_PERCENTILE, TAU_SPAN_HIGH_PERCENTILE],
    }


def cliffs_delta(a: Sequence[float], b: Sequence[float]) -> float:
    """P(a > b) - P(a < b), in [-1, 1]. Rank-based, so no distributional assumption.

    Computed by sorting rather than by the naive O(n*m) double loop: the panel compares
    every cognate dG against every non-cognate dG, which is ~190,000 comparisons here and
    grows quadratically with the species cap.
    """
    x = np.sort(np.asarray(a, dtype=float))
    y = np.asarray(b, dtype=float)
    if x.size == 0 or y.size == 0:
        raise ValueError("cliffs_delta needs two non-empty samples")
    # For each y, how many x are strictly below / at-or-below it.
    below = np.searchsorted(x, y, side="left")
    at_or_below = np.searchsorted(x, y, side="right")
    greater = x.size - at_or_below          # count of x > y
    less = below                            # count of x < y
    return float((greater.sum() - less.sum()) / (x.size * y.size))


def distinct_cognate_dgs(observations: Sequence[Dict[str, Any]]) -> List[float]:
    """The cognate dGs, deduplicated.

    Each cognate fold is reused across every foreign partner of its focal isolate, so the
    observation list holds each cognate dG `foreign_partners` times. Summarising the raw
    column would report an n inflated by that factor and an SD narrowed by it.
    """
    return [dg for _key, _five, _three, dg in sorted(
        {(o["focal_isolate_key"], o["five_prime_segment"],
          o["three_prime_segment"], o["dg_cognate"]) for o in observations})]


def assess_discrimination(observations: Sequence[Dict[str, Any]],
                          rescaling: Dict[str, Any]) -> Dict[str, Any]:
    """Does ddG actually separate cognate from non-cognate, and what does M_pack do?

    THIS IS THE EVIDENCE WS3.3 NEEDS AND THE PLAN DID NOT ASK FOR. Deriving tau establishes
    a scale; it does not establish that there is anything to scale. The RPI's central
    finding (PROJECT-STATE section 3) is that only M_pack and M_pol can see which
    particular S+M+L assembly is being scored, so if ddG does not separate here, M_pack
    will not discriminate cross-species shuffles no matter how it is rescaled -- and that
    has to be visible before WS3.2 is built on top of it, not discovered in WS4.

    `fraction_mpack_saturated` is the number that matters most: ddG <= 0 means the
    non-cognate duplex is at least as strong as the cognate one, M_pack clamps to exactly
    1.0, and that observation contributes no signal at all.
    """
    ddg = np.array([o["ddg"] for o in observations], dtype=float)
    non_cognate = np.array([o["dg_non_cognate"] for o in observations], dtype=float)
    cognate_dg = np.array(distinct_cognate_dgs(observations), dtype=float)

    if rescaling["form"] == "tau_scaled_boltzmann":
        tau = rescaling["parameters"]["tau"]
        m_pack = np.array([rd.boltzmann_compatibility(float(d), 0.0, rt=tau) for d in ddg])
    else:
        p = rescaling["parameters"]
        m_pack = 1.0 / (1.0 + np.exp((ddg - p["ddg_50"]) / p["s"]))

    return {
        "fraction_ddg_positive": float((ddg > 0).mean()),
        "fraction_mpack_saturated": float((ddg <= 0).mean()),
        "saturated_meaning": "ddG <= 0: the non-cognate duplex is at least as strong as "
                             "the cognate one, so M_pack clamps to 1.0 and the "
                             "observation carries no packaging signal",
        "cliffs_delta_noncognate_vs_cognate": cliffs_delta(non_cognate, cognate_dg),
        "cliffs_delta_reading": "P(non-cognate dG > cognate dG) - P(<). Positive means "
                                "non-cognate duplexes are WEAKER (less negative), the "
                                "direction M_pack assumes. |d| < 0.147 negligible, "
                                "< 0.33 small, < 0.474 medium (Romano et al. 2006).",
        "m_pack_under_chosen_form": _summarise(m_pack),
    }


def derive_logistic(assessment: Dict[str, Any]) -> Dict[str, float]:
    """M_pack = 1 / (1 + exp((ddG - ddG_50) / s)), parameterised from the measured gap.

    ddG_50 sits at the density valley -- the threshold the data actually shows, which is
    the whole reason Option C would be preferred. s is set so the logistic's 10%-90%
    transition spans the gap between the component means: that width is 2*ln(9)*s.
    """
    valley = assessment.get("density_valley_ddg")
    two = assessment["two_component"]
    ddg_50 = float(valley) if valley is not None else 0.5 * (two["mean_low"] + two["mean_high"])
    gap = abs(two["mean_high"] - two["mean_low"])
    return {"ddg_50": ddg_50, "s": gap / (2.0 * math.log(9.0)), "gap_width": gap}


def choose_rescaling(values: Sequence[float],
                     assessment: Dict[str, Any]) -> Dict[str, Any]:
    if assessment["bimodal_with_clean_gap"]:
        return {
            "form": "logistic",
            "option": "C",
            "formula": "M_pack = 1 / (1 + exp((ddG - ddg_50) / s))",
            "parameters": derive_logistic(assessment),
            "reason": "the pre-registered D1 rule passed on all four checks: the ddG "
                      "distribution is bimodal with a clean gap, so the sigmoid describes "
                      "a real threshold rather than imposing one",
        }
    failed = sorted(k for k, ok in assessment["checks"].items() if not ok)
    return {
        "form": "tau_scaled_boltzmann",
        "option": "A",
        "formula": "M_pack = exp(-ddG / tau)",
        "parameters": derive_tau(values),
        "reason": "the pre-registered D1 rule failed on: " + ", ".join(failed)
                  + " -- the ddG distribution is not bimodal with a clean gap, so a "
                    "logistic would impose a threshold the data does not show",
    }


# ---------------------------------------------------------------------------
# Artifact
# ---------------------------------------------------------------------------

RPI_CALIBRATION_SCHEMA_VERSION = 1


def calibrate(reference: "rd.ReferenceData",
              window_nt: int = rd.RPI_PANHANDLE_WINDOW_NT,
              seed: int = RPI_CALIBRATION_SEED,
              max_per_species: int = CALIBRATION_MAX_ISOLATES_PER_SPECIES,
              foreign_partners: int = CALIBRATION_FOREIGN_PARTNERS_PER_FOCAL,
              ) -> Dict[str, Any]:
    """Build the panel, measure the ddG distribution, and apply the D1 branch."""
    observations, panel_provenance = build_panel(
        reference, window_nt=window_nt, seed=seed,
        max_per_species=max_per_species, foreign_partners=foreign_partners)

    ddg = [o["ddg"] for o in observations]
    assessment = assess_bimodality(ddg)
    rescaling = choose_rescaling(ddg, assessment)

    # The isolates this panel consumed. WS4's b1/b2/random group builders must exclude
    # these, which is the half of the disjointness guarantee that cannot be enforced from
    # here -- those groups do not exist yet.
    consumed = sorted({o["focal_isolate_key"] for o in observations}
                      | {o["partner_isolate_key"] for o in observations})

    return {
        "schema_version": RPI_CALIBRATION_SCHEMA_VERSION,
        "decision": "D1 -- M_pack rescaling, branched on the measured ddG distribution",
        "seed": seed,
        "window_nt": window_nt,
        "fold_provenance": rd.viennarna_provenance(),
        "panel": panel_provenance,
        "dg_cognate": _summarise(distinct_cognate_dgs(observations)),
        "dg_non_cognate": _summarise([o["dg_non_cognate"] for o in observations]),
        "ddg": _summarise(ddg),
        "bimodality": assessment,
        "rescaling": rescaling,
        "discrimination": assess_discrimination(observations, rescaling),
        "reserved_for_calibration": {
            "note": "WS4's b1/b2/random group builders MUST exclude these isolate keys. "
                    "tau is derived from this panel's ddG span, so any overlap would tune "
                    "the scale parameter on the validation set.",
            "isolate_keys": consumed,
            "count": len(consumed),
        },
        "provisional": {
            "tau_4_86_superseded": "D13's provisional tau ~= 4.86 was derived from a "
                                   "ddG span of 22.4 measured before the N = 50 window "
                                   "was chosen. This artifact supersedes it.",
        },
    }


def write_artifact(result: Dict[str, Any], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(rd._json_safe(result), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    return path


def _report(result: Dict[str, Any]) -> str:
    lines = []
    panel, ddg, bim, res = (result["panel"], result["ddg"],
                            result["bimodality"], result["rescaling"])
    lines.append("=== WS3.1 calibration panel ===")
    lines.append(f"  mechanism      : {result['fold_provenance']['mechanism']} "
                 f"{result['fold_provenance']['version']}")
    lines.append(f"  window         : {result['window_nt']} nt    seed: {result['seed']}")
    lines.append(f"  pool           : {panel['wild_full_length_terminally_complete']} "
                 f"terminally complete -> {panel['eligible_after_ws4_exclusion']} after "
                 f"WS4 exclusion ({panel['excluded_as_ws4_reserved']} reserved)")
    lines.append(f"  selected       : {panel['selected']} isolates / "
                 f"{panel['species_selected']} species (cap "
                 f"{panel['max_isolates_per_species']}/species)")
    lines.append(f"  folds          : {panel['total_folds']} "
                 f"({panel['cognate_folds']} cognate, {panel['non_cognate_folds']} non-cognate)")
    for name in ("dg_cognate", "dg_non_cognate", "ddg"):
        s = result[name]
        lines.append(f"  {name:<14}: n={s['n']} mean={s['mean']:.2f} sd={s['sd']:.2f} "
                     f"min={s['min']:.2f} p5={s['p5']:.2f} med={s['median']:.2f} "
                     f"p95={s['p95']:.2f} max={s['max']:.2f}")
    lines.append("")
    lines.append("=== D1 bimodality rule ===")
    for check, ok in sorted(bim["checks"].items()):
        lines.append(f"  [{'PASS' if ok else 'FAIL'}] {check}")
    lines.append(f"  BIC 1-comp={bim['one_component']['bic']:.1f}  "
                 f"2-comp={bim['two_component']['bic']:.1f}  "
                 f"delta={bim['bic_delta_two_minus_one']:+.1f}")
    lines.append(f"  components: mu={bim['two_component']['mean_low']:.2f}"
                 f"/{bim['two_component']['mean_high']:.2f}  "
                 f"sd={bim['two_component']['sd_low']:.2f}"
                 f"/{bim['two_component']['sd_high']:.2f}  "
                 f"w={bim['two_component']['weight_low']:.3f}"
                 f"/{bim['two_component']['weight_high']:.3f}")
    lines.append(f"  separation={bim['mean_separation_sd']:.2f} SD "
                 f"(need >= {BIMODALITY_MIN_MEAN_SEPARATION_SD})   "
                 f"density modes={bim['density_local_maxima']}   "
                 f"BC={bim['bimodality_coefficient']:.3f} "
                 f"(flag > {bim['bimodality_coefficient_threshold']:.3f})")
    lines.append("")
    lines.append(f"=== BRANCH: Option {res['option']} -- {res['form']} ===")
    lines.append(f"  {res['formula']}")
    for key, value in sorted(res["parameters"].items()):
        lines.append(f"    {key} = {value}" if not isinstance(value, float)
                     else f"    {key} = {value:.4f}")
    lines.append(f"  reason: {res['reason']}")
    disc = result.get("discrimination")
    if disc:
        lines.append("")
        lines.append("=== Does ddG carry signal? (evidence for WS3.3) ===")
        lines.append(f"  ddG > 0 (non-cognate weaker, expected direction): "
                     f"{disc['fraction_ddg_positive']:.1%}")
        lines.append(f"  M_pack saturated at 1.0 (ddG <= 0, no signal)   : "
                     f"{disc['fraction_mpack_saturated']:.1%}")
        lines.append(f"  Cliff's delta, non-cognate vs cognate dG        : "
                     f"{disc['cliffs_delta_noncognate_vs_cognate']:+.3f}")
        m = disc["m_pack_under_chosen_form"]
        lines.append(f"  M_pack under the chosen form: mean={m['mean']:.3f} "
                     f"med={m['median']:.3f} p5={m['p5']:.3f} p95={m['p95']:.3f} "
                     f"min={m['min']:.4f}")
    return "\n".join(lines)


STABILITY_SEEDS = (RPI_CALIBRATION_SEED, 1, 999)
STABILITY_WINDOWS = (25, 50, 75)


def stability_sweep(reference: "rd.ReferenceData",
                    seeds: Sequence[int] = STABILITY_SEEDS,
                    windows: Sequence[int] = STABILITY_WINDOWS,
                    ) -> List[Dict[str, Any]]:
    """Re-run the branch across seeds and the D2 window sweep.

    A branch that flips with the RNG seed is a sampling artifact, not a measurement, and
    tau is only quotable if it survives the panel being resampled. Kept as an on-demand
    CLI mode rather than baked into the artifact: the artifact must describe ONE pinned
    panel, and a stability table stored beside it would go stale the moment either the
    dataset or the selection logic changed.
    """
    rows = []
    for window in windows:
        for seed in seeds:
            result = calibrate(reference, window_nt=window, seed=seed)
            bimodality = result["bimodality"]
            parameters = result["rescaling"]["parameters"]
            rows.append({
                "window_nt": window, "seed": seed,
                "option": result["rescaling"]["option"],
                "bimodal": bimodality["bimodal_with_clean_gap"],
                "mean_separation_sd": bimodality["mean_separation_sd"],
                "bic_delta": bimodality["bic_delta_two_minus_one"],
                "tau": parameters.get("tau"),
                "mean_ddg": result["ddg"]["mean"],
                "fraction_mpack_saturated":
                    result["discrimination"]["fraction_mpack_saturated"],
            })
    return rows


def _stability_report(rows: Sequence[Dict[str, Any]]) -> str:
    lines = ["=== Stability of the D1 branch across seeds and the D2 window sweep ===",
             f"  {'N':>4} {'seed':>9} {'opt':>4} {'bimodal':>8} {'sep_SD':>7} "
             f"{'dBIC':>8} {'tau':>7} {'mean_ddG':>9} {'sat':>6}"]
    for row in rows:
        tau = f"{row['tau']:.3f}" if row["tau"] is not None else "n/a"
        lines.append(f"  {row['window_nt']:>4} {row['seed']:>9} {row['option']:>4} "
                     f"{str(row['bimodal']):>8} {row['mean_separation_sd']:>7.2f} "
                     f"{row['bic_delta']:>+8.1f} {tau:>7} {row['mean_ddg']:>9.2f} "
                     f"{row['fraction_mpack_saturated']:>5.0%}")
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--root", type=Path, default=Path("all_hanta"),
                        help="dataset root passed to load_reference_data")
    parser.add_argument("--out", type=Path, default=DEFAULT_ARTIFACT_PATH)
    parser.add_argument("--window-nt", type=int, default=rd.RPI_PANHANDLE_WINDOW_NT,
                        help="panhandle window; D2 sweeps 25/50/75")
    parser.add_argument("--seed", type=int, default=RPI_CALIBRATION_SEED)
    parser.add_argument("--max-per-species", type=int,
                        default=CALIBRATION_MAX_ISOLATES_PER_SPECIES)
    parser.add_argument("--stability", action="store_true",
                        help="re-run the branch across seeds and windows; writes nothing")
    args = parser.parse_args(argv)

    reference = rd.load_reference_data(args.root)
    if args.stability:
        print(_stability_report(stability_sweep(reference)))
        return 0
    result = calibrate(reference, window_nt=args.window_nt, seed=args.seed,
                       max_per_species=args.max_per_species)
    print(_report(result))
    path = write_artifact(result, args.out)
    print(f"\nwrote {path}")
    return 0


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TerminalCompletenessTest(unittest.TestCase):
    def test_accepts_a_sequence_carrying_both_conserved_tags(self):
        self.assertTrue(has_terminal_panhandle(
            HANTAVIRUS_TERMINAL_5P_TAG + "ACGT" * 10 + HANTAVIRUS_TERMINAL_3P_TAG))

    def test_rejects_a_sequence_truncated_at_either_end(self):
        body = "ACGT" * 10
        self.assertFalse(has_terminal_panhandle(body + HANTAVIRUS_TERMINAL_3P_TAG))
        self.assertFalse(has_terminal_panhandle(HANTAVIRUS_TERMINAL_5P_TAG + body))
        self.assertFalse(has_terminal_panhandle(body))

    def test_the_two_tags_are_reverse_complements(self):
        """The filter is only meaningful if the tags describe the two ends of ONE
        panhandle. If these ever drift apart the filter would silently select for
        something that is not a panhandle at all."""
        self.assertEqual(rd._reverse_complement(HANTAVIRUS_TERMINAL_5P_TAG),
                         HANTAVIRUS_TERMINAL_3P_TAG)


class CofoldOrientationTest(unittest.TestCase):
    """The plan's Task 8 Step 1 specified the 3' window reverse-complemented. It is not,
    and these are the measurements that settle it. If someone reintroduces the
    reverse-complement, these fail."""

    def test_cofold_pairs_the_two_strands_antiparallel(self):
        a = "AUGGCAUCGAUCGAUUAGCA"
        antiparallel = rd.cofold_dg(a, rd._reverse_complement(a))
        parallel = rd.cofold_dg(a, rd._reverse_complement(a)[::-1])
        self.assertLess(antiparallel, parallel)

    def test_a_synthetic_panhandle_folds_better_as_is_than_reverse_complemented(self):
        """A genomic 5' terminus and its own 3' terminus, both written 5'->3' as a FASTA
        stores them, are already reverse complements of each other -- so they cofold
        as-is. Reverse-complementing the 3' window destroys the duplex."""
        five = "TAGTAGTAGACTCCCTAAAGAGCTACTAGAACAACGATGGCAACTATGGA"
        three = rd._reverse_complement(five)
        self.assertLess(rd.cofold_dg(five, three),
                        rd.cofold_dg(five, rd._reverse_complement(three)))


class BimodalityDecisionRuleTest(unittest.TestCase):
    """The D1 rule against distributions whose answer is known by construction. This is
    what makes the rule pre-registered rather than fitted: it is pinned to synthetic data
    here, independent of whatever the real panel turns out to look like."""

    def test_two_well_separated_clusters_are_called_bimodal(self):
        rng = np.random.default_rng(7)
        values = np.concatenate([rng.normal(0.0, 1.0, 400), rng.normal(10.0, 1.0, 400)])
        assessment = assess_bimodality(values)
        self.assertTrue(assessment["bimodal_with_clean_gap"], assessment["checks"])
        self.assertGreaterEqual(assessment["density_local_maxima"], 2)

    def test_a_single_gaussian_is_not_called_bimodal(self):
        rng = np.random.default_rng(8)
        assessment = assess_bimodality(rng.normal(0.0, 1.0, 800))
        self.assertFalse(assessment["bimodal_with_clean_gap"])

    def test_a_skewed_unimodal_distribution_is_not_called_bimodal(self):
        """The case Sarle's coefficient gets wrong. A long right tail pulls BC above 5/9
        while the density still has one peak, so the rule must not rest on BC -- and
        this test is the guard that it does not."""
        rng = np.random.default_rng(9)
        values = rng.exponential(1.0, 1200)
        assessment = assess_bimodality(values)
        self.assertGreater(assessment["bimodality_coefficient"],
                           assessment["bimodality_coefficient_threshold"])
        self.assertFalse(assessment["bimodal_with_clean_gap"])

    def test_two_clusters_closer_than_the_threshold_are_not_called_bimodal(self):
        """Guard-the-guard on the separation clause: genuinely two components, but 1.2 SD
        apart, so the density has no valley. Without this the first test would pass for
        any two-component data whatsoever."""
        rng = np.random.default_rng(10)
        values = np.concatenate([rng.normal(0.0, 1.0, 400), rng.normal(1.2, 1.0, 400)])
        assessment = assess_bimodality(values)
        self.assertTrue(assessment["checks"]["bic_improves"]
                        or assessment["checks"]["both_components_substantial"])
        self.assertFalse(assessment["bimodal_with_clean_gap"])

    def test_a_negligible_second_cluster_is_not_called_bimodal(self):
        """Guard-the-guard on the weight clause: 2% of points sitting far away are
        outliers, not a mode."""
        rng = np.random.default_rng(11)
        values = np.concatenate([rng.normal(0.0, 1.0, 980), rng.normal(12.0, 0.5, 20)])
        assessment = assess_bimodality(values)
        self.assertFalse(assessment["checks"]["both_components_substantial"])
        self.assertFalse(assessment["bimodal_with_clean_gap"])


class MixtureFitTest(unittest.TestCase):
    def test_recovers_known_component_parameters(self):
        rng = np.random.default_rng(12)
        values = np.concatenate([rng.normal(-3.0, 0.5, 2000), rng.normal(4.0, 1.0, 2000)])
        fit = fit_two_component_gaussian(values)
        self.assertAlmostEqual(fit["mean_low"], -3.0, delta=0.1)
        self.assertAlmostEqual(fit["mean_high"], 4.0, delta=0.1)
        self.assertAlmostEqual(fit["weight_low"], 0.5, delta=0.05)

    def test_components_are_returned_low_mean_first(self):
        rng = np.random.default_rng(13)
        fit = fit_two_component_gaussian(
            np.concatenate([rng.normal(5.0, 1.0, 300), rng.normal(-5.0, 1.0, 300)]))
        self.assertLess(fit["mean_low"], fit["mean_high"])

    def test_does_not_diverge_when_a_component_collapses_onto_duplicates(self):
        """Real ddG panels contain exact duplicate values -- ViennaRNA reports energies
        to 0.01 kcal/mol, so ties are common. Without the variance floor a component
        collapses onto a tie, its sd goes to 0 and the log-likelihood to +inf, which then
        wins on BIC while describing nothing."""
        values = np.array([1.0] * 200 + list(np.linspace(5.0, 9.0, 200)))
        fit = fit_two_component_gaussian(values)
        self.assertTrue(math.isfinite(fit["loglik"]))
        self.assertTrue(math.isfinite(fit["bic"]))
        self.assertGreater(fit["sd_low"], 0.0)


class RescalingBranchTest(unittest.TestCase):
    def test_bimodal_input_selects_the_logistic(self):
        rng = np.random.default_rng(14)
        values = np.concatenate([rng.normal(0.0, 1.0, 400), rng.normal(12.0, 1.0, 400)])
        choice = choose_rescaling(values, assess_bimodality(values))
        self.assertEqual(choice["form"], "logistic")
        self.assertEqual(choice["option"], "C")
        self.assertGreater(choice["parameters"]["s"], 0.0)
        self.assertLess(choice["parameters"]["ddg_50"], 12.0)
        self.assertGreater(choice["parameters"]["ddg_50"], 0.0)

    def test_unimodal_input_selects_the_tau_scaled_boltzmann(self):
        rng = np.random.default_rng(15)
        values = rng.normal(5.0, 2.0, 800)
        choice = choose_rescaling(values, assess_bimodality(values))
        self.assertEqual(choice["form"], "tau_scaled_boltzmann")
        self.assertEqual(choice["option"], "A")
        self.assertGreater(choice["parameters"]["tau"], 0.0)

    def test_tau_maps_the_measured_span_onto_the_target_dynamic_range(self):
        """The defining property of Option A: exp(-span/tau) must land on 1/D. Asserted
        against boltzmann_compatibility itself, not against the arithmetic, so a change to
        that function's clamping cannot silently break the calibration."""
        values = np.linspace(0.0, 20.0, 1001)
        tau = derive_tau(values)
        span = tau["ddg_span_robust_p5_p95"]
        self.assertAlmostEqual(
            rd.boltzmann_compatibility(span, 0.0, rt=tau["tau"]),
            1.0 / RPI_MPACK_TARGET_DYNAMIC_RANGE, places=6)

    def test_tau_reports_both_the_robust_and_the_raw_span(self):
        """derive_tau deviates from D1's literal min-max reading. Both numbers must be in
        the artifact so the deviation is visible rather than buried."""
        values = np.concatenate([np.linspace(0.0, 10.0, 999), [1000.0]])
        tau = derive_tau(values)
        self.assertLess(tau["ddg_span_robust_p5_p95"], tau["ddg_span_raw_min_max"])
        self.assertLess(tau["tau"], tau["tau_raw_span"])


class BimodalityCoefficientTest(unittest.TestCase):
    def test_a_normal_sample_sits_below_the_uniform_threshold(self):
        rng = np.random.default_rng(16)
        self.assertLess(bimodality_coefficient(rng.normal(0, 1, 5000)), 5.0 / 9.0)

    def test_a_two_cluster_sample_sits_above_it(self):
        rng = np.random.default_rng(17)
        values = np.concatenate([rng.normal(0, 1, 2500), rng.normal(8, 1, 2500)])
        self.assertGreater(bimodality_coefficient(values), 5.0 / 9.0)

    def test_zero_variance_raises_rather_than_dividing_by_zero(self):
        with self.assertRaises(ValueError):
            bimodality_coefficient([3.0] * 50)


class ArtifactShapeTest(unittest.TestCase):
    """Asserts on the SHIPPED JSON, not on calibrate()'s return value. Three reviews of
    this codebase have caught the same defect -- an invariant asserted on an intermediate
    while the value that ships violates it -- so the artifact is round-tripped through
    json and checked on the way back."""

    def _artifact(self) -> Dict[str, Any]:
        import tempfile
        rng = np.random.default_rng(18)
        values = rng.normal(6.0, 2.0, 500)
        assessment = assess_bimodality(values)
        result = {
            "schema_version": RPI_CALIBRATION_SCHEMA_VERSION,
            "seed": RPI_CALIBRATION_SEED,
            "window_nt": rd.RPI_PANHANDLE_WINDOW_NT,
            "fold_provenance": rd.viennarna_provenance(),
            "ddg": _summarise(values),
            "bimodality": assessment,
            "rescaling": choose_rescaling(values, assessment),
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = write_artifact(result, Path(tmp) / "sub" / "artifact.json")
            return json.loads(path.read_text(encoding="utf-8"))

    def test_shipped_artifact_records_mechanism_and_version(self):
        """D2b: a bare version would not distinguish these energies from CLI-produced
        ones, and the whole tau derivation depends on which produced them."""
        shipped = self._artifact()
        self.assertEqual(shipped["fold_provenance"]["mechanism"],
                         "ViennaRNA Python bindings")
        self.assertTrue(shipped["fold_provenance"]["version"])

    def test_shipped_artifact_records_the_window(self):
        self.assertEqual(self._artifact()["window_nt"], 50)

    def test_shipped_numbers_are_json_numbers_not_strings(self):
        """N5's failure mode: numpy scalars serialised as JSON strings. Every statistic
        here comes out of numpy, so the round trip is checked rather than assumed."""
        shipped = self._artifact()
        for key, value in shipped["ddg"].items():
            self.assertIsInstance(value, (int, float), f"ddg.{key} is {type(value)}")
        self.assertIsInstance(shipped["bimodality"]["mean_separation_sd"], float)
        for value in shipped["rescaling"]["parameters"].values():
            self.assertIsInstance(value, (int, float, list))

    def test_shipped_artifact_names_a_branch_and_its_parameters(self):
        shipped = self._artifact()
        self.assertIn(shipped["rescaling"]["form"],
                      ("logistic", "tau_scaled_boltzmann"))
        self.assertTrue(shipped["rescaling"]["parameters"])
        self.assertTrue(shipped["rescaling"]["reason"])

    def test_checks_survive_the_round_trip_as_booleans(self):
        """np.bool_ is not a bool and json's default=str would stringify it, turning
        `false` into the truthy string "False" for any downstream reader."""
        for value in self._artifact()["bimodality"]["checks"].values():
            self.assertIsInstance(value, bool)


class CliffsDeltaTest(unittest.TestCase):
    def test_is_plus_one_when_every_a_exceeds_every_b(self):
        self.assertAlmostEqual(cliffs_delta([10.0, 11.0, 12.0], [1.0, 2.0]), 1.0)

    def test_is_minus_one_when_every_a_is_below_every_b(self):
        self.assertAlmostEqual(cliffs_delta([1.0, 2.0], [10.0, 11.0, 12.0]), -1.0)

    def test_is_zero_for_identical_samples(self):
        """Ties must count as neither greater nor less. A searchsorted implementation that
        used one side for both counts would report a spurious +/-1 here."""
        values = [1.0, 2.0, 3.0, 4.0]
        self.assertAlmostEqual(cliffs_delta(values, list(values)), 0.0)

    def test_matches_the_naive_double_loop(self):
        """Guard on the sorted implementation: it replaces an O(n*m) loop and must agree
        with it exactly, ties included."""
        rng = np.random.default_rng(21)
        a = np.round(rng.normal(0, 2, 60), 1)
        b = np.round(rng.normal(1, 2, 40), 1)
        naive = sum((1 if x > y else -1 if x < y else 0) for x in a for y in b)
        self.assertAlmostEqual(cliffs_delta(a, b), naive / (a.size * b.size))

    def test_empty_sample_raises(self):
        with self.assertRaises(ValueError):
            cliffs_delta([], [1.0])


class DiscriminationTest(unittest.TestCase):
    @staticmethod
    def _observations(ddgs: Sequence[float]) -> List[Dict[str, Any]]:
        return [{"focal_isolate_key": f"focal{i}", "five_prime_segment": "S",
                 "three_prime_segment": "M", "dg_cognate": -30.0,
                 "dg_non_cognate": -30.0 + d, "ddg": d} for i, d in enumerate(ddgs)]

    def test_saturated_fraction_counts_the_non_positive_ddgs(self):
        observations = self._observations([-2.0, -1.0, 0.0, 1.0, 2.0])
        rescaling = {"form": "tau_scaled_boltzmann", "parameters": {"tau": 4.0}}
        result = assess_discrimination(observations, rescaling)
        self.assertAlmostEqual(result["fraction_mpack_saturated"], 0.6)
        self.assertAlmostEqual(result["fraction_ddg_positive"], 0.4)

    def test_the_two_fractions_partition_the_panel(self):
        """Guard-the-guard: ddG == 0.0 must land in exactly one of the two buckets. If
        both used a strict comparison the fractions would silently fail to sum to 1."""
        rng = np.random.default_rng(22)
        ddgs = list(np.round(rng.normal(0.0, 3.0, 400), 1))
        result = assess_discrimination(
            self._observations(ddgs),
            {"form": "tau_scaled_boltzmann", "parameters": {"tau": 4.0}})
        self.assertAlmostEqual(result["fraction_ddg_positive"]
                               + result["fraction_mpack_saturated"], 1.0)

    def test_m_pack_is_bounded_and_never_exactly_zero(self):
        """The property WS3.2 must preserve, asserted here on the shipped summary rather
        than on boltzmann_compatibility in isolation."""
        observations = self._observations([-50.0, -5.0, 0.0, 5.0, 50.0, 500.0])
        result = assess_discrimination(
            observations, {"form": "tau_scaled_boltzmann", "parameters": {"tau": 4.0}})
        summary = result["m_pack_under_chosen_form"]
        self.assertLessEqual(summary["max"], 1.0)
        self.assertGreater(summary["min"], 0.0)

    def test_m_pack_is_monotonic_decreasing_in_ddg(self):
        for form, parameters in (("tau_scaled_boltzmann", {"tau": 4.0}),
                                 ("logistic", {"ddg_50": 5.0, "s": 1.5})):
            result = assess_discrimination(
                self._observations([1.0, 2.0, 3.0]),
                {"form": form, "parameters": parameters})
            self.assertLessEqual(result["m_pack_under_chosen_form"]["max"], 1.0)
            ascending = assess_discrimination(
                self._observations([1.0]), {"form": form, "parameters": parameters})
            descending = assess_discrimination(
                self._observations([9.0]), {"form": form, "parameters": parameters})
            self.assertGreater(ascending["m_pack_under_chosen_form"]["mean"],
                               descending["m_pack_under_chosen_form"]["mean"], form)


class DistinctCognateDgsTest(unittest.TestCase):
    def test_deduplicates_the_cognate_fold_reused_across_partners(self):
        """Each cognate fold appears once per foreign partner in the observation list.
        Summarising the raw column would report an n inflated by that factor -- so the
        shipped dg_cognate.n must be the number of FOLDS, not of observations."""
        observations = [
            {"focal_isolate_key": "a", "five_prime_segment": "S",
             "three_prime_segment": "M", "dg_cognate": -30.0},
            {"focal_isolate_key": "a", "five_prime_segment": "S",
             "three_prime_segment": "M", "dg_cognate": -30.0},
            {"focal_isolate_key": "a", "five_prime_segment": "S",
             "three_prime_segment": "L", "dg_cognate": -25.0},
            {"focal_isolate_key": "b", "five_prime_segment": "S",
             "three_prime_segment": "M", "dg_cognate": -30.0},
        ]
        # Ordered by fold identity (isolate, 5' segment, 3' segment), so a/S/L precedes
        # a/S/M. Asserting the exact order, not just the multiset: the artifact's
        # dg_cognate summary must be reproducible run to run.
        self.assertEqual(distinct_cognate_dgs(observations), [-25.0, -30.0, -30.0])

    def test_keeps_equal_dgs_from_different_folds(self):
        """Guard-the-guard: deduplication is by FOLD IDENTITY, not by value. Two
        different segment pairs that happen to fold to the same energy are two
        measurements and must both count."""
        observations = [
            {"focal_isolate_key": "a", "five_prime_segment": "S",
             "three_prime_segment": "M", "dg_cognate": -30.0},
            {"focal_isolate_key": "a", "five_prime_segment": "M",
             "three_prime_segment": "L", "dg_cognate": -30.0},
        ]
        self.assertEqual(len(distinct_cognate_dgs(observations)), 2)


class SegmentDuplexTest(unittest.TestCase):
    def test_there_are_six_ordered_inter_segment_duplexes(self):
        self.assertEqual(len(SEGMENT_DUPLEXES), 6)
        self.assertEqual(len(set(SEGMENT_DUPLEXES)), 6)

    def test_no_segment_is_paired_with_itself(self):
        """Intra-segment duplexes are the classical intramolecular panhandle, not the
        inter-segment quantity M_pack scores. A reassortant never swaps a segment against
        itself, so including them would dilute the contrast with a term that cannot
        differ between cognate and non-cognate."""
        for five, three in SEGMENT_DUPLEXES:
            self.assertNotEqual(five, three)

    def test_both_orientations_of_each_pair_are_present(self):
        for pair in (("S", "M"), ("S", "L"), ("M", "L")):
            self.assertIn(pair, SEGMENT_DUPLEXES)
            self.assertIn(pair[::-1], SEGMENT_DUPLEXES)


if __name__ == "__main__":
    sys.exit(main())
