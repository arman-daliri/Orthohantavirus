# rpi_discrimination_stats.py
"""WS4 Task 14 -- the pass criterion over rpi_discrimination.py's scored panel.

Deliberately a SEPARATE module from rpi_discrimination.py (see that module's own
docstring): panel construction is subject to real-data availability (how many sympatric
pairs exist, which isolates qualify), while this module is subject to methodology (which
test, which effect size, what "equivalence" is allowed to mean). Keeping them apart lets
either change without the other's tests needing to move.

THE CRITERION (spec/plan, D7): `reals ~ b1 > b2 > random`, reported as three pairwise
comparisons on `rpi_generation`:
  1. b2 vs random   -- expect b2 SIGNIFICANTLY HIGHER. The weak floor: even the hardest
                        cross-species group should beat pure noise.
  2. b1 vs b2        -- expect b1 SIGNIFICANTLY HIGHER. THE REAL TEST: b2 cannot be
                        rejected by Factor E (it is ecologically sympatric by
                        construction), so any gap here is M_pack + M_pol earning their
                        keep.
  3. reals vs b1     -- expect NO significant difference -- see EQUIVALENCE below.

NO SCIPY IN THIS ENVIRONMENT (requirements-pip.txt / environment.yml reconciliation,
PROJECT-STATE.md) -- the Mann-Whitney U p-value is computed by hand via the normal
approximation with a tie correction, the same "hand-roll rather than add a dependency"
convention rpi_calibration.py already uses for its own EM mixture fit and Cliff's delta.
Adequate here: every group in this panel has n >= 11 (see rpi_discrimination.py's
sample-size comments), well inside the range the normal approximation is conventionally
trusted for (n1, n2 >= 8 with a tie correction).

EFFECT SIZE is Cliff's delta, REUSED from rpi_calibration.py (not reimplemented) --
P(x > y) - P(x < y), already the effect size WS3.1's own M_pack finding is reported in
(delta = +0.222), so a reader comparing this report to that one is comparing the same
statistic. Magnitude bands follow Romano, Kromrey, Coraggio & Skowronek (2006): negligible
< 0.147, small < 0.33, medium < 0.474, large >= 0.474 -- the same convention already
implicit in WS3.1's own "small" label for 0.222.

EQUIVALENCE, HONESTLY (D7 / Task 14's own binding text): "reals ~ b1 asserted as an
equivalence interval, never as a separation -- a non-significant U-test is absence of
evidence, not evidence of the tie being claimed." This module does NOT run a formal TOST
equivalence test (that needs a pre-registered margin in RPI's own units, which no prior
decision in this project has set, and inventing one now would be exactly the kind of
after-the-fact threshold-picking D1's bimodality rule was written to prevent). Instead,
`compare_groups` reports BOTH: the U-test's significance (or lack of it) AND Cliff's
delta's magnitude, and `evaluate_pass_criterion` labels a comparison
"plausible_tie_not_proven_equivalent" only when EITHER of the two would prompt someone to
call it a tie -- p > 0.05 for the location-shift test, OR delta magnitude "negligible" by
the band above -- WITH THAT EXACT CAVEAT STRING attached, never a bare "equivalent".

Run:
    python -m unittest rpi_discrimination_stats    # unit tests
"""

import math
import unittest
from typing import Any, Dict, List, Optional, Sequence, Tuple

import rpi_calibration as cal

ALPHA = 0.05

# Romano, Kromrey, Coraggio & Skowronek (2006) magnitude bands for Cliff's delta.
CLIFFS_DELTA_BANDS: Tuple[Tuple[float, str], ...] = (
    (0.147, "negligible"), (0.33, "small"), (0.474, "medium"), (float("inf"), "large"),
)


def cliffs_delta_magnitude(delta: float) -> str:
    magnitude = abs(delta)
    for threshold, label in CLIFFS_DELTA_BANDS:
        if magnitude < threshold:
            return label
    return "large"  # unreachable given the inf sentinel above; kept for clarity


def _normal_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def mann_whitney_u(a: Sequence[float], b: Sequence[float]) -> Dict[str, float]:
    """Two-sided Mann-Whitney U test via the normal approximation, WITH a tie
    correction to the variance (Mann & Whitney 1947 / Lehmann 1975's standard formula).

    Returns U (for group `a`), z, and the two-sided p-value. Raises ValueError on an
    empty sample -- an empty group is a panel-construction problem to surface loudly,
    not a comparison to silently skip.
    """
    if not a or not b:
        raise ValueError("mann_whitney_u needs two non-empty samples")
    n1, n2 = len(a), len(b)
    combined = sorted((value, 0) for value in a) + sorted((value, 1) for value in b)
    combined.sort(key=lambda pair: pair[0])

    # Average ranks for ties (1-indexed), then sum the ranks belonging to group `a`.
    ranks = [0.0] * len(combined)
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        average_rank = (i + 1 + j) / 2.0  # ranks i+1..j, 1-indexed, averaged
        for k in range(i, j):
            ranks[k] = average_rank
        i = j

    rank_sum_a = sum(rank for rank, (_, group) in zip(ranks, combined) if group == 0)
    u_a = rank_sum_a - n1 * (n1 + 1) / 2.0

    # Tie correction: sum over tied groups of (t^3 - t), t = group size.
    tie_term = 0.0
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        t = j - i
        if t > 1:
            tie_term += t ** 3 - t
        i = j

    n = n1 + n2
    mean_u = n1 * n2 / 2.0
    variance_u = (n1 * n2 / 12.0) * ((n + 1) - tie_term / (n * (n - 1))) if n > 1 else 0.0
    if variance_u <= 0.0:
        # Every value tied across both groups -- no information to test.
        return {"u": u_a, "z": 0.0, "p_two_sided": 1.0, "n1": n1, "n2": n2}
    sd_u = math.sqrt(variance_u)
    # Continuity correction: shrink |u_a - mean_u| by 0.5 toward zero.
    diff = u_a - mean_u
    diff = diff - 0.5 if diff > 0 else (diff + 0.5 if diff < 0 else 0.0)
    z = diff / sd_u
    p_two_sided = 2.0 * (1.0 - _normal_cdf(abs(z)))
    return {"u": u_a, "z": z, "p_two_sided": min(1.0, p_two_sided), "n1": n1, "n2": n2}


def compare_groups(
        a_values: Sequence[float], b_values: Sequence[float],
        a_label: str, b_label: str, expect: str,
) -> Dict[str, Any]:
    """One pairwise comparison: Mann-Whitney U (a vs b) + Cliff's delta (a vs b).

    `expect` is a FREE-TEXT description of what this comparison is FOR (e.g. "a > b:
    the real test of M_pack + M_pol"), recorded in the output so a reader does not have
    to cross-reference this module's docstring to know why the comparison exists. It is
    never used to alter the computed statistics -- see `_directional_result` for the
    only place the expectation and the measurement are compared against each other.
    """
    u_result = mann_whitney_u(a_values, b_values)
    delta = cal.cliffs_delta(a_values, b_values)
    return {
        "a_label": a_label, "b_label": b_label, "expect": expect,
        "n_a": len(a_values), "n_b": len(b_values),
        "a_median": _median(a_values), "b_median": _median(b_values),
        "u_statistic": u_result["u"], "z": u_result["z"],
        "p_two_sided": u_result["p_two_sided"],
        "significant_at_alpha_0_05": u_result["p_two_sided"] < ALPHA,
        "cliffs_delta": delta,
        "cliffs_delta_magnitude": cliffs_delta_magnitude(delta),
    }


def _median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2 == 1:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def _values(scored_group: Sequence[Dict[str, Any]], key: str = "rpi_generation") -> List[float]:
    return [row[key] for row in scored_group if row.get(key) is not None]


def evaluate_pass_criterion(
        scored_groups: Dict[str, Sequence[Dict[str, Any]]],
) -> Dict[str, Any]:
    """The three pairwise comparisons the D7 ranking `reals ~ b1 > b2 > random` reduces
    to (see module docstring). Requires the four rpi_discrimination.py group keys.

    HONEST-FAILURE CLAUSE (Task 14): this function reports whatever the numbers say. A
    comparison failing its `expect` direction is recorded as `"met_expectation": False`
    and is NOT hidden, retried with different sampling, or used to suggest retuning tau
    -- see reassortment_detector.py's RPI_WEIGHTS comment on `m_pack_tau` and this
    project's repeated 2026-08-29 decisions not to retune against the validation set.
    """
    random_values = _values(scored_groups["a_random_fakes"])
    b1_values = _values(scored_groups["b1_within_species_shuffle"])
    b2_values = _values(scored_groups["b2_cross_species_sympatric"])
    real_values = _values(scored_groups["c_real_positive"])

    b2_vs_random = compare_groups(
        b2_values, random_values, "b2_cross_species_sympatric", "a_random_fakes",
        expect="b2 > random: the weak floor -- even ecologically sympatric cross-species "
               "fakes should beat isolates assembled with no ecological or molecular "
               "plausibility at all.")
    b2_vs_random["met_expectation"] = (
        b2_vs_random["significant_at_alpha_0_05"] and b2_vs_random["cliffs_delta"] > 0)

    b1_vs_b2 = compare_groups(
        b1_values, b2_values, "b1_within_species_shuffle", "b2_cross_species_sympatric",
        expect="b1 > b2: THE REAL TEST. Factor E cannot separate these (b2 is "
               "ecologically sympatric by construction), so a gap here is M_pack + "
               "M_pol's own discrimination, not ecology's.")
    b1_vs_b2["met_expectation"] = (
        b1_vs_b2["significant_at_alpha_0_05"] and b1_vs_b2["cliffs_delta"] > 0)

    reals_vs_b1 = compare_groups(
        real_values, b1_values, "c_real_positive", "b1_within_species_shuffle",
        expect="reals ~ b1: NOT a claim of separation. RPI measures plausibility, not "
               "occurrence (decision D7) -- real reassortants and plausible "
               "within-species shuffles are expected to be indistinguishable to this "
               "index; separating them is gate 4's job.")
    reals_vs_b1["plausible_tie_not_proven_equivalent"] = (
        not reals_vs_b1["significant_at_alpha_0_05"]
        or reals_vs_b1["cliffs_delta_magnitude"] == "negligible")
    reals_vs_b1["equivalence_caveat"] = (
        "A non-significant U-test is ABSENCE OF EVIDENCE of a difference, not PROOF of "
        "equivalence -- no formal equivalence margin is defined for RPI's own units. "
        "'plausible_tie_not_proven_equivalent' names exactly what is and is not being "
        "claimed here.")

    return {
        "b2_vs_random": b2_vs_random,
        "b1_vs_b2": b1_vs_b2,
        "reals_vs_b1": reals_vs_b1,
        "overall_ranking_supported": (
            b2_vs_random["met_expectation"] and b1_vs_b2["met_expectation"]),
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class MannWhitneyUTest(unittest.TestCase):
    def test_identical_distributions_are_not_significant(self):
        a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        b = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]
        result = mann_whitney_u(a, b)
        self.assertGreater(result["p_two_sided"], 0.05)

    def test_completely_separated_groups_are_significant(self):
        a = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
        b = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        result = mann_whitney_u(a, b)
        self.assertLess(result["p_two_sided"], 0.01)
        self.assertGreater(result["u"], result["n1"] * result["n2"] / 2.0)

    def test_symmetric_in_which_group_is_a_up_to_sign(self):
        a = [10.0, 11.0, 12.0, 13.0, 14.0]
        b = [1.0, 2.0, 3.0, 4.0, 5.0]
        forward = mann_whitney_u(a, b)
        backward = mann_whitney_u(b, a)
        self.assertAlmostEqual(forward["p_two_sided"], backward["p_two_sided"], places=9)
        self.assertAlmostEqual(forward["u"] + backward["u"],
                               forward["n1"] * forward["n2"])

    def test_handles_ties_without_a_zero_or_negative_variance(self):
        a = [5.0] * 10
        b = [5.0] * 10
        result = mann_whitney_u(a, b)
        self.assertEqual(result["p_two_sided"], 1.0)

    def test_raises_on_an_empty_group(self):
        with self.assertRaises(ValueError):
            mann_whitney_u([], [1.0, 2.0])


class CliffsDeltaMagnitudeTest(unittest.TestCase):
    def test_bands_match_romano_et_al(self):
        self.assertEqual(cliffs_delta_magnitude(0.05), "negligible")
        self.assertEqual(cliffs_delta_magnitude(0.20), "small")
        self.assertEqual(cliffs_delta_magnitude(0.40), "medium")
        self.assertEqual(cliffs_delta_magnitude(0.60), "large")

    def test_sign_does_not_affect_magnitude(self):
        self.assertEqual(cliffs_delta_magnitude(-0.20), cliffs_delta_magnitude(0.20))

    def test_ws3_1_m_pack_finding_bands_as_small(self):
        """Regression guard: WS3.1 reported Cliff's delta = +0.222 for M_pack and
        called it "small" in every doc this project shipped. If the bands ever drift,
        this is the test that would silently relabel a number that has already been
        quoted publicly."""
        self.assertEqual(cliffs_delta_magnitude(0.222), "small")


class CompareGroupsTest(unittest.TestCase):
    def test_reports_direction_consistent_with_cliffs_delta_sign(self):
        higher = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
        lower = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        result = compare_groups(higher, lower, "higher", "lower", expect="higher > lower")
        self.assertGreater(result["cliffs_delta"], 0)
        self.assertTrue(result["significant_at_alpha_0_05"])
        self.assertEqual(result["expect"], "higher > lower")


class EvaluatePassCriterionTest(unittest.TestCase):
    def _rows(self, values):
        return [{"rpi_generation": v} for v in values]

    def test_a_clean_pass_reports_met_expectation_true_throughout(self):
        groups = {
            "a_random_fakes": self._rows([0.001, 0.002, 0.001, 0.003, 0.002, 0.001,
                                          0.002, 0.001, 0.002, 0.001, 0.002]),
            "b2_cross_species_sympatric": self._rows([0.05, 0.06, 0.05, 0.07, 0.06, 0.05,
                                                      0.06, 0.05, 0.06, 0.05, 0.06]),
            "b1_within_species_shuffle": self._rows([0.8, 0.85, 0.9, 0.82, 0.88, 0.79,
                                                      0.91, 0.83, 0.86, 0.84, 0.87]),
            "c_real_positive": self._rows([0.79, 0.86, 0.91, 0.81, 0.89, 0.78, 0.92,
                                           0.82, 0.87, 0.85, 0.88]),
        }
        result = evaluate_pass_criterion(groups)
        self.assertTrue(result["b2_vs_random"]["met_expectation"])
        self.assertTrue(result["b1_vs_b2"]["met_expectation"])
        self.assertTrue(result["reals_vs_b1"]["plausible_tie_not_proven_equivalent"])
        self.assertTrue(result["overall_ranking_supported"])

    def test_reals_vs_b1_never_claims_bare_equivalence(self):
        """The exact wording this test guards: 'plausible_tie_not_proven_equivalent',
        never a bare 'equivalent' key, and the caveat string is always present
        regardless of which way the comparison comes out."""
        groups = {
            "a_random_fakes": self._rows([0.001] * 11),
            "b2_cross_species_sympatric": self._rows([0.05] * 11),
            "b1_within_species_shuffle": self._rows([0.8] * 11),
            "c_real_positive": self._rows([0.8] * 11),
        }
        result = evaluate_pass_criterion(groups)
        self.assertIn("plausible_tie_not_proven_equivalent", result["reals_vs_b1"])
        self.assertNotIn("equivalent", result["reals_vs_b1"])
        self.assertIn("ABSENCE OF EVIDENCE", result["reals_vs_b1"]["equivalence_caveat"])

    def test_a_failed_b2_vs_random_separation_is_reported_not_hidden(self):
        """Honest-failure clause: if b2 does not separate from random, this function
        must say so plainly (met_expectation False), not omit the comparison or
        silently pass it."""
        groups = {
            "a_random_fakes": self._rows([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                          0.5, 0.5]),
            "b2_cross_species_sympatric": self._rows([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                                      0.5, 0.5, 0.5, 0.5]),
            "b1_within_species_shuffle": self._rows([0.8] * 11),
            "c_real_positive": self._rows([0.8] * 11),
        }
        result = evaluate_pass_criterion(groups)
        self.assertFalse(result["b2_vs_random"]["met_expectation"])
        self.assertFalse(result["overall_ranking_supported"])


if __name__ == "__main__":
    unittest.main()
