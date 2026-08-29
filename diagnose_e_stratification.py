"""Standalone diagnostic (WS4 Task 13/14, PROJECT-STATE.md section 7c).

Loads an already-generated rpi_discrimination.py panel JSON (does NOT regenerate it),
stratifies groups.c_real_positive by rpi_e_total, and re-runs
rpi_discrimination_stats.compare_groups for reals(E~1) vs b1 only -- to check whether
the reals-vs-b1 tie reported by evaluate_pass_criterion is driven by the published
inter-species lab crosses (which carry a near-zero Factor E ecological prior) diluting
an otherwise-intra-species real-positive panel.

Pure read-and-recompute: does not import or modify rpi_discrimination.py or
rpi_discrimination_stats.py beyond calling their public compare_groups().

Usage:
    python diagnose_e_stratification.py output/discrimination_scratch/discrimination_panel.json
"""
import json
import sys
from collections import Counter

from rpi_discrimination_stats import compare_groups

E_THRESHOLD = 0.5  # chosen after inspecting the distribution below -- see printed histogram


def main() -> None:
    path = sys.argv[1] if len(sys.argv) > 1 else "output/discrimination_scratch/discrimination_panel.json"
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    groups = data["groups"]
    c_real = groups["c_real_positive"]
    b1 = groups["b1_within_species_shuffle"]

    scored_real = [row for row in c_real if row.get("rpi_generation") is not None]
    unscored_real = len(c_real) - len(scored_real)

    e_counts = Counter(row["rpi_e_total"] for row in scored_real)
    print(f"Loaded: {path}")
    print(f"c_real_positive: n={len(c_real)}, scored={len(scored_real)}, "
          f"unscored/junk_sequence={unscored_real}")
    print()
    print("rpi_e_total distribution across scored c_real_positive rows:")
    for value, count in sorted(e_counts.items(), key=lambda kv: kv[0]):
        print(f"  rpi_e_total = {value!r}: {count} case(s)")
    print()
    print(f"Chosen threshold: rpi_e_total >= {E_THRESHOLD} "
          f"(distribution is bimodal -- 1.0 vs 0.01 -- this threshold falls cleanly "
          f"in the gap and does not need tuning)")
    print()

    reals_e1 = [row for row in scored_real if row["rpi_e_total"] >= E_THRESHOLD]
    reals_e_low = [row for row in scored_real if row["rpi_e_total"] < E_THRESHOLD]
    print(f"reals(E~1) stratum: n={len(reals_e1)} -- case_ids: "
          f"{[row['case_id'] for row in reals_e1]}")
    print(f"reals(E-low) stratum (excluded from this comparison): n={len(reals_e_low)} -- "
          f"case_ids: {[row['case_id'] for row in reals_e_low]}")
    print()

    b1_values = [row["rpi_generation"] for row in b1 if row.get("rpi_generation") is not None]
    reals_e1_values = [row["rpi_generation"] for row in reals_e1]

    result = compare_groups(
        reals_e1_values, b1_values,
        "c_real_positive(E>=0.5)", "b1_within_species_shuffle",
        expect="reals(E~1) ~ b1: same D7 tie-expectation as the unstratified "
               "reals_vs_b1 comparison, restricted to the intra-species subset of "
               "reals to check whether Factor E composition (not M_pack/M_pol) drove "
               "the original result.")
    result["plausible_tie_not_proven_equivalent"] = (
        not result["significant_at_alpha_0_05"]
        or result["cliffs_delta_magnitude"] == "negligible")

    print("compare_groups(reals(E>=0.5), b1) result:")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
