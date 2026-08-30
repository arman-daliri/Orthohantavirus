"""rpi_discrimination_stats.py -- thin shim.

Re-exports this module's original public surface from rpi_core.py, the
single consolidated module the RPI toolchain now lives in (branch
05-module-consolidation). See notes_consolidation.md for exactly what
moved and why; see rpi_core.py's own banner comments for this module's
section within it.
"""
import unittest

from rpi_core import (
    ALPHA, CLIFFS_DELTA_BANDS, cliffs_delta_magnitude, _normal_cdf,
    mann_whitney_u, compare_groups, _median, _values, evaluate_pass_criterion,
    MannWhitneyUTest, CliffsDeltaMagnitudeTest, CompareGroupsTest,
    EvaluatePassCriterionTest
)


if __name__ == "__main__":
    unittest.main()
