"""rpi_calibration.py -- thin shim.

Re-exports this module's original public surface from rpi_core.py, the
single consolidated module the RPI toolchain now lives in (branch
05-module-consolidation). See notes_consolidation.md for exactly what
moved and why; see rpi_core.py's own banner comments for this module's
section within it.
"""
import sys

from rpi_core import (
    RPI_CALIBRATION_SEED, DEFAULT_ARTIFACT_PATH, HANTAVIRUS_TERMINAL_5P_TAG,
    HANTAVIRUS_TERMINAL_3P_TAG, CALIBRATION_MAX_ISOLATES_PER_SPECIES,
    CALIBRATION_FOREIGN_PARTNERS_PER_FOCAL, UNCLASSIFIED_SPECIES_LABEL,
    SEGMENTS, SEGMENT_DUPLEXES, BIMODALITY_REQUIRES_BIC_IMPROVEMENT,
    BIMODALITY_MIN_COMPONENT_WEIGHT, BIMODALITY_MIN_MEAN_SEPARATION_SD,
    BIMODALITY_DENSITY_GRID_POINTS, RPI_MPACK_TARGET_DYNAMIC_RANGE,
    TAU_SPAN_LOW_PERCENTILE, TAU_SPAN_HIGH_PERCENTILE, segment_accession,
    genomic_sequence, has_terminal_panhandle, _is_wild, _segment_is_full_length,
    terminally_complete_isolates, ws4_reserved_accessions, _isolate_touches,
    select_calibration_isolates, duplex_dg, _summarise, bimodality_coefficient,
    _normal_logpdf, fit_two_component_gaussian, fit_one_component_gaussian,
    _mixture_density_modes, assess_bimodality, derive_tau, cliffs_delta,
    distinct_cognate_dgs, assess_discrimination, derive_logistic,
    choose_rescaling, RPI_CALIBRATION_SCHEMA_VERSION, calibrate, write_artifact,
    _report, STABILITY_SEEDS, STABILITY_WINDOWS, stability_sweep,
    _stability_report, TerminalCompletenessTest, CofoldOrientationTest,
    BimodalityDecisionRuleTest, MixtureFitTest, RescalingBranchTest,
    BimodalityCoefficientTest, ArtifactShapeTest, CliffsDeltaTest,
    DiscriminationTest, DistinctCognateDgsTest, SegmentDuplexTest
)
from rpi_core import _calibration_build_panel as build_panel
from rpi_core import _calibration_main as main


if __name__ == "__main__":
    sys.exit(main())
