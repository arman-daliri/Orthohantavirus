"""pipeline.py -- thin shim.

Re-exports this module's original public surface from rpi_core.py, the
single consolidated module the RPI toolchain now lives in (branch
05-module-consolidation). See notes_consolidation.md for exactly what
moved and why; see rpi_core.py's own banner comments for this module's
section within it.
"""
import sys

from rpi_core import (
    COMMANDS, PipelineStageError, _run_stage, _stage_calibrate,
    _stage_discriminate, _stage_stats, _default_output_path, _stage_report,
    _format_summary, run, _parse_run_args, StageCalibrateTest,
    RunHaltsOnStageFailureTest, StageReportTest, MainDispatchTest
)
from rpi_core import _pipeline_main as main


if __name__ == "__main__":
    sys.exit(main())
