"""Single orchestrating entry point over the RPI fitness-function toolchain.

Does NOT merge reassortment_detector.py, rpi_calibration.py, rpi_discrimination.py,
rpi_discrimination_stats.py, validate_reassortment_detector.py,
diagnose_e_stratification.py or diagnose_gate0_junk.py -- those stay separate files for
the reasons documented in their own docstrings (calibration must be a pinned artifact,
not recomputed live; panel construction and pass-criterion statistics change for
different reasons and so live in different modules; the detector-validation harness
checks a different system than RPI; the diagnostics are read-only one-offs). This module
only calls their existing public functions in one place, so "what can I run and in what
order" has one answer.

`run` (the default command) executes the intended WS3/WS4 order end to end:
    calibrate (skipped if data/rpi_calibration.json already exists, unless --recalibrate)
    -> discriminate (score the four-group panel)
    -> stats (the reals ~ b1 > b2 > random pass criterion)
    -> report (one combined JSON + console summary)
A stage that raises halts the run there -- each stage's input is the previous stage's
output, so continuing into a stage whose input doesn't exist would just fail again less
legibly.

`validate`, `diagnose-e` and `diagnose-gate0` are separate subcommands, not part of
`run`: validate checks a different system (the gate 0-4 binary detector, not RPI), and
the diagnostics are one-off investigations, not a step anyone runs every time.

Run:
    python pipeline.py                 # same as `run`
    python pipeline.py run --recalibrate
    python pipeline.py validate --json output/validation/validation_report.json
    python pipeline.py diagnose-e
    python pipeline.py diagnose-gate0
    python -m unittest pipeline
"""
import argparse
import json
import sys
import unittest
import unittest.mock as mock
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import diagnose_e_stratification as diag_e
import diagnose_gate0_junk as diag_gate0
import reassortment_detector as rd
import rpi_calibration as cal
import rpi_discrimination as disc
import rpi_discrimination_stats as stats
import validate_reassortment_detector as vd

COMMANDS = ("run", "validate", "diagnose-e", "diagnose-gate0")


class PipelineStageError(RuntimeError):
    """Raised when a `run` stage fails; carries which stage and the original cause."""

    def __init__(self, stage: str, cause: BaseException):
        super().__init__(f"pipeline halted at stage '{stage}': {cause}")
        self.stage = stage
        self.cause = cause


def _run_stage(name: str, fn):
    try:
        return fn()
    except Exception as exc:  # noqa: BLE001 -- deliberately broad, re-raised typed
        raise PipelineStageError(name, exc) from exc


def _stage_calibrate(
        reference: "rd.ReferenceData", recalibrate: bool,
        path: Path = cal.DEFAULT_ARTIFACT_PATH) -> Dict[str, Any]:
    """Reuses the pinned artifact unless it is missing or --recalibrate was passed --
    calibration is a distinct, auditable step (rpi_calibration.py's own docstring),
    not something a pipeline run should silently redo."""
    if recalibrate or not path.exists():
        result = cal.calibrate(reference)
        cal.write_artifact(result, path)
    return json.loads(path.read_text(encoding="utf-8"))


def _stage_discriminate(
        reference: "rd.ReferenceData", seed: int) -> Dict[str, List[Dict[str, Any]]]:
    precomputed = rd._build_precomputed(reference)
    panel = disc.build_discrimination_panel(reference, seed=seed)
    return {name: disc.score_cases(reference, precomputed, cases)
            for name, cases in panel.items()}


def _stage_stats(scored_groups: Dict[str, List[Dict[str, Any]]]) -> Dict[str, Any]:
    return stats.evaluate_pass_criterion(scored_groups)


def _default_output_path() -> Path:
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return Path("output") / "pipeline" / f"run_{ts}.json"


def _stage_report(
        calibration: Dict[str, Any], scored_groups: Dict[str, List[Dict[str, Any]]],
        stats_result: Dict[str, Any], output_path: Path) -> Path:
    payload = {
        "tau": calibration["rescaling"]["parameters"]["tau"],
        "group_summary": disc.summarise_groups(scored_groups),
        "groups": scored_groups,
        "pass_criterion": stats_result,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(rd._json_safe(payload), indent=2, default=str, allow_nan=False),
        encoding="utf-8")
    return output_path


def _format_summary(
        calibration: Dict[str, Any], scored_groups: Dict[str, List[Dict[str, Any]]],
        stats_result: Dict[str, Any]) -> str:
    lines = [f"tau = {calibration['rescaling']['parameters']['tau']:.3f}"]
    for name, rows in scored_groups.items():
        lines.append(f"  {name}: n={len(rows)}")
    lines.append(f"overall_ranking_supported = {stats_result.get('overall_ranking_supported')}")
    return "\n".join(lines)


def run(args: argparse.Namespace) -> int:
    reference = rd.load_reference_data(Path(args.reference_input) / "all_hanta")
    try:
        calibration = _run_stage(
            "calibrate", lambda: _stage_calibrate(reference, args.recalibrate))
        scored_groups = _run_stage(
            "discriminate", lambda: _stage_discriminate(reference, args.seed))
        stats_result = _run_stage("stats", lambda: _stage_stats(scored_groups))
    except PipelineStageError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    output_path = args.output or _default_output_path()
    report_path = _stage_report(calibration, scored_groups, stats_result, output_path)
    print(_format_summary(calibration, scored_groups, stats_result))
    print(f"\nwrote {report_path}")
    return 0


def _parse_run_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="pipeline.py run")
    parser.add_argument("--recalibrate", action="store_true",
                        help="rebuild data/rpi_calibration.json even if it exists")
    parser.add_argument("--reference-input", default=".",
                        help="package root containing all_hanta/ncbi_dataset/data")
    parser.add_argument("--output", type=Path, default=None,
                        help="combined report path (default output/pipeline/run_<ts>.json)")
    parser.add_argument("--seed", type=int, default=disc.RPI_DISCRIMINATION_SEED)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in COMMANDS:
        command, rest = argv[0], argv[1:]
    else:
        command, rest = "run", argv

    if command == "validate":
        return vd.main(rest)
    if command == "diagnose-e":
        diag_e.main()
        return 0
    if command == "diagnose-gate0":
        diag_gate0.main()
        return 0
    return run(_parse_run_args(rest))


# ---------------------------------------------------------------------------
# Tests
#
# The real stages are slow (rpi_discrimination.py's panel alone takes ~15 min) and
# rpi_calibration.py's ddG folding takes real ViennaRNA time -- these tests stub the
# stage functions rather than calling them live. What is under test is THIS module's
# own orchestration logic (skip-vs-recalibrate, halt-on-failure, command dispatch,
# report shape), not the modules it calls.
# ---------------------------------------------------------------------------


class StageCalibrateTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(self._testMethodName + "_calib_scratch.json")
        self.addCleanup(lambda: self.tmp.unlink(missing_ok=True))

    def test_skips_calibration_when_artifact_already_exists(self):
        self.tmp.write_text(json.dumps({"rescaling": {"parameters": {"tau": 1.0}}}),
                            encoding="utf-8")
        with mock.patch.object(cal, "calibrate") as fake_calibrate:
            result = _stage_calibrate(reference=mock.Mock(), recalibrate=False, path=self.tmp)
        fake_calibrate.assert_not_called()
        self.assertEqual(result["rescaling"]["parameters"]["tau"], 1.0)

    def test_calibrates_when_artifact_missing(self):
        self.assertFalse(self.tmp.exists())
        fake_result = {"rescaling": {"parameters": {"tau": 2.0}}}
        with mock.patch.object(cal, "calibrate", return_value=fake_result) as fake_calibrate, \
             mock.patch.object(cal, "write_artifact",
                               side_effect=lambda result, path: path.write_text(
                                   json.dumps(result), encoding="utf-8")) as fake_write:
            result = _stage_calibrate(reference=mock.Mock(), recalibrate=False, path=self.tmp)
        fake_calibrate.assert_called_once()
        fake_write.assert_called_once()
        self.assertEqual(result["rescaling"]["parameters"]["tau"], 2.0)

    def test_recalibrate_flag_forces_rebuild_even_when_artifact_exists(self):
        self.tmp.write_text(json.dumps({"rescaling": {"parameters": {"tau": 1.0}}}),
                            encoding="utf-8")
        fake_result = {"rescaling": {"parameters": {"tau": 9.0}}}
        with mock.patch.object(cal, "calibrate", return_value=fake_result) as fake_calibrate, \
             mock.patch.object(cal, "write_artifact",
                               side_effect=lambda result, path: path.write_text(
                                   json.dumps(result), encoding="utf-8")):
            result = _stage_calibrate(reference=mock.Mock(), recalibrate=True, path=self.tmp)
        fake_calibrate.assert_called_once()
        self.assertEqual(result["rescaling"]["parameters"]["tau"], 9.0)


class RunHaltsOnStageFailureTest(unittest.TestCase):
    def _args(self, **overrides):
        args = argparse.Namespace(
            recalibrate=False, reference_input=".", output=None,
            seed=disc.RPI_DISCRIMINATION_SEED)
        for key, value in overrides.items():
            setattr(args, key, value)
        return args

    def test_discriminate_failure_halts_before_stats_or_report(self):
        with mock.patch.object(rd, "load_reference_data", return_value=mock.Mock()), \
             mock.patch("pipeline._stage_calibrate",
                        return_value={"rescaling": {"parameters": {"tau": 1.0}}}), \
             mock.patch("pipeline._stage_discriminate",
                        side_effect=RuntimeError("boom")), \
             mock.patch("pipeline._stage_stats") as fake_stats, \
             mock.patch("pipeline._stage_report") as fake_report:
            exit_code = run(self._args())
        self.assertEqual(exit_code, 1)
        fake_stats.assert_not_called()
        fake_report.assert_not_called()

    def test_calibrate_failure_halts_before_discriminate(self):
        with mock.patch.object(rd, "load_reference_data", return_value=mock.Mock()), \
             mock.patch("pipeline._stage_calibrate", side_effect=RuntimeError("boom")), \
             mock.patch("pipeline._stage_discriminate") as fake_discriminate:
            exit_code = run(self._args())
        self.assertEqual(exit_code, 1)
        fake_discriminate.assert_not_called()

    def test_successful_run_writes_report_and_returns_zero(self):
        calibration = {"rescaling": {"parameters": {"tau": 3.914}}}
        scored_groups = {"a_random_fakes": [{"rpi_generation": 0.1}]}
        stats_result = {"overall_ranking_supported": False}
        with mock.patch.object(rd, "load_reference_data", return_value=mock.Mock()), \
             mock.patch("pipeline._stage_calibrate", return_value=calibration), \
             mock.patch("pipeline._stage_discriminate", return_value=scored_groups), \
             mock.patch("pipeline._stage_stats", return_value=stats_result), \
             mock.patch("pipeline._stage_report",
                        return_value=Path("output/pipeline/fake.json")) as fake_report:
            exit_code = run(self._args())
        self.assertEqual(exit_code, 0)
        fake_report.assert_called_once_with(
            calibration, scored_groups, stats_result, mock.ANY)


class StageReportTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(self._testMethodName + "_report_scratch.json")
        self.addCleanup(lambda: self.tmp.unlink(missing_ok=True))

    def test_writes_valid_json_with_tau_summary_and_pass_criterion(self):
        calibration = {"rescaling": {"parameters": {"tau": 3.914}}}
        scored_groups = {
            "a_random_fakes": [{"rpi_generation": 0.1}, {"rpi_generation": 0.2}],
        }
        stats_result = {"overall_ranking_supported": False}
        path = _stage_report(calibration, scored_groups, stats_result, self.tmp)
        self.assertEqual(path, self.tmp)
        written = json.loads(self.tmp.read_text(encoding="utf-8"))
        self.assertAlmostEqual(written["tau"], 3.914)
        self.assertIn("a_random_fakes", written["group_summary"])
        self.assertEqual(written["groups"], scored_groups)
        self.assertEqual(written["pass_criterion"], stats_result)


class MainDispatchTest(unittest.TestCase):
    def test_no_args_dispatches_to_run(self):
        with mock.patch("pipeline.run", return_value=0) as fake_run:
            exit_code = main([])
        fake_run.assert_called_once()
        self.assertEqual(exit_code, 0)

    def test_explicit_run_command_dispatches_to_run(self):
        with mock.patch("pipeline.run", return_value=0) as fake_run:
            main(["run", "--recalibrate"])
        args = fake_run.call_args[0][0]
        self.assertTrue(args.recalibrate)

    def test_validate_command_passes_through_remaining_args(self):
        with mock.patch.object(vd, "main", return_value=0) as fake_validate:
            exit_code = main(["validate", "--json", "out.json"])
        fake_validate.assert_called_once_with(["--json", "out.json"])
        self.assertEqual(exit_code, 0)

    def test_diagnose_e_command_dispatches_with_no_args(self):
        with mock.patch.object(diag_e, "main") as fake_diag:
            exit_code = main(["diagnose-e"])
        fake_diag.assert_called_once_with()
        self.assertEqual(exit_code, 0)

    def test_diagnose_gate0_command_dispatches_with_no_args(self):
        with mock.patch.object(diag_gate0, "main") as fake_diag:
            exit_code = main(["diagnose-gate0"])
        fake_diag.assert_called_once_with()
        self.assertEqual(exit_code, 0)


if __name__ == "__main__":
    sys.exit(main())
