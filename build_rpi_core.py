"""Phase 1 mechanical merge script: builds rpi_core.py from the six RPI toolchain
modules. Pure code motion + the minimum set of surgical reference-fixups required for
a single shared namespace to behave identically to six separate modules -- every
fixup documented here AND in the Phase 1 commit message. Re-run to regenerate
identically from the untouched six source files (Phase 2 replaces them with shims).
"""
import ast
import re

SECTIONS = [
    ("reassortment_detector.py", "DETECTOR"),
    ("rpi_calibration.py", "CALIBRATION"),
    ("rpi_discrimination_stats.py", "DISCRIMINATION_STATS"),
    ("rpi_discrimination.py", "DISCRIMINATION"),
    ("validate_reassortment_detector.py", "VALIDATE"),
    ("pipeline.py", "PIPELINE"),
]

BANNERS = {
    "reassortment_detector.py": "reassortment_detector.py -- reference data, reference loading, Gates 0-4, the verdict rubric, and RPI scoring.",
    "rpi_calibration.py": "rpi_calibration.py -- tau calibration against a held-out isolate panel.",
    "rpi_discrimination_stats.py": "rpi_discrimination_stats.py -- the reals~b1>b2>random pass criterion (Mann-Whitney U + Cliff's delta).",
    "rpi_discrimination.py": "rpi_discrimination.py -- the four-group discrimination panel builder and scorer.",
    "validate_reassortment_detector.py": "validate_reassortment_detector.py -- confusion-matrix validation against published reassortants.",
    "pipeline.py": "pipeline.py -- single orchestrating entry point over the whole toolchain.",
}


def read(fname):
    with open(fname, encoding="utf-8") as f:
        return f.read()


def extract_and_strip(fname):
    src = read(fname)
    lines = src.splitlines(keepends=True)
    tree = ast.parse(src)

    docstring = None
    doc_range = None
    import_ranges = []
    guard_range = None

    for i, node in enumerate(tree.body):
        if i == 0 and isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) \
                and isinstance(node.value.value, str):
            docstring = node.value.value
            doc_range = (node.lineno, node.end_lineno)
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            import_ranges.append((node.lineno, node.end_lineno))
        if isinstance(node, ast.If):
            test = node.test
            if (isinstance(test, ast.Compare) and isinstance(test.left, ast.Name)
                    and test.left.id == "__name__"):
                guard_range = (node.lineno, node.end_lineno)

    remove_ranges = list(import_ranges)
    if doc_range is not None:
        remove_ranges.append(doc_range)
    if guard_range is not None:
        remove_ranges.append(guard_range)

    remove_lines = set()
    for start, end in remove_ranges:
        for ln in range(start, end + 1):
            remove_lines.add(ln)

    import_text_lines = []
    for start, end in import_ranges:
        for ln in range(start, end + 1):
            import_text_lines.append(lines[ln - 1])

    kept = [line for idx, line in enumerate(lines, start=1) if idx not in remove_lines]
    body = "".join(kept)
    return docstring, import_text_lines, body


def sub_exact(text, pattern, replacement, expected_count, label):
    # Replacement passed as a callable so re.sub inserts it VERBATIM -- no
    # backslash-escape template processing (which would otherwise turn a
    # literal `\n` in the replacement text into an actual newline character).
    new_text, n = re.subn(pattern, lambda m: replacement, text, count=0, flags=re.MULTILINE)
    if n != expected_count:
        raise SystemExit(
            f"FIXUP MISMATCH [{label}]: expected {expected_count} substitutions, got {n}\n"
            f"pattern={pattern!r}")
    return new_text


# ---------------------------------------------------------------------------
# Extract all six files
# ---------------------------------------------------------------------------
results = {}
all_import_lines = []
for fname, prefix in SECTIONS:
    doc, imports, body = extract_and_strip(fname)
    results[fname] = {"doc": doc, "body": body, "prefix": prefix}
    all_import_lines.extend(imports)

# ---------------------------------------------------------------------------
# Deduplicated import block (order-preserving, first-seen wins), EXCLUDING:
#
# (a) each module's imports of ONE ANOTHER (e.g. "import
#     validate_reassortment_detector as vd") -- these five are exactly what
#     the self-alias footer replaces. Kept inside rpi_core.py, they cause a
#     circular import once Phase 2's shims exist (rpi_core imports a shim,
#     which imports names FROM rpi_core, which hasn't finished executing yet).
#
# (b) pipeline.py's two EXTERNAL imports (diagnose_e_stratification,
#     diagnose_gate0_junk) -- these are genuinely separate modules, not part
#     of this merge, but diagnose_e_stratification.py itself does `from
#     rpi_discrimination_stats import compare_groups`, reaching back into the
#     shim ecosystem. That is harmless for pipeline.py alone (a standalone
#     leaf importer), but becomes a genuine CYCLE once pipeline's code shares
#     one module with rpi_discrimination_stats's: importing the
#     rpi_discrimination_stats SHIM first pulls in all of rpi_core (pipeline's
#     own section included), which would need diagnose_e_stratification fully
#     loaded, which needs the rpi_discrimination_stats shim fully loaded --
#     the same still-in-progress import that started the chain. No physical
#     placement inside rpi_core.py fixes this (it isn't a definition-order
#     problem, unlike (a) or the def-time-default cases -- it's a cycle
#     between two shims). Fixed instead by making both imports LAZY: moved
#     from module top-level into _pipeline_main's own body (and the two
#     dispatch tests that patch them), at the exact point each is used -- see
#     the pipeline-specific fixups below.
#
# Neither exclusion was assumed safe -- both were caught by actually running
# the shims' own test suites/imports, not by reasoning alone.
# ---------------------------------------------------------------------------
SIBLING_IMPORTS = {
    "import reassortment_detector as rd",
    "import validate_reassortment_detector as vd",
    "import rpi_calibration as cal",
    "import rpi_discrimination as disc",
    "import rpi_discrimination_stats as stats",
}
DEFERRED_EXTERNAL_IMPORTS = [
    "import diagnose_e_stratification as diag_e",
    "import diagnose_gate0_junk as diag_gate0",
]
seen = set()
dedup_imports = []
for line in all_import_lines:
    key = line.rstrip("\n")
    if key in DEFERRED_EXTERNAL_IMPORTS:
        continue
    if key and key not in seen and key not in SIBLING_IMPORTS:
        seen.add(key)
        dedup_imports.append(key)

# ---------------------------------------------------------------------------
# Per-file fixups (main/build_panel collision renames, __doc__ -> named
# constant, and the specific cross-module/reference-tracking call sites that
# collision renames force). Every substitution is count-checked -- a mismatch
# halts the script rather than silently applying the wrong number of edits.
# ---------------------------------------------------------------------------

detector_body = results["reassortment_detector.py"]["body"]
detector_body = sub_exact(
    detector_body, r"^def main\(argv=None\) -> int:$", "def _detector_main(argv=None) -> int:",
    1, "detector: def main rename")
detector_body = sub_exact(
    detector_body, re.escape("description=__doc__)"), "description=_DETECTOR_MODULE_DOC)",
    1, "detector: description=__doc__")
results["reassortment_detector.py"]["body"] = detector_body

cal_body = results["rpi_calibration.py"]["body"]
cal_body = sub_exact(
    cal_body, r"^def main\(argv: Optional\[Sequence\[str\]\] = None\) -> int:$",
    "def _calibration_main(argv: Optional[Sequence[str]] = None) -> int:",
    1, "calibration: def main rename")
cal_body = sub_exact(
    cal_body, r"^def build_panel\($", "def _calibration_build_panel(",
    1, "calibration: def build_panel rename")
cal_body = sub_exact(
    cal_body, re.escape("panel, _skipped = vd.build_panel(builder)"),
    "panel, _skipped = _validate_build_panel(builder)",
    1, "calibration: vd.build_panel call site")
cal_body = sub_exact(
    cal_body, re.escape("observations, panel_provenance = build_panel("),
    "observations, panel_provenance = _calibration_build_panel(",
    1, "calibration: own build_panel call site")
cal_body = sub_exact(
    cal_body, re.escape("description=__doc__.splitlines()[1])"),
    "description=_CALIBRATION_MODULE_DOC.splitlines()[1])",
    1, "calibration: description=__doc__")
results["rpi_calibration.py"]["body"] = cal_body

disc_body = results["rpi_discrimination.py"]["body"]
disc_body = sub_exact(
    disc_body, r"^def main\(argv: Optional\[Sequence\[str\]\] = None\) -> int:$",
    "def _discrimination_main(argv: Optional[Sequence[str]] = None) -> int:",
    1, "discrimination: def main rename")
disc_body = sub_exact(
    disc_body, re.escape("description=__doc__)"), "description=_DISCRIMINATION_MODULE_DOC)",
    1, "discrimination: description=__doc__")
results["rpi_discrimination.py"]["body"] = disc_body

validate_body = results["validate_reassortment_detector.py"]["body"]
validate_body = sub_exact(
    validate_body, r"^def main\(argv=None\) -> int:$", "def _validate_main(argv=None) -> int:",
    1, "validate: def main rename")
validate_body = sub_exact(
    validate_body, re.escape("def build_panel(builder: PanelBuilder,"),
    "def _validate_build_panel(builder: PanelBuilder,",
    1, "validate: def build_panel rename")
validate_body = sub_exact(
    validate_body, re.escape("panel, skipped = build_panel(builder)"),
    "panel, skipped = _validate_build_panel(builder)",
    1, "validate: own build_panel call site")
validate_body = sub_exact(
    validate_body, re.escape('description=__doc__.split("\\n")[0])'),
    'description=_VALIDATE_MODULE_DOC.split("\\n")[0])',
    1, "validate: description=__doc__")
results["validate_reassortment_detector.py"]["body"] = validate_body

pipeline_body = results["pipeline.py"]["body"]
pipeline_body = sub_exact(
    pipeline_body, r"^def main\(argv: Optional\[Sequence\[str\]\] = None\) -> int:$",
    "def _pipeline_main(argv: Optional[Sequence[str]] = None) -> int:",
    1, "pipeline: def main rename")
pipeline_body = sub_exact(
    pipeline_body, re.escape("return vd.main(rest)"), "return _validate_main(rest)",
    1, "pipeline: vd.main call site")
# The 5 test-suite self-calls to pipeline's own main([...]).
pipeline_body = sub_exact(
    pipeline_body, r"(?<![_.\w])main\(\[", "_pipeline_main([",
    5, "pipeline: test self-calls to main([...])")
pipeline_body = sub_exact(
    pipeline_body, re.escape('mock.patch.object(vd, "main", return_value=0) as fake_validate'),
    'mock.patch.object(vd, "_validate_main", return_value=0) as fake_validate',
    1, "pipeline: mock.patch.object(vd, 'main', ...)")
# String-path mock.patch("pipeline.X", ...) -> mock.patch(__name__ + ".X", ...):
# the code these patch now lives in rpi_core.py, not the pipeline.py shim, so the
# patch target must resolve dynamically to wherever this test code actually runs.
for name, count in [
    ("_stage_calibrate", 3), ("_stage_discriminate", 3),
    ("_stage_stats", 2), ("_stage_report", 2), ("run", 2),
]:
    pipeline_body = sub_exact(
        pipeline_body, re.escape(f'"pipeline.{name}"'), f'__name__ + ".{name}"',
        count, f"pipeline: mock.patch string path for {name}")
# diag_e/diag_gate0: originally module-top imports of two genuinely SEPARATE
# (non-merged) modules. diagnose_e_stratification.py itself reaches back into
# the shim ecosystem (`from rpi_discrimination_stats import compare_groups`),
# which is fine when pipeline.py is its own standalone module (no dependency
# in the other direction) but becomes a real circular import once pipeline's
# code and rpi_discrimination_stats's code share ONE module: importing the
# rpi_discrimination_stats shim first pulls in all of rpi_core (pipeline's
# section included), which would need diagnose_e_stratification fully loaded,
# which needs rpi_discrimination_stats's shim fully loaded -- but that's the
# still-in-progress import that started this chain. No section reordering
# fixes this (it isn't a def-time-vs-call-time timing issue, it's a genuine
# cycle across two shims). Fixed by making the import LAZY -- local to
# _pipeline_main's body, at the exact dispatch branch that needs it, so
# `import rpi_core` (via any shim) never touches diagnose_e_stratification /
# diagnose_gate0_junk at module-load time at all. The two dispatch tests that
# patch `diag_e.main`/`diag_gate0.main` via a pre-existing module-level name
# get their own local import instead (safe at TEST-EXECUTION time -- by then
# rpi_core/all shims have finished loading; diagnose_e_stratification.py
# itself has no dependency back on rpi_core at module level, only on
# rpi_discrimination_stats, which is fully loaded by the time tests run).
pipeline_body = sub_exact(
    pipeline_body, re.escape('    if command == "diagnose-e":\n        diag_e.main()\n'),
    '    if command == "diagnose-e":\n'
    '        import diagnose_e_stratification as diag_e\n'
    '        diag_e.main()\n',
    1, "pipeline: lazy diag_e import in _pipeline_main")
pipeline_body = sub_exact(
    pipeline_body, re.escape('    if command == "diagnose-gate0":\n        diag_gate0.main()\n'),
    '    if command == "diagnose-gate0":\n'
    '        import diagnose_gate0_junk as diag_gate0\n'
    '        diag_gate0.main()\n',
    1, "pipeline: lazy diag_gate0 import in _pipeline_main")
pipeline_body = sub_exact(
    pipeline_body, re.escape('        with mock.patch.object(diag_e, "main") as fake_diag:'),
    '        import diagnose_e_stratification as diag_e\n'
    '        with mock.patch.object(diag_e, "main") as fake_diag:',
    1, "pipeline: lazy diag_e import in its dispatch test")
pipeline_body = sub_exact(
    pipeline_body, re.escape('        with mock.patch.object(diag_gate0, "main") as fake_diag:'),
    '        import diagnose_gate0_junk as diag_gate0\n'
    '        with mock.patch.object(diag_gate0, "main") as fake_diag:',
    1, "pipeline: lazy diag_gate0 import in its dispatch test")
results["pipeline.py"]["body"] = pipeline_body

# ---------------------------------------------------------------------------
# Module-doc constants (named, so the 4 argparse `description=__doc__` call
# sites -- now renamed above to reference these constants -- keep showing
# each ORIGINAL module's own docstring text, not rpi_core.py's).
# ---------------------------------------------------------------------------
doc_constants = []
DOC_CONST_NAMES = {
    "reassortment_detector.py": "_DETECTOR_MODULE_DOC",
    "rpi_calibration.py": "_CALIBRATION_MODULE_DOC",
    "rpi_discrimination.py": "_DISCRIMINATION_MODULE_DOC",
    "validate_reassortment_detector.py": "_VALIDATE_MODULE_DOC",
}
for fname, const_name in DOC_CONST_NAMES.items():
    doc_text = results[fname]["doc"]
    doc_constants.append(f"{const_name} = {doc_text!r}\n")

# ---------------------------------------------------------------------------
# Assemble rpi_core.py
# ---------------------------------------------------------------------------
out = []
out.append('"""rpi_core.py -- consolidated RPI fitness-validation toolchain.\n\n')
out.append("Single real module behind six thin shim files (reassortment_detector.py,\n")
out.append("rpi_calibration.py, rpi_discrimination.py, rpi_discrimination_stats.py,\n")
out.append("validate_reassortment_detector.py, pipeline.py), each re-exporting its own\n")
out.append("original public names from here. Built by mechanical merge (Phase 1 of\n")
out.append("branch 05-module-consolidation) -- see notes_consolidation.md for the exact\n")
out.append("renames and reference-fixups this required, and each shim file for what it\n")
out.append('re-exports.\n"""\n\n')

out.append("# --- Consolidated imports (deduplicated, order of first appearance across\n")
out.append("# the six source modules: detector, calibration, discrimination_stats,\n")
out.append("# discrimination, validate, pipeline) ---\n")
for line in dedup_imports:
    out.append(line + "\n")
out.append("\n")

out.append("# --- Original per-module docstrings, preserved as named constants so each\n")
out.append("# module's own argparse `description=` text is unchanged (see below) ---\n")
for line in doc_constants:
    out.append(line)
out.append("\n")

for i, (fname, prefix) in enumerate(SECTIONS):
    banner = BANNERS[fname]
    out.append("#" * 79 + "\n")
    out.append(f"# {banner}\n")
    out.append("#" * 79 + "\n\n")
    out.append(results[fname]["body"])
    out.append("\n\n")

    if fname == "reassortment_detector.py":
        # Cross-module aliases placed HERE, right after the base (detector)
        # section, not at the end of the file. Detector itself uses none of
        # rd/cal/disc/stats/vd (it has no cross-module references at all), so
        # this is the earliest point every later section's def-time-evaluated
        # cross-references (rpi_calibration.py's `window_nt: int =
        # rd.RPI_PANHANDLE_WINDOW_NT` parameter defaults, evaluated
        # immediately when that `def` executes, NOT deferred into the
        # function body) can already find `rd` bound AND
        # `RPI_PANHANDLE_WINDOW_NT` already defined on it. A footer placed
        # after ALL SIX sections (the first version of this script) leaves
        # `rd` completely undefined at the exact line those defaults need it
        # -- caught by actually importing the shims, not assumed safe. Every
        # LATER section's own def-time cross-references (rpi_discrimination.py
        # and pipeline.py's `cal.DEFAULT_ARTIFACT_PATH`) are satisfied the
        # same way, since calibration's section runs to completion before
        # either of those.
        out.append("#" * 79 + "\n")
        out.append("# Cross-module aliases: every original module's own `import X as alias`\n")
        out.append("# (rd/cal/disc/stats/vd) resolves here to THIS module, since every name any\n")
        out.append("# of the six files reached via such an alias is defined either directly\n")
        out.append("# above (detector) or in a LATER section below -- safe for later sections'\n")
        out.append("# own def-time-evaluated defaults (e.g. rpi_calibration.py's `window_nt:\n")
        out.append("# int = rd.RPI_PANHANDLE_WINDOW_NT`) because `rd` is bound here, before any\n")
        out.append("# of them execute, and safe for every body-deferred (call-time) reference\n")
        out.append("# regardless of section order, verified with an AST scan before this merge.\n")
        out.append("# The two exceptions (main, build_panel colliding as bare top-level names)\n")
        out.append("# were fixed by renaming instead, since a self-alias can't disambiguate two\n")
        out.append("# different things sharing one attribute name on the same object.\n")
        out.append("#" * 79 + "\n")
        out.append("rd = cal = disc = stats = vd = sys.modules[__name__]\n\n")

final = "".join(out)
with open("rpi_core.py", "w", encoding="utf-8") as f:
    f.write(final)

print(f"wrote rpi_core.py, {len(final.splitlines())} lines")
