# Module consolidation — implementation record (branch 05-module-consolidation)

Consolidates the six RPI toolchain modules (`reassortment_detector.py`,
`rpi_calibration.py`, `rpi_discrimination.py`, `rpi_discrimination_stats.py`,
`validate_reassortment_detector.py`, `pipeline.py`) into one real module,
`rpi_core.py`, with the six original filenames becoming thin re-export shims.
**Absolute rule honored throughout: zero behavior change.** Every fixup below is
wiring required for six separate namespaces to coexist correctly in one shared
namespace — no line of RPI logic (`RPI_WEIGHTS`, tau, gate checks, seeds, output
columns) was touched. Verified, not assumed: see Phase 3 below.

## Phase 0 — golden baseline

`output/consolidation_baseline/`: full test suite (322 tests, green,
`test_log.txt`), the seeded discrimination panel at seed 20260829
(`discrimination_panel_golden.json`), and `python pipeline.py --help` /
`python pipeline.py run --help` output — the reference every later phase diffs
against. Committed `468e90f` before any merge code was written.

## Phase 1 — mechanical merge

Built by `build_rpi_core.py` (kept in the repo for reproducibility/review — it
regenerates `rpi_core.py` byte-identically from the six untouched source files;
its exact fixup counts are enforced in the script itself and recorded here)
that:

1. Parses each of the six files' AST, strips (a) its own top-of-file import
   statements, (b) its own leading module docstring, (c) its own
   `if __name__ == "__main__":` guard block — everything else is copied
   **verbatim, unedited**.
2. Concatenates the six bodies behind one deduplicated import block, in
   dependency order.
3. Applies the minimum set of identifier renames + reference fixups a shared
   namespace requires (full list below).
4. Appends a self-aliasing footer so every original `import X as alias`
   qualified reference (`rd.foo`, `cal.foo`, `vd.foo`, `disc.foo`, `stats.foo`)
   keeps resolving correctly.

### Section order — deviates from the literal instruction, and why

Instructed order was "stats, calibration, detector, discrimination, validate,
pipeline". Actual order used: **detector, calibration, discrimination_stats,
discrimination, validate, pipeline**.

This is not a logic change — it is where each unedited block of code physically
sits in the file — but it deviates from the literal word order given, so it is
flagged explicitly per the task's own "if a line of logic seems to need
changing, STOP and report" spirit, applied here to structure rather than logic.

Reason: an AST scan (checking `FunctionDef.args.defaults`/`kw_defaults`
specifically, since Python evaluates parameter defaults in the *enclosing*
scope at `def`-time, not deferred into the function body — a distinction the
first, naive version of this scan missed) found 4 places where a cross-module
reference is evaluated **immediately at module-load time**, not deferred to a
function call:

| file | line | expression |
|---|---|---|
| `rpi_calibration.py` | 337 | `window_nt: int = rd.RPI_PANHANDLE_WINDOW_NT` |
| `rpi_calibration.py` | 756 | `window_nt: int = rd.RPI_PANHANDLE_WINDOW_NT` |
| `rpi_discrimination.py` | 119 | `artifact_path: Path = cal.DEFAULT_ARTIFACT_PATH` |
| `pipeline.py` | 72 | `path: Path = cal.DEFAULT_ARTIFACT_PATH` |

These force **detector before calibration** and **calibration before
discrimination and pipeline** — the self-aliasing trick (below) only helps
*deferred* (function-body) cross-references; it cannot rescue a name that must
already exist at the exact line Python is executing during module load.
`validate` and `discrimination_stats` have no such constraint (confirmed: an
AST scan for ALL cross-module `Attribute` nodes outside function/lambda bodies,
including class bodies, decorators, base classes, and parameter defaults, found
exactly these 4 and no others) so they were placed to track the instructed
order as closely as the hard constraint allows.

### Self-aliasing (handles every remaining cross-reference, unedited)

```python
rd = cal = disc = stats = vd = sys.modules[__name__]
```

Every other qualified cross-module reference in the six files (confirmed by
the same AST scan: every one outside the 4 rows above) is evaluated inside a
function or method body — deferred until the function is actually *called*,
by which point the whole module (including this line) has finished loading.
Patching via `mock.patch.object(rd, "load_reference_data", ...)`-style tests
also keeps working: `rd` now literally *is* the module `load_reference_data`
is defined in, so patching the attribute through the alias patches the exact
same dict entry any caller (qualified or not) resolves against.

### Renames — required because Python has one namespace, not five

Five modules each define their own `main()`; two define their own
`build_panel()`. These are **public** name collisions (not private helpers,
which is what the task's own Phase 1 text anticipated) — merging six modules
into one Python namespace makes this collision unavoidable, not optional: a
module cannot hold five different objects under one bare name `main`. Handled
by generalizing the task's own private-collision convention to these forced
public cases: module-prefixed private names inside `rpi_core.py`, with each
**shim** re-exporting its own one under the original public name `main` /
`build_panel` (Phase 2) — so from every external caller's point of view,
`reassortment_detector.main`, `pipeline.main`, `rpi_calibration.build_panel`,
etc. are exactly as they were.

| original | new name in `rpi_core.py` |
|---|---|
| `reassortment_detector.main` | `_detector_main` |
| `rpi_calibration.main` | `_calibration_main` |
| `rpi_discrimination.main` | `_discrimination_main` |
| `validate_reassortment_detector.main` | `_validate_main` |
| `pipeline.main` | `_pipeline_main` |
| `rpi_calibration.build_panel` | `_calibration_build_panel` |
| `validate_reassortment_detector.build_panel` | `_validate_build_panel` |

`rpi_discrimination_stats.py`'s own `if __name__ == "__main__": unittest.main()`
calls the **standard library's** `unittest.main`, not a project-defined one —
no rename needed there.

### Reference fixups these renames force (self-aliasing cannot cover these)

The self-alias trick only helps when the *calling* code keeps its original
`vd.foo(...)` qualified form. Two call sites had to stop being qualified
(since the target became a renamed, module-scoped private name), and their
own module's self-calls to the pre-rename bare name had to follow:

| file : line (pre-merge) | before | after |
|---|---|---|
| `rpi_calibration.py:233` | `vd.build_panel(builder)` | `_validate_build_panel(builder)` |
| `rpi_calibration.py:762` | `build_panel(...)` (own call) | `_calibration_build_panel(...)` |
| `validate_reassortment_detector.py:674` | `build_panel(builder)` (own call) | `_validate_build_panel(builder)` |
| `pipeline.py:164` | `return vd.main(rest)` | `return _validate_main(rest)` |

### `pipeline.py`'s own test suite — mock-patching had to follow its code

`pipeline.py`'s tests patch internal stage functions via **string path**
(`mock.patch("pipeline._stage_calibrate", ...)`) or via `mock.patch.object`
on the `vd` alias. Once the actual `run()`/`main()`/`_stage_*` functions
physically live in `rpi_core.py`, an unqualified call inside them resolves
against `rpi_core`'s own globals — patching the STRING path `"pipeline.X"`
would silently patch the *shim's* copy instead, which nothing calls, so the
mock would never intercept and the test would either silently exercise real
(slow, un-stubbed) code or fail its `assert_called`-style assertions.

Fixed by switching the 12 string-path patches to the same `__name__ + ".X"`
dynamic-resolution idiom `reassortment_detector.py`'s own test suite already
uses in exactly this situation (`BoltzmannCompatibilityTest`, patching
`__name__ + ".boltzmann_compatibility"`) — not a new pattern, borrowed
precedent already in this codebase:

| target name | occurrences | file:lines (pre-merge) |
|---|---|---|
| `_stage_calibrate` | 3 | `pipeline.py:234,247,258` |
| `_stage_discriminate` | 3 | `pipeline.py:236,248,259` |
| `_stage_stats` | 2 | `pipeline.py:238,260` |
| `_stage_report` | 2 | `pipeline.py:239,261` |
| `run` | 2 | `pipeline.py:291,297` |

Plus the 5 in-file self-calls to `main([...])` in `pipeline.py`'s
`MainDispatchTest` → `_pipeline_main([...])`, and the one
`mock.patch.object(vd, "main", ...)` → `mock.patch.object(vd, "_validate_main", ...)`
(since the call site it exercises, above, no longer goes through the `vd.`
qualified form at all).

Every one of these substitutions was applied via a script that **count-checked**
each fixup (halts on any mismatch rather than silently applying the wrong
number) — none were eyeballed.

### `__doc__` — four CLIs build their `argparse` description from it

`main.py`-style `argparse.ArgumentParser(description=__doc__)` reads the
*calling function's own module's* `__doc__` — inside `rpi_core.py`, that would
be `rpi_core.py`'s own docstring for **all five** merged CLIs, not each one's
original description text. Fixed by capturing each of the four affected
files' original leading docstring **verbatim** (via `repr()` of the AST-parsed
string, not retyped) into a named constant, and updating the 4
`description=__doc__`-style call sites to reference it:

| file | constant | call site (pre-merge) |
|---|---|---|
| `reassortment_detector.py` | `_DETECTOR_MODULE_DOC` | `reassortment_detector.py:3120` |
| `rpi_calibration.py` | `_CALIBRATION_MODULE_DOC` | `rpi_calibration.py:925` |
| `rpi_discrimination.py` | `_DISCRIMINATION_MODULE_DOC` | `rpi_discrimination.py:534` |
| `validate_reassortment_detector.py` | `_VALIDATE_MODULE_DOC` | `validate_reassortment_detector.py:656` |

Verified byte-identical against the originals (`m._DETECTOR_MODULE_DOC ==
reassortment_detector.__doc__`, etc. — all 4 assert equal) and the resulting
`--help` output for all four CLIs is **byte-identical** to `python
<original file> --help` once `sys.argv[0]` is the shim's own name (which it
will be under real invocation — `argparse`'s `prog` defaults from
`sys.argv[0]`, never from `__file__`, so this was already safe by construction
and only needed verifying, not fixing).

`pipeline.py`'s own `--help` path never reads `__doc__` at all (its `run`
subparser uses a hardcoded `prog="pipeline.py run"`, no `description=`) —
confirmed unaffected, matches Phase 0's captured golden output.

### `Path(__file__).parent` — five uses in `reassortment_detector.py`, all safe

All five resolve a *directory*, not a filename, and `rpi_core.py` lives in the
same directory as the shims — `Path(__file__).parent` is identical either way.
No test asserts on the filename itself; confirmed via `grep` (no
`__module__`/`.__name__`/`__qualname__` introspection anywhere in the six
files).

## Verification performed before this phase's commit

- `python -m py_compile rpi_core.py` — clean.
- `import rpi_core` — clean; `rpi_core.rd is rpi_core` (self-alias confirmed
  live).
- All 4 renamed-`__doc__` constants byte-identical to the originals'
  `__doc__`.
- All 4 affected CLIs' `--help` output byte-identical to the originals
  (`sys.argv[0]` corrected to simulate real shim invocation).
- All 7 renamed symbols present, `main`/`build_panel` bare names absent, all 5
  `main` functions distinct objects.
- Manual smoke test: `mock.patch("rpi_core._stage_calibrate", ...)` etc.
  correctly intercepts `rpi_core.run()`'s internal calls (proves the
  mock-string-path fix works, not just compiles).
- **`python -m unittest rpi_core` — 322 tests, green — same count as Phase 0's
  golden baseline**, run directly against the merged module before any shim
  exists.

## What Phase 2 does (shims)

Each of the six original files is replaced with a thin shim that imports, from
`rpi_core`, **exactly the set of top-level names originally defined in that
specific file** (public and private, tests included) — not `from rpi_core
import *`. This is deliberate and necessary: `rpi_core.py`'s namespace holds
every name from all six modules, so an unscoped wildcard import would make
`python -m unittest reassortment_detector` (for example) discover and run
every OTHER module's test classes too, inflating every single-module test
count to all 322 instead of that module's own documented count (e.g.
`rpi_calibration`'s documented 39). Each shim's name list is computed the same
way the collision scan was — walking that ORIGINAL file's own AST — so it is
exactly the original public surface, private helpers other files depend on
included (e.g. `pipeline.py`'s shim must still expose `_stage_discriminate`,
since `capture_raw_scored_groups.py` — outside the six-module set entirely —
calls it by that name).

## Files outside the six that depend on them (checked, unaffected by design)

- `capture_raw_scored_groups.py`: `import pipeline`, `import
  reassortment_detector as rd`; calls `pipeline._stage_discriminate(...)`,
  `rd.load_reference_data(...)`, `rd._json_safe(...)`.
- `diagnose_e_stratification.py`: `from rpi_discrimination_stats import
  compare_groups`.
- `diagnose_gate0_junk.py`: `import reassortment_detector as rd`, `import
  rpi_calibration as cal`; calls `rd._build_precomputed`,
  `cal.genomic_sequence`, `cal.segment_accession`, `rd.CandidateInput`,
  `rd.gate0_check`.

All of these resolve correctly through the shims described above (each name
used here is in its origin file's own top-level name set, so each shim's
explicit re-export list covers it).
