# RPI Fitness Function & Validation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the Reassortment Plausibility Index work end-to-end and demonstrate that it ranks
reals ≈ within-species shuffles > cross-species sympatric shuffles > random shuffles, reproducibly
on Windows and Linux.

**Architecture:** Three bug fixes land first and are independent of everything else. Reproducibility
pins the toolchain. M_pack and M_pol are then turned on for real — sourced from the implicated
parent/donor consensus, never self-comparison — which is what makes the index sequence-dependent.
Validation scores four groups and reports both an RPI ranking and a gate-4 tree separation.

**Tech Stack:** Python 3.14 (system interpreter), numpy 2.5.1, pandas 3.0.5, Biopython, dendropy,
ViennaRNA 2.7.2 (RNAcofold), MAFFT 7.525, IQ-TREE 3.1.3, unittest, matplotlib.

**Spec:** `docs/superpowers/specs/2026-08-29-rpi-fitness-function-design.md`

## Global Constraints

- **Never change existing gate logic, verdicts, or output columns.** `combine_verdict` and every
  `gate*_` column are frozen. RPI stays a parallel column set. `VerdictUnchangedRegressionTest`
  (`reassortment_detector.py:4022`) must stay green at every commit.
- **One commit per workstream.** Never rewrite history; never force-push.
- **If a step fails, log it and continue.** No aborting a run on one bad case.
- M_pack must be **bounded in (0, 1]**, **monotonic** in ddG, **no exact zeros**, **sweepable ±50%**.
- **D1:** calibrate first — bimodal ddG with a clean gap → logistic; otherwise τ-scaled Boltzmann
  with `τ = ddG_span / ln(100)` ≈ 4.86.
- **D2:** RNAcofold window **N = 50 nt**, swept 25 / 50 / 75.
- **D3:** M_pack and M_pol computed from the candidate's terminal windows against the **implicated
  parent/donor consensus** — never self-comparison.
- **D4:** Gate 4 fixture runs in a **scratch dir only**, never touching `output/`.
- **D6:** byte-identity for MAFFT and RPI outputs; **topology + support agreement** for IQ-TREE.
- **D7:** `reals ≈ b1` is an **expected tie**, asserted as equivalence, never as a separation.
- Suite baseline: **169 tests, ~14 min**. Every task ends with the suite green.
- `RPI_WEIGHTS_VERSION` (`reassortment_detector.py:1467`) is bumped in the *same commit* as any
  change to `RPI_WEIGHTS`.

---

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `reassortment_detector.py` | Modify | N4/N5/N6 fixes; τ/window constants; real M_pack + M_pol call sites |
| `rpi_calibration.py` | Create | Pinned ddG calibration panel; derives τ; writes `data/rpi_calibration.json` |
| `rpi_discrimination.py` | Create | Four-group discrimination harness, chimera ladder, sweep driver |
| `validate_reassortment_detector.py` | Modify | SNV × BCCV cross spec; expose `PanelBuilder` |
| `requirements-pip.txt` | Create | System-Python path, reconciled against `environment.yml` |
| `RATIONALE.md` | Create | Design decisions and their reasoning — the defence document |
| `docs/methods-paragraph.md` | Create | Corrected Methods text |

---

# WS1 — Bug fixes (N4, N5, N6)

No dependencies. Independent of D1/D2. Start here.

### Task 1: N4 — non-finite ddG must raise, not silently become maximal penalty

**Files:**
- Modify: `reassortment_detector.py:1265-1301` (`boltzmann_compatibility`)
- Test: `reassortment_detector.py` — `BoltzmannCompatibilityTest`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `boltzmann_compatibility(dg_non_cognate, dg_cognate, rt=None) -> float`, now raising
  `ValueError` when `dg_non_cognate - dg_cognate` is not finite. Task 8 (WS3.3) catches this at
  the call site and records `rpi_m_pack_status = "undefined_energy"`.

**Background:** `max(-700.0, nan)` returns `-700.0` (Python's `max` returns the first argument on
an unordered comparison), so `boltzmann_compatibility(float('nan'), 0.0)` currently returns
`9.86e-305` — maximal penalty. Before the underflow clamp it returned `1.0` — no penalty. Both are
wrong: a nan ddG means RNAcofold returned garbage, and this module's stated rule is that a
silently-defaulted energy must never propagate.

- [ ] **Step 1: Write the failing test**

```python
class BoltzmannNonFiniteTest(unittest.TestCase):
    def test_nan_ddg_raises_rather_than_silently_scoring(self):
        """max(-700.0, nan) is -700.0, so a nan ddG silently became MAXIMAL penalty
        (9.86e-305); before the underflow clamp the same input silently became NO
        penalty (1.0). A nan here means RNAcofold returned garbage -- neither answer is
        defensible, and this module's no-silent-sentinel rule says it must surface."""
        with self.assertRaises(ValueError):
            boltzmann_compatibility(float("nan"), 0.0)
        with self.assertRaises(ValueError):
            boltzmann_compatibility(0.0, float("nan"))
        with self.assertRaises(ValueError):
            boltzmann_compatibility(float("inf"), 0.0)

    def test_finite_inputs_are_unaffected(self):
        """Guard the guard: the raise must not have narrowed the working range."""
        self.assertAlmostEqual(boltzmann_compatibility(0.0, 0.0), 1.0)
        self.assertLess(boltzmann_compatibility(20.0, 0.0), 1.0)
        self.assertGreater(boltzmann_compatibility(20.0, 0.0), 0.0)
        # A far-stronger non-cognate duplex still clamps at 1.0 rather than overflowing.
        self.assertEqual(boltzmann_compatibility(-5000.0, 0.0), 1.0)
        # A far-weaker one still floors above zero rather than underflowing to 0.0.
        self.assertGreater(boltzmann_compatibility(5000.0, 0.0), 0.0)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m unittest reassortment_detector.BoltzmannNonFiniteTest -v`
Expected: FAIL — `test_nan_ddg_raises_rather_than_silently_scoring` reports
`ValueError not raised`.

- [ ] **Step 3: Write minimal implementation**

In `boltzmann_compatibility`, immediately after `ddg = dg_non_cognate - dg_cognate`:

```python
    ddg = dg_non_cognate - dg_cognate
    # A non-finite ddG means RNAcofold returned something unusable. It must NOT be
    # scored. Both available silent answers are wrong and they are opposites:
    # max(-700.0, nan) is -700.0, so nan scored as MAXIMAL penalty (9.86e-305), while
    # before the underflow clamp min(0.0, nan) was 0.0 and it scored as NO penalty
    # (1.0). Raising is the only honest option; the caller in score_candidate catches
    # this and records rpi_m_pack_status = "undefined_energy" rather than aborting.
    if not math.isfinite(ddg):
        raise ValueError(
            f"non-finite ddG ({ddg!r}) from dg_non_cognate={dg_non_cognate!r}, "
            f"dg_cognate={dg_cognate!r} -- RNAcofold returned an unusable energy"
        )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m unittest reassortment_detector.BoltzmannNonFiniteTest -v`
Expected: PASS, 2 tests.

- [ ] **Step 5: Update the docstring's clamp rationale**

The clamp comment currently justifies −700 by reference to "whole-segment duplexes whose MFEs run
to hundreds of kcal/mol" — false per D2 (the window is 50 nt). Replace that clause with:

```
    # -700 keeps exp() in normal range (its smallest normal result is ~1e-308) while
    # still being far below any physically meaningful ddg. The clamp is retained under
    # the 50-nt panhandle window of D2 -- observed ddg runs ~0-22 kcal/mol there, so the
    # clamp is never reached in practice, but it is cheap and it is the guarantee that
    # keeps the no-literal-zeros rule true for adversarial input.
```

- [ ] **Step 6: Run the RPI test classes to confirm nothing regressed**

Run: `python -m unittest reassortment_detector.RpiAssemblyTest reassortment_detector.BoltzmannNonFiniteTest -v`
Expected: PASS.

---

### Task 2: N5 — `_json_safe` must handle numpy scalars

**Files:**
- Modify: `reassortment_detector.py:2500-2513` (`_json_safe`)
- Test: `reassortment_detector.py` — near `WriteRpiReportTest` (`:4219`)

**Interfaces:**
- Consumes: nothing.
- Produces: `_json_safe(value) -> Any`, now normalising `np.floating` and `np.integer` to Python
  scalars and mapping non-finite values to `None`.

**Background — the bug is wider than recorded.** `np.float32('nan')` is not `isinstance(x, float)`
(only `np.float64` subclasses it), so it survives sanitisation; `default=str` then intercepts before
`allow_nan=False` can fire. Verified in this environment:

| input | current output |
|---|---|
| `np.float32('nan')` | `"nan"` — the JSON **string** |
| `np.float32(1.5)` | `"1.5"` — a finite value, also a **string** |
| `np.int64(3)` | `"3"` — also a **string** |

So finite numpy scalars are silently stringified too, not just nan. WS3 introduces numpy-typed
energies into this report, which is exactly what makes all three reachable.

- [ ] **Step 1: Write the failing test**

```python
class JsonSafeNumpyScalarTest(unittest.TestCase):
    def test_numpy_float32_nan_becomes_null_not_the_string_nan(self):
        """np.float32 is not a float subclass (only np.float64 is), so it survived
        _json_safe; default=str then intercepted it before allow_nan=False could fire
        and it serialised as the JSON STRING "nan" -- a value where null was meant."""
        payload = _json_safe({"auc": np.float32("nan")})
        self.assertIsNone(payload["auc"])
        self.assertEqual(json.dumps(payload, allow_nan=False), '{"auc": null}')

    def test_finite_numpy_scalars_survive_as_numbers_not_strings(self):
        """The same type hole stringified FINITE numpy scalars: np.float32(1.5)
        serialised as "1.5" and np.int64(3) as "3"."""
        payload = _json_safe({"f32": np.float32(1.5), "f64": np.float64(2.5),
                              "i64": np.int64(3)})
        self.assertEqual(payload["f32"], 1.5)
        self.assertEqual(payload["f64"], 2.5)
        self.assertEqual(payload["i64"], 3)
        self.assertEqual(json.dumps(payload, sort_keys=True, allow_nan=False),
                         '{"f32": 1.5, "f64": 2.5, "i64": 3}')

    def test_python_types_and_bool_are_unchanged(self):
        """Guard the guard: bool is an int subclass, so a careless numbers.Real check
        would turn True into 1 and change the report's schema."""
        payload = _json_safe({"b": True, "i": 3, "s": "x", "n": None,
                              "f": 1.5, "nan": float("nan")})
        self.assertIs(payload["b"], True)
        self.assertEqual(payload["i"], 3)
        self.assertEqual(payload["s"], "x")
        self.assertIsNone(payload["n"])
        self.assertEqual(payload["f"], 1.5)
        self.assertIsNone(payload["nan"])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m unittest reassortment_detector.JsonSafeNumpyScalarTest -v`
Expected: FAIL — `test_numpy_float32_nan_becomes_null_not_the_string_nan` reports the value is
`np.float32('nan')`, not `None`.

- [ ] **Step 3: Write minimal implementation**

Replace `_json_safe`'s first branch:

```python
    # np.floating BEFORE the plain-float branch, and np.integer alongside it. np.float32
    # is NOT a float subclass (only np.float64 is), so it fell through every check here,
    # and json's default=str then intercepted it before allow_nan=False could fire --
    # np.float32('nan') serialised as the JSON *string* "nan", and finite np.float32 /
    # np.int64 stringified too ("1.5", "3"). Converting to a Python scalar fixes both the
    # non-finite case (-> null) and the finite one (-> a number).
    #
    # bool is checked FIRST because it is an int subclass: without this guard True would
    # be normalised to 1 and the report's schema would change.
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        value = float(value)
        return value if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m unittest reassortment_detector.JsonSafeNumpyScalarTest -v`
Expected: PASS, 3 tests.

- [ ] **Step 5: Assert on the shipped artifact, not the intermediate**

Per `CLAUDE.md`'s recurring-defect note, add to `WriteRpiReportTest` a check that reads the written
JSON back off disk with a parser that rejects the bare tokens Python's `json.loads` accepts:

```python
    def test_written_report_contains_no_stringified_numerics(self):
        """Assert on the row that ships. json.loads accepts bare NaN/Infinity, which is
        exactly why the round-trip test never noticed the stringified-numeric hole; this
        checks the raw file text instead."""
        with tempfile.TemporaryDirectory() as tmp:
            path = write_rpi_report([], Path(tmp))
            text = path.read_text(encoding="utf-8")
        self.assertNotIn('"nan"', text)
        self.assertNotIn('"inf"', text)
        self.assertNotIn("NaN", text)
        self.assertNotIn("Infinity", text)
        json.loads(text, parse_constant=_reject_json_constant)
```

with a module-level helper beside `_json_safe`:

```python
def _reject_json_constant(token: str) -> Any:
    """Make json.loads refuse the bare NaN/Infinity tokens it accepts by default.

    Python's json.loads accepts `NaN`, `Infinity` and `-Infinity`, none of which are
    legal JSON -- jq and JavaScript's JSON.parse both reject the file. That leniency is
    why the existing round-trip test could not see a non-conforming report.
    """
    raise AssertionError(f"report contains non-JSON constant {token!r}")
```

- [ ] **Step 6: Run the report tests**

Run: `python -m unittest reassortment_detector.WriteRpiReportTest reassortment_detector.JsonSafeNumpyScalarTest -v`
Expected: PASS.

---

### Task 3: N6 — prove the retention/generation ordering guard is live

**Files:**
- Modify: `reassortment_detector.py:3947-3961`
  (`test_rounding_never_reports_retention_above_generation`)
- Test: same file

**Interfaces:**
- Consumes: nothing.
- Produces: no production change — this task hardens a test only.

**Background:** the existing test passes `floor**5` to *both* columns at the call site, so it cannot
fail via any change to `_rpi_decimals`, and with `odd_segment="M", intact=True` the two columns are
trivially equal (F = 1.0). Unlike its N1 sibling
`test_rounding_never_reports_e_total_above_its_own_components` (`:3920`, which has a `seen_fractional`
guard-the-guard) it has no assertion that the fixture is still exercising anything.

The fix has two halves: a guard-the-guard on the fixture, and a companion test that demonstrates
the split-precision failure mode directly — which is the only way to prove the equal-precision
choice in `compute_rpi` is what prevents it.

- [ ] **Step 1: Write the failing test**

Replace the existing test body and add the companion:

```python
    def test_rounding_never_reports_retention_above_generation(self):
        """F <= 1 makes rpi_retention <= rpi_generation an invariant of the model, but
        rounding only preserves order at EQUAL precision. Reporting retention finer than
        generation put 1.57233e-05 above 1.572e-05 -- an inversion created purely by
        formatting."""
        far_a = _eco_isolate("SpA", "USA", "2001", host_order="Rodentia",
                             host_family="Cricetidae", host_name="Peromyscus maniculatus")
        far_b = _eco_isolate("SpB", "Chile", "2004", host_order="Soricomorpha",
                             host_family="Soricidae", host_name="Sorex araneus")
        seen_strict = False
        for odd_segment in ("S", "M", "L"):
            for intact in (True, False):
                row = compute_rpi(["SpA", "SpB"], [far_a, far_b], odd_segment, 0.9,
                                  1.0, 1.0, "fixed_constant", intact)
                if row["rpi_retention"] < row["rpi_generation"]:
                    seen_strict = True
                self.assertLessEqual(row["rpi_retention"], row["rpi_generation"],
                                     f"{odd_segment}/{intact}")
        # Guard the guard: with F = 1.0 the two columns are trivially equal and the
        # assertion above measures nothing. At least one combination must have F < 1 --
        # i.e. must actually exercise the scaling this invariant is about.
        self.assertTrue(seen_strict,
                        "fixture stopped producing a row where F < 1.0")

    def test_split_precision_is_what_would_invert_the_pair(self):
        """The guard above cannot fail via _rpi_decimals: compute_rpi passes floor**5 to
        BOTH columns, so no change to that function can split them. This companion proves
        the invariant it protects is real rather than vacuous, by reproducing the exact
        inversion the equal-precision choice prevents.

        1.57233e-05 rounded at 10 dp is 1.57233e-05; the same value rounded at 8 dp is
        1.572e-05. Reporting retention at the finer precision put it ABOVE a generation
        it can never actually exceed."""
        raw = 1.572330e-05
        generation_at_8 = round(raw, 8)
        retention_at_10 = round(raw, 10)
        self.assertGreater(retention_at_10, generation_at_8)
        # And at equal precision -- what compute_rpi actually does -- it does not invert.
        self.assertLessEqual(round(raw, 10), round(raw, 10))
        self.assertEqual(_rpi_decimals(RPI_WEIGHTS["eco_floor"] ** 5), 10)
```

- [ ] **Step 2: Run test to verify the new assertions are meaningful**

Run: `python -m unittest reassortment_detector.RpiAssemblyTest -v`
Expected: PASS. Then verify the guard-the-guard is live by temporarily changing the loop to
`for odd_segment in ("M",):` and `for intact in (True,):` — the test must now FAIL with
"fixture stopped producing a row where F < 1.0". **Revert the temporary change.**

- [ ] **Step 3: Confirm no production code changed**

Run: `git diff --stat reassortment_detector.py`
Expected: only test-region lines plus Tasks 1–2's changes; no change inside `compute_rpi`.

---

### Task 4: WS1 verification and commit

- [ ] **Step 1: Run the full suite**

Run: `python -m unittest reassortment_detector`
Expected: PASS. Baseline was 169 tests; expect **175** (Task 1 adds 2, Task 2 adds 4 including the
artifact check, Task 3 adds 1 net). Record the actual count — do not assert completion until this
run is observed green.

- [ ] **Step 2: Run the validation harness to confirm no behavioural drift**

Run: `python validate_reassortment_detector.py --json output/validation/validation_report.json`
Expected: sensitivity 18/18, specificity 27/27, blind spot 5/5, adversarial 4/4 — unchanged.

- [ ] **Step 3: Commit**

```bash
git add reassortment_detector.py
git commit -m "Fix N4 nan-ddG flip, N5 numpy JSON hole, N6 vacuous ordering guard"
```

---

# WS2 — Reproducibility

Independent of WS1; may run in parallel. Requires D4 (approved).

### Task 5: `requirements-pip.txt` and the environment reconciliation

**Files:**
- Create: `requirements-pip.txt`
- Modify: `environment.yml`
- Test: `reassortment_detector.py` — new `EnvironmentReconciliationTest`

- [ ] **Step 1: Record the actual installed versions**

```bash
python -c "import sys, numpy, pandas, dendropy, Bio; print(sys.version.split()[0], numpy.__version__, pandas.__version__, dendropy.__version__, Bio.__version__)"
python -c "import RNA; print('viennarna', RNA.__version__)"
```

- [ ] **Step 2: Write `requirements-pip.txt`**

Header must state which interpreter it describes and why it exists, so it does not read as a rival
to `environment.yml`:

```
# System-Python path (Python 3.14, Windows), NOT the hanta_phylo conda env.
#
# environment.yml describes the conda environment. This file describes the interpreter
# the project is ACTUALLY run on in the authoring environment, where conda/mamba are not
# on PATH. The two are reconciled here rather than left to contradict each other
# silently: every pin below must be version-compatible with environment.yml's, and
# EnvironmentReconciliationTest asserts that.
#
# mafft and iqtree are NOT installable from PyPI. `pymafft` is an unrelated package --
# do not install it. Those two come from conda or from vendored binaries; see WS2/Task 6.
```

- [ ] **Step 3: Write the reconciliation test**

```python
class EnvironmentReconciliationTest(unittest.TestCase):
    def test_pip_and_conda_pins_do_not_contradict(self):
        """environment.yml and requirements-pip.txt describe two different interpreters.
        They are allowed to differ in scope; they are NOT allowed to pin incompatible
        versions of the same package. Nothing enforced that before this test."""
        ...
```

- [ ] **Step 4: Run it, then commit** — deferred to Task 7's shared commit.

---

### Task 6: Install and pin MAFFT + IQ-TREE (D4)

**Files:**
- Modify: `environment.yml`, `requirements-pip.txt` (comment block naming the binary source)

- [ ] **Step 1: Install** — conda if available, otherwise vendored Windows binaries. Pin
  `mafft=7.525`, `iqtree=3.1.3` to match `environment.yml`.
- [ ] **Step 2: Verify** — `mafft --version`, `iqtree --version`; record exact strings.
- [ ] **Step 3: Gate 4 fixture run in a scratch dir only.** Per D4 this must not touch `output/`.
  Assert it, do not merely intend it:

```python
    def test_gate4_fixture_writes_only_inside_the_scratch_dir(self):
        """D4: the fixture run must never touch output/. Asserted by snapshotting the
        real output tree before and after."""
        before = {p: p.stat().st_mtime_ns for p in Path("output").rglob("*") if p.is_file()}
        ...
        after = {p: p.stat().st_mtime_ns for p in Path("output").rglob("*") if p.is_file()}
        self.assertEqual(before, after)
```

- [ ] **Step 4: Remove the "never reconciled" caveat** from `environment.yml`'s viennarna comment
  and the `_parse_rnacofold_dg` "ASSUMED FORMAT, NOT VALIDATED" docstring — both were resolved
  last session and are now stale.

---

### Task 7: Cross-OS reproducibility fixture (D6)

**Files:**
- Create: `tests/fixtures/repro/` (small: 4 taxa, one segment)
- Test: `reassortment_detector.py` — `ReproducibilityFixtureTest`

- [ ] **Step 1: Seed audit.** `random.Random(seed)` at `validate_reassortment_detector.py:396` is
  seeded; `--phylo-seed` defaults to 42. Every new RNG in WS3/WS4 takes an explicit seed, and the
  seed is written into every artifact.
- [ ] **Step 2: Byte-hash the deterministic artifacts** — MAFFT alignment output and the RPI CSV.
- [ ] **Step 3: Topology + support comparison for IQ-TREE, NOT bytes.** Per D6: IQ-TREE likelihood
  values are not guaranteed bit-identical across OS, compiler, or thread count even with a fixed
  seed. Compare via Robinson–Foulds distance == 0 and per-node support within tolerance. A byte
  hash here would be a promise the pipeline cannot keep.
- [ ] **Step 4: Commit WS2.**

```bash
git add requirements-pip.txt environment.yml tests/fixtures/repro reassortment_detector.py
git commit -m "Pin and verify toolchain; add cross-OS reproducibility fixture"
```

---

# WS3 — M_pack and M_pol scored for real

Gated on WS1. Implements D1, D2, D3.

### Task 8: `rpi_calibration.py` — the pinned ddG panel (D1 branch point)

**Files:**
- Create: `rpi_calibration.py`, `data/rpi_calibration.json`

**Interfaces:**
- Produces: `panhandle_window(genomic_seq, segment, end, n=50) -> str`;
  `calibrate(reference, n=50, seed=...) -> Dict[str, Any]` writing the artifact.

- [ ] **Step 1: Implement the 50-nt panhandle window (D2)** — terminal 50 nt of the 5′ end of one
  segment against the terminal 50 nt of the 3′ end of the other, reverse-complemented.
  `_reverse_complement` already exists at `reassortment_detector.py:415`.
- [ ] **Step 2: Build cognate and non-cognate pair sets** from the reference index. **This panel
  must be disjoint from the WS4 validation groups** — otherwise τ is tuned on the validation set
  and the whole result is circular.
- [ ] **Step 3: Run RNAcofold over the panel; record ViennaRNA version in the artifact.**
- [ ] **Step 4: Test bimodality and branch (D1).** Report the ddG distribution, a dip-statistic or
  two-component Gaussian-mixture BIC comparison, and the gap width. **Bimodal with a clean gap →
  logistic (Option C). Otherwise → τ-scaled Boltzmann with `τ = ddG_span / ln(100)`.**
  Write the chosen form and its parameters into `data/rpi_calibration.json`.
- [ ] **Step 5: Report the branch outcome before implementing Task 9.**

### Task 9: Implement the chosen rescaling

- [ ] Add τ (or ddG₅₀/s) and window N to `RPI_WEIGHTS`; **bump `RPI_WEIGHTS_VERSION`** in the same
  commit.
- [ ] Assert the four required properties as tests: bounded (0,1], monotonic decreasing in ddG,
  no exact zeros at extremes, and ±50% sweepable.

### Task 10: Score M_pack for real (D3)

- [ ] Replace `m_pack=1.0` (`reassortment_detector.py:2456`) with a cognate-vs-non-cognate
  comparison against the **implicated parent/donor consensus** sourced from
  `reference.genomic_index`. **Never self-comparison.**
- [ ] Add `rpi_m_pack_status`: `scored` / `fixed_constant` / `undefined_energy` /
  `no_parental_sequence`. `undefined_energy` is where Task 1's `ValueError` is caught — log and
  continue, per the global constraint.

### Task 11: Score M_pol for real (D3)

- [ ] Stop passing `terminal_windows` as both arguments to `polymerase_compatibility`
  (`reassortment_detector.py:2431`). Source real parental termini.
- [ ] Delete the forced `m_pol_status="fixed_constant"` override so `m_pol_variance_check` can do
  its job.
- [ ] **`VerdictUnchangedRegressionTest` must stay green** — RPI must not feed back into verdicts.

---

# WS4 — Validation (core deliverable)

Gated on WS3.

### Task 12: SNV × BCCV cross spec (D5)

- [ ] Resolve BCCV's ICTV species name against the reference index (it is absent from
  `PUBLISHED_CROSSES`; only counts exist in `RPI_INVITRO_BENCHMARK`).
- [ ] Add the cross spec. If no complete triplet exists, **log a skip — do not silently drop.**

### Task 13: `rpi_discrimination.py` — the four groups (D7)

- [ ] **(a) random fakes** — different species *and* different continents. Weak test; E alone
  separates it. Report it as such.
- [ ] **(b1) within-species shuffles** — **plausibility-ceiling control, expected HIGH.**
- [ ] **(b2) cross-species sympatric** — same host order, country, overlapping years.
  **The real test of M_pack + M_pol.**
- [ ] **(c) real positives** — 8 supported candidates + published lab pairs.

### Task 14: Pass criterion, both halves

- [ ] **RPI ranking:** `reals ≈ b1 > b2 > random`. Separations reported with Mann–Whitney U +
  effect size. **`reals ≈ b1` asserted as an equivalence interval, never as a separation** — a
  non-significant U-test is absence of evidence, not evidence of the tie being claimed.
- [ ] **Detector separation:** gate 4 tree calls separate reals from b1.
- [ ] **Honest-failure clause:** if b2 does not separate, report the measurement and the
  limitation. **Do not retune τ against the validation set.**

### Task 15: Chimera ladder and ±50% sweep

- [ ] Planted chimeras across the observed divergence range; report the detection limit.
- [ ] Extend `sensitivity_sweep` (`reassortment_detector.py:1932`) to cover τ and N (25/50/75).

---

# WS5 — Outputs

- [ ] **Task 16:** results tables (CSV) + plots (per-group distributions, chimera ladder, sweep
  tornado).
- [ ] **Task 17:** `docs/methods-paragraph.md` — corrected Methods text.
- [ ] **Task 18:** `RATIONALE.md` — why τ is not a temperature; why min–max was rejected; why the
  window is 50 nt; why IQ-TREE is compared on topology not bytes; **what b1 measures and why RPI
  scoring it high is correct, not a false positive.**
- [ ] **Task 19:** update `CLAUDE.md`'s RESUME section and the README's RPI limitations list —
  items about M_pack/M_pol being hardcoded at 1.0 are now stale.
