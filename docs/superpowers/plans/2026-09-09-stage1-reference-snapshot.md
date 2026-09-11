# Stage 1 Reference Snapshot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `prepare_dataset.py` (Stage 1) the only file that reads `all_hanta/ncbi_dataset/`; every downstream file (`phylogeny.py`, `rpi_core.py`, `Algorithm.py`, `figures.py`) consumes a JSON reference snapshot Stage 1 writes, instead of independently re-parsing the raw NCBI package.

**Architecture:** Phase 1 is a pure mechanical rename (`01-prepare_dataset.py` → `prepare_dataset.py`, `import all_hanta as ah` → `import prepare_dataset as prep`, every `ah.*` call site → `prep.*`) with zero behavior change — every file still reads `all_hanta/` directly afterward. Phase 2 adds three new JSON outputs to Stage 1 (`ortho_reference_isolates.json`, `ortho_reference_genomic.json`, `ortho_reference_cds.json`) plus a `prep.load_reference_data(output_dir)` reader, then rewrites each downstream file's ingest to call that reader instead of touching `all_hanta/`.

**Tech Stack:** Python 3, stdlib `unittest` (one `TestCase`-per-concern block embedded at the bottom of each module, run via `python -m unittest <module> -v`), `pandas` for CSV writing, plain `json` for the new snapshot format.

**Spec:** `docs/superpowers/specs/2026-09-09-stage1-reference-snapshot-design.md`

## Global Constraints

- Alias for the renamed module is `prep`, never `pd` (collides with `import pandas as pd`, already used by `prepare_dataset.py`) and never `ah` (ties it to the dead module name).
- Renames use word-boundary substitution only (e.g. `sed 's/\bah\./prep./g'`), never a blind replace of the letters "ah" — they appear inside unrelated identifiers and prose.
- The raw data folder `all_hanta/` (no hyphen, no `.py`) keeps its name and location. Never rename it. `README.md`'s download/unzip commands referencing it are untouched.
- Only `prepare_dataset.py` may read `all_hanta/ncbi_dataset/data` after Phase 2. Every other file gets its isolate/sequence data from `output/ortho_reference_*.json`.
- Snapshot format is plain JSON via the existing `safe_write_json` helper (atomic write + sha256 for the checksum manifest) — no pickle, no parquet, no new I/O helper.
- Stage 1 stays filter-free: it writes every isolate/accession, unfiltered. The 206-coding-complete, 272-terminally-complete, and strict-merged-set filters are computed downstream (Stage 2 / Stage 3), never in `prepare_dataset.py`.
- `parse_cds_fasta` (the function) lives in `prepare_dataset.py` after Phase 2 — Stage 1 must not import from Stage 2 (`phylogeny.py`) to read `cds.fna`.
- Comments/docstrings that name the module `all_hanta.py` get reworded to `prepare_dataset.py`. Comments/docstrings that mean the data folder `all_hanta/` are left alone.
- Every task ends with `python -m unittest <module> -v` passing (326 tests for `rpi_core`, per current README) before moving to the next task.

---

## Phase 1: Rename + alias fix (behavior-preserving)

### Task 1: Rename `01-prepare_dataset.py` → `prepare_dataset.py`

**Files:**
- Modify (rename): `01-prepare_dataset.py` → `prepare_dataset.py`
- Modify: `prepare_dataset.py` (self-referential docstrings only)

**Interfaces:**
- Produces: an importable module `prepare_dataset` with the exact same public symbols `01-prepare_dataset.py` had (no code changes in this task beyond the rename and doc text below).

- [ ] **Step 1: Rename the file with git so history is preserved**

```bash
git mv 01-prepare_dataset.py prepare_dataset.py
```

- [ ] **Step 2: Fix the module's own stale run-instructions comment**

In `prepare_dataset.py`, find (near the bottom, in the `# Tests` banner comment):

```python
# Merged from test_all_hanta.py. Run with: python3 -m unittest all_hanta -v
```

Replace with:

```python
# Merged from test_all_hanta.py. Run with: python3 -m unittest prepare_dataset -v
```

(Leave `test_all_hanta.py` as-is — it names a historical source file this test code was merged from, not the current module.)

- [ ] **Step 3: Confirm the module imports cleanly under its new name**

```bash
python -c "import prepare_dataset"
```

Expected: no output, exit code 0.

- [ ] **Step 4: Run its test suite**

```bash
python -m unittest prepare_dataset -v
```

Expected: all tests pass, same count as before the rename (this task changed no logic).

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "$(cat <<'EOF'
Rename 01-prepare_dataset.py to prepare_dataset.py

01-prepare_dataset.py was not a legal Python module name (leading
digit, hyphen), which silently broke every downstream import all_hanta
as ah statement. This rename alone fixes nothing downstream yet --
Task 2-4 update the import sites.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 2: `phylogeny.py` — import + alias rename (`ah` → `prep`)

**Files:**
- Modify: `phylogeny.py` (48 call sites as of this plan's writing, plus the `import` line and prose)

**Interfaces:**
- Consumes: `prepare_dataset` module (Task 1), unchanged public symbols.
- Produces: no change to `phylogeny.py`'s own public API — this task is a pure identifier rename.

- [ ] **Step 1: Confirm `prep` is not already used in this file**

```bash
grep -n "\bprep\b" phylogeny.py
```

Expected: no output. If there is output, STOP and report the collision instead of proceeding.

- [ ] **Step 2: Change the import line**

In `phylogeny.py`, find:

```python
import all_hanta as ah
```

Replace with:

```python
import prepare_dataset as prep
```

- [ ] **Step 3: Rename every `ah.` call site to `prep.`, word-boundary only**

```bash
sed -i 's/\bah\./prep./g' phylogeny.py
```

This also fixes occurrences inside string type-hints (e.g. `"ah.IsolateRecord"` → `"prep.IsolateRecord"`) and docstrings that literally read `ah.PHYLOGROUP` etc., since those already use the `ah.` token this substitution targets.

- [ ] **Step 4: Reword prose that names the module `all_hanta.py`**

Each of these lines currently reads "all_hanta.py" meaning the *module* (not the data folder) — reword each to "prepare_dataset.py". Locations (line numbers as of this plan; re-find by content if they've shifted):

- Module docstring (`phylogeny.py:7,9,10,13`): "divergence diagnostics all_hanta.py (Stage 1)...", "Independent of all_hanta.py's output files...", "the same way all_hanta.py's own run_pipeline did...", "Reuses all_hanta.py's ingest primitives via `import all_hanta as ah`." — rewrite this whole paragraph to describe the *new* behavior instead of the old one (this module no longer re-derives from raw input or reuses ingest primitives against a raw path — that language becomes wrong in Phase 2, but for this task just fix the module name references; Task 9 rewrites this docstring's substance).
- `phylogeny.py:111`: "all_hanta.py's k-mer diagnostics..." → "prepare_dataset.py's k-mer diagnostics..."
- `phylogeny.py:395`: "not all_hanta.py's ortho_pipeline.log" → "not prepare_dataset.py's ortho_pipeline.log"
- `phylogeny.py:926`: "this fed all_hanta.py's surveillance_priority(...)" → "this fed prepare_dataset.py's surveillance_priority(...)"
- `phylogeny.py:958`: "Mirrors all_hanta.py's run_pipeline()..." → "Mirrors prepare_dataset.py's run_pipeline()..."

Do not touch `phylogeny.py:1009`'s `"Package root containing all_hanta/ncbi_dataset/data"` in this task — that names the *data folder*, not the module; it is handled in Task 9 when the CLI itself changes.

- [ ] **Step 5: Verify zero `ah.` leftovers**

```bash
grep -n "\bah\.\|import all_hanta" phylogeny.py
```

Expected: no output.

- [ ] **Step 6: Run the test suite**

```bash
python -m unittest phylogeny -v
```

Expected: all tests pass (this task is behavior-preserving — `prep` is `ah` under a new name and `all_hanta/` is still read exactly as before).

- [ ] **Step 7: Commit**

```bash
git add phylogeny.py
git commit -m "$(cat <<'EOF'
phylogeny.py: rename all_hanta import alias ah -> prep

Mechanical rename following prepare_dataset.py's own rename (previous
commit). No behavior change -- phylogeny.py still reads all_hanta/
directly; Task 9 removes that.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 3: `rpi_core.py` — import + alias rename (`ah` → `prep`)

**Files:**
- Modify: `rpi_core.py` (87 call sites as of this plan's writing, plus the `import` line and prose)

**Interfaces:**
- Consumes: `prepare_dataset` module (Task 1).
- Produces: no API change — pure rename.

- [ ] **Step 1: Confirm `prep` is not already used in this file**

```bash
grep -n "\bprep\b" rpi_core.py
```

Expected: no output. STOP and report if there is.

- [ ] **Step 2: Change the import line**

```python
import all_hanta as ah
```
→
```python
import prepare_dataset as prep
```

(Leave the adjacent `import phylogeny as ph` line untouched — that alias is unaffected by this rename.)

- [ ] **Step 3: Rename every `ah.` call site to `prep.`**

```bash
sed -i 's/\bah\./prep./g' rpi_core.py
```

- [ ] **Step 4: Reword prose naming the module `all_hanta.py`**

`rpi_core.py` has ~15 comment/docstring lines saying "all_hanta.py" meaning the module (e.g. `_DETECTOR_MODULE_DOC`'s "Imports all_hanta.py's ingest and phylogenetics machinery", the `load_reference_data` comment "replicates the ingest portion of all_hanta.py's run_pipeline", the gate-comments citing `all_hanta.py:86`/`all_hanta.py:90`/`all_hanta.py:1560-1565`, `SEGMENT_LENGTH_ABSURD_MAX stays a hard ceiling ... already works in all_hanta.py`, `classify_reassortment_call` reference, `SPECIES_INFO` reference, and the `all_hanta.py genome-wide phylogenetic surveillance` source string at line ~8435). For each: reword "all_hanta.py" → "prepare_dataset.py" where it names the module or one of its functions/constants.

Do **not** reword the ~20 lines that mean the *data folder* `all_hanta/` (e.g. "real sequences in `all_hanta/`", "REAL all_hanta/ reference set", "package root containing all_hanta/ncbi_dataset/data", "test_all_hanta.py's _make_isolate" — that's a historical test-file name reference, leave it) — those are handled by Task 10, not this one. When in doubt: does the sentence make sense with "prepare_dataset.py" substituted? If it's talking about code/functions, yes, reword it. If it's talking about where files live on disk, no, leave it.

- [ ] **Step 5: Verify zero `ah.` leftovers**

```bash
grep -n "\bah\.\|import all_hanta" rpi_core.py
```

Expected: no output.

- [ ] **Step 6: Run the test suite**

```bash
python -m unittest rpi_core -v
```

Expected: all 326 tests pass.

- [ ] **Step 7: Commit**

```bash
git add rpi_core.py
git commit -m "$(cat <<'EOF'
rpi_core.py: rename all_hanta import alias ah -> prep

Mechanical rename, no behavior change -- rpi_core.py still reads
all_hanta/ directly via load_reference_data; Task 10 removes that.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 4: `figures.py` — import + alias rename (`ah` → `prep`)

**Files:**
- Modify: `figures.py` (16 call sites)

**Interfaces:**
- Consumes: `prepare_dataset` module (Task 1).
- Produces: no API change — pure rename.

- [ ] **Step 1: Confirm `prep` is not already used**

```bash
grep -n "\bprep\b" figures.py
```

Expected: no output.

- [ ] **Step 2: Change the import line**

```python
import all_hanta as ah
```
→
```python
import prepare_dataset as prep
```

(Leave `import phylogeny as ph` untouched.)

- [ ] **Step 3: Rename every `ah.` call site to `prep.`**

```bash
sed -i 's/\bah\./prep./g' figures.py
```

- [ ] **Step 4: Reword the one module-naming prose line**

`figures.py:18`: "`output/ortho_reassortment_table.csv` (all_hanta.py's own genome-wide phylogenetic surveillance)" → "...(prepare_dataset.py's own genome-wide phylogenetic surveillance)". (This is actually describing `phylogeny.py`'s Stage 2 output historically attributed to the old monolith — leave the attribution as "prepare_dataset.py" per the mechanical rule; Task 12 can further correct attribution if needed, but that's a documentation nuance outside this task's scope.)

Leave `figures.py:558,564,565` alone in this task — those name the data folder `all_hanta/` and its smoke-test path; Task 12 changes them when the raw-access removal happens.

- [ ] **Step 5: Verify zero `ah.` leftovers**

```bash
grep -n "\bah\.\|import all_hanta" figures.py
```

Expected: no output.

- [ ] **Step 6: Run the test suite**

```bash
python -m unittest figures -v
```

Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add figures.py
git commit -m "$(cat <<'EOF'
figures.py: rename all_hanta import alias ah -> prep

Mechanical rename, no behavior change.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 5: `Algorithm.py` — confirm no `ah` usage, run Phase 1 acceptance check

**Files:**
- Read-only check: `Algorithm.py` (uses `rc.load_reference_data`, never imports `all_hanta` directly — confirmed via `grep -c "\bah\." Algorithm.py` returning 0 during design)

**Interfaces:**
- Consumes: nothing new.
- Produces: Phase 1's completion signal.

- [ ] **Step 1: Confirm `Algorithm.py` has no `ah.` usage to rename**

```bash
grep -n "\bah\.\|import all_hanta" Algorithm.py
```

Expected: no output (it only ever went through `rpi_core.load_reference_data`, never imported `all_hanta` itself).

- [ ] **Step 2: Run Phase 1's full acceptance grep**

```bash
grep -rn "\bah\.\|import all_hanta" --include="*.py" .
```

Expected: no output across the whole repo. If anything remains, find and fix it (repeat Steps 2-5 of the relevant Task 2/3/4 pattern for that file) before proceeding.

- [ ] **Step 3: Run every module's test suite one more time, together**

```bash
python -m unittest prepare_dataset -v
python -m unittest phylogeny -v
python -m unittest rpi_core -v
python -m unittest figures -v
python -m unittest Algorithm -v
```

Expected: all green. This is the Phase 1 exit gate — every file still reads `all_hanta/` directly (unchanged behavior), but nothing references the dead `all_hanta` module name anymore.

- [ ] **Step 4: Commit (only if Step 1-2 required fixes; otherwise this task produced no diff — skip the commit)**

If fixes were needed:

```bash
git add -A
git commit -m "$(cat <<'EOF'
Phase 1 acceptance: zero remaining all_hanta module references

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

## Phase 2: Snapshot artifacts + raw-access removal

### Task 6: `prepare_dataset.py` — relocate `parse_cds_fasta` from `phylogeny.py`

**Files:**
- Modify: `prepare_dataset.py` (add the function + its test)
- Modify: `phylogeny.py` (remove the function + its test, update call site)

**Interfaces:**
- Produces: `prep.parse_cds_fasta(path: Optional[Path]) -> Dict[str, List[Tuple[str, str]]]` — same signature and behavior as `phylogeny.py`'s current one.
- Consumes (in `phylogeny.py` after this task): `prep.parse_cds_fasta` instead of a locally-defined function.

Only `parse_cds_fasta` moves — `select_primary_cds` stays in `phylogeny.py` (it doesn't read raw files; it operates on already-parsed `(product, sequence)` tuples) and `rpi_core.py` keeps calling `ph.select_primary_cds`.

- [ ] **Step 1: Copy `parse_cds_fasta` into `prepare_dataset.py`**

In `phylogeny.py`, the function (as of this plan, lines 176-206) reads:

```python
def parse_cds_fasta(path: Optional[Path]) -> Dict[str, List[Tuple[str, str]]]:
    """Index cds.fna by NUCLEOTIDE accession -- not protein accession.

    Unlike protein.faa (keyed by protein accessionVersion), cds.fna's defline carries the
    parent nucleotide accession plus a coordinate suffix (e.g. 'PX492530.1:1-430'),
    because a CDS is a coordinate range on its parent nucleotide record. A single
    nucleotide accession can carry more than one CDS -- 134 in this package, mostly
    ambisense S segments encoding both N and NSs on overlapping reading frames, plus a
    few M segments annotated with separate G1/G2 CDS alongside a combined precursor
    entry. parse_fasta's setdefault-first-wins index would silently keep whichever CDS
    happened to be read first and drop the other, so this keeps every distinct
    (product, sequence) pair per accession and leaves selection to select_primary_cds.
    """
    index: DefaultDict[str, List[Tuple[str, str]]] = defaultdict(list)
    seen: set = set()
    for desc, seq in ah._iter_fasta_records(path):
        if not desc:
            continue
        tokens = desc.split(None, 1)
        acc = tokens[0].split(":")[0]
        base = ah.accession_base(acc)
        product = tokens[1].split("[")[0].strip() if len(tokens) > 1 else ""
        seq = seq.upper()
        if not seq or not base:
            continue
        dedupe_key = (base, product, seq)
        if dedupe_key in seen:
            continue
        seen.add(dedupe_key)
        index[base].append((product, seq))
    return dict(index)
```

(Note: after Task 2's rename, `ah._iter_fasta_records`/`ah.accession_base` in this function read as `prep._iter_fasta_records`/`prep.accession_base` inside `phylogeny.py` — but since this function is moving *into* `prepare_dataset.py`, those become plain unqualified `_iter_fasta_records`/`accession_base` calls, since both are already defined in that same module.)

Paste this function into `prepare_dataset.py` immediately after `parse_fasta` (so it sits next to its sibling FASTA-parsing function), with the two `prep.`-qualified calls made unqualified:

```python
def parse_cds_fasta(path: Optional[Path]) -> Dict[str, List[Tuple[str, str]]]:
    """Index cds.fna by NUCLEOTIDE accession -- not protein accession.

    Unlike protein.faa (keyed by protein accessionVersion), cds.fna's defline carries the
    parent nucleotide accession plus a coordinate suffix (e.g. 'PX492530.1:1-430'),
    because a CDS is a coordinate range on its parent nucleotide record. A single
    nucleotide accession can carry more than one CDS -- 134 in this package, mostly
    ambisense S segments encoding both N and NSs on overlapping reading frames, plus a
    few M segments annotated with separate G1/G2 CDS alongside a combined precursor
    entry. parse_fasta's setdefault-first-wins index would silently keep whichever CDS
    happened to be read first and drop the other, so this keeps every distinct
    (product, sequence) pair per accession and leaves selection to select_primary_cds
    (phylogeny.py).
    """
    index: DefaultDict[str, List[Tuple[str, str]]] = defaultdict(list)
    seen: set = set()
    for desc, seq in _iter_fasta_records(path):
        if not desc:
            continue
        tokens = desc.split(None, 1)
        acc = tokens[0].split(":")[0]
        base = accession_base(acc)
        product = tokens[1].split("[")[0].strip() if len(tokens) > 1 else ""
        seq = seq.upper()
        if not seq or not base:
            continue
        dedupe_key = (base, product, seq)
        if dedupe_key in seen:
            continue
        seen.add(dedupe_key)
        index[base].append((product, seq))
    return dict(index)
```

- [ ] **Step 2: Delete the function from `phylogeny.py`**

Remove the copy of `parse_cds_fasta` from `phylogeny.py` (lines 176-206 as of this plan).

- [ ] **Step 3: Find and move its test**

```bash
grep -n "class.*ParseCds\|def test.*cds_fasta\|parse_cds_fasta" phylogeny.py
```

Move any `unittest.TestCase` class that tests `parse_cds_fasta` directly (not `select_primary_cds`) from `phylogeny.py`'s test section into `prepare_dataset.py`'s test section, changing its calls from `parse_cds_fasta(...)` (unqualified, since it was in the same file) to still `parse_cds_fasta(...)` (also unqualified, now in `prepare_dataset.py` where the function now lives — no change needed to the test body itself, only which file it lives in).

- [ ] **Step 4: `phylogeny.py`'s call site**

Update the one call site (line ~994 as of this plan, `cds_index = parse_cds_fasta(paths.cds_fasta)`) to `cds_index = prep.parse_cds_fasta(paths.cds_fasta)`. (This is a stepping-stone — Task 9 replaces this whole ingest block with a snapshot load, at which point this call site disappears entirely. It's updated here first so the file keeps working after this task in isolation.)

- [ ] **Step 5: Run both modules' test suites**

```bash
python -m unittest prepare_dataset -v
python -m unittest phylogeny -v
```

Expected: all pass — `prepare_dataset.py` gains the moved test, `phylogeny.py` loses it, functionality identical.

- [ ] **Step 6: Commit**

```bash
git add prepare_dataset.py phylogeny.py
git commit -m "$(cat <<'EOF'
Move parse_cds_fasta from phylogeny.py to prepare_dataset.py

Only Stage 1 (prepare_dataset.py) reads raw NCBI files after this
refactor; phylogeny.py must not import Stage-1-owned raw-file parsing
from Stage 2 duplication or vice versa. This relocates the one
raw-file reader phylogeny.py owned into Stage 1, ahead of removing
phylogeny.py's raw ingest entirely (next task).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 7: `prepare_dataset.py` — write the 3 snapshot files + provenance block

**Files:**
- Modify: `prepare_dataset.py` (`run_pipeline`, `summarize`)

**Interfaces:**
- Produces: `output/ortho_reference_isolates.json`, `output/ortho_reference_genomic.json`, `output/ortho_reference_cds.json`, all added to `ortho_output_checksums.txt`; `ortho_summary.json` gains a `reference_snapshot` key.
- Consumes: `isolates: List[IsolateRecord]`, `genomic_index: Dict[str, FastaRecord]`, `cds_index: Dict[str, List[Tuple[str,str]]]` (all already built earlier in `run_pipeline`; `cds_index` newly built via Task 6's relocated `parse_cds_fasta`).

- [ ] **Step 1: Write a failing test for the new writes**

Add to `prepare_dataset.py`'s test section:

```python
class ReferenceSnapshotWriteTest(unittest.TestCase):
    """run_pipeline must write the 3 reference snapshot files, unfiltered."""

    def test_writes_three_snapshot_files_with_expected_shapes(self):
        isolates = [_make_isolate("iso1", "Orthohantavirus x", "S1.1", "M1.1", "L1.1")]
        genomic_index = {
            "S1.1": FastaRecord(accession="S1.1", accession_base="S1", header="S1.1 h",
                                sequence="ACGT", length=4),
        }
        cds_index = {"S1": [("nucleocapsid protein", "ACGTACGT")]}

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            isolates_digest = safe_write_json(
                out / "ortho_reference_isolates.json", [asdict(i) for i in isolates])
            genomic_digest = safe_write_json(
                out / "ortho_reference_genomic.json",
                [asdict(r) for r in genomic_index.values()])
            cds_digest = safe_write_json(out / "ortho_reference_cds.json", cds_index)

            self.assertTrue((out / "ortho_reference_isolates.json").exists())
            self.assertTrue((out / "ortho_reference_genomic.json").exists())
            self.assertTrue((out / "ortho_reference_cds.json").exists())
            self.assertEqual(len(isolates_digest), 64)  # sha256 hex digest length
            self.assertEqual(len(genomic_digest), 64)
            self.assertEqual(len(cds_digest), 64)

            loaded_isolates = json.loads((out / "ortho_reference_isolates.json").read_text())
            self.assertEqual(loaded_isolates[0]["isolate_key"], "iso1")
            loaded_cds = json.loads((out / "ortho_reference_cds.json").read_text())
            self.assertEqual(loaded_cds["S1"], [["nucleocapsid protein", "ACGTACGT"]])
```

(This test exercises the exact write calls `run_pipeline` will make — it's written against `safe_write_json` directly rather than invoking the whole `run_pipeline` CLI, since `run_pipeline` needs a real/fixture NCBI package on disk to run end-to-end; that end-to-end path is already covered by this module's existing fixture-based tests, e.g. `RealPackageIntegrationTest`-style classes, which Step 4 below re-runs.)

Check `_make_isolate` exists already in this file's test helpers (it's referenced at `rpi_core.py:1116` as a pattern other files copy) — if the exact helper signature differs, adapt the test to whatever isolate-construction helper `prepare_dataset.py`'s own tests already use; do not invent a second one.

- [ ] **Step 2: Run it to confirm it fails**

```bash
python -m unittest prepare_dataset.ReferenceSnapshotWriteTest -v
```

Expected: PASS actually, since this test only exercises `safe_write_json` (already implemented) — it is a specification test for the shape, not a TDD-red test for new code. That's fine: it locks in the exact shape `run_pipeline`'s new calls (Step 3) must produce. Confirm it passes before continuing.

- [ ] **Step 3: Add the writes to `run_pipeline`**

In `prepare_dataset.py`'s `run_pipeline` (around the existing `write_csv`/`write_json` calls, currently ending with `write_csv("ortho_risk_table.csv", ...)`), add, right after the CDS index is built (see Task 6 — if Task 6 already added a `parse_cds_fasta` call to `run_pipeline`'s own ingest, reuse that variable; if not, add it here):

```python
    logging.info("Loading CDS index")
    cds_index = parse_cds_fasta(paths.cds_fasta)
    logging.info("  cds: %s accessions with at least one CDS", len(cds_index))

    logging.info("Writing reference snapshot")
    write_json("ortho_reference_isolates.json", [asdict(i) for i in isolates])
    write_json("ortho_reference_genomic.json", [asdict(r) for r in genomic_index.values()])
    write_json("ortho_reference_cds.json", cds_index)
```

Place this after `isolates = build_isolate_records(accessions)` and before `write_output_checksums(...)` (so the 3 new files' digests land in `output_digests` and get included in `ortho_output_checksums.txt`, same as every other output).

- [ ] **Step 4: Add the provenance block to `summarize()`**

In `summarize()`'s return dict, add a new top-level key (place it after `"counts"`):

```python
        "reference_snapshot": {
            "isolate_count": len(isolates),
            "genomic_accession_count": len(genomic_index),
            "cds_indexed_accession_count": len(cds_index),
            "files": {
                "ortho_reference_isolates.json": output_digests.get("ortho_reference_isolates.json"),
                "ortho_reference_genomic.json": output_digests.get("ortho_reference_genomic.json"),
                "ortho_reference_cds.json": output_digests.get("ortho_reference_cds.json"),
            },
        },
```

This requires threading `genomic_index`, `cds_index`, and `output_digests` into `summarize()`'s parameter list (it currently takes `accessions, isolates, divergence, qc_rows, checksums, curation_issues, paths, cluster_threshold, min_species_size, data_records`). Add three parameters: `genomic_index: Dict[str, FastaRecord]`, `cds_index: Dict[str, List[Tuple[str, str]]]`, `output_digests: Dict[str, str]`, and update `run_pipeline`'s call to `summarize(...)` to pass them.

Note the `sha256` values referenced here are computed by the `write_json` calls in Step 3, which populate `output_digests` — `summarize()` is called after those writes in `run_pipeline`, so the dict will already contain them; confirm the call order in `run_pipeline` is: writes (Step 3) → `summarize(...)` → `write_json("ortho_summary.json", summary)`. This matches `run_pipeline`'s existing order (summary is written after every other output, so it can report on them).

- [ ] **Step 5: Run the module's test suite**

```bash
python -m unittest prepare_dataset -v
```

Expected: all pass, including the new `ReferenceSnapshotWriteTest` and any existing end-to-end fixture test that runs the full `run_pipeline` (which will now additionally assert, or at minimum not break on, the 3 new output files — check for an existing test asserting the exact list of files `run_pipeline` writes, e.g. one that checks `ortho_output_checksums.txt`'s contents or iterates `output/*.csv`/`*.json`; if such a test hardcodes the file list, update it to include the 3 new names).

- [ ] **Step 6: Commit**

```bash
git add prepare_dataset.py
git commit -m "$(cat <<'EOF'
prepare_dataset.py: write the reference snapshot (isolates/genomic/cds)

Three new unfiltered JSON outputs -- ortho_reference_isolates.json,
ortho_reference_genomic.json, ortho_reference_cds.json -- so
downstream stages can stop re-parsing all_hanta/ themselves. Stage 1
stays filter-free: every isolate/accession is written, not just
coding-complete or terminally-complete subsets (those filters are
Stage 2/3's job). ortho_summary.json gains a reference_snapshot
provenance block.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 8: `prepare_dataset.py` — `ReferenceData` + `load_reference_data` reader, with the keying-equality test

**Files:**
- Modify: `prepare_dataset.py` (add `ReferenceData` dataclass, `load_reference_data`, and its tests)

**Interfaces:**
- Produces:
  - `@dataclass class ReferenceData: isolates: List[IsolateRecord]; genomic_index: Dict[str, FastaRecord]; cds_index: Dict[str, List[Tuple[str, str]]]`
  - `def load_reference_data(output_dir: Path) -> ReferenceData`
- Consumes: the 3 files Task 7 writes.

- [ ] **Step 1: Write the failing keying-equality test first**

This is the assertion the spec review required: `load_reference_data(...).genomic_index` must be key-for-key identical to `parse_fasta(...)`'s own in-memory output on the same input, not just "roughly the same size."

```python
class ReferenceDataRoundTripTest(unittest.TestCase):
    """load_reference_data must reconstruct exactly what run_pipeline wrote --
    key-for-key on genomic_index, not just matching counts."""

    def test_genomic_index_keying_matches_parse_fasta_exactly(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta_path = Path(tmp) / "genomic.fna"
            fasta_path.write_text(
                ">PX111111.1 Orthohantavirus x segment S\nACGTACGTACGT\n"
                ">PX222222.1 Orthohantavirus x segment M\nTTTTGGGGCCCC\n",
                encoding="utf-8")
            original_index = parse_fasta(fasta_path)

            out = Path(tmp) / "output"
            out.mkdir()
            safe_write_json(out / "ortho_reference_isolates.json", [])
            safe_write_json(
                out / "ortho_reference_genomic.json",
                [asdict(r) for r in original_index.values()])
            safe_write_json(out / "ortho_reference_cds.json", {})

            reference = load_reference_data(out)

            self.assertEqual(set(reference.genomic_index.keys()), set(original_index.keys()))
            for key, original_record in original_index.items():
                self.assertEqual(reference.genomic_index[key].sequence, original_record.sequence)
                self.assertEqual(reference.genomic_index[key].accession, original_record.accession)

    def test_cds_index_round_trips_as_list_of_pairs(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            safe_write_json(out / "ortho_reference_isolates.json", [])
            safe_write_json(out / "ortho_reference_genomic.json", [])
            safe_write_json(out / "ortho_reference_cds.json",
                            {"PX111111": [("nucleocapsid protein", "ACGT")]})

            reference = load_reference_data(out)

            self.assertEqual(reference.cds_index["PX111111"], [("nucleocapsid protein", "ACGT")])
            self.assertIsInstance(reference.cds_index["PX111111"][0], tuple)

    def test_isolates_reconstruct_as_isolate_record_instances(self):
        isolate = _make_isolate("iso1", "Orthohantavirus x", "S1.1", "M1.1", "L1.1")
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            safe_write_json(out / "ortho_reference_isolates.json", [asdict(isolate)])
            safe_write_json(out / "ortho_reference_genomic.json", [])
            safe_write_json(out / "ortho_reference_cds.json", {})

            reference = load_reference_data(out)

            self.assertEqual(len(reference.isolates), 1)
            self.assertIsInstance(reference.isolates[0], IsolateRecord)
            self.assertEqual(reference.isolates[0].isolate_key, "iso1")
```

- [ ] **Step 2: Run to confirm these fail**

```bash
python -m unittest prepare_dataset.ReferenceDataRoundTripTest -v
```

Expected: FAIL with `NameError: name 'load_reference_data' is not defined` (and `ReferenceData` undefined).

- [ ] **Step 3: Implement `ReferenceData` and `load_reference_data`**

Add near `parse_cds_fasta` (after Task 6) in `prepare_dataset.py`:

```python
@dataclass
class ReferenceData:
    """Everything a downstream stage needs, reconstructed from Stage 1's snapshot --
    never from a raw all_hanta/ read. isolates is the full, UNFILTERED set; Stage 2's
    206-coding-complete and Stage 3's 272-terminally-complete/strict-merged filters are
    applied by their own callers, not here.
    """
    isolates: List[IsolateRecord]
    genomic_index: Dict[str, FastaRecord]
    cds_index: Dict[str, List[Tuple[str, str]]]


def load_reference_data(output_dir: Path) -> ReferenceData:
    """Reads the 3 snapshot files Stage 1's run_pipeline writes. Raises FileNotFoundError
    (not a silent empty result) if Stage 1 hasn't been run against this output_dir yet --
    a missing snapshot is a setup error, not an empty-dataset condition.
    """
    isolates_path = output_dir / "ortho_reference_isolates.json"
    genomic_path = output_dir / "ortho_reference_genomic.json"
    cds_path = output_dir / "ortho_reference_cds.json"
    for path in (isolates_path, genomic_path, cds_path):
        if not path.exists():
            raise FileNotFoundError(
                f"Reference snapshot file missing: {path}. Run prepare_dataset.py "
                f"(Stage 1) against this output directory first.")

    isolates = [IsolateRecord(**row) for row in json.loads(isolates_path.read_text(encoding="utf-8"))]
    genomic_index = {
        row["accession"]: FastaRecord(**row)
        for row in json.loads(genomic_path.read_text(encoding="utf-8"))
    }
    raw_cds = json.loads(cds_path.read_text(encoding="utf-8"))
    cds_index = {k: [tuple(pair) for pair in v] for k, v in raw_cds.items()}

    return ReferenceData(isolates=isolates, genomic_index=genomic_index, cds_index=cds_index)
```

- [ ] **Step 4: Run to confirm the new tests pass**

```bash
python -m unittest prepare_dataset.ReferenceDataRoundTripTest -v
```

Expected: all 3 pass.

- [ ] **Step 5: Run the full module suite**

```bash
python -m unittest prepare_dataset -v
```

Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add prepare_dataset.py
git commit -m "$(cat <<'EOF'
prepare_dataset.py: add ReferenceData + load_reference_data reader

The snapshot reader lives in Stage 1 alongside the writer, since all
three consumers (phylogeny.py, rpi_core.py, figures.py) need it, not
just one. Includes a key-for-key equality test against parse_fasta's
own in-memory keying, per code review: matching counts alone would not
have caught a keying scheme mismatch (e.g. base-accession vs
full-accession, or last-wins vs first-wins on duplicates).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 9: `phylogeny.py` — replace raw ingest with the snapshot load

**Files:**
- Modify: `phylogeny.py` (`run_phylogeny`, `parse_args`, `main`, module docstring)

**Interfaces:**
- Consumes: `prep.load_reference_data(output_dir) -> ReferenceData` (Task 8), `prep.assign_kmer_clusters`, `prep.compute_segment_divergence` (pure functions, unchanged, already renamed in Task 2).
- Produces: `run_phylogeny(output_dir: Path, args: argparse.Namespace) -> None` (signature changes from `(paths: "prep.InputPaths", args)`).

- [ ] **Step 1: Rewrite the module docstring**

Replace the current docstring (lines 1-14) with:

```python
"""Orthohantavirus reassortment phylogenetics stage (MAFFT + IQ-TREE).

Stage 2 of the three-stage pipeline documented in README.md's "Pipeline stages" section.
Analyzes the reference isolate set from Stage 1's snapshot for real, tree-based
reassortment evidence -- the rigorous replacement for the k-mer divergence diagnostics
prepare_dataset.py (Stage 1) still emits as informational-only columns.

Reads ONLY output/ortho_reference_isolates.json and output/ortho_reference_cds.json (via
prep.load_reference_data) -- never all_hanta/ncbi_dataset/ directly. Stage 1 is the sole
reader of the raw NCBI package; this stage re-applies its own k-mer clustering and
divergence diagnostics (assign_kmer_clusters / compute_segment_divergence) to the loaded
isolates, since those need this stage's own --cluster-threshold / --min-species-size
values, which can differ from Stage 1's defaults.
"""
```

- [ ] **Step 2: Rewrite `run_phylogeny`**

Replace the whole function (lines 957-1000 as of this plan) with:

```python
def run_phylogeny(output_dir: Path, args: argparse.Namespace) -> None:
    logging.info("Loading Stage 1 reference snapshot from %s", output_dir)
    reference = prep.load_reference_data(output_dir)
    isolates = reference.isolates
    logging.info("  %s isolates, %s complete triplets",
                 len(isolates), sum(1 for i in isolates if i.has_complete_triplet))

    prep.assign_kmer_clusters(isolates, args.cluster_threshold)
    divergence = prep.compute_segment_divergence(isolates, min_species_size=args.min_species_size)
    divergence_by_key = {d.isolate_key: d for d in divergence}

    cds_index = reference.cds_index
    logging.info("  cds: %s accessions with at least one CDS", len(cds_index))

    logging.info("Running phylogenetics stage")
    phylo_rows = run_phylogenetics_stage(isolates, divergence_by_key, cds_index, output_dir, args)
    prep.safe_write_csv(output_dir / "ortho_reassortment_table.csv", phylo_rows)
    logging.info("Done. %s candidate rows written to ortho_reassortment_table.csv", len(phylo_rows))
```

- [ ] **Step 3: Simplify the CLI — drop `--input`, keep `--output`**

In `parse_args`, replace:

```python
    parser.add_argument("--input", default=".",
                        help="Package root containing all_hanta/ncbi_dataset/data")
    parser.add_argument("--output", help="Output directory; defaults to <input>/output")
```

with:

```python
    parser.add_argument("--output", default="output",
                        help="Stage 1's output directory (contains the reference snapshot "
                             "this stage reads, and where it writes its own results)")
```

- [ ] **Step 4: Update `main`**

Replace:

```python
def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    paths = ah.resolve_inputs(Path(args.input), output=Path(args.output) if args.output else None)
    ah.configure_logging(paths.output_dir / "ortho_phylogeny.log", verbose=args.verbose)
    logging.info("Orthohantavirus phylogenetics stage (Stage 2)")
    logging.info("Package: %s", paths.data_dir)
    logging.info("Output : %s", paths.output_dir)
    run_phylogeny(paths, args)
    return 0
```

with:

```python
def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    prep.configure_logging(output_dir / "ortho_phylogeny.log", verbose=args.verbose)
    logging.info("Orthohantavirus phylogenetics stage (Stage 2)")
    logging.info("Output: %s", output_dir)
    run_phylogeny(output_dir, args)
    return 0
```

(Note: after Task 2's rename this file already says `prep.` for these calls — this step is presented with `prep.` throughout since Task 2 landed first.)

- [ ] **Step 4b: Persist the tip-label map as `labels.json`**

In `_analyze_group` (see the code around `label_map: Dict[str, PhyloCandidate] = dict(zip(labels, members))`), add immediately after that line:

```python
    prep.safe_write_json(group_dir / "labels.json",
                         {label: candidate.isolate.isolate_key for label, candidate in label_map.items()})
```

This writes one `labels.json` per group (shared across its S/M/L segment subdirectories, since the label assignment is the same for all three), so `figures.py` (Task 12) can relabel tree tips without re-deriving the grouping via `collect_phylo_candidates`/`group_candidates` and risking a `--cluster-threshold`/`--min-species-size` mismatch against whatever values this Stage 2 run actually used.

- [ ] **Step 4c: Test that `labels.json` is written and matches the label map**

```python
class GroupLabelsPersistedTest(unittest.TestCase):
    """_analyze_group must persist label -> isolate_key so figures.py never has to
    reconstruct group membership to relabel a committed tree."""

    def test_labels_json_matches_label_assignment(self):
        # Reuses whatever fixture members/mafft-iqtree-faking pattern this file's
        # existing _analyze_group tests already use (see the PhylogeneticsScopeEndToEndTest
        # -style fake-binary tests) -- do not invent a second fixture pattern. After
        # calling _analyze_group with a fake mafft/iqtree pair against a temp group_dir_root:
        group_dir = group_dir_root / prep.normalize_attr_name(group_key)
        labels_path = group_dir / "labels.json"
        self.assertTrue(labels_path.exists())
        written = json.loads(labels_path.read_text(encoding="utf-8"))
        expected = {f"I{i:04d}": m.isolate.isolate_key for i, m in enumerate(members)}
        self.assertEqual(written, expected)
```

Adapt this to whichever existing fake-mafft/fake-iqtree test scaffolding `phylogeny.py` already has for `_analyze_group` (do not build a second, parallel fixture mechanism) — the assertion that matters is `labels.json`'s content equals `{f"I{i:04d}": m.isolate.isolate_key ...}` for the same `members` ordering `_analyze_group` received.

- [ ] **Step 5: Fix the `--phylo-outgroup`/`ah.PHYLOGROUP` help-text reference and the `--input` help removal's downstream mention**

Search for any remaining `all_hanta/ncbi_dataset/data` mention in `phylogeny.py`'s CLI help text (the `--input` help string just removed in Step 3 was the only one) — confirm none remain:

```bash
grep -n "all_hanta" phylogeny.py
```

Expected: no output (the module docstring rewrite in Step 1 already removed the last prose references).

- [ ] **Step 6: Update tests that constructed `ah.InputPaths`/called `run_phylogeny(paths, args)` directly**

```bash
grep -n "InputPaths\|run_phylogeny(" phylogeny.py
```

For each test calling `run_phylogeny(paths, ...)` with a constructed `prep.InputPaths`, change it to call `run_phylogeny(output_dir, ...)` with a plain `Path` to a temp directory, after first writing a small fixture snapshot into that directory (using `prep.safe_write_json` directly, same pattern as Task 8's round-trip tests) — for example:

```python
def _write_fixture_snapshot(output_dir: Path, isolates, cds_index=None) -> None:
    prep.safe_write_json(output_dir / "ortho_reference_isolates.json",
                         [asdict(i) for i in isolates])
    prep.safe_write_json(output_dir / "ortho_reference_genomic.json", [])
    prep.safe_write_json(output_dir / "ortho_reference_cds.json", cds_index or {})
```

Add this helper near the top of `phylogeny.py`'s test section and use it in place of whatever raw-package fixture setup those tests currently do.

- [ ] **Step 7: Run the test suite**

```bash
python -m unittest phylogeny -v
```

Expected: all pass. If a test still references `all_hanta/` as a real on-disk fixture (e.g. a skip-if-missing integration test), convert it to build a fixture snapshot via `_write_fixture_snapshot` instead, per Step 6.

- [ ] **Step 8: Confirm no raw access remains**

```bash
grep -n "all_hanta\|resolve_inputs\|iter_jsonl\|transform_one_record\|build_isolate_records" phylogeny.py
```

Expected: no output.

- [ ] **Step 9: Commit**

```bash
git add phylogeny.py
git commit -m "$(cat <<'EOF'
phylogeny.py: stop reading all_hanta/, load Stage 1's snapshot instead

Deletes the raw re-derivation (resolve_inputs / iter_jsonl / parse_fasta
/ transform_one_record / build_isolate_records) that duplicated Stage
1's ingest. run_phylogeny now takes an output_dir and calls
prep.load_reference_data(output_dir); the CLI drops --input entirely
(there is no raw package to point it at anymore) and --output becomes
both the read and write location.

Also persists each tree group's synthetic tip-label -> isolate_key
map as labels.json, so figures.py (next task) can relabel committed
trees without re-deriving group membership itself.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 10: `rpi_core.py` — replace `ReferenceData`/`load_reference_data` with Stage 1's, fix all CLI defaults

**Files:**
- Modify: `rpi_core.py` (delete local `ReferenceData`/`load_reference_data`; update 4 CLI call sites: `score`, `calibrate`, `validate`, `run`; drop `ph.parse_cds_fasta`)

**Interfaces:**
- Consumes: `prep.ReferenceData`, `prep.load_reference_data` (Task 8).
- Produces: `rpi_core.py` no longer defines its own `ReferenceData`/`load_reference_data` — every internal reference (`rd.load_reference_data`, `ReferenceData(...)` constructor calls in tests) becomes `prep.load_reference_data`/`prep.ReferenceData`.

- [ ] **Step 1: Delete the local definitions**

Remove `rpi_core.py`'s own `ReferenceData` dataclass and `load_reference_data` function (lines 388-414 as of this plan — the block starting `@dataclass\nclass ReferenceData:` through the `return ReferenceData(...)` in `load_reference_data`).

- [ ] **Step 2: Rename every in-file reference**

```bash
grep -n "\bReferenceData\b\|\bload_reference_data\b" rpi_core.py | wc -l
```

Then, word-boundary rename every bare `ReferenceData` → `prep.ReferenceData` and every bare `load_reference_data` → `prep.load_reference_data` EXCEPT where already qualified as `rd.load_reference_data` (since `rd = sys.modules[__name__]` self-aliases this module — `rd.load_reference_data` meant "this module's load_reference_data"; those call sites become `prep.load_reference_data` directly, not `rd.prep.load_reference_data`):

```bash
sed -i 's/\brd\.load_reference_data\b/prep.load_reference_data/g; s/\bload_reference_data(/prep.load_reference_data(/g; s/\bReferenceData(/prep.ReferenceData(/g' rpi_core.py
```

Then fix the double-qualification the second substitution just created on `load_reference_data`'s own (now-deleted) definition line — since Step 1 already removed that definition, there is nothing to fix there, but check for any doubled `prep.prep.` artifacts from running both substitutions:

```bash
grep -n "prep\.prep\.\|prep\.load_reference_data(" rpi_core.py
```

Fix any `prep.prep.` down to `prep.` by hand.

- [ ] **Step 3: `score` subcommand (line ~3219-3220, ~3334)**

Update the help text (no default to change — `required=True`):

```python
    parser.add_argument("--reference-input", required=True,
                        help="Stage 1 output directory containing the reference snapshot")
```

Confirm the call site now reads `reference = prep.load_reference_data(Path(args.reference_input))` (from Step 2's rename).

- [ ] **Step 4: `calibrate` subcommand (line ~7395)**

```python
    parser.add_argument("--root", type=Path, default=Path("all_hanta"),
                        help="dataset root passed to load_reference_data")
```
→
```python
    parser.add_argument("--root", type=Path, default=Path("output"),
                        help="Stage 1 output directory passed to load_reference_data")
```

- [ ] **Step 5: `validate` subcommand (lines ~9494, ~9501-9507)**

```python
    parser.add_argument("--reference-input", default="all_hanta")
    ...
    root = Path(args.reference_input)
    if not (root / "ncbi_dataset" / "data" / "data_report.jsonl").exists():
        print(f"ERROR: {root}/ NCBI package not found.", file=sys.stderr)
        return 2

    print(f"Loading reference from {root}/ ...", flush=True)
    reference = prep.load_reference_data(root)
```

becomes:

```python
    parser.add_argument("--reference-input", default="output")
    ...
    root = Path(args.reference_input)
    if not (root / "ortho_reference_isolates.json").exists():
        print(f"ERROR: {root}/ reference snapshot not found "
              f"(run prepare_dataset.py first).", file=sys.stderr)
        return 2

    print(f"Loading reference from {root}/ ...", flush=True)
    reference = prep.load_reference_data(root)
```

- [ ] **Step 6: `run` subcommand (lines ~9704, ~9726-9727)**

```python
def run(args: argparse.Namespace) -> int:
    reference = rd.load_reference_data(Path(args.reference_input) / "all_hanta")
```
→
```python
def run(args: argparse.Namespace) -> int:
    reference = prep.load_reference_data(Path(args.reference_input) / "output")
```

```python
    parser.add_argument("--reference-input", default=".",
                        help="package root containing all_hanta/ncbi_dataset/data")
```
→
```python
    parser.add_argument("--reference-input", default=".",
                        help="project root containing the output/ directory "
                             "(Stage 1's reference snapshot)")
```

- [ ] **Step 7: Drop `ph.parse_cds_fasta`**

```bash
grep -n "ph\.parse_cds_fasta" rpi_core.py
```

That call site (previously inside the deleted `load_reference_data`, Step 1) is already gone. Confirm no other call site references it:

```bash
grep -n "\bph\." rpi_core.py
```

Expected: only `PhyloCandidate`, `select_primary_cds`, `resolve_binary`, `_analyze_group` remain (these are legitimate shared pure-code reuse, unchanged).

- [ ] **Step 8: Update test fixtures that require the real `all_hanta/` package**

```bash
grep -n 'Path("all_hanta")\|Path(__file__).parent / "all_hanta"\|"--reference-input", "all_hanta"' rpi_core.py
```

For each hit (skip-if-missing integration tests around lines 4209-4302, 5431-5436, 6315-6318, 8611-8626, 9013-9144, 9704-area CLI-arg tests as of this plan):

- Tests that currently do `root = Path(__file__).parent / "all_hanta"` then `self.skipTest("all_hanta/ NCBI package not present...")` if missing: change to `root = Path(__file__).parent / "output"` and `self.skipTest("output/ reference snapshot not present -- run prepare_dataset.py first")` if `(root / "ortho_reference_isolates.json")` doesn't exist. These tests exercise the REAL committed snapshot (once Stage 1 has been run in this checkout), same spirit as before but pointed at the new artifact.
- Tests that build a CLI argv list with `"--reference-input", "all_hanta"`: change the string to `"output"`.
- Tests that call `rc.load_reference_data(Path("all_hanta"))` or `rd.load_reference_data(Path("all_hanta"))` directly: change to `prep.load_reference_data(Path("output"))` (already covered by Step 2's `sed`, but re-verify the `Path("all_hanta")` argument itself is updated — the `sed` in Step 2 only renamed the function/class names, not this string argument).

```bash
sed -i 's/Path("all_hanta")/Path("output")/g' rpi_core.py
```

Then re-check:

```bash
grep -n '"all_hanta"' rpi_core.py
```

Expected: no output (every remaining `"all_hanta"` string, if any, should be a comment/docstring prose reference to the data folder — leave those; there should be no executable-code string literal left).

- [ ] **Step 9: Run the full test suite**

```bash
python -m unittest rpi_core -v
```

Expected: all 326 tests pass (or more, if any got added/split during fixture conversion — none should have been dropped; if a previously-`skipTest`-guarded test now runs because `output/ortho_reference_isolates.json` exists in this checkout, that's expected and fine as long as it passes).

- [ ] **Step 10: Commit**

```bash
git add rpi_core.py
git commit -m "$(cat <<'EOF'
rpi_core.py: read the Stage 1 snapshot instead of all_hanta/ directly

Removes rpi_core's own ReferenceData/load_reference_data in favor of
prepare_dataset.py's (Stage 1 owns both the write and read side of the
format it defines). Updates all 4 CLI entry points (score, calibrate,
validate, run) to default to output/ instead of all_hanta/, and
validate's existence check now looks for the snapshot file instead of
a raw NCBI report. Drops the ph.parse_cds_fasta borrow from
phylogeny.py -- Stage 1 builds cds_index itself now.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 11: `Algorithm.py` — point at `output/` instead of `all_hanta/`

**Files:**
- Modify: `Algorithm.py` (lines ~290, ~412)

**Interfaces:**
- Consumes: `rc.load_reference_data` (which, after Task 10, resolves to `prep.load_reference_data` inside `rpi_core.py` — `Algorithm.py`'s own call site is unaffected by that internal rename since it calls `rc.load_reference_data(...)`, and `rc` is `rpi_core`, which still exposes `load_reference_data` as a name via Task 10's `sed` — confirm this by checking `rpi_core.py` still has a top-level callable named `load_reference_data`, not only `prep.load_reference_data`. If Task 10's rename replaced ALL bare `load_reference_data` occurrences including a re-export, `Algorithm.py` needs `rc.prep.load_reference_data(...)` instead — check which is actually true before editing.)

- [ ] **Step 1: Confirm how `rpi_core.py` exposes the loader after Task 10**

```bash
python -c "import rpi_core as rc; print(rc.load_reference_data)"
```

If this prints a function (because `rpi_core.py` still has its own module-level name `load_reference_data` bound to `prep.load_reference_data`, e.g. via `from prepare_dataset import load_reference_data` or an explicit re-export), use `rc.load_reference_data(...)` below. If it raises `AttributeError`, use `rc.prep.load_reference_data(...)` instead. Use whichever actually works — do not guess.

- [ ] **Step 2: Update the default argument**

```python
    parser.add_argument("--reference-input", default="all_hanta")
```
→
```python
    parser.add_argument("--reference-input", default="output")
```

- [ ] **Step 3: Update the call site (line ~308)**

Confirmed working form from Step 1, e.g.:

```python
    reference = rc.load_reference_data(Path(args.reference_input))
```

(unchanged code — only the CLI default in Step 2 changes the effective path; `Path(args.reference_input)` was always passed straight through).

- [ ] **Step 4: Update the test fixture at line ~412**

```python
    def test_pool_matches_usable_pool(self):
        ...
        from pathlib import Path
        reference = rc.load_reference_data(Path("all_hanta"))
```
→
```python
    def test_pool_matches_usable_pool(self):
        ...
        from pathlib import Path
        reference = rc.load_reference_data(Path("output"))
```

(Or `rc.prep.load_reference_data(...)`, per Step 1's finding.) Check whether this test has a skip-guard for the fixture's existence; if not and it currently relies on `all_hanta/` always being present in CI/dev checkouts, add one mirroring Task 10 Step 8's pattern:

```python
        if not Path("output/ortho_reference_isolates.json").exists():
            self.skipTest("output/ reference snapshot not present -- run prepare_dataset.py first")
```

- [ ] **Step 5: Run the test suite**

```bash
python -m unittest Algorithm -v
```

Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add Algorithm.py
git commit -m "$(cat <<'EOF'
Algorithm.py: --reference-input defaults to output/, not all_hanta/

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 12: `figures.py` — drop group reconstruction, read `labels.json` + tree paths directly

**Files:**
- Modify: `figures.py` (`load_ingest`, `group_for_scope`, `label_map_for_members`, `find_group_for_isolate`, `load_tree_relabeled`, `render_before_figures`, `render_after_figures`, `main`, smoke test)

**Interfaces:**
- Consumes: `prep.load_reference_data` (Task 8); `labels.json` + `s_tree_path`/`m_tree_path`/`l_tree_path` columns in `ortho_reassortment_table.csv`, both written by `phylogeny.py` (Task 9, Step 4b).
- Produces: `load_ingest(output_dir: str = "output") -> List[IsolateRecord]` (shrinks from a 4-tuple to a plain list — grouping/divergence/cds_index are no longer needed here at all). `load_tree_relabeled(tree_path: Path) -> Bio.Phylo.BaseTree.Tree` (drops its `group_dir_root`/`group_key`/`segment`/`label_map` parameters in favor of one path). `group_for_scope`, `label_map_for_members`, `find_group_for_isolate` are deleted.

- [ ] **Step 1: Delete `group_for_scope`, `label_map_for_members`, `find_group_for_isolate`**

Remove these three functions (lines 91-119 as of this plan) entirely — nothing in `figures.py` reconstructs tree groups after this task.

- [ ] **Step 2: Simplify `load_ingest`**

Replace (lines 60-88 as of this plan):

```python
def load_ingest(
    input_dir: str = ".", cluster_threshold: float = 0.02, min_species_size: int = 8,
) -> Tuple[List["ah.IsolateRecord"], Dict[str, "ah.SegmentDivergenceProfile"],
          Dict[str, List[Tuple[str, str]]], "ah.InputPaths"]:
    """Mirrors phylogeny.run_phylogeny's ingest chain through cds_index -- the same
    ~15-line duplication the three-stage split spec already accepts (§3.2), here for a
    fourth reason: figures need the same isolate population Stage 2 scored, not a
    filtered output table.
    """
    paths = ah.resolve_inputs(Path(input_dir))
    data_records = list(ah.iter_jsonl(paths.data_report))
    annotation_index = ah.build_annotation_index(ah.iter_jsonl(paths.annotation_report))
    biosample_index = ah.build_biosample_index(ah.iter_jsonl(paths.biosample_report))
    genomic_index = ah.parse_fasta(paths.genomic_fasta)
    protein_index = ah.parse_fasta(paths.protein_fasta)

    accessions = []
    for record in data_records:
        transformed, _qc_row = ah.transform_one_record(
            record, annotation_index, biosample_index, genomic_index, protein_index)
        if transformed is not None:
            accessions.append(transformed)

    isolates = ah.build_isolate_records(accessions)
    ah.assign_kmer_clusters(isolates, cluster_threshold)
    divergence = ah.compute_segment_divergence(isolates, min_species_size=min_species_size)
    divergence_by_key = {d.isolate_key: d for d in divergence}
    cds_index = ph.parse_cds_fasta(paths.cds_fasta)
    return isolates, divergence_by_key, cds_index, paths
```

with:

```python
def load_ingest(output_dir: str = "output") -> List["prep.IsolateRecord"]:
    """Isolate metadata only (species/country/host, for the AFTER panel's parent-metadata
    table) -- loaded from Stage 1's reference snapshot, never all_hanta/ncbi_dataset/
    directly. Grouping/divergence/cds_index are not needed here: phylogeny.py now
    persists tip-label -> isolate_key maps as labels.json per tree group, and tree paths
    are read straight off ortho_reassortment_table.csv, so this module never reconstructs
    group membership itself (see docs/superpowers/specs/2026-09-09-stage1-reference-
    snapshot-design.md's amendment).
    """
    return prep.load_reference_data(Path(output_dir)).isolates
```

Note the module-level `import phylogeny as ph` may now be unused — check with `grep -n "\bph\." figures.py` after this task's other steps land; if nothing references `ph.` anymore, remove the import.

- [ ] **Step 3: Simplify `load_tree_relabeled`**

Replace (lines 121-132 as of this plan):

```python
def load_tree_relabeled(group_dir_root: Path, group_key: str, segment: str,
                        label_map: Dict[str, "ph.PhyloCandidate"]) -> "Bio.Phylo.BaseTree.Tree":
    from Bio import Phylo
    path = group_dir_root / ah.normalize_attr_name(group_key) / segment / "tree.treefile"
    if not path.exists():
        raise FileNotFoundError(f"committed tree not found: {path}")
    tree = Phylo.read(str(path), "newick")
    for tip in tree.get_terminals():
        candidate = label_map.get(tip.name)
        if candidate is not None:
            tip.name = candidate.isolate.isolate_key
    return tree
```

with:

```python
def load_tree_relabeled(tree_path: Path) -> "Bio.Phylo.BaseTree.Tree":
    """Loads a committed tree.treefile and relabels its synthetic tip names (I0000, ...)
    back to isolate_key, using the labels.json phylogeny.py wrote alongside the group's
    S/M/L segment directories (tree_path.parent is the segment dir, tree_path.parent.parent
    is the group dir labels.json lives in).
    """
    from Bio import Phylo
    if not tree_path.exists():
        raise FileNotFoundError(f"committed tree not found: {tree_path}")
    labels_path = tree_path.parent.parent / "labels.json"
    if not labels_path.exists():
        raise FileNotFoundError(f"labels.json not found next to {tree_path} "
                                f"(expected at {labels_path})")
    label_to_isolate_key = json.loads(labels_path.read_text(encoding="utf-8"))
    tree = Phylo.read(str(tree_path), "newick")
    for tip in tree.get_terminals():
        isolate_key = label_to_isolate_key.get(tip.name)
        if isolate_key is not None:
            tip.name = isolate_key
    return tree
```

- [ ] **Step 4: Rewrite `render_before_figures` to use tree paths from the CSV row**

Replace the ingest/grouping preamble and the per-segment loop body (lines 246-283 as of this plan):

```python
def render_before_figures(output_dir: Path, input_dir: str = ".",
                          reassortment_table: Optional[Path] = None,
                          min_tips: int = DEFAULT_MIN_TIPS) -> List[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import ConnectionPatch

    isolates, divergence_by_key, cds_index, paths = load_ingest(input_dir)
    by_group = {
        "species": group_for_scope(isolates, divergence_by_key, cds_index, "species"),
        "phylogroup": group_for_scope(isolates, divergence_by_key, cds_index, "phylogroup"),
    }
    group_dir_root = {
        "species": paths.output_dir / "phylogenetics",
        "phylogroup": paths.output_dir / "phylogenetics_phylogroup",
    }

    table_path = reassortment_table or (paths.output_dir / "ortho_reassortment_table.csv")
    rows = supported_reassortment_rows(table_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for row in rows:
        isolate_key = row["isolate_key"]
        scope = row["analysis_scope"]
        discordant = row.get("discordant_segment") or ""
        fig, axes = plt.subplots(1, 3, figsize=(13, 5))
        tip_pos: Dict[str, Tuple[float, float]] = {}
        for ax, segment in zip(axes, SEGMENTS):
            group_key, members = find_group_for_isolate(by_group[scope], isolate_key)
            if group_key is None:
                ax.set_axis_off()
                ax.text(0.5, 0.5, f"{segment}: not in a\n{scope}-scope tree",
                        ha="center", va="center", wrap=True, transform=ax.transAxes)
                continue
            label_map = label_map_for_members(members)
            tree = load_tree_relabeled(group_dir_root[scope], group_key, segment, label_map)
```

with:

```python
def render_before_figures(output_dir: Path, data_dir: str = "output",
                          reassortment_table: Optional[Path] = None,
                          min_tips: int = DEFAULT_MIN_TIPS) -> List[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import ConnectionPatch

    stage_dir = Path(data_dir)
    table_path = reassortment_table or (stage_dir / "ortho_reassortment_table.csv")
    rows = supported_reassortment_rows(table_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for row in rows:
        isolate_key = row["isolate_key"]
        discordant = row.get("discordant_segment") or ""
        fig, axes = plt.subplots(1, 3, figsize=(13, 5))
        tip_pos: Dict[str, Tuple[float, float]] = {}
        for ax, segment in zip(axes, SEGMENTS):
            tree_path_str = row.get(f"{segment.lower()}_tree_path") or ""
            if not tree_path_str:
                ax.set_axis_off()
                ax.text(0.5, 0.5, f"{segment}: not in a committed tree",
                        ha="center", va="center", wrap=True, transform=ax.transAxes)
                continue
            tree = load_tree_relabeled(stage_dir / tree_path_str)
```

The rest of the loop body (`clade = local_neighborhood(...)` through `written.append(out_path)`) is unchanged, EXCEPT the caption text at the bottom, which referenced the removed `scope` variable:

```python
        fig.text(0.5, 0.01,
                 "Real, unmodified S/M/L sequences from a single deposited isolate; "
                 "trees pruned to a local neighborhood for legibility (full tree in "
                 f"output/phylogenetics{'_phylogroup' if scope == 'phylogroup' else ''}/).",
                 ha="center", fontsize=7, color="#555555")
```

becomes (reading the scope straight from the row instead of a loop-local variable):

```python
        row_scope = row["analysis_scope"]
        fig.text(0.5, 0.01,
                 "Real, unmodified S/M/L sequences from a single deposited isolate; "
                 "trees pruned to a local neighborhood for legibility (full tree in "
                 f"output/phylogenetics{'_phylogroup' if row_scope == 'phylogroup' else ''}/).",
                 ha="center", fontsize=7, color="#555555")
```

Also update `ah.normalize_attr_name` → `prep.normalize_attr_name` at the `out_path = output_dir / f"{...}.png"` line if Task 4's rename hasn't already covered it (it should have — confirm with `grep -n "ah\." figures.py` returning nothing at this point).

- [ ] **Step 5: Rewrite `render_after_figures` similarly**

Replace the ingest/grouping preamble and per-segment loop body (lines 331-372 as of this plan):

```python
def render_after_figures(ga_search_json: Path, output_dir: Path, input_dir: str = ".",
                         min_tips: int = DEFAULT_MIN_TIPS) -> Tuple[List[Path], Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    payload = json.loads(ga_search_json.read_text(encoding="utf-8"))
    winners = payload["winners"]

    isolates, divergence_by_key, cds_index, paths = load_ingest(input_dir)
    by_key = {i.isolate_key: i for i in isolates}
    by_group = {
        "species": group_for_scope(isolates, divergence_by_key, cds_index, "species"),
        "phylogroup": group_for_scope(isolates, divergence_by_key, cds_index, "phylogroup"),
    }
    group_dir_root = {
        "species": paths.output_dir / "phylogenetics",
        "phylogroup": paths.output_dir / "phylogenetics_phylogroup",
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    table_rows: List[Dict[str, str]] = []

    for idx, winner in enumerate(winners):
        triplet = winner["triplet"]
        candidate_id = f"ga-winner-{idx:02d}"
        fig, axes = plt.subplots(1, 3, figsize=(13, 5))
        for ax, segment, donor_key in zip(axes, SEGMENTS, triplet):
            group_key = members = scope_used = None
            for scope in ("species", "phylogroup"):
                group_key, members = find_group_for_isolate(by_group[scope], donor_key)
                if group_key is not None:
                    scope_used = scope
                    break
            if group_key is None:
                ax.set_axis_off()
                ax.text(0.5, 0.5, f"{segment} donor:\nnot in a committed tree",
                        ha="center", va="center", wrap=True, transform=ax.transAxes)
                continue
            label_map = label_map_for_members(members)
            tree = load_tree_relabeled(group_dir_root[scope_used], group_key, segment, label_map)
```

with:

```python
def render_after_figures(ga_search_json: Path, output_dir: Path, data_dir: str = "output",
                         min_tips: int = DEFAULT_MIN_TIPS) -> Tuple[List[Path], Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    payload = json.loads(ga_search_json.read_text(encoding="utf-8"))
    winners = payload["winners"]

    stage_dir = Path(data_dir)
    isolates = load_ingest(data_dir)
    by_key = {i.isolate_key: i for i in isolates}
    tree_paths_by_isolate = tree_paths_by_isolate_from_table(stage_dir / "ortho_reassortment_table.csv")

    output_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    table_rows: List[Dict[str, str]] = []

    for idx, winner in enumerate(winners):
        triplet = winner["triplet"]
        candidate_id = f"ga-winner-{idx:02d}"
        fig, axes = plt.subplots(1, 3, figsize=(13, 5))
        for ax, segment, donor_key in zip(axes, SEGMENTS, triplet):
            tree_path_str = tree_paths_by_isolate.get(donor_key, {}).get(segment, "")
            if not tree_path_str:
                ax.set_axis_off()
                ax.text(0.5, 0.5, f"{segment} donor:\nnot in a committed tree",
                        ha="center", va="center", wrap=True, transform=ax.transAxes)
                continue
            tree = load_tree_relabeled(stage_dir / tree_path_str)
```

The rest of the loop body is unchanged.

- [ ] **Step 6: Add the `tree_paths_by_isolate_from_table` helper**

`render_after_figures` needs one isolate's tree paths regardless of `reassortment_call` (a GA winner need not be a `supported_reassortment_candidate`) — unlike `supported_reassortment_rows`, which filters to that call. Add, near `supported_reassortment_rows`:

```python
def tree_paths_by_isolate_from_table(reassortment_table: Path) -> Dict[str, Dict[str, str]]:
    """isolate_key -> {segment: tree_path}, for every analyzed isolate regardless of its
    reassortment_call -- a GA winner need not be a supported_reassortment_candidate to
    have a committed tree. Species-scope preferred over phylogroup-scope per isolate,
    same precedence supported_reassortment_rows uses, and only rows with at least one
    non-empty tree path are kept (phylo_status='species_below_min_taxa'/'failed' rows
    have none).
    """
    with reassortment_table.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    by_isolate: Dict[str, Dict[str, str]] = {}
    scope_used: Dict[str, str] = {}
    for row in rows:
        paths = {seg: row.get(f"{seg.lower()}_tree_path") or "" for seg in SEGMENTS}
        if not any(paths.values()):
            continue
        key = row["isolate_key"]
        scope = row["analysis_scope"]
        if key not in by_isolate or (scope_used[key] == "phylogroup" and scope == "species"):
            by_isolate[key] = paths
            scope_used[key] = scope
    return by_isolate
```

- [ ] **Step 7: Update `main`'s CLI wiring**

```python
    before = sub.add_parser("before")
    before.add_argument("--input", default=".")
    before.add_argument("--output", default="output/figures/before")
    before.add_argument("--min-tips", type=int, default=DEFAULT_MIN_TIPS)

    after = sub.add_parser("after")
    after.add_argument("--input", default=".")
    after.add_argument("--ga-search-json", required=True)
    after.add_argument("--output", default="output/figures/after")
    after.add_argument("--min-tips", type=int, default=DEFAULT_MIN_TIPS)

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if args.command == "before":
        written = render_before_figures(Path(args.output), args.input, min_tips=args.min_tips)
        print(f"wrote {len(written)} figure(s) to {args.output}")
    else:
        written, table_path = render_after_figures(
            Path(args.ga_search_json), Path(args.output), args.input, min_tips=args.min_tips)
        print(f"wrote {len(written)} figure(s) and {table_path}")
```

becomes (renaming `--input`'s meaning from "project root containing all_hanta/" to "Stage 1/2 output directory", default changes from `"."` to `"output"`; `--output` is unchanged — it's always been figures' own render-output directory, a different thing):

```python
    before = sub.add_parser("before")
    before.add_argument("--input", default="output",
                        help="Stage 1/2 output directory (reference snapshot + "
                             "ortho_reassortment_table.csv + phylogenetics*/)")
    before.add_argument("--output", default="output/figures/before")
    before.add_argument("--min-tips", type=int, default=DEFAULT_MIN_TIPS)

    after = sub.add_parser("after")
    after.add_argument("--input", default="output",
                        help="Stage 1/2 output directory (reference snapshot + "
                             "ortho_reassortment_table.csv + phylogenetics*/)")
    after.add_argument("--ga-search-json", required=True)
    after.add_argument("--output", default="output/figures/after")
    after.add_argument("--min-tips", type=int, default=DEFAULT_MIN_TIPS)

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if args.command == "before":
        written = render_before_figures(Path(args.output), args.input, min_tips=args.min_tips)
        print(f"wrote {len(written)} figure(s) to {args.output}")
    else:
        written, table_path = render_after_figures(
            Path(args.ga_search_json), Path(args.output), args.input, min_tips=args.min_tips)
        print(f"wrote {len(written)} figure(s) and {table_path}")
```

(The positional argument `args.input` is passed to `render_before_figures`/`render_after_figures`'s `data_dir` parameter — the parameter rename in Steps 4-5 is internal; the CLI flag name `--input` is kept as-is to minimize churn, only its default and help text change.)

- [ ] **Step 8: Fix the module docstring and smoke test**

Module docstring line 18 was already fixed to say "prepare_dataset.py's own..." in Task 4 — leave as-is.

Smoke test (lines 558, 564-565 as of this plan):

```python
    """Smoke test against the real committed all_hanta/ package and output/ -- confirms
    ...
        if not Path("all_hanta/ncbi_dataset/data").exists():
            self.skipTest("all_hanta/ package not present in this environment")
```

→

```python
    """Smoke test against the real committed reference snapshot and phylogenetics output
    under output/ -- confirms
    ...
        if not Path("output/ortho_reference_isolates.json").exists():
            self.skipTest("output/ reference snapshot not present -- run prepare_dataset.py first")
```

- [ ] **Step 9: Update any test that constructed `by_group`/called the deleted functions directly**

```bash
grep -n "group_for_scope\|label_map_for_members\|find_group_for_isolate\|load_tree_relabeled(" figures.py
```

The test section has direct unit tests for `label_map_for_members` and `find_group_for_isolate` (referenced around lines 456, 472, 479-485 in the original file) — delete those test cases along with the functions (Step 1). Any test exercising `load_tree_relabeled` needs updating to call the new one-argument form, constructing a temp directory with a `tree.treefile` and sibling (parent-of-parent) `labels.json` fixture instead of a `label_map` of fake `PhyloCandidate` objects.

- [ ] **Step 10: Confirm no raw access or deleted-function references remain**

```bash
grep -n "all_hanta\|resolve_inputs\|iter_jsonl\|transform_one_record\|build_isolate_records\|group_for_scope\|label_map_for_members\|find_group_for_isolate" figures.py
```

Expected: no output.

- [ ] **Step 11: Run the test suite**

```bash
python -m unittest figures -v
```

Expected: all pass.

- [ ] **Step 12: Commit**

```bash
git add figures.py
git commit -m "$(cat <<'EOF'
figures.py: read labels.json + CSV tree paths, delete group reconstruction

load_ingest shrinks to isolate metadata only (species/country/host) --
grouping, divergence, and cds_index were only ever needed to feed the
now-deleted group_for_scope/label_map_for_members/find_group_for_isolate,
which re-ran phylogeny.py's grouping logic and risked a
--cluster-threshold/--min-species-size mismatch against whatever Stage
2 was actually run with. Tree paths now come straight from
ortho_reassortment_table.csv's s/m/l_tree_path columns, and tip
relabeling reads phylogeny.py's newly-persisted labels.json instead of
reconstructing the label assignment.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 13: `rpi_core.py` — strict-merged set + overlap reporting

**Files:**
- Modify: `rpi_core.py` (the fitness-function reference-loading call site — the same location `run()`/whichever function currently calls `load_reference_data` to build the pool the GA scores against; confirm exact function name via grep before editing)

**Interfaces:**
- Consumes: `terminally_complete_isolates(reference)` (existing, unchanged), `reference.isolates[*].has_complete_triplet` (existing field).
- Produces: a `merged` isolate list and a printed one-line overlap summary at that call site.

- [ ] **Step 1: Find the call site that currently builds the fitness function's isolate pool**

```bash
grep -n "terminally_complete_isolates(" rpi_core.py
```

Identify the function that calls this AND is the entry point for the GA fitness function's pool-building (likely `usable_pool` or similar, given `Algorithm.py:416`'s reference to `rc.usable_pool(reference, reserved)`).

- [ ] **Step 2: Add the strict-merged computation and overlap line**

At that call site, add (adapting variable names to match what's actually there):

```python
    coding_complete = {i.isolate_key for i in reference.isolates if i.has_complete_triplet}
    terminal_complete_isolates = terminally_complete_isolates(reference)
    terminal_complete = {i.isolate_key for i in terminal_complete_isolates}
    merged = [i for i in terminal_complete_isolates if i.has_complete_triplet]

    both = coding_complete & terminal_complete
    only_coding = coding_complete - terminal_complete
    only_terminal = terminal_complete - coding_complete
    print(f"completeness overlap: in both {len(both)}, only coding-complete (206-set) "
          f"{len(only_coding)}, only terminally-complete (272-set) {len(only_terminal)}, "
          f"strict merged set {len(merged)}")
```

Place this where the pool is actually consumed by the fitness function / GA search entry point (do not duplicate it at every `terminally_complete_isolates` call site — only the one that's the real Stage 3 entry point for scoring).

- [ ] **Step 3: Use `merged` as the fitness function's input pool**

Wire `merged` (not `terminal_complete_isolates`, not `reference.isolates`) into whatever variable currently feeds the GA/fitness function at this call site.

- [ ] **Step 4: Write a test for the overlap arithmetic**

```python
class CompletenessOverlapTest(unittest.TestCase):
    def test_overlap_counts_and_merged_set(self):
        coding_only = _make_isolate("c1", "Orthohantavirus x", "S1", "M1", "L1")
        coding_only.has_complete_triplet = True
        terminal_only = _make_isolate("t1", "Orthohantavirus x", "S2", "M2", "L2")
        terminal_only.has_complete_triplet = False
        both = _make_isolate("b1", "Orthohantavirus x", "S3", "M3", "L3")
        both.has_complete_triplet = True

        reference = prep.ReferenceData(
            isolates=[coding_only, terminal_only, both],
            genomic_index={}, cds_index={})

        # terminally_complete_isolates needs genomic sequences with the terminal tags;
        # rather than fabricating those here, monkeypatch it for this unit test to
        # isolate the overlap arithmetic from the terminal-completeness filter itself
        # (which already has its own tests -- TerminalCompletenessTest).
        with unittest.mock.patch(
                "rpi_core.terminally_complete_isolates",
                return_value=[terminal_only, both]):
            coding_complete = {i.isolate_key for i in reference.isolates if i.has_complete_triplet}
            terminal_complete_isolates = terminally_complete_isolates(reference)
            terminal_complete = {i.isolate_key for i in terminal_complete_isolates}
            merged = [i for i in terminal_complete_isolates if i.has_complete_triplet]

        self.assertEqual(coding_complete, {"c1", "b1"})
        self.assertEqual(terminal_complete, {"t1", "b1"})
        self.assertEqual([i.isolate_key for i in merged], ["b1"])
```

- [ ] **Step 5: Run it**

```bash
python -m unittest rpi_core.CompletenessOverlapTest -v
```

Expected: pass.

- [ ] **Step 6: Run the full suite**

```bash
python -m unittest rpi_core -v
```

Expected: all pass, including whatever existing test covers the pool-building call site you modified in Steps 2-3 (its assertions may need updating if it hardcoded the previous pool size — check for a test asserting an exact isolate count from this call site and update the expected number if the strict-merged filter now changes what's returned there; note "changes what's returned" here means "makes explicit a filter that either already existed implicitly or is newly introduced by this task" — do not silently change scoring behavior beyond what Task 13 is scoped to do. If the existing pool-building function already applied exactly this merged filter under a different name, Task 13 is a refactor (rename/relocate the logic, same result); if it did not, flag this to the user before proceeding, since changing what isolates feed the fitness function is a scoring-behavior change outside a "wire the plumbing" refactor).

- [ ] **Step 7: Commit**

```bash
git add rpi_core.py
git commit -m "$(cat <<'EOF'
rpi_core.py: strict-merged completeness set + one-line overlap report

merged = terminally-complete AND coding-complete, computed once at the
fitness function's pool-building call site. Prints the coding-complete
(206-set) / terminally-complete (272-set) / merged / overlap counts so
the two rulers' relationship is visible on every run, not just implied.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 14: README — update CLI examples and document the three completeness rulers

**Files:**
- Modify: `README.md`

**Interfaces:**
- None (documentation only).

- [ ] **Step 1: Update the Stage 1 rename**

Find any `README.md` reference to `all_hanta.py`/`01-prepare_dataset.py` as a command (e.g. `python all_hanta.py` or `python 01-prepare_dataset.py`) and update to `python prepare_dataset.py`:

```bash
grep -n "all_hanta\.py\|01-prepare_dataset" README.md
```

- [ ] **Step 2: Update CLI examples with the new defaults**

```bash
grep -n "reference-input\|--input " README.md
```

For each example using `--reference-input .`  or similar against the old `all_hanta/`-relative semantics, update to match Task 9-11's new defaults (`output` for most, per that task's exact default chosen). Cross-check against the actual `parse_args`/`add_argument` defaults in each file after Tasks 9-11 land, rather than assuming — the exact default strings are the source of truth.

- [ ] **Step 3: Document the three completeness definitions**

Add a short section (near wherever the pipeline stages are already documented, per the "Pipeline stages" section referenced in `phylogeny.py`'s docstring) defining:

- **Coding-complete (206 set)**: isolates where `has_complete_triplet` is true (all three segments coding-complete) — used by Stage 2 (phylogeny.py) for tree building.
- **Terminally-complete (272 set)**: isolates where all three segments are both full-length AND carry genuine 5'/3' termini (`terminally_complete_isolates` in rpi_core.py) — the broader search/calibration pool for Stage 3.
- **Strict merged set**: the intersection (coding-complete AND terminally-complete) — the actual fitness-function input pool, expected to be smaller than 206, since packaging/reassortment signal lives in the termini that a coding-complete-but-trimmed deposit lacks.

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "$(cat <<'EOF'
README: update CLI examples for the Stage 1 snapshot, document the
three completeness definitions (coding-complete / terminally-complete
/ strict merged)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013CpAeTwfkmTkgJUnk2Uwm1
EOF
)"
```

---

### Task 15: Final acceptance

**Files:** none modified — verification only.

- [ ] **Step 1: Run both acceptance greps from the spec**

```bash
grep -rn "\bah\.\|import all_hanta" --include="*.py" .
```
Expected: empty.

```bash
grep -rn '"all_hanta"' --include="*.py" .
```
Expected: empty except inside `prepare_dataset.py` and any test fixture explicitly documented in Tasks 9-12 as still legitimately pointing at the raw package (there should be none by this point — every test fixture was converted to point at `output/` — if this grep finds anything outside `prepare_dataset.py`, treat it as a bug and fix it, don't just note it).

- [ ] **Step 2: Run every module's full test suite**

```bash
python -m unittest prepare_dataset -v
python -m unittest phylogeny -v
python -m unittest rpi_core -v
python -m unittest Algorithm -v
python -m unittest figures -v
```

Expected: all green.

- [ ] **Step 3: Run Stage 1 end-to-end and confirm the 3 new files appear**

```bash
python prepare_dataset.py --input . --output /tmp/ortho_snapshot_check
ls /tmp/ortho_snapshot_check/ortho_reference_isolates.json /tmp/ortho_snapshot_check/ortho_reference_genomic.json /tmp/ortho_snapshot_check/ortho_reference_cds.json
sha256sum -c /tmp/ortho_snapshot_check/ortho_output_checksums.txt
```

Expected: all 3 files exist, checksum verification passes.

- [ ] **Step 4: Run Stage 2 and Stage 3 against that snapshot with no filesystem access to `all_hanta/`**

```bash
python phylogeny.py --output /tmp/ortho_snapshot_check
python rpi_core.py validate --reference-input /tmp/ortho_snapshot_check --json /tmp/ortho_snapshot_check/validation_report.json
```

Expected: both run without any `FileNotFoundError` pointing at `all_hanta/ncbi_dataset`, confirming the raw package is genuinely no longer touched by these stages (you can additionally strace/monitor filesystem access if available, but a clean run against a snapshot directory that has no `all_hanta/` sibling is sufficient confirmation).

- [ ] **Step 5: Report completion**

No commit for this task (verification only) — report the final acceptance grep outputs and test results to confirm the refactor is complete.

---

## Self-Review Notes

- **Spec coverage:** Rules 1 (rename) → Task 1. Rule 2 (alias rename + safety steps) → Tasks 2-5. Rule 3 (raw-data ban + per-file disposition, including the figures.py correction) → Tasks 6, 9, 10, 11, 12. Rule 4 (find-then-classify) → the grep-driven steps throughout Tasks 2-4 and 9-12. Rule 5 (leave README downloads alone) → Task 14 does not touch download/unzip commands. Rule 6 (approved snapshot design) → Tasks 7-8 (writer + reader), Task 13 (merged set + overlap, kept out of Stage 1). Rule 7 (acceptance) → Task 15, using the amended `"all_hanta"` string-literal grep. The spec's Task-6-adjacent amendment (parse_cds_fasta relocation) is its own task (Task 6). The keying-equality assertion is Task 8, Step 1. The FINAL spec amendment (labels.json persistence, figures.py's group-reconstruction removal) is Task 9 Steps 4b-4c (writer side) and the rewritten Task 12 (reader side).
- **Type consistency:** `ReferenceData`/`load_reference_data` are defined once, in `prepare_dataset.py` (Task 8), and every other task (9, 10, 11, 12, 13) consumes that exact name via `prep.` — no task redefines it locally.
- **Ambiguity flagged inline:** Task 11 Step 1 and Task 12 Step 1 both instruct the executor to check actual behavior (`python -c "import ..."`, grep for existing usage) before committing to one of two possible edits, since this plan cannot know rpi_core.py's exact re-export shape or figures.py's exact call-site unpacking without running Task 10/earlier tasks first. Task 13 Step 6 explicitly requires flagging to the user if the merged-set filter changes actual scoring behavior beyond a plumbing refactor, rather than silently proceeding — this is a deliberate stop-and-ask condition, not a gap.
