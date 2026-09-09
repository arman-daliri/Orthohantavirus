# Stage 1 reference snapshot: decoupling downstream pipelines from all_hanta/

Date: 2026-09-09

## Problem

`all_hanta.py` (Stage 1 ingest) was renamed to `01-prepare_dataset.py`. That
name is not importable (`import 01-prepare_dataset` is a syntax error — leading
digit, hyphen), so every downstream file's `import all_hanta as ah` is
currently broken:

- `phylogeny.py` (Stage 2) does `import all_hanta as ah` and independently
  re-derives isolates, the genomic/protein FASTA index, and the CDS index
  straight from `--input` (the raw `all_hanta/ncbi_dataset/` package), by
  design, to stay "independent of all_hanta.py's output files."
- `rpi_core.py` (Stage 3)'s `load_reference_data(root)` does the same
  re-derivation via `ah.resolve_inputs`/`ah.iter_jsonl`/`ah.parse_fasta`/
  `ah.transform_one_record`/`ah.build_isolate_records`, plus
  `ph.parse_cds_fasta` (borrowed from phylogeny.py) for the CDS index.
- `Algorithm.py` calls `rpi_core.load_reference_data(Path("all_hanta"))`.
- `figures.py` does `import all_hanta as ah` and its smoke test checks for
  `all_hanta/ncbi_dataset/data`.

Three files independently re-parsing the same ~15 MB raw NCBI package is both
the proximate cause of the broken import and the underlying defect: Stage 1
is supposed to be the single accurate ingest, and Stages 2/3 are supposed to
consume its output, not repeat its work.

## Goal

`prepare_dataset.py` (renamed from `01-prepare_dataset.py`) becomes the
**only** file that reads `all_hanta/ncbi_dataset/data`. It ingests once and
writes a reference snapshot; every downstream pipeline eats that snapshot.
Pure code (dataclasses, constants, stateless functions) is still shared by
import — only *raw file access* is restricted to Stage 1.

## Non-goals

- No change to the raw data folder: `all_hanta/` keeps its name and location.
  It is NCBI's data package, not code.
- No change to the two-ruler design (206 coding-complete for Stage 2, 272
  terminally-complete pool / strict-merged set for Stage 3) — this spec only
  changes *how the underlying isolate/sequence data reaches each stage*, not
  which isolates each stage selects.
- No change to `rpi_core.py`'s legitimate reuse of `phylogeny.py`'s pure code
  (`PhyloCandidate`, `select_primary_cds`, `resolve_binary`, `_analyze_group`).

## Design

### 1. Rename: `01-prepare_dataset.py` → `prepare_dataset.py`

Fixes the broken import at the root. No behavior change.

### 2. Import + alias rename: `ah` → `prep`

Every file currently doing `import all_hanta as ah` switches to:

```python
import prepare_dataset as prep
```

and every `ah.<symbol>` call site in that file is renamed to `prep.<symbol>`.
This is a mechanical rename across `phylogeny.py`, `rpi_core.py`, and
`figures.py` (and any other file with the import) — not a partial fix.

**Alias choice:** `prep`, not `pd` (collides with `import pandas as pd`,
which `prepare_dataset.py` itself uses) and not `ah` (still tied to the old
module name).

**Safety steps for the rename** (mechanical, but real files are large —
`rpi_core.py` alone has ~150 `ah.*` call sites):

1. Before touching anything: `grep -rn "\bprep\b" --include="*.py" .` — confirm
   `prep` is not already bound to something else in these files. If it is,
   stop and report the collision instead of shadowing it.
2. Rename with word-boundary substitution only (e.g. `sed 's/\bah\./prep./g'`
   per file) — never a blind replace of the letters "ah", which appear inside
   unrelated identifiers and prose.
3. Update comments/docstrings that name the *module* `all_hanta.py` (e.g.
   "Reuses all_hanta.py's ingest primitives") to say `prepare_dataset.py`.
   Do **not** touch text that means the *data folder* `all_hanta/` — that
   name doesn't change (see Non-goals).
4. After editing, `grep -rn "\bah\.\|import all_hanta" --include="*.py" .`
   must return nothing, then run the full test suite.

### 3. Raw-data access rule

- `all_hanta/` (the folder) keeps its name; nothing here renames it.
- `prepare_dataset.py` is the only file allowed to read
  `all_hanta/ncbi_dataset/data`. That is its job (ingesting raw data).
- Downstream files never open that folder again. Per-file disposition:

  | File | Change |
  |---|---|
  | `phylogeny.py` | Delete its raw re-derivation (`resolve_inputs`, `iter_jsonl`, `parse_fasta`, `transform_one_record`, `build_isolate_records`, `assign_kmer_clusters`, `compute_segment_divergence`). Load the Stage 1 snapshot instead. |
  | `rpi_core.py` | Rewrite `load_reference_data(Path("all_hanta"))` → `load_reference_data(Path("output"))`, reading the 3 JSON snapshot files. Drop `ph.parse_cds_fasta`. |
  | `Algorithm.py` | `--reference-input` default `"all_hanta"` → `"output"`; the one `Path("all_hanta")` call site becomes `Path("output")`. |
  | `figures.py` | Keep its `prep` import (pure code only — no raw access to begin with). Move its smoke test's existence-check from `all_hanta/ncbi_dataset/data` to `output/ortho_reference_isolates.json`. |

- `README.md`'s download commands (`datasets download ... --filename
  all_hanta.zip`, `unzip -d all_hanta`) are untouched — that's just where the
  raw data lives, not a code path this spec changes.

### 4. Snapshot artifacts (Stage 1 output, format: JSON)

Three new files, written with the existing `safe_write_json` (atomic write,
sha256 returned for the checksum manifest) — no new I/O helper needed:

- **`output/ortho_reference_isolates.json`** — JSON array, one object per
  `IsolateRecord`, every field, for **every** isolate (not filtered to 206 or
  272 — filtering is a downstream concern, not Stage 1's).
- **`output/ortho_reference_genomic.json`** — JSON array, one object per
  genomic accession: `{accession, accession_base, header, sequence, length}`
  (the `FastaRecord` fields), for every accession in the genomic FASTA index,
  not just complete-triplet isolates.
- **`output/ortho_reference_cds.json`** — JSON object keyed by base
  accession → list of `[product, sequence]` pairs, the same shape as today's
  in-memory `cds_index` (built today by `phylogeny.py`'s `parse_cds_fasta`,
  moved into Stage 1's own ingest run).

All three are added to `ortho_output_checksums.txt` like every other Stage 1
output.

**Why JSON, not pickle or parquet:** pickle is unreadable to a reviewer,
couples the artifact to the exact dataclass shape that wrote it (a field
rename breaks loading silently, no clear error), and is a code-execution
vector on load — wrong for a repo whose outputs are judged by whether a
reviewer would accept them. Parquet adds a new dependency (pyarrow) for no
benefit at this row count (thousands, not millions) and nothing else in the
repo uses it. JSON reuses `safe_write_json`'s existing atomic-write + sha256
path unchanged, so provenance is uniform across every Stage 1 output, not a
special case for this one.

### 5. Stage 1 stays filter-free

`prepare_dataset.py` writes the full, unfiltered isolate/genomic/CDS data. It
does not know about "206 coding-complete" or "272 terminally-complete" — those
are downstream concepts:

- Stage 2 (`phylogeny.py`) keeps computing and printing its own 206-count
  from `isolate.has_complete_triplet`, now read off the loaded snapshot
  instead of a live re-derivation.
- Stage 3 (`rpi_core.py`) keeps `terminally_complete_isolates(reference)`
  unchanged (it already operates purely on `reference.isolates` +
  `reference.genomic_index`, which the snapshot now supplies). The strict
  merged set for the fitness function becomes:

  ```python
  merged = [i for i in terminally_complete_isolates(reference)
            if i.has_complete_triplet]
  ```

  computed once, in `rpi_core.py`, at the same call site that currently loads
  the reference data.

### 6. Overlap reporting

The same `rpi_core.py` call site prints the one-line overlap using set
arithmetic on `isolate_key`:

```python
coding_complete = {i.isolate_key for i in reference.isolates if i.has_complete_triplet}
terminal_complete = {i.isolate_key for i in terminally_complete_isolates(reference)}
# "in both: N, only 206: N, only 272: N"
```

No new data is needed for this — it's a Stage 3-only change. The three
completeness definitions (coding-complete / terminally-complete / strict
merged) get documented in `README.md`.

### 7. Provenance

`ortho_summary.json` (already written by Stage 1) gets a new
`reference_snapshot` block: isolate count, genomic-accession count,
cds-indexed-accession count, and the three new files' sha256s — chained to
the same `verify_inputs` md5sum digests every other Stage 1 output already
carries. The 206/272/merged counts themselves are NOT recorded here (Stage 1
doesn't compute them); they're printed and logged by Stage 2 and Stage 3
respectively.

### 8. Interface after this change

```python
# Stage 1 (prepare_dataset.py) -- unchanged ingest, +3 new writes:
write_json("ortho_reference_isolates.json", [asdict(i) for i in isolates])
write_json("ortho_reference_genomic.json", [asdict(r) for r in genomic_index.values()])
write_json("ortho_reference_cds.json", cds_index)  # already {accession: [(product, seq), ...]}

# Stage 2 (phylogeny.py) / Stage 3 (rpi_core.py):
def load_reference_data(output_dir: Path) -> ReferenceData:
    isolates = [prep.IsolateRecord(**row) for row in load_json(output_dir / "ortho_reference_isolates.json")]
    genomic_index = {r["accession"]: prep.FastaRecord(**r)
                      for r in load_json(output_dir / "ortho_reference_genomic.json")}
    cds_index = {k: [tuple(pair) for pair in v]
                 for k, v in load_json(output_dir / "ortho_reference_cds.json").items()}
    return ReferenceData(isolates=isolates, genomic_index=genomic_index, cds_index=cds_index)
```

No file after Stage 1 calls `resolve_inputs`, `iter_jsonl`, `parse_fasta`,
`transform_one_record`, or `build_isolate_records` against a raw path again.

## Existing call sites needing disposition (inventory, not exhaustive editing here)

`grep -rn "all_hanta" --include="*.py" .` currently returns ~70 hits outside
`prepare_dataset.py`, spanning `phylogeny.py`, `rpi_core.py`, `Algorithm.py`,
`figures.py`. Each falls into one of three buckets — the implementation plan
must classify and handle every one, not just the illustrative examples above:

1. **Module import / `ah.*` call site** → rule 2 (rename to `prep`).
2. **Raw data path in a downstream file** (e.g. `Path("all_hanta")` default
   args, `load_reference_data(Path("all_hanta"))` call sites, test fixtures
   that build `Path(__file__).parent / "all_hanta"`) → rule 3 (replace with
   `Path("output")` / snapshot-file fixtures).
3. **Prose referencing `all_hanta.py` the module, or `all_hanta/` the data
   folder, in a comment/docstring** → reword the module references only
   (rule 2c); leave data-folder references alone.

Test fixtures in `rpi_core.py` (`test_all_hanta.py`-pattern classes around
lines 4209-4302, 5431-5436, 6315-6318, 8611-8626, 9013-9144, 9704) that
currently skip unless `all_hanta/` exists on disk need converting to build or
load a small snapshot fixture instead, so they no longer depend on the raw
package being present in the checkout.

## Acceptance criteria

```
grep -rn "\bah\.\|import all_hanta" --include="*.py" .   # must be empty
grep -rn "Path(\"all_hanta\")" --include="*.py" .        # must be empty except
                                                          # inside prepare_dataset.py
```

Full test suite green. `python prepare_dataset.py` (Stage 1) still produces
every existing output plus the 3 new snapshot files; `phylogeny.py` and
`rpi_core.py`/`Algorithm.py` run against `output/` only, with no filesystem
access to `all_hanta/ncbi_dataset/` from those files.

## Open questions for the implementation plan

- Exact ordering: rename+alias-fix (rules 1-2) can land as a mechanical,
  behavior-preserving commit before the snapshot/raw-access work (rules 3-6),
  which changes actual data flow. Recommend splitting into two plan phases
  for reviewability, each with its own test run.
- `IsolateRecord`/`FastaRecord` field renames in the future would silently
  break `prep.IsolateRecord(**row)` reconstruction from old snapshots — not
  a concern for this change, but worth a one-line note in
  `prepare_dataset.py` near the dataclass definitions once this lands.
