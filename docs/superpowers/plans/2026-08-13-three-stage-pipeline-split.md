# Three-Stage Pipeline Split Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split `all_hanta.py`'s inline `--phylogenetics` MAFFT+IQ-TREE analysis into a standalone `phylogeny.py` (Stage 2), fix the resulting import surface in the three Stage-3 files, prune two redundant output files, and document all three stages in the README.

**Architecture:** Mechanical code-motion refactor, not new functionality. A verified, exhaustive line-range inventory (below) is extracted from `all_hanta.py` into `phylogeny.py` via a Python script (not manual retyping, to eliminate transcription risk), given its own CLI, then the source ranges are deleted from `all_hanta.py`. Tests move in lockstep. The suite must stay green after every task — this is a refactor's only meaningful "test."

**Tech Stack:** Python 3, stdlib `unittest`, `numpy`, `pandas`, Biopython + `dendropy` (Stage 2 only, lazily imported), `mafft`/`iqtree` binaries (Stage 2/3 only).

**Spec:** `docs/superpowers/specs/2026-08-13-three-stage-pipeline-split-design.md`

## Global Constraints

- No behavior change to any output file's column schema — only which files exist (drop `ortho_segment_divergence.csv`, `ortho_ml_segments.csv`).
- `ortho_risk_table.csv` becomes unconditionally phylo-blind (matches today's no-`--phylogenetics` default exactly).
- Every task ends with the full test suite green. A refactor task that leaves the suite red is not done.
- No new shared utility module — `all_hanta.py` keeps ingest primitives, `phylogeny.py` keeps tree-analysis primitives (spec §8).
- `output/` is git-tracked by existing convention; regenerate and re-commit it, don't stop tracking it.

---

## Verified Extraction Inventory (read this before Task 1)

Four exact line ranges in `all_hanta.py` (1-indexed, inclusive) move to `phylogeny.py`, in this order:

| # | Range | Contents |
|---|---|---|
| 1 | 451–456 | `HOST_ORDER_PHYLOGROUP_FALLBACK` dict |
| 2 | 732–745 | `@dataclass class PhyloCandidate` |
| 3 | 2028–2066 | `_segment_rank_metrics` |
| 4 | 2295–3115 | Section banner comment through `_phylo_call_by_isolate` (24 defs: `HUMAN_AFFECTING_CLASSES`, `MIN_CDS_LENGTH`, `is_human_affecting_or_likely`, `parse_cds_fasta`, `select_primary_cds`, `collect_phylo_candidates`, `_phylogroup_for_isolate`, `_phylogroup_for_candidate`, `_assignment_basis`, `_candidate_group_key`, `group_candidates`, `group_candidates_by_species`, `resolve_binary`, `write_segment_fasta`, `run_mafft`, `run_iqtree`, `_patristic_distance_matrix`, `tree_divergence_for_species`, `_nearest_supported_descendant`, `bootstrap_support_for_pair`, `classify_reassortment_call`, `_phylo_row`, `_resolve_outgroup_label`, `_analyze_group`, `_run_phylogenetics_scope`, `run_phylogenetics_stage`, `_REASSORTMENT_CALL_RANK`, `_phylo_call_by_isolate`) |

Range 4's true start (2295, not 2329 as the design spec's first draft said) was found by checking every symbol used-but-not-called in the block, not just function calls — `HUMAN_AFFECTING_CLASSES` (line 2325) sits inside a section banner comment (2295–2323) that starts 34 lines before `MIN_CDS_LENGTH`. Missing this would have split a defined-and-used-together pair across two files.

**Ten bare names inside these ranges need an `ah.` prefix once moved** (they stay in `all_hanta.py`; verified by exhaustive grep of every call/subscript in the extracted text, not just function-call syntax, which is what caught the four listed after `IsolateRecord` on a second pass):

`_rank_positions` (1 site), `accession_base` (4 sites), `species_risk_profile` (1 site), `classify_protein_role` (1 site), `_iter_fasta_records` (1 site), `IsolateRecord` (3 type-hint sites), `PHYLOGROUP` (6 sites, one inside a comment), `SPECIES_RISK` (1 site), `SegmentDivergenceProfile` (2 sites), `SEGMENT_PRIMARY_PRODUCT` (1 site).

**Two names defined inside the moved ranges need NO prefix** (they move with their own definition): `HUMAN_AFFECTING_CLASSES`, `MIN_CDS_LENGTH`.

**One dependency runs the other direction and must NOT be prefixed away**: `_rank_positions` is defined at line 1997 (outside every moved range, stays in `all_hanta.py`) because it's also used by Stage 1's own k-mer path at line 2153. This is exactly why it's in the ah.-prefix list above, not a symbol that moves.

**Stage 3's import-surface fix — six symbols, verified by grepping every `ah.` call site in all three files, not skimming:**

| File | Line | Current | Becomes |
|---|---|---|---|
| `reassortment_detector.py` | 151–152 | `ah.resolve_binary(...)` ×2 | `ph.resolve_binary(...)` |
| `detector_gates.py` | 513, 518, 532, 572 | `ah.PhyloCandidate` | `ph.PhyloCandidate` |
| `detector_gates.py` | 528 | `ah.select_primary_cds(...)` | `ph.select_primary_cds(...)` |
| `detector_gates.py` | 529 | `ah.MIN_CDS_LENGTH` | `ph.MIN_CDS_LENGTH` |
| `detector_gates.py` | 581 | `ah._analyze_group(...)` | `ph._analyze_group(...)` |
| `detector_reference.py` | 40 | `ah.parse_cds_fasta(...)` | `ph.parse_cds_fasta(...)` |

---

### Task 1: Extract `phylogeny.py`'s core code from `all_hanta.py`

Purely additive — `all_hanta.py` is not modified in this task, so nothing can regress yet.

**Files:**
- Create: `phylogeny.py`
- Create (temporary, delete after use): `/tmp/extract_phylogeny.py`

**Interfaces:**
- Produces: `phylogeny.py` importable as a module exposing all 26 moved names (see inventory table) plus `HOST_ORDER_PHYLOGROUP_FALLBACK`, `PhyloCandidate`, `_segment_rank_metrics`. Depends on `import all_hanta as ah` for the ten cross-referenced names above.

- [ ] **Step 1: Write the extraction script**

```python
# /tmp/extract_phylogeny.py
import re
from pathlib import Path

src = Path("all_hanta.py").read_text().splitlines(keepends=True)


def extract(start, end):
    return "".join(src[start - 1:end])


host_fallback = extract(451, 456)
phylo_candidate = extract(732, 745)
segment_rank_metrics = extract(2028, 2066)
main_block = extract(2295, 3115)

body = "\n\n".join([host_fallback.rstrip("\n"), phylo_candidate.rstrip("\n"),
                    segment_rank_metrics.rstrip("\n"), main_block.rstrip("\n")]) + "\n"

# Ten names that stay in all_hanta.py but are referenced bare inside the moved text.
for name in ("_rank_positions", "accession_base", "species_risk_profile",
             "classify_protein_role", "_iter_fasta_records", "IsolateRecord",
             "PHYLOGROUP", "SPECIES_RISK", "SegmentDivergenceProfile",
             "SEGMENT_PRIMARY_PRODUCT"):
    body = re.sub(rf"(?<!ah\.)\b{name}\b", f"ah.{name}", body)

# The moved section-banner comment's setup instructions predate environment.yml
# (2026-08-13) and recommend the exact unpinned `-c bioconda -c conda-forge` ordering
# that pins Biopython to bioconda's ancient 1.70 build -- see environment.yml's commit
# message. Replace with the current guidance rather than carrying stale advice forward.
body = body.replace(
    "# Setup:\n"
    "#   conda create -n hanta_phylo -c bioconda -c conda-forge mafft iqtree biopython dendropy\n"
    "#   conda activate hanta_phylo   # or pass --mafft-bin / --iqtree-bin with full paths\n",
    "# Setup:\n"
    "#   conda env create -f environment.yml\n"
    "#   conda activate hanta_phylo   # or pass --mafft-bin / --iqtree-bin with full paths\n",
)
assert "conda create -n hanta_phylo -c bioconda" not in body, "stale setup comment not replaced"

header = '''#!/usr/bin/env python3
"""Orthohantavirus reassortment phylogenetics stage (MAFFT + IQ-TREE).

Stage 2 of the three-stage pipeline documented in README.md's "Pipeline stages" section.
Analyzes the reference isolate set ingested from an NCBI Datasets Orthohantavirus package
for real, tree-based reassortment evidence -- the rigorous replacement for the k-mer
divergence diagnostics all_hanta.py (Stage 1) still emits as informational-only columns.

Independent of all_hanta.py's output files: re-derives isolates, the genomic/protein FASTA
index, and the CDS index from --input itself, the same way all_hanta.py's own run_pipeline
did before this code was split out (docs/superpowers/specs/2026-08-13-three-stage-pipeline-
split-design.md has the split rationale and the exact symbol inventory this file was
extracted from). Reuses all_hanta.py's ingest primitives via `import all_hanta as ah`.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, DefaultDict, Dict, List, Optional, Sequence, Tuple

import numpy as np

import all_hanta as ah


'''

Path("phylogeny.py").write_text(header + body)
print("Wrote phylogeny.py:", len((header + body).splitlines()), "lines")
```

- [ ] **Step 2: Run the extraction script**

```bash
python3 /tmp/extract_phylogeny.py
```

Expected: prints `Wrote phylogeny.py: <N> lines` (roughly 900), no assertion error.

- [ ] **Step 3: Verify the module imports cleanly**

```bash
python3 -c "import phylogeny as ph; print('ok')"
```

Expected: `ok`, no `NameError`/`ImportError`/`SyntaxError`. If this fails with a `NameError` on
some identifier, that identifier was missed by the ten-name fixup list above — add it, re-run
Step 1's script (it's idempotent, safe to re-run against the untouched `all_hanta.py`), and
retry this step.

- [ ] **Step 4: Spot-check the moved symbols are present and correctly typed**

```bash
python3 -c "
import phylogeny as ph
import inspect
names = ['HOST_ORDER_PHYLOGROUP_FALLBACK', 'PhyloCandidate', '_segment_rank_metrics',
         'MIN_CDS_LENGTH', 'HUMAN_AFFECTING_CLASSES', 'parse_cds_fasta', 'select_primary_cds',
         'resolve_binary', 'run_mafft', 'run_iqtree', 'classify_reassortment_call',
         '_analyze_group', 'run_phylogenetics_stage', '_phylo_call_by_isolate']
missing = [n for n in names if not hasattr(ph, n)]
assert not missing, f'missing from phylogeny.py: {missing}'
assert callable(ph.resolve_binary) and callable(ph._analyze_group)
print('all', len(names), 'symbols present')
"
```

Expected: `all 14 symbols present`.

- [ ] **Step 5: Delete the temporary extraction script and commit**

```bash
rm /tmp/extract_phylogeny.py
git add phylogeny.py
git commit -m "$(cat <<'EOF'
Extract phylogenetics code from all_hanta.py into phylogeny.py

Additive only -- all_hanta.py is untouched in this commit, so nothing
regresses. The four source line ranges and the ten-symbol ah.-prefix
fixup list are the verified inventory from
docs/superpowers/specs/2026-08-13-three-stage-pipeline-split-design.md.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Give `phylogeny.py` its own CLI and top-level run function

**Files:**
- Modify: `phylogeny.py` (append to end of file)

**Interfaces:**
- Consumes: everything from Task 1 (`ah.resolve_inputs`, `ah.verify_inputs`, `ah.report_checksums`, `ah.configure_logging`, `ah.iter_jsonl`, `ah.build_annotation_index`, `ah.build_biosample_index`, `ah.parse_fasta`, `ah.transform_one_record`, `ah.build_isolate_records`, `ah.assign_kmer_clusters`, `ah.compute_segment_divergence`, `ah.safe_write_csv` — all unchanged Stage-1 functions; `parse_cds_fasta`, `collect_phylo_candidates`, `run_phylogenetics_stage` — Task 1's moved functions, called bare since they're now in this module).
- Produces: `phylogeny.parse_args(argv)`, `phylogeny.main(argv)`, `phylogeny.run_phylogeny(paths, args)` — the entry points Task 6's manual smoke test and any future CLI use call.

- [ ] **Step 1: Append the run function and CLI to `phylogeny.py`**

```python


def run_phylogeny(paths: "ah.InputPaths", args: argparse.Namespace) -> None:
    """Mirrors all_hanta.py's run_pipeline() ingest chain exactly through
    compute_segment_divergence (needed for divergence_by_key, which
    collect_phylo_candidates uses to tag divergence_candidate isolates) -- see spec
    section 3.2 for why this ~duplication is a deliberate choice, not an oversight.
    """
    logging.info("Verifying package integrity against NCBI md5sum.txt")
    checksums = ah.verify_inputs(paths)
    ah.report_checksums(checksums, allow_mismatch=args.allow_checksum_mismatch)

    logging.info("Loading reports")
    data_records = list(ah.iter_jsonl(paths.data_report))
    annotation_index = ah.build_annotation_index(ah.iter_jsonl(paths.annotation_report))
    biosample_index = ah.build_biosample_index(ah.iter_jsonl(paths.biosample_report))

    logging.info("Parsing FASTA")
    genomic_index = ah.parse_fasta(paths.genomic_fasta)
    protein_index = ah.parse_fasta(paths.protein_fasta)

    logging.info("Transforming accession records")
    accessions = []
    for record in data_records:
        transformed, _qc_row = ah.transform_one_record(
            record, annotation_index, biosample_index, genomic_index, protein_index)
        if transformed is not None:
            accessions.append(transformed)

    logging.info("Reconstructing isolates")
    isolates = ah.build_isolate_records(accessions)
    logging.info("  %s isolates, %s complete triplets",
                 len(isolates), sum(1 for i in isolates if i.has_complete_triplet))

    ah.assign_kmer_clusters(isolates, args.cluster_threshold)
    divergence = ah.compute_segment_divergence(isolates, min_species_size=args.min_species_size)
    divergence_by_key = {d.isolate_key: d for d in divergence}

    logging.info("Loading CDS")
    cds_index = parse_cds_fasta(paths.cds_fasta)
    logging.info("  cds: %s accessions with at least one CDS", len(cds_index))

    logging.info("Running phylogenetics stage")
    phylo_rows = run_phylogenetics_stage(isolates, divergence_by_key, cds_index, paths.output_dir, args)
    ah.safe_write_csv(paths.output_dir / "ortho_reassortment_table.csv", phylo_rows)
    logging.info("Done. %s candidate rows written to ortho_reassortment_table.csv", len(phylo_rows))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Orthohantavirus phylogenetics stage: MAFFT + IQ-TREE reassortment analysis "
                    "over the reference isolate set (Stage 2 of the three-stage pipeline)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", default=".",
                        help="Package root containing all_hanta/ncbi_dataset/data")
    parser.add_argument("--output", help="Output directory; defaults to <input>/output")
    parser.add_argument("--allow-checksum-mismatch", action="store_true",
                        help="Continue even if inputs do not match NCBI's md5sum.txt")
    parser.add_argument("--cluster-threshold", type=float, default=0.02,
                        help="Cosine k-mer distance for single-linkage clustering within a species")
    parser.add_argument("--min-species-size", type=int, default=8,
                        help="Minimum assessable isolates in a species before divergence "
                             "diagnostics (and divergence_candidate tagging) are computed")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    parser.add_argument("--mafft-bin", default=None, help="Path to mafft (default: PATH lookup)")
    parser.add_argument("--iqtree-bin", default=None, help="Path to iqtree (default: PATH lookup)")
    parser.add_argument("--phylo-min-taxa", type=int, default=4,
                        help="Minimum candidate isolates in a species before building segment trees")
    parser.add_argument("--bootstrap", type=int, default=1000,
                        help="IQ-TREE ultrafast bootstrap replicates (>=1000)")
    parser.add_argument("--min-bootstrap-support", type=float, default=95.0,
                        help="Minimum ultrafast bootstrap support to call a topology conflict supported")
    parser.add_argument("--phylo-seed", type=int, default=42,
                        help="Fixed IQ-TREE random seed, for reproducible trees")
    parser.add_argument("--phylo-threads", type=int, default=max(1, (os.cpu_count() or 4) - 1),
                        help="Threads for mafft and iqtree")
    parser.add_argument("--phylo-scope", choices=("species", "phylogroup", "both"), default="both",
                        help="species: within-species trees only. phylogroup: trees spanning "
                             "every species in a reservoir-host phylogroup (see ah.PHYLOGROUP), "
                             "for interspecies signal a species-only tree structurally cannot "
                             "show. both: run each independently.")
    parser.add_argument("--phylo-outgroup", default=None,
                        help="Accession or isolate_key to root trees on, via IQ-TREE's -o.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    paths = ah.resolve_inputs(Path(args.input), output=Path(args.output) if args.output else None)
    ah.configure_logging(paths.log_file, verbose=args.verbose)
    logging.info("Orthohantavirus phylogenetics stage (Stage 2)")
    logging.info("Package: %s", paths.data_dir)
    logging.info("Output : %s", paths.output_dir)
    run_phylogeny(paths, args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Verify the CLI parses and the module still imports**

```bash
python3 -c "import phylogeny as ph; ph.parse_args(['--input', '.']); print('parse_args ok')"
python3 phylogeny.py --help > /dev/null && echo "CLI help ok"
```

Expected: `parse_args ok` and `CLI help ok`. No automated test covers `main`/`parse_args` directly
here, matching `all_hanta.py`'s own CLI (which also has no `ParseArgsTest`) — the logic these
wrap (`run_phylogenetics_stage`, `collect_phylo_candidates`, etc.) is fully covered by the tests
Task 3 moves.

- [ ] **Step 3: Commit**

```bash
git add phylogeny.py
git commit -m "$(cat <<'EOF'
Give phylogeny.py its own CLI and top-level run function

run_phylogeny() mirrors all_hanta.py's ingest chain through
compute_segment_divergence (needed for divergence_by_key), then calls
the moved run_phylogenetics_stage and writes ortho_reassortment_table.csv
directly, independent of any of all_hanta.py's output files.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Split `test_all_hanta.py` into `test_phylogeny.py` + trimmed `test_all_hanta.py`

At this point `all_hanta.py` still has its own (soon-to-be-deleted) copy of every moved
function, so both the new `test_phylogeny.py` (calling `ph.X`) and the old `test_all_hanta.py`
(still calling `ah.X` for anything not yet trimmed) can be verified green independently before
Task 5 removes the duplication. This ordering means the suite is never red mid-plan.

**Files:**
- Create: `test_phylogeny.py`
- Modify: `test_all_hanta.py`
- Create (temporary, delete after use): `/tmp/split_tests.py`

**Interfaces:**
- Consumes: `phylogeny as ph` (Task 1/2), `all_hanta as ah` (unchanged).
- Produces: `test_phylogeny.py` and `test_all_hanta.py`, both passing under `python3 -m unittest`.

**Verified move/stay/delete inventory (by exact line range in the current file):**

| Range | Name | Disposition |
|---|---|---|
| 37–43 | `try: import Bio.Phylo... _HAVE_PHYLO_DEPS` block | moves |
| 46–62 | `FakeConfidence`, `FakeTree` | moves |
| 64–90 | `SegmentRankMetricsBackwardCompatTest` | **stays**, one internal fix (see Step 3) |
| 91–155 | `ClassifyReassortmentCallTest` | moves |
| 156–188 | `BootstrapSupportRootFallbackTest` | moves |
| 189–212 | `_make_isolate` helper | **duplicated** into both files (used by moving classes at 227/278/322 and staying classes at 502/503/564) |
| 215–230 | `_make_candidates` helper | moves (uses `ah.PhyloCandidate` → becomes `ph.PhyloCandidate`) |
| 233–248 | `_make_accession` helper | stays (only used by `SpeciesEffortStatsTest`) |
| 251–266 | `_cds_index_for` helper | moves |
| 269–331 | `CollectPhyloCandidatesScopeTest` | moves |
| 332–363 | `GroupCandidatesTest` | moves |
| 364–483 | `PhylogeneticsScopeEndToEndTest` | moves |
| 484–526 | `SpeciesEffortStatsTest` | stays unchanged |
| 527–555 | `PhyloCallByIsolateTest` | moves |
| 556–589 | `SurveillancePriorityReassortmentTermTest` | stays unchanged |
| 590–618 | `BuildRiskRowsPhyloIntegrationTest` | **deleted** — tests `build_risk_rows`'s `phylo_call_by_key` param, which Task 5 removes entirely (spec §5/§6 correction) |

Import trim for `test_all_hanta.py`: drop `subprocess` (already unused, 0 call sites — pre-existing
dead import, fine to drop while touching this file), `tempfile`, `sys`, `StringIO`, the
`Bio.Phylo`/`dendropy` try/except block — none are used by the three staying classes (verified:
`tempfile.` only at line 420, `sys.` only at 375/396, `StringIO(` only at 164/172/179/184, all
inside moving classes). Keep `random`? No — `random.` usage (216, 256, 400–414) is entirely
inside moving code; drop it too. Keep `numpy as np` — used at lines 72/79/80, inside the staying
`SegmentRankMetricsBackwardCompatTest`, nowhere else. Add `import phylogeny as ph` (needed for
Step 3's fix).

- [ ] **Step 1: Write the split script**

```python
# /tmp/split_tests.py
import re
from pathlib import Path

src = Path("test_all_hanta.py").read_text().splitlines(keepends=True)


def extract(start, end):
    return "".join(src[start - 1:end])


# ---- Pieces moving to test_phylogeny.py, in original file order ----
phylo_deps_block = extract(37, 43)
fake_classes = extract(46, 62)
classify_test = extract(91, 155)
bootstrap_test = extract(156, 188)
make_isolate = extract(189, 212)          # duplicated -- also copied into test_all_hanta.py
make_candidates = extract(215, 230)
cds_index_for = extract(251, 266)
collect_test = extract(269, 331)
group_test = extract(332, 363)
e2e_test = extract(364, 483)
phylo_call_test = extract(527, 555)

phylo_body = "\n".join(part.rstrip("\n") for part in [
    fake_classes, classify_test, bootstrap_test, make_isolate, make_candidates,
    cds_index_for, collect_test, group_test, e2e_test, phylo_call_test,
]) + "\n"

# make_candidates constructs ah.PhyloCandidate -- that symbol moved to phylogeny.py.
phylo_body = re.sub(r"(?<!ph\.)\bah\.PhyloCandidate\b", "ph.PhyloCandidate", phylo_body)
# Every remaining ah.<moved-symbol> call in the moving classes must become ph.<symbol>.
for name in ("classify_reassortment_call", "collect_phylo_candidates", "_assignment_basis",
             "_candidate_group_key", "group_candidates", "group_candidates_by_species",
             "_run_phylogenetics_scope", "_phylo_call_by_isolate"):
    phylo_body = re.sub(rf"(?<!ph\.)\bah\.{name}\b", f"ph.{name}", phylo_body)

phylo_header = '''#!/usr/bin/env python3
"""Regression tests for phylogeny.py -- moved from test_all_hanta.py on 2026-08-13 when the
--phylogenetics code was split out of all_hanta.py into its own module. See
docs/superpowers/specs/2026-08-13-three-stage-pipeline-split-design.md for why.

stdlib unittest, matching the repo's house style. Run with:

    python3 -m unittest test_phylogeny -v

Fakes Bio.Phylo tree objects directly for the bootstrap/classify tests; fakes mafft and iqtree
as tiny throwaway scripts for the end-to-end test, since neither binary is installed in this
environment -- skipped automatically if Bio.Phylo/dendropy aren't importable.
"""

import random
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path

import all_hanta as ah
import phylogeny as ph

'''
phylo_header += phylo_deps_block + "\n\n"
phylo_footer = '''

if __name__ == "__main__":
    unittest.main()
'''

Path("test_phylogeny.py").write_text(phylo_header + phylo_body + phylo_footer)

# ---- test_all_hanta.py: keep only the staying pieces, in original order ----
# Written directly rather than via extract-and-replace on lines 1-43: the try/except
# Bio.Phylo block spans lines 37-43, one past extract(1, 36)'s range, so a replace()
# targeting text that includes it would silently no-op against a too-short extract.
original_docstring = extract(1, 23)  # through the closing triple-quote, unchanged
staying_header = original_docstring + '''
import unittest

import numpy as np

import all_hanta as ah
import phylogeny as ph

'''
segment_rank_test = extract(64, 90)
make_accession = extract(233, 248)
species_effort_test = extract(484, 526)
surveillance_test = extract(556, 589)
footer = extract(618, 619)  # if __name__ == "__main__":\n    unittest.main()

staying_body = "\n".join(part.rstrip("\n") for part in [
    segment_rank_test, make_isolate, make_accession, species_effort_test, surveillance_test,
]) + "\n\n\n" + footer

# The one internal fix: _segment_rank_metrics moved to phylogeny.py.
staying_body = staying_body.replace("ah._segment_rank_metrics(", "ph._segment_rank_metrics(")
assert "ah._segment_rank_metrics(" not in staying_body

Path("test_all_hanta.py").write_text(staying_header + staying_body)
print("Wrote test_phylogeny.py and test_all_hanta.py")
```

- [ ] **Step 2: Run the split script**

```bash
python3 /tmp/split_tests.py
```

Expected: `Wrote test_phylogeny.py and test_all_hanta.py`, no assertion error.

- [ ] **Step 3: Verify the assertion in Step 1 caught the internal fix, then run both test files**

```bash
python3 -m unittest test_all_hanta -v 2>&1 | tail -20
python3 -m unittest test_phylogeny -v 2>&1 | tail -20
```

Expected: both `OK`. `test_all_hanta` should report exactly the 3 staying classes'-worth of
tests (8 test methods: 1 + 2 + 5 across `SegmentRankMetricsBackwardCompatTest`,
`SpeciesEffortStatsTest`, `SurveillancePriorityReassortmentTermTest`). `test_phylogeny` should
report the rest. If either fails on a `NameError`/`AttributeError` for `ah.<X>` or `ph.<X>`, that
call site was missed by Step 1's regex list — add it and re-run from Step 1 (both scripts are
idempotent against the original committed `test_all_hanta.py`... **except this is the second
run of a script that already overwrote `test_all_hanta.py`** — if you need to retry, first
restore with `git checkout -- test_all_hanta.py`, since Step 1's script reads the file it also
overwrites.

- [ ] **Step 4: Run the full suite to confirm nothing else regressed**

```bash
python3 -m unittest discover -s . -p "test_*.py" 2>&1 | tail -15
```

Expected: `Ran <N> tests ... OK` where `<N>` matches the pre-split count (96) minus the 2 deleted
`BuildRiskRowsPhyloIntegrationTest` methods, i.e. 94. (`all_hanta.py` itself is untouched in this
task, so its old phylogenetics code still backs any stray `ah.` reference that Step 1 missed —
this task cannot silently break Stage 1's own behavior, only fail loudly on a bad test split.)

- [ ] **Step 5: Delete the temporary script and commit**

```bash
rm /tmp/split_tests.py
git add test_phylogeny.py test_all_hanta.py
git commit -m "$(cat <<'EOF'
Split phylogenetics tests into test_phylogeny.py

Moves the 6 phylogenetics-dependent test classes and their shared
fixtures to test_phylogeny.py (calling ph.X instead of ah.X), following
the repo's one-test-file-per-module convention. Deletes
BuildRiskRowsPhyloIntegrationTest, which tested build_risk_rows's
phylo_call_by_key parameter -- Task 5 removes that parameter entirely,
so the behavior it tested no longer exists in either module.
_make_isolate is duplicated into both files since both staying and
moving classes use it.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Fix the Stage-3 import surface

**Files:**
- Modify: `reassortment_detector.py:15` (import), `:151-152` (calls)
- Modify: `detector_gates.py:14` (import), `:513,518,528,529,532,572,581` (calls)
- Modify: `detector_reference.py:13` (import), `:40` (call)

**Interfaces:**
- Consumes: `phylogeny as ph` (Task 1/2).
- Produces: no change to any function signature in these three files — call-site rewiring only.

- [ ] **Step 1: Fix `reassortment_detector.py`**

```python
# Edit old_string / new_string:
```
old:
```python
import all_hanta as ah
```
new:
```python
import all_hanta as ah
import phylogeny as ph
```
old:
```python
    mafft_bin = ah.resolve_binary(args.mafft_bin, "mafft")
    iqtree_bin = ah.resolve_binary(args.iqtree_bin, "iqtree")
```
new:
```python
    mafft_bin = ph.resolve_binary(args.mafft_bin, "mafft")
    iqtree_bin = ph.resolve_binary(args.iqtree_bin, "iqtree")
```

- [ ] **Step 2: Fix `detector_reference.py`**

old:
```python
import all_hanta as ah
```
new:
```python
import all_hanta as ah
import phylogeny as ph
```
old:
```python
    cds_index = ah.parse_cds_fasta(paths.cds_fasta)
```
new:
```python
    cds_index = ph.parse_cds_fasta(paths.cds_fasta)
```

- [ ] **Step 3: Fix `detector_gates.py`**

old:
```python
import all_hanta as ah
```
new:
```python
import all_hanta as ah
import phylogeny as ph
```

Then, five call-site fixes (use `replace_all` where the same substitution recurs):

- `List["ah.PhyloCandidate"], str]` → `List["ph.PhyloCandidate"], str]`
- `members: List["ah.PhyloCandidate"] = []` → `members: List["ph.PhyloCandidate"] = []`
- `seq = ah.select_primary_cds(records, segment)` → `seq = ph.select_primary_cds(records, segment)`
- `if len(seq) >= ah.MIN_CDS_LENGTH:` → `if len(seq) >= ph.MIN_CDS_LENGTH:`
- `members.append(ah.PhyloCandidate(isolate=isolate, tags=["comparator"], sequences=sequences))` → `members.append(ph.PhyloCandidate(isolate=isolate, tags=["comparator"], sequences=sequences))`
- `candidate_member = ah.PhyloCandidate(` → `candidate_member = ph.PhyloCandidate(`
- `rows = ah._analyze_group(` → `rows = ph._analyze_group(`

- [ ] **Step 4: Verify all Stage-3 tests still pass**

```bash
python3 -m unittest test_detector_gates test_detector_reference test_reassortment_detector -v 2>&1 | tail -20
```

Expected: `OK`. This exercises `resolve_binary`, `PhyloCandidate`, `select_primary_cds`,
`MIN_CDS_LENGTH`, `_analyze_group`, and `parse_cds_fasta` through real (unmocked) code paths —
if any call site was missed, this fails with `AttributeError` immediately, not silently.

- [ ] **Step 5: Commit**

```bash
git add reassortment_detector.py detector_gates.py detector_reference.py
git commit -m "$(cat <<'EOF'
Repoint Stage 3's phylogenetics calls from all_hanta to phylogeny

Six call sites across reassortment_detector.py, detector_gates.py, and
detector_reference.py switch from ah.X to ph.X for the six symbols that
moved in the pipeline split: resolve_binary, PhyloCandidate,
select_primary_cds, MIN_CDS_LENGTH, _analyze_group, parse_cds_fasta.
No logic change, verified by the existing test suite passing unmocked.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: Trim `all_hanta.py`

Everything referencing the moved code has been repointed (Tasks 3–4), so it's now safe to delete
the duplicated originals.

**Files:**
- Modify: `all_hanta.py`
- Create (temporary, delete after use): `/tmp/trim_all_hanta.py`

**Interfaces:**
- Produces: `all_hanta.py` with none of the 26 moved names defined, `build_risk_rows` with its
  `phylo_call_by_key` parameter removed, no `--phylogenetics`/`--mafft-bin`/etc CLI args, no
  `ortho_segment_divergence.csv`/`ortho_ml_segments.csv` writes.

- [ ] **Step 1: Simplify `build_risk_rows`** (do this by hand first, it's small and precise —
      not part of the bulk-deletion script)

old:
```python
def build_risk_rows(
    isolates: Sequence[IsolateRecord],
    divergence: Dict[str, SegmentDivergenceProfile],
    phylo_call_by_key: Optional[Dict[str, str]] = None,
) -> List[Dict[str, Any]]:
    """phylo_call_by_key is None (or empty) whenever --phylogenetics wasn't run, or ran
    but failed for this isolate's tree -- surveillance_priority already treats a missing
    call as "no reassortment term", so this table is byte-identical to a phylogenetics-
    free run except for the always-present, always-zero priority_reassortment_term column.
    """
    phylo_call_by_key = phylo_call_by_key or {}
    rows = []
    for isolate in isolates:
        risk = species_risk_profile(isolate.species_name)
        profile = divergence.get(isolate.isolate_key)
        phylo_call = phylo_call_by_key.get(isolate.isolate_key)
        terms = surveillance_priority(isolate, risk, profile, phylo_call)
```
new:
```python
def build_risk_rows(
    isolates: Sequence[IsolateRecord],
    divergence: Dict[str, SegmentDivergenceProfile],
) -> List[Dict[str, Any]]:
    """ortho_risk_table.csv is always phylo-blind now (Stage 2's phylogenetics analysis is a
    separate script/run; combine on isolate_key -- see README's Outputs table). This is
    byte-identical to the pre-split phylogenetics-free run's output.
    """
    rows = []
    for isolate in isolates:
        risk = species_risk_profile(isolate.species_name)
        profile = divergence.get(isolate.isolate_key)
        terms = surveillance_priority(isolate, risk, profile, None)
```

- [ ] **Step 2: Verify `build_risk_rows`'s only caller compiles** (its call site is fixed in
      Step 4 below — this step is a syntax sanity check only)

```bash
python3 -c "import ast; ast.parse(open('all_hanta.py').read())" && echo "syntax ok"
```

Expected: `syntax ok`.

- [ ] **Step 3: Write and run the bulk-deletion script**

```python
# /tmp/trim_all_hanta.py
from pathlib import Path

lines = Path("all_hanta.py").read_text().splitlines(keepends=True)

# Delete ranges in reverse order so earlier line numbers stay valid as later ones are cut.
# (1-indexed, inclusive; matches the plan's Verified Extraction Inventory table exactly.)
for start, end in [(2295, 3115), (2028, 2066), (732, 745), (451, 456)]:
    del lines[start - 1:end]

Path("all_hanta.py").write_text("".join(lines))
print("Trimmed all_hanta.py to", len(lines), "lines")
```

```bash
python3 /tmp/trim_all_hanta.py
```

Expected: `Trimmed all_hanta.py to <N> lines` (roughly 3270, down from 4169).

- [ ] **Step 4: Remove the CLI phylo arg group and the run_pipeline phylo block**

old (in `parse_args`):
```python
    # --concordance-tolerance / --discordance-threshold were removed along with the
    # withdrawn reassortment classifier; no threshold made its ranking agree with the
    # literature, so exposing one invited false confidence.
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")

    phylo = parser.add_argument_group(
        "phylogenetics (opt-in; requires mafft + iqtree -- see README)")
    phylo.add_argument("--phylogenetics", action="store_true",
                       help="Run MAFFT + IQ-TREE reassortment analysis over human-affecting, "
                            "divergence-candidate and complete-triplet isolates")
    phylo.add_argument("--mafft-bin", default=None, help="Path to mafft (default: PATH lookup)")
    phylo.add_argument("--iqtree-bin", default=None, help="Path to iqtree (default: PATH lookup)")
    phylo.add_argument("--phylo-min-taxa", type=int, default=4,
                       help="Minimum candidate isolates in a species before building segment trees")
    phylo.add_argument("--bootstrap", type=int, default=1000,
                       help="IQ-TREE ultrafast bootstrap replicates (>=1000)")
    phylo.add_argument("--min-bootstrap-support", type=float, default=95.0,
                       help="Minimum ultrafast bootstrap support to call a topology conflict supported")
    phylo.add_argument("--phylo-seed", type=int, default=42,
                       help="Fixed IQ-TREE random seed, for reproducible trees")
    phylo.add_argument("--phylo-threads", type=int, default=max(1, (os.cpu_count() or 4) - 1),
                       help="Threads for mafft and iqtree")
    phylo.add_argument("--phylo-scope", choices=("species", "phylogroup", "both"), default="both",
                       help="species: within-species trees only -- can only ever show "
                            "reassortment WITHIN a species. phylogroup: trees spanning every "
                            "species in a reservoir-host phylogroup (see PHYLOGROUP), for "
                            "interspecies signal a species-only tree structurally cannot show. "
                            "both: run each independently; a candidate isolate can then appear "
                            "in ortho_reassortment_table.csv up to twice, once per scope "
                            "(see the analysis_scope column)")
    phylo.add_argument("--phylo-outgroup", default=None,
                       help="Accession or isolate_key to root trees on, via IQ-TREE's -o. Only "
                            "applied to a given tree if that isolate is one of its taxa. Without "
                            "this, phylogroup trees are unrooted: they can show that segments "
                            "from different species cluster together, but not which species was "
                            "donor and which was recipient. No default is assumed.")
    return parser.parse_args(argv)
```
new:
```python
    # --concordance-tolerance / --discordance-threshold were removed along with the
    # withdrawn reassortment classifier; no threshold made its ranking agree with the
    # literature, so exposing one invited false confidence.
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    return parser.parse_args(argv)
```

old (in `run_pipeline`, immediately after `write_json("ortho_summary.json", summary)`):
```python
    # ortho_risk_table.csv is deliberately the one core table NOT written above, and is
    # built after the phylogenetics block below rather than before it (Task 5): a
    # bootstrap-supported reassortment call from --phylogenetics now feeds
    # surveillance_priority's score, which means the risk table needs phylo_rows to exist
    # first. Every OTHER core output is already safely on disk by this point, so nothing
    # about the optional phylogenetics stage -- including this reordering -- can cost the
    # run any of them (all_hanta.py had the whole phylogenetics block *before* "Writing
    # outputs" until a code review caught it; it was moved after core outputs, and stays
    # exception-isolated from them here for the same reason).
    #
    # If --phylogenetics isn't passed, or the stage fails, phylo_call_by_key stays empty
    # and ortho_risk_table.csv is written exactly as it always was -- a phylogenetics
    # failure costs the run reassortment evidence in the risk score, never the table itself.
    phylo_call_by_key: Dict[str, str] = {}
    if args.phylogenetics:
        logging.info("Phylogenetics stage requested (--phylogenetics)")
        try:
            cds_index = parse_cds_fasta(paths.cds_fasta)
            logging.info("  cds: %s accessions with at least one CDS", len(cds_index))
            phylo_rows = run_phylogenetics_stage(isolates, divergence_by_key, cds_index, out, args)
            write_csv("ortho_reassortment_table.csv", phylo_rows)
            phylo_call_by_key = _phylo_call_by_isolate(phylo_rows)
        except (Exception, SystemExit) as exc:  # noqa: BLE001 - risk table below is unaffected
            logging.error(
                "Phylogenetics stage failed; ortho_risk_table.csv will be written without "
                "reassortment evidence, same as a run without --phylogenetics: %s", exc)
    else:
        logging.info(
            "Skipping phylogenetics stage (pass --phylogenetics to run MAFFT + IQ-TREE "
            "reassortment analysis; requires both on PATH, see README)")

    write_csv("ortho_risk_table.csv", build_risk_rows(isolates, divergence_by_key, phylo_call_by_key))
```
new:
```python
    # ortho_risk_table.csv is always phylo-blind now -- Stage 2 (phylogeny.py) is a separate
    # script/run; combine its ortho_reassortment_table.csv on isolate_key if you need both
    # (see README's Outputs table). This is byte-identical to the pre-split behavior when
    # --phylogenetics wasn't passed, which was already the common case.
    write_csv("ortho_risk_table.csv", build_risk_rows(isolates, divergence_by_key))
```

old (still in `run_pipeline`'s output-writing block):
```python
    write_csv("ortho_isolate_table.csv", build_isolate_rows(isolates, divergence_by_key))
    write_csv("ortho_host_table.csv", build_host_rows(accessions, isolates))
    write_csv("ortho_segment_divergence.csv", build_segment_divergence_rows(divergence))
    write_csv("ortho_ml_feature_matrix.csv", ml_rows)
    write_csv("ortho_ml_segments.csv", build_ml_segment_rows(accessions))
    write_csv("ortho_ml_sequence_table.csv", build_ml_sequence_rows(isolates, genomic_index))
```
new:
```python
    write_csv("ortho_isolate_table.csv", build_isolate_rows(isolates, divergence_by_key))
    write_csv("ortho_host_table.csv", build_host_rows(accessions, isolates))
    write_csv("ortho_ml_feature_matrix.csv", ml_rows)
    write_csv("ortho_ml_sequence_table.csv", build_ml_sequence_rows(isolates, genomic_index))
```

(Leave `build_segment_divergence_rows` and `build_ml_segment_rows` function definitions in place
for now — deleting now-unused functions is a separate, mechanical dead-code pass; flag them but
don't remove in this task, to keep this task's diff focused on the phylogenetics split. If a
linter/review flags them as dead in Task 5's own review, remove then.)

- [ ] **Step 5: Verify `all_hanta.py` imports and the automated symbol-boundary check passes**

```bash
python3 -c "import all_hanta as ah; print('ok')"
python3 -c "
import all_hanta as ah
import phylogeny as ph
moved = ['HOST_ORDER_PHYLOGROUP_FALLBACK', 'PhyloCandidate', '_segment_rank_metrics',
         'MIN_CDS_LENGTH', 'HUMAN_AFFECTING_CLASSES', 'is_human_affecting_or_likely',
         'parse_cds_fasta', 'select_primary_cds', 'collect_phylo_candidates',
         '_phylogroup_for_isolate', '_phylogroup_for_candidate', '_assignment_basis',
         '_candidate_group_key', 'group_candidates', 'group_candidates_by_species',
         'resolve_binary', 'write_segment_fasta', 'run_mafft', 'run_iqtree',
         '_patristic_distance_matrix', 'tree_divergence_for_species',
         '_nearest_supported_descendant', 'bootstrap_support_for_pair',
         'classify_reassortment_call', '_phylo_row', '_resolve_outgroup_label',
         '_analyze_group', '_run_phylogenetics_scope', 'run_phylogenetics_stage',
         '_REASSORTMENT_CALL_RANK', '_phylo_call_by_isolate']
still_in_ah = [n for n in moved if hasattr(ah, n)]
missing_from_ph = [n for n in moved if not hasattr(ph, n)]
assert not still_in_ah, f'still defined in all_hanta.py: {still_in_ah}'
assert not missing_from_ph, f'missing from phylogeny.py: {missing_from_ph}'
print('boundary check passed:', len(moved), 'symbols confirmed moved cleanly')
"
```

Expected: `ok`, then `boundary check passed: 31 symbols confirmed moved cleanly`. This is the
automated verification promised during design — not a human re-reading the diff, a script that
fails loudly if either half of the move is incomplete.

- [ ] **Step 6: Run the full suite**

```bash
python3 -m unittest discover -s . -p "test_*.py" -v 2>&1 | tail -30
```

Expected: `OK`, same 94-test count as Task 3's checkpoint (this task removes code, not tests).

- [ ] **Step 7: Delete the temporary script and commit**

```bash
rm /tmp/trim_all_hanta.py
git add all_hanta.py
git commit -m "$(cat <<'EOF'
Remove phylogenetics code from all_hanta.py, now live in phylogeny.py

Deletes the four extracted line ranges, the --phylogenetics/--mafft-bin/
etc CLI args, and the run_pipeline phylogenetics block. build_risk_rows
loses its phylo_call_by_key parameter -- ortho_risk_table.csv is now
unconditionally phylo-blind, byte-identical to a pre-split run without
--phylogenetics. Verified by an automated boundary check (none of the
31 moved names remain in all_hanta, all are present in phylogeny) plus
the full test suite.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Regenerate `output/` and drop the two pruned files

**Files:**
- Modify: `output/*.csv`, `output/*.json`, `output/ortho_output_checksums.txt`, `output/ortho_pipeline.log`
- Delete: `output/ortho_segment_divergence.csv`, `output/ortho_ml_segments.csv`

- [ ] **Step 1: Run Stage 1 fresh**

```bash
python3 all_hanta.py --input . --output output --verbose 2>&1 | tail -20
```

Expected: completes (~16s per README), no traceback, final log line reports accession/species/
isolate/triplet counts matching the current committed `output/ortho_summary.json`'s counts
(spot-check: `grep -A5 '"counts"' output/ortho_summary.json` before and after should show
identical `accession_records`/`canonical_species`/`isolates`/`complete_triplets` values, since
this is a re-run of the same command against the same input, not a data change).

- [ ] **Step 2: Verify only the two intended files disappeared**

```bash
git status --porcelain output/
```

Expected: `ortho_segment_divergence.csv` and `ortho_ml_segments.csv` show as deleted (`D`); every
other tracked file shows modified (`M`) only if its content actually changed (`ortho_pipeline.log`
always changes, it's timestamped; `ortho_output_checksums.txt` always changes, it hashes the
other files). No file should show as newly added — Stage 1 doesn't create anything not already
tracked.

- [ ] **Step 3: Spot-check no other file's row count shifted**

```bash
git diff --stat output/ | grep -v "ortho_pipeline.log\|ortho_output_checksums.txt\|ortho_segment_divergence.csv\|ortho_ml_segments.csv"
```

Expected: any files listed here should show only whitespace-level or zero diffs — if
`ortho_accession_table.csv` or similar shows hundreds of changed lines, something besides the
intended removal changed and needs investigation before committing (per spec §5: "no other
content should change").

- [ ] **Step 4: Commit**

```bash
git add output/
git commit -m "$(cat <<'EOF'
Regenerate output/ without the two pruned files

ortho_segment_divergence.csv and ortho_ml_segments.csv are dropped
(duplicated data, see the pipeline-split spec section 5). Everything
else is a same-command re-run against unchanged input -- only
ortho_pipeline.log and ortho_output_checksums.txt change beyond that.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: Update `README.md`

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Replace the `--phylogenetics` subsection with a Pipeline stages section**

old:
```markdown
### Optional: `--phylogenetics` (MAFFT + IQ-TREE reassortment analysis)

| Argument | Default | Purpose |
|---|---|---|
| `--phylogenetics` | off | Run the MAFFT + IQ-TREE stage described below |
| `--mafft-bin` | PATH lookup | Path to `mafft` |
| `--iqtree-bin` | PATH lookup | Path to `iqtree` |
| `--phylo-min-taxa` | `4` | Minimum candidate isolates in a species before building segment trees |
| `--bootstrap` | `1000` | IQ-TREE ultrafast bootstrap replicates |
| `--min-bootstrap-support` | `95.0` | Minimum ufboot support to call a topology conflict *supported* |
| `--phylo-seed` | `42` | Fixed IQ-TREE random seed, for reproducible trees |
| `--phylo-threads` | cores − 1 | Threads for `mafft` and `iqtree` |

Requires `mafft` and `iqtree` (IQ-TREE 3.x; `iqtree2` builds work too, pass `--iqtree-bin`) plus
Biopython and `dendropy` -- none of which are needed for a default run. Set up once with the
pinned environment file rather than an ad hoc `conda create`:

```bash
conda config --set channel_priority strict
conda env create -f environment.yml
```

Channel order matters here, not just as style: `-c bioconda -c conda-forge` (bioconda first)
resolves Biopython from bioconda's own channel, whose newest build is `1.70` and drags in a
`Python 3.6` / `numpy 1.12` toolchain to satisfy it -- silently, with no error. `environment.yml`
lists `conda-forge` first and `bioconda` second (the order bioconda's own docs recommend) so
Biopython, pandas, numpy and pyarrow resolve from conda-forge's current builds, while `mafft` and
`iqtree` -- which only exist in bioconda -- still resolve there. `mafft`/`iqtree` are pinned to
exact versions (`7.525` / `3.1.3`); the Python-side packages are floors, not patches, since exact
patch availability on conda-forge isn't something this repo can assert.

Then either `conda activate hanta_phylo` before running this script, or leave your normal
environment active and pass `--mafft-bin`/`--iqtree-bin` pointing at that env's `bin/`. **Don't
just prepend that env's `bin/` to `PATH`** without activating fully -- an unpinned env resolved
this way has been observed to install a Python 3.6 interpreter as `python3`, which would shadow
the one actually running this script.

Runtime is dominated by IQ-TREE's ModelFinder + ultrafast bootstrap and scales with both the
number of candidate isolates per species and CDS length (L segments are the slowest, ~6.4 kb).
A handful of taxa finishes in seconds; species with hundreds of qualifying isolates (PUUV, TULV,
HTNV, SEOV in this package) can take much longer. There is no built-in subsampling -- every
qualifying isolate is analyzed.

---
```
new:
```markdown
## Pipeline stages

The pipeline is three independent scripts, each runnable any time after the raw NCBI package
exists (none requires another to have already run):

| Stage | Script | Needs | Produces |
|---|---|---|---|
| 1 | `all_hanta.py` | `pandas`/`pyarrow`/`numpy` only | The core tables above, including `ortho_ml_sequence_table.csv` |
| 2 | `phylogeny.py` | `mafft`/`iqtree` on PATH, Biopython, `dendropy` | `ortho_reassortment_table.csv`, `phylogenetics/<species>/<segment>/` |
| 3 | `reassortment_detector.py` | `mafft`/`iqtree` on PATH, Biopython, `dendropy` | Score a new, not-yet-published candidate genome against the 5-gate framework |

`ortho_risk_table.csv` (Stage 1) is always phylo-blind. Combine it with Stage 2's evidence by
joining `ortho_reassortment_table.csv` on `isolate_key` — there's no automatic recombination
between the two scripts by design (see `docs/superpowers/specs/2026-08-13-three-stage-pipeline-
split-design.md` §4 for why).

### Stage 2: `phylogeny.py` (MAFFT + IQ-TREE reassortment analysis)

```bash
python3 phylogeny.py --input . --output output
```

| Argument | Default | Purpose |
|---|---|---|
| `--mafft-bin` | PATH lookup | Path to `mafft` |
| `--iqtree-bin` | PATH lookup | Path to `iqtree` |
| `--phylo-min-taxa` | `4` | Minimum candidate isolates in a species before building segment trees |
| `--bootstrap` | `1000` | IQ-TREE ultrafast bootstrap replicates |
| `--min-bootstrap-support` | `95.0` | Minimum ufboot support to call a topology conflict *supported* |
| `--phylo-seed` | `42` | Fixed IQ-TREE random seed, for reproducible trees |
| `--phylo-threads` | cores − 1 | Threads for `mafft` and `iqtree` |

Requires `mafft` and `iqtree` (IQ-TREE 3.x; `iqtree2` builds work too, pass `--iqtree-bin`) plus
Biopython and `dendropy`. Set up once with the pinned environment file rather than an ad hoc
`conda create`:

```bash
conda config --set channel_priority strict
conda env create -f environment.yml
```

Channel order matters here, not just as style: `-c bioconda -c conda-forge` (bioconda first)
resolves Biopython from bioconda's own channel, whose newest build is `1.70` and drags in a
`Python 3.6` / `numpy 1.12` toolchain to satisfy it -- silently, with no error. `environment.yml`
lists `conda-forge` first and `bioconda` second (the order bioconda's own docs recommend) so
Biopython, pandas, numpy and pyarrow resolve from conda-forge's current builds, while `mafft` and
`iqtree` -- which only exist in bioconda -- still resolve there. `mafft`/`iqtree` are pinned to
exact versions (`7.525` / `3.1.3`); the Python-side packages are floors, not patches, since exact
patch availability on conda-forge isn't something this repo can assert.

Then either `conda activate hanta_phylo` before running this script, or leave your normal
environment active and pass `--mafft-bin`/`--iqtree-bin` pointing at that env's `bin/`. **Don't
just prepend that env's `bin/` to `PATH`** without activating fully -- an unpinned env resolved
this way has been observed to install a Python 3.6 interpreter as `python3`, which would shadow
the one actually running this script.

Runtime is dominated by IQ-TREE's ModelFinder + ultrafast bootstrap and scales with both the
number of candidate isolates per species and CDS length (L segments are the slowest, ~6.4 kb).
A handful of taxa finishes in seconds; species with hundreds of qualifying isolates (PUUV, TULV,
HTNV, SEOV in this package) can take much longer. There is no built-in subsampling -- every
qualifying isolate is analyzed.

---
```

- [ ] **Step 2: Update the Outputs table**

old:
```markdown
| `ortho_segment_divergence.csv` | 9,852 × 17 | Per-segment nearest neighbours and rank separations — **diagnostics, not reassortment calls** |
| `ortho_ml_feature_matrix.csv` | 7,651 × 236 | Raw features + labels, with `group_id` and `year` for leakage-free CV |
| `ortho_ml_segments.csv` | 15,781 × 19 | Segment-level table for per-segment models |
| `ortho_ml_sequence_table.csv` | 206 × 6 | One row per isolate with a coding-complete S+M+L triplet (lab-derived and vaccine-strain isolates excluded): species name, isolate name, raw S/M/L genomic sequences, CFR range |
| `ortho_suppressed.csv` | 1,032 × 16 | Non-representative accessions with the reason |
| `ortho_qc_report.csv` | 15,991 × 11 | Per-accession QC including failed transforms with traceback |
| `ortho_feature_scaling.json` | — | Per-feature mean/std, to be fitted **inside** your CV fold |
| `ortho_summary.json` | — | Provenance, integrity results, counts, methods and stated limitations |
| `ortho_reassortment_table.csv` *(`--phylogenetics` only)* | one row per candidate isolate | Real per-segment tree diagnostics, RF distances, bootstrap-gated `reassortment_call` |
| `phylogenetics/<species>/<segment>/` *(`--phylogenetics` only)* | — | MAFFT alignment + IQ-TREE outputs per species/segment, for inspection or reuse |
```
new:
```markdown
| `ortho_ml_feature_matrix.csv` | 7,651 × 236 | Raw features + labels, with `group_id` and `year` for leakage-free CV |
| `ortho_ml_sequence_table.csv` | 206 × 6 | One row per isolate with a coding-complete S+M+L triplet (lab-derived and vaccine-strain isolates excluded): species name, isolate name, raw S/M/L genomic sequences, CFR range |
| `ortho_suppressed.csv` | 1,032 × 16 | Non-representative accessions with the reason |
| `ortho_qc_report.csv` | 15,991 × 11 | Per-accession QC including failed transforms with traceback |
| `ortho_feature_scaling.json` | — | Per-feature mean/std, to be fitted **inside** your CV fold |
| `ortho_summary.json` | — | Provenance, integrity results, counts, methods and stated limitations |
| `ortho_reassortment_table.csv` *(Stage 2, `phylogeny.py`)* | one row per candidate isolate | Real per-segment tree diagnostics, RF distances, bootstrap-gated `reassortment_call` — join on `isolate_key` with `ortho_risk_table.csv`, which is always phylo-blind |
| `phylogenetics/<species>/<segment>/` *(Stage 2, `phylogeny.py`)* | — | MAFFT alignment + IQ-TREE outputs per species/segment, for inspection or reuse |
```

Also remove the `ortho_ml_segments.csv` and `ortho_segment_divergence.csv` `<details>` blocks
from the Column reference section (search for `<summary><code>ortho_ml_segments.csv</code>` and
`<summary><code>ortho_segment_divergence.csv</code>` — delete each full `<details>...</details>`
block).

- [ ] **Step 3: Update any remaining prose reference to `--phylogenetics`**

```bash
grep -n "\-\-phylogenetics" README.md
```

For each remaining hit (expected: the "Known limitations" section and possibly the reassortment
methodology section), replace `--phylogenetics` with "Stage 2 (`phylogeny.py`)" and re-read the
surrounding sentence to confirm it still reads correctly as prose, not just a mechanical
find-replace — these are explanatory sentences, not tables, so check each one by eye.

- [ ] **Step 4: Verify no more references to the removed flag or files remain**

```bash
grep -n "\-\-phylogenetics\|ortho_segment_divergence\|ortho_ml_segments" README.md
```

Expected: no output (empty).

- [ ] **Step 5: Commit**

```bash
git add README.md
git commit -m "$(cat <<'EOF'
Document the three-stage pipeline in README

Replaces the --phylogenetics subsection (the flag no longer exists)
with a Pipeline stages overview covering all_hanta.py, phylogeny.py,
and reassortment_detector.py, and updates the Outputs table for the
two removed files and phylogeny.py's new ownership of the
phylogenetics-only outputs.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Final Verification

- [ ] Run the complete suite one more time end to end:

```bash
python3 -m unittest discover -s . -p "test_*.py" -v 2>&1 | tail -15
```

Expected: `Ran 94 tests ... OK`.

- [ ] Confirm `git log --oneline -7` shows the seven task commits in order, and `git status` is clean.
