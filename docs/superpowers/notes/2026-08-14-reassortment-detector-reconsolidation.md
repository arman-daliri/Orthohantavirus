# Reassortment detector: reconsolidation into one file

**Date:** 2026-08-14
**Status:** approved, reverses part of `docs/superpowers/specs/2026-08-13-three-stage-pipeline-split-design.md`

## What this reverses

The 2026-08-13 design split Stage 3 into four files: `reassortment_detector.py`
(thin CLI), `detector_gates.py` (the 5-gate scoring logic, 676 lines),
`detector_reference.py` (reference-isolate loading, 80 lines), and
`reassortment_reference_data.py` (curated structural/literature reference
data, 172 lines). That spec's §8 ("Out of scope") explicitly ruled out
further restructuring beyond the `ah.`/`ph.` import fix.

This note records a direct, explicit decision to reverse that split: the
four files are merged back into a single `reassortment_detector.py`
(1,101 lines of production code combined). `phylogeny.py` was already a
single file and needs no change.

## Why

Requested directly, after being told the request conflicted with the
prior day's approved split and being offered the option to keep the
split instead. The four-file layout was one design's judgment call about
navigability versus file count; consolidating back to one file is a
legitimate opposite judgment call, not a correctness fix.

## What does not change

- No behavior, function names, or public call signatures change — this
  is a file-layout change only.
- Test files stay as separate files (`test_detector_gates.py`,
  `test_reassortment_detector.py`, `test_detector_reference.py`,
  `test_reassortment_reference_data.py`); only their imports move to
  point at the single merged module.
- `all_hanta.py` and `phylogeny.py` are untouched.
