"""Scratch, read-only acceptance-test helper: captures pipeline._stage_discriminate's
raw per-row scored output (same seed as output/discrimination_scratch/discrimination_panel.json)
so it can be diffed against that reference file. Not part of the pipeline itself."""
import json
from pathlib import Path

import pipeline
import reassortment_detector as rd

reference = rd.load_reference_data(Path("all_hanta"))
scored = pipeline._stage_discriminate(reference, seed=20260829)
payload = {"seed": 20260829, "groups": scored}
with open("output/pipeline_acceptance/raw_scored_groups.json", "w", encoding="utf-8") as f:
    f.write(json.dumps(rd._json_safe(payload), indent=2, default=str, allow_nan=False))
print("done")
