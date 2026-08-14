#!/usr/bin/env python3
"""Stage 3 entry point -- thin numbered wrapper for reassortment_detector.py.

Run this (not reassortment_detector.py directly) so `ls` sorts the three pipeline
stages in run order: 01_all_hanta.py -> 02_phylogeny.py -> 03_reassortment_detector.py.
See README.md's "Pipeline stages" section for what each stage needs and produces.
"""

import reassortment_detector

if __name__ == "__main__":
    raise SystemExit(reassortment_detector.main())
