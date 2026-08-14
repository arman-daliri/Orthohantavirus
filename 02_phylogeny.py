#!/usr/bin/env python3
"""Stage 2 entry point -- thin numbered wrapper for phylogeny.py.

Run this (not phylogeny.py directly) so `ls` sorts the three pipeline stages in run
order: 01_all_hanta.py -> 02_phylogeny.py -> 03_reassortment_detector.py. See
README.md's "Pipeline stages" section for what each stage needs and produces.
"""

import phylogeny

if __name__ == "__main__":
    raise SystemExit(phylogeny.main())
