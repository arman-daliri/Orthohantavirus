"""Re-run gate4 winner phylogeny using already-completed GA results.

Usage (from project root in Git Bash):
    python gate4_replay.py \
        --mafft-bin  tools/mafft-win/mafft.bat \
        --iqtree-bin tools/iqtree-3.1.3-Windows/bin/iqtree3.exe

Reads output/ga_search/<pool>/seed_20260904.json, calls confirm_winners with
the correct binary paths, then patches confirmation and gate4_phylogeny back
into the JSON in place.  The GA is NOT re-run -- only the phylogeny step.
"""

import argparse
import json
import sys
from pathlib import Path

import rpi_core as rc
from Algorithm import Oracle, confirm_winners, load_reserved_isolate_keys

SEED = 20260904
POOLS = ["terminal_complete", "strict_both"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference-input", default="output")
    parser.add_argument("--ga-search-root",  default="output/ga_search")
    parser.add_argument("--mafft-bin",  required=True)
    parser.add_argument("--iqtree-bin", required=True)
    args = parser.parse_args(argv)

    # Resolve binary paths to absolute so they survive subprocess CWD changes.
    mafft_bin  = str(Path(args.mafft_bin).resolve())
    iqtree_bin = str(Path(args.iqtree_bin).resolve())

    ref_dir = Path(args.reference_input)
    print("Loading reference data ...")
    reference = rc.prep.load_reference_data(ref_dir)
    precomputed = rc._build_precomputed(reference)
    reserved = load_reserved_isolate_keys()

    for pool_name in POOLS:
        pool_out_dir = Path(args.ga_search_root) / pool_name
        seed_json = pool_out_dir / f"seed_{SEED}.json"
        if not seed_json.exists():
            print(f"[{pool_name}] SKIP: {seed_json} not found")
            continue

        payload = json.loads(seed_json.read_text(encoding="utf-8"))
        winners = payload.get("winners", [])
        if not winners:
            print(f"[{pool_name}] SKIP: no winners in {seed_json}")
            continue

        print(f"[{pool_name}] Running gate4 for {len(winners)} winner(s) ...")
        fast_args = rc._score_args(str(pool_out_dir))
        oracle = Oracle(reference, precomputed, fast_args)

        confirmation = confirm_winners(
            oracle, winners, pool_out_dir, with_phylogeny=True,
            mafft_bin=mafft_bin, iqtree_bin=iqtree_bin)

        payload["confirmation"] = confirmation
        payload["gate4_phylogeny"] = [
            {k: w.get(k) for k in ("triplet", "fitness",
                                   "gate4_phylo_status", "gate4_reassortment_call",
                                   "gate4_error")}
            for w in winners
        ]
        seed_json.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
        print(f"[{pool_name}] confirmation: {confirmation}")
        print(f"[{pool_name}] wrote {seed_json}")


if __name__ == "__main__":
    sys.exit(main())
