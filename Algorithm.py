"""population-based combinatorial search over segment-donor assignments.

Spec: docs/superpowers/specs/2026-09-04-combinatorial-search-design.md

Purpose: retrospective, in-silico analysis of public NCBI orthohantavirus data.
An individual is a triplet of EXISTING reference isolates assigned to the three
segment slots (S, M, L). The already-validated scorer in rpi_core.py ranks how
plausible each combination is. No sequence is created, edited, or synthesized
anywhere in this module; the search only recombines accessioned reference
material. rpi_core is treated as a frozen library: imported, never modified.
"""

import argparse
import csv
import json
import random
import statistics
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple
from unittest import mock

import phylogeny as ph
import rpi_core as rc

Triplet = Tuple[str, str, str]

DEFAULT_SEED = 20260904
CALIBRATION_ARTIFACT = Path("data/rpi_calibration.json")


def _field(obj: Any, *names: str) -> Any:
    for name in names:
        if hasattr(obj, name):
            return getattr(obj, name)
    raise AttributeError(f"{type(obj).__name__} has none of {names}")


def load_reserved_isolate_keys(path: Path = CALIBRATION_ARTIFACT) -> set:
    if not path.exists():
        return set()
    artifact = json.loads(path.read_text(encoding="utf-8"))
    block = artifact.get("reserved_for_calibration", {}) or {}
    keys = block.get("isolate_keys") or block.get("isolatekeys") or []
    return set(keys)


def isolate_key_of(isolate: Any) -> str:
    return _field(isolate, "isolate_key", "isolatekey")


def completeness_overlap(reference: Any, reserved_keys: set) -> Dict[str, int]:
    """Set-arithmetic comparison between the two completeness rulers: coding_complete
    (has_complete_triplet) vs terminal_complete (rc.usable_pool's terminally-complete-
    and-not-reserved filter). Both rulers are evaluated over the SAME universe --
    rc.wild_classified_pool (wild, classified, not-reserved) -- so only_coding/
    only_terminal report genuine ruler disagreement, not isolates that were never
    eligible for either pool (non-wild, unclassified, or reserved-for-calibration).
    Reported once per run in the ga_pool_comparison.json payload."""
    eligible = rc.wild_classified_pool(reference, reserved_keys)
    coding_complete = {isolate_key_of(i) for i in eligible if i.has_complete_triplet}
    terminal_complete = {isolate_key_of(i) for i in rc.usable_pool(reference, reserved_keys)}
    both = coding_complete & terminal_complete
    only_coding = coding_complete - terminal_complete
    only_terminal = terminal_complete - coding_complete
    return {"both": len(both), "only_coding": len(only_coding), "only_terminal": len(only_terminal)}


def default_candidate_factory(reference: Any) -> Callable[[str, Dict[str, Any]], Any]:
    factory = getattr(rc, "_candidate_from_isolates", None)
    if factory is None:
        raise SystemExit(
            "rpi_core._candidate_from_isolates not found; grep rpi_core.py for "
            "'def candidate' and wire the correct factory name here")
    return lambda candidate_id, sources: factory(candidate_id, reference, sources)


def neutral_factors(row: Dict[str, Any]) -> str:
    """Which RPI factors in this row are a neutral fallback rather than a measured value,
    as ';'-joined 'Factor:reason' entries, or '' when E, M_pack and M_pol were all
    measured. Derived only from fields rpi_core.py already computed -- descriptive, never
    read by scoring or selection.

    A bare 1.0 does not say which it is: M_pack and M_pol return 1.0 both when scored
    (M_pack saturates at exactly 1.0 for ~38% of calibration observations) and when they
    fall back (no_parental_sequence, no_signal_beyond_tag, undefined_energy). Only the
    status distinguishes them. K is not listed: it is a prior, never a neutral 1, and
    every row carries rpi_k_stratum_provisional instead.
    """
    parts = []
    unknown_e = row.get("rpi_e_unknown_components")
    if unknown_e:
        parts.append("E:" + str(unknown_e).replace(",", "+"))
    for name, key in (("M_pack", "rpi_m_pack_status"), ("M_pol", "rpi_m_pol_status")):
        status = row.get(key)
        if status and status != "scored":
            parts.append(f"{name}:{status}")
    return ";".join(parts)


class Oracle:
    def __init__(self, reference: Any, precomputed: Any, args: Any,
                 candidate_factory: Optional[Callable[[str, Dict[str, Any]], Any]] = None):
        self.reference = reference
        self.precomputed = precomputed
        self.args = args
        self.by_key = {isolate_key_of(i): i for i in reference.isolates}
        self.candidate_factory = candidate_factory or default_candidate_factory(reference)
        self.cache: Dict[Triplet, float] = {}
        self.rows: Dict[Triplet, Dict[str, Any]] = {}
        # (triplet, repr(exception)) for every candidate score_candidate raised on. An
        # unattended multi-hour search must not lose all progress to one bad triplet --
        # scored as 0.0 (same treatment as junk_sequence) and recorded here so the run
        # report can say HOW MANY were skipped and why, rather than silently pretending
        # every candidate scored cleanly.
        self.errors: List[Tuple[Triplet, str]] = []

    def score(self, triplet: Triplet) -> float:
        if triplet in self.cache:
            return self.cache[triplet]
        sources = {"S": self.by_key[triplet[0]], "M": self.by_key[triplet[1]],
                   "L": self.by_key[triplet[2]]}
        try:
            candidate = self.candidate_factory("ga", sources)
            row = rc.score_candidate(candidate, self.reference, self.precomputed, self.args)
        except Exception as exc:  # noqa: BLE001 -- deliberately broad, see errors docstring above
            self.errors.append((triplet, repr(exc)))
            value = 0.0
            self.cache[triplet] = value
            self.rows[triplet] = {"verdict": "scoring_error", "scoring_error": repr(exc)}
            return value
        value = 0.0 if row.get("verdict") == "junk_sequence" else float(row.get("rpi_generation") or 0.0)
        self.cache[triplet] = value
        self.rows[triplet] = {k: row.get(k) for k in (
            "verdict", "rpi_generation", "rpi_retention", "rpi_e_total",
            "rpi_e_isolate_basis", "rpi_e_unknown_components", "rpi_m_pack",
            "rpi_m_pack_status", "rpi_m_pol", "rpi_m_pol_status", "rpi_k_stratum",
            "rpi_k_stratum_provisional", "gate1_pattern", "gate1_discordant_segment",
            "gate1_implicated_species")}
        self.rows[triplet]["rpi_neutral_factors"] = neutral_factors(row)
        return value


def tournament(scored: List[Tuple[Triplet, float]], rng: random.Random, size: int) -> Triplet:
    contenders = rng.sample(scored, min(size, len(scored)))
    return max(contenders, key=lambda item: item[1])[0]


def crossover(a: Triplet, b: Triplet, rng: random.Random) -> Tuple[Triplet, Triplet]:
    point = rng.choice([1, 2])
    return a[:point] + b[point:], b[:point] + a[point:]


def mutate(triplet: Triplet, pool: Sequence[str], rng: random.Random, rate: float) -> Triplet:
    genes = list(triplet)
    for i in range(3):
        if rng.random() < rate:
            genes[i] = rng.choice(pool)
    return tuple(genes)


def run_search(pool: Sequence[str], score_fn: Callable[[Triplet], float], seed: int,
               population_size: int, generations: int, mutation_rate: float,
               tournament_size: int = 3, elitism: int = 2, plateau: int = 20,
               diversity: bool = False) -> Dict[str, Any]:
    rng = random.Random(seed)
    population = [tuple(rng.choice(pool) for _ in range(3)) for _ in range(population_size)]
    log = []
    best_score = -1.0
    stall = 0
    final_scored: List[Tuple[Triplet, float]] = []
    generations_run = 0
    for generation in range(generations):
        generations_run = generation + 1
        scored = [(t, score_fn(t)) for t in population]
        scored.sort(key=lambda item: item[1], reverse=True)
        final_scored = scored
        values = [v for _, v in scored]
        log.append({"generation": generation, "best": values[0],
                    "mean": statistics.fmean(values), "median": statistics.median(values),
                    "distinct_triplets": len({t for t, _ in scored})})
        if values[0] > best_score + 1e-12:
            best_score = values[0]
            stall = 0
        else:
            stall += 1
        if stall >= plateau:
            break
        children = [t for t, _ in scored[:elitism]]
        while len(children) < population_size:
            parent_a = tournament(scored, rng, tournament_size)
            parent_b = tournament(scored, rng, tournament_size)
            child_a, child_b = crossover(parent_a, parent_b, rng)
            children.append(mutate(child_a, pool, rng, mutation_rate))
            if len(children) < population_size:
                children.append(mutate(child_b, pool, rng, mutation_rate))
        if diversity:
            seen, unique = set(), []
            for t in children:
                if t not in seen:
                    unique.append(t)
                    seen.add(t)
            while len(unique) < population_size:
                t = tuple(rng.choice(pool) for _ in range(3))
                if t not in seen:
                    unique.append(t)
                    seen.add(t)
            children = unique
        population = children
    return {"generations_run": generations_run, "log": log,
            "final_population": final_scored, "stopped_on_plateau": stall >= plateau}


def confirm_winners(oracle: Oracle, winners: List[Dict[str, Any]], output_dir: Path,
                    with_phylogeny: bool, checkpoint_path: Optional[Path] = None) -> str:
    """checkpoint_path, if given, appends each winner's result to a JSONL file as soon
    as it finishes -- not only to the in-memory `winners` list this function already
    mutates. A real gate4 run over hundreds-of-taxa comparator trees can run for hours;
    without this, a crash or kill at hour N loses every winner confirmed before it,
    since the only other record (the caller's final JSON payload) is written once, at
    the very end. Purely additive: omitting checkpoint_path reproduces the exact prior
    behavior (see test_no_checkpoint_path_behaves_exactly_as_before)."""
    if not with_phylogeny:
        return ("deferred: full-gate confirmation not requested "
                "(re-run with --with-phylogeny on a machine with mafft/iqtree)")
    try:
        mafft_bin = ph.resolve_binary(None, "mafft")
        iqtree_bin = ph.resolve_binary(None, "iqtree3")
    except SystemExit as exc:
        return f"failed: could not resolve mafft/iqtree3 -- {exc}"
    args = argparse.Namespace(output=str(output_dir), mafft_bin=mafft_bin, iqtree_bin=iqtree_bin,
                              phylo_min_taxa=4, bootstrap=1000,
                              min_bootstrap_support=95.0, phylo_seed=42, phylo_threads=1)
    error_count = 0
    checkpoint = checkpoint_path.open("a", encoding="utf-8") if checkpoint_path else None
    try:
        for idx, winner in enumerate(winners):
            triplet = tuple(winner["triplet"])
            sources = {"S": oracle.by_key[triplet[0]], "M": oracle.by_key[triplet[1]],
                       "L": oracle.by_key[triplet[2]]}
            try:
                row = rc.score_candidate(oracle.candidate_factory(f"ga-confirm-{idx}", sources),
                                         oracle.reference, oracle.precomputed, args)
                winner["gate4_phylo_status"] = row.get("gate4_phylo_status")
                winner["gate4_reassortment_call"] = row.get("gate4_reassortment_call")
            except Exception as exc:  # binaries missing etc. -- recorded, never hidden
                winner["gate4_error"] = str(exc)
                error_count += 1
            if checkpoint is not None:
                checkpoint.write(json.dumps(winner) + "\n")
                checkpoint.flush()
    finally:
        if checkpoint is not None:
            checkpoint.close()
    if winners and error_count == len(winners):
        return "failed: gate4 errored for every winner (see per-winner gate4_error fields)"
    return "attempted: per-winner gate4 fields recorded (errors included, if any)"


def _check_no_total_scoring_failure(oracle: "Oracle") -> None:
    """Raise if every triplet the search actually scored errored -- a systemic failure
    (bad --reference-input, a factory KeyError on every candidate) must not be reported
    as a successful run with best fitness 0.0, indistinguishable from a genuine flat
    landscape. See §7e Finding 10. Never fires on a healthy run, which always has at
    least one non-errored score."""
    if oracle.cache and len(oracle.errors) == len(oracle.cache):
        raise SystemExit(
            f"GA search aborted: all {len(oracle.cache)} scored triplets errored "
            "(see oracle.errors) -- refusing to report a flat 0.0 landscape as a "
            "successful run")


def dump_scored_table(oracle: "Oracle", output_dir: Path, seed: int) -> Tuple[Path, Path, int]:
    """Write every unique triplet the search actually scored (not just the top --top
    winners) to all_scored_seed_<seed>.jsonl and .csv -- the full RPI row (E, M_pack,
    M_pol, K, rpi_generation, status flags) oracle.score() already computed and cached,
    for candidates that never made the final population. Row order follows
    oracle.cache's insertion order (first-scored first), which is deterministic given a
    fixed seed since the search itself is (see GaSearchTest.test_same_seed_produces_
    identical_run). Does not read or write anything scoring/selection depends on --
    purely an additional view onto oracle.cache/oracle.rows after the search is done."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records: List[Dict[str, Any]] = []
    for triplet, fitness in oracle.cache.items():
        row: Dict[str, Any] = {"triplet_s": triplet[0], "triplet_m": triplet[1],
                               "triplet_l": triplet[2], "fitness": fitness}
        row.update(oracle.rows.get(triplet, {}))
        records.append(row)

    jsonl_path = output_dir / f"all_scored_seed_{seed}.jsonl"
    with jsonl_path.open("w", encoding="utf-8") as f:
        for row in records:
            f.write(json.dumps(row, default=str) + "\n")

    fieldnames: List[str] = []
    seen_fields = set()
    for row in records:
        for key in row:
            if key not in seen_fields:
                seen_fields.add(key)
                fieldnames.append(key)
    csv_path = output_dir / f"all_scored_seed_{seed}.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, restval="")
        writer.writeheader()
        for row in records:
            writer.writerow(row)

    return jsonl_path, csv_path, len(records)


def _resolve_isolate_segments(isolate: Any, reference: Any) -> List[Tuple[str, str, str]]:
    """Yields (segment, accession, sequence) for each of an isolate's S/M/L segments that
    resolves to a genomic sequence via reference.genomic_index. A segment with no
    accession, or an accession absent from genomic_index (e.g. CDS-only coverage), is
    silently omitted here -- callers count the gap themselves, matching confirm_winners'
    "don't let one missing lookup kill the whole run" posture.

    Falls back to the version-stripped base accession, same as every other genomic_index
    lookup in the codebase (rpi_core.py's genomic_sequences_by_segment and three other
    call sites) -- an isolate's s/m/l_accession field can carry a different NCBI version
    suffix than the FASTA record genomic_index indexed it under.
    """
    resolved = []
    for segment, accession_field in (
            ("S", "s_accession"), ("M", "m_accession"), ("L", "l_accession")):
        accession = getattr(isolate, accession_field, "")
        fasta_record = (reference.genomic_index.get(accession)
                        or reference.genomic_index.get(rc.prep.accession_base(accession))
                        ) if accession else None
        if fasta_record is not None:
            resolved.append((segment, accession, fasta_record.sequence))
    return resolved


def _write_fasta(records: List[Tuple[str, str]], output_dir: Path, filename: str) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = output_dir / filename
    with fasta_path.open("w", encoding="utf-8") as f:
        for header, sequence in records:
            f.write(f">{header}\n{sequence}\n")
    return fasta_path


def write_winner_sequences(reference: Any, winners: List[Dict[str, Any]], output_dir: Path,
                           seed: int) -> Tuple[Path, int]:
    """Write each winner's actual S/M/L segment sequences to a single multi-FASTA,
    winner_sequences_seed_<seed>.fasta.

    Winners (and every row in all_scored_seed_<seed>) carry isolate_keys, not sequence
    data -- see dump_scored_table's docstring for why: the same few hundred reference
    isolates get reused across thousands of scored triplets, so inlining sequence data
    into every scored row would mean re-serializing the same segment over and over. This
    is one of two places actual sequence data joins the GA's output (see also
    write_pool_sequences), and it only resolves the (small) winners list.

    Looks up each winner's isolate_key -> accession (via reference.isolates) -> sequence
    (via reference.genomic_index), matching the same two-hop lookup phylogeny.py uses to
    build per-group sequences.fasta.
    """
    by_key = {isolate.isolate_key: isolate for isolate in reference.isolates}
    resolved_cache: Dict[str, Dict[str, Tuple[str, str]]] = {}
    records: List[Tuple[str, str]] = []
    missing = 0
    for idx, winner in enumerate(winners):
        triplet = winner["triplet"]
        for segment, isolate_key in zip(("S", "M", "L"), triplet):
            if isolate_key not in resolved_cache:
                isolate = by_key.get(isolate_key)
                resolved_cache[isolate_key] = {s: (acc, seq) for s, acc, seq in
                    _resolve_isolate_segments(isolate, reference)} if isolate else {}
            resolved = resolved_cache[isolate_key]
            if segment not in resolved:
                missing += 1
                continue
            accession, sequence = resolved[segment]
            species = by_key[isolate_key].species_name
            header = f"winner{idx:02d}_{segment}|{isolate_key}|{accession}|{species}"
            records.append((header, sequence))

    fasta_path = _write_fasta(records, output_dir, f"winner_sequences_seed_{seed}.fasta")
    return fasta_path, missing


def write_pool_sequences(reference: Any, pool: List[str], output_dir: Path,
                         seed: int) -> Tuple[Path, int]:
    """Write every isolate in the search pool's S/M/L segment sequences to a single
    multi-FASTA, pool_sequences_seed_<seed>.fasta -- one record per isolate per segment,
    covering every isolate that could appear in any row of all_scored_seed_<seed>.csv.

    The GA draws every scored triplet from this same fixed pool (a few hundred isolates,
    not the thousands of triplet rows built from it -- see dump_scored_table), so this is
    the full set of sequences any scored candidate, not just the top --top winners, could
    be made of. One record per isolate per resolved segment (not per triplet appearance),
    so unlike winner_sequences_seed_<seed>.fasta a given isolate's sequence appears at
    most once per segment regardless of how many scored triplets used it.
    """
    by_key = {isolate.isolate_key: isolate for isolate in reference.isolates}
    records: List[Tuple[str, str]] = []
    missing = 0
    for isolate_key in pool:
        isolate = by_key.get(isolate_key)
        if isolate is None:
            missing += 3
            continue
        resolved = _resolve_isolate_segments(isolate, reference)
        missing += 3 - len(resolved)
        for segment, accession, sequence in resolved:
            header = f"{isolate_key}_{segment}|{accession}|{isolate.species_name}"
            records.append((header, sequence))

    fasta_path = _write_fasta(records, output_dir, f"pool_sequences_seed_{seed}.fasta")
    return fasta_path, missing


# Fixed declaration order -- the pool file is the contract the GA consumes.
# The reserved-calibration exclusion is applied at load time (same semantics as today).
_POOL_SPECS: List[Tuple[str, str]] = [
    ("terminal_complete", "ga_pool_terminal_complete.csv"),
    ("strict_both", "ga_pool_strict_both.csv"),
]


def load_pool_from_file(pool_file: Path, reserved_keys: set) -> List[str]:
    """Read isolate_keys from a Stage-1 pool CSV and apply reserved-calibration exclusion.

    The pool CSV is the contract -- the GA consumes it rather than recomputing the ruler.
    Reserved-key removal is applied here, matching the semantics of usable_pool
    (which also excludes reserved keys at call time).
    Raises FileNotFoundError if the pool file is missing (Stage 1 must run first).
    """
    if not pool_file.exists():
        raise FileNotFoundError(
            f"GA pool file missing: {pool_file}. Run prepare_dataset.py (Stage 1) first.")
    with pool_file.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    if rows and "isolate_key" not in rows[0]:
        raise KeyError(
            f"pool CSV {pool_file} lacks 'isolate_key' column; "
            f"found columns: {list(rows[0].keys())}")
    all_keys = [row["isolate_key"] for row in rows]
    pool = sorted(k for k in all_keys if k not in reserved_keys)
    if not pool:
        raise SystemExit(f"pool {pool_file.name} is empty after reserved-key exclusion")
    return pool


def _run_one_pool(pool_name: str, pool: List[str], reference: Any, precomputed: Any,
                  pool_out_dir: Path, args: Any, reserved: set,
                  overlap: Dict[str, int]) -> Dict[str, Any]:
    """Run the full GA search for one pool and write per-pool artifacts.
    Returns a comparison-summary dict. All GA hyperparameters are taken from args;
    the pool is the only variable across calls. overlap is pre-computed once by main().
    """
    import time
    pool_out_dir.mkdir(parents=True, exist_ok=True)
    fast_args = rc._score_args(str(pool_out_dir))
    oracle = Oracle(reference, precomputed, fast_args)

    t0 = time.monotonic()
    result = run_search(pool, oracle.score, args.seed, args.population,
                        args.generations, args.mutation_rate, args.tournament_size,
                        args.elitism, args.plateau, args.diversity)
    wall_seconds = time.monotonic() - t0
    _check_no_total_scoring_failure(oracle)

    jsonl_path, csv_path, n_scored = dump_scored_table(oracle, pool_out_dir, args.seed)
    print(f"[{pool_name}] wrote {n_scored} scored triplet(s) to {jsonl_path} and {csv_path}")

    seen, winners = set(), []
    for triplet, value in result["final_population"]:
        if triplet in seen:
            continue
        seen.add(triplet)
        winners.append({"triplet": list(triplet), "fitness": value,
                        **oracle.rows.get(triplet, {})})
        if len(winners) >= args.top:
            break
    confirmation = confirm_winners(oracle, winners, pool_out_dir, args.with_phylogeny)

    fasta_path, missing_sequences = write_winner_sequences(
        reference, winners, pool_out_dir, args.seed)
    print(f"[{pool_name}] wrote {len(winners)} winner(s)' sequences to {fasta_path}"
          + (f" ({missing_sequences} segment(s) had no sequence available)"
             if missing_sequences else ""))

    pool_fasta_path, pool_missing = write_pool_sequences(reference, pool, pool_out_dir, args.seed)
    print(f"[{pool_name}] wrote {len(pool)} pool isolate(s)' sequences to {pool_fasta_path}"
          + (f" ({pool_missing} segment(s) had no sequence available)" if pool_missing else ""))

    payload = {"config": vars(args), "pool_name": pool_name,
               "calibration_reserved_excluded": len(reserved),
               "pool_ruler": pool_name, "pool_size": len(pool),
               "completeness_overlap": overlap,
               "generations_run": result["generations_run"],
               "stopped_on_plateau": result["stopped_on_plateau"],
               "log": result["log"], "winners": winners, "confirmation": confirmation,
               "scoring_errors": [{"triplet": list(t), "error": e} for t, e in oracle.errors]}
    out_path = pool_out_dir / f"seed_{args.seed}.json"
    out_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")

    if result["log"]:
        first, last = result["log"][0], result["log"][-1]
        print(f"[{pool_name}] generations: {result['generations_run']} "
              f"(plateau stop: {result['stopped_on_plateau']})")
        print(f"[{pool_name}] best fitness: first {first['best']:.4f} -> final {last['best']:.4f}")
        print(f"[{pool_name}] winners written: {len(winners)}; confirmation: {confirmation}")
    if oracle.errors:
        print(f"[{pool_name}] scoring errors (scored 0.0): {len(oracle.errors)}", file=sys.stderr)
    print(f"[{pool_name}] wrote {out_path}")

    return {
        "pool_name": pool_name,
        "pool_size": len(pool),
        "generations_run": result["generations_run"],
        "stopped_on_plateau": result["stopped_on_plateau"],
        "best_fitness": result["log"][-1]["best"] if result["log"] else None,
        "winner_fitnesses": [w["fitness"] for w in winners],
        "scoring_error_count": len(oracle.errors),
        "wall_seconds": round(wall_seconds, 1),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference-input", default="output")
    parser.add_argument("--output", default="output/ga_search")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--population", type=int, default=150)
    parser.add_argument("--generations", type=int, default=100)
    parser.add_argument("--mutation-rate", type=float, default=0.02)
    parser.add_argument("--tournament-size", type=int, default=3)
    parser.add_argument("--elitism", type=int, default=2)
    parser.add_argument("--plateau", type=int, default=20)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--diversity", action="store_true")
    parser.add_argument("--with-phylogeny", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    ref_dir = Path(args.reference_input)
    reference = rc.prep.load_reference_data(ref_dir)
    precomputed = rc._build_precomputed(reference)
    reserved = load_reserved_isolate_keys()
    out_root = Path(args.output)
    overlap = completeness_overlap(reference, reserved)

    comparison_rows: List[Dict[str, Any]] = []
    for pool_name, pool_file_name in _POOL_SPECS:
        pool_file = ref_dir / pool_file_name
        pool = load_pool_from_file(pool_file, reserved)
        pool_out_dir = out_root / pool_name
        summary = _run_one_pool(pool_name, pool, reference, precomputed,
                                pool_out_dir, args, reserved, overlap)
        comparison_rows.append(summary)

    # Combined comparison
    comparison = {
        "seed": args.seed,
        "hyperparameters": {
            "population": args.population, "generations": args.generations,
            "mutation_rate": args.mutation_rate, "tournament_size": args.tournament_size,
            "elitism": args.elitism, "plateau": args.plateau,
        },
        "runs": comparison_rows,
    }
    out_root.mkdir(parents=True, exist_ok=True)
    comparison_path = out_root / "ga_pool_comparison.json"
    comparison_path.write_text(json.dumps(comparison, indent=2, default=str), encoding="utf-8")
    print("\n=== GA pool comparison ===")
    for row in comparison_rows:
        best_str = f"{row['best_fitness']:.4f}" if row['best_fitness'] is not None else "N/A"
        print(f"  {row['pool_name']}: pool={row['pool_size']}, "
              f"gens={row['generations_run']} (plateau={row['stopped_on_plateau']}), "
              f"best={best_str}, "
              f"errors={row['scoring_error_count']}, wall={row['wall_seconds']}s")
    print(f"wrote {comparison_path}")
    return 0


class GaSearchTest(unittest.TestCase):
    POOL = [f"iso{i}" for i in range(8)]

    def fake_score(self, triplet: Triplet) -> float:
        return (hash(triplet) % 1000) / 1000.0

    def test_same_seed_produces_identical_run(self):
        a = run_search(self.POOL, self.fake_score, seed=7, population_size=20,
                       generations=15, mutation_rate=0.05)
        b = run_search(self.POOL, self.fake_score, seed=7, population_size=20,
                       generations=15, mutation_rate=0.05)
        self.assertEqual(a, b)

    def test_different_seed_usually_differs(self):
        a = run_search(self.POOL, self.fake_score, seed=7, population_size=20,
                       generations=15, mutation_rate=0.05)
        b = run_search(self.POOL, self.fake_score, seed=8, population_size=20,
                       generations=15, mutation_rate=0.05)
        self.assertNotEqual(a["log"], b["log"])

    def test_plateau_stop_ends_early(self):
        result = run_search(self.POOL, lambda t: 0.5, seed=7, population_size=20,
                            generations=100, mutation_rate=0.05, plateau=5)
        self.assertTrue(result["stopped_on_plateau"])
        self.assertLess(result["generations_run"], 100)

    def test_crossover_swaps_whole_slots_only(self):
        rng = random.Random(1)
        a, b = ("a1", "a2", "a3"), ("b1", "b2", "b3")
        for _ in range(20):
            c1, c2 = crossover(a, b, rng)
            self.assertEqual(len(c1), 3)
            self.assertEqual(len(c2), 3)
            for i in range(3):
                self.assertIn(c1[i], (a[i], b[i]))
                self.assertIn(c2[i], (a[i], b[i]))

    def test_mutation_only_replaces_with_pool_donors(self):
        rng = random.Random(1)
        for _ in range(50):
            mutated = mutate(("a", "b", "c"), self.POOL, rng, rate=0.9)
            self.assertEqual(len(mutated), 3)
            for gene in mutated:
                self.assertIn(gene, list(self.POOL) + ["a", "b", "c"])

    def test_tournament_returns_a_population_member(self):
        rng = random.Random(1)
        scored = [(("a", "b", "c"), 0.1), (("d", "e", "f"), 0.9)]
        winner = tournament(scored, rng, size=2)
        self.assertEqual(winner, ("d", "e", "f"))

    def test_pool_matches_strict_merged_pool(self):
        """The strict_both pool file must contain the same isolate keys as
        rc.strict_merged_pool (before reserved-key exclusion, since both pool files
        are written without exclusion by Stage 1)."""
        strict_file = Path("output/ga_pool_strict_both.csv")
        if not Path("output/phylo_rpi_ga_isolates.json").exists() or not strict_file.exists():
            self.skipTest("output/ snapshot or pool files not present -- run prepare_dataset.py first")
        reference = rc.prep.load_reference_data(Path("output"))
        reserved: set = set()  # no exclusion in pool files
        file_pool = set(load_pool_from_file(strict_file, reserved))
        canonical = set(i.isolate_key for i in rc.strict_merged_pool(reference, reserved))
        self.assertEqual(file_pool, canonical)

    def test_strict_pool_is_subset_of_terminal_pool(self):
        """The strict pool (terminal AND coding-complete) must be a subset of the
        terminal pool (terminal-complete only). Pool files are written by Stage 1."""
        terminal_file = Path("output/ga_pool_terminal_complete.csv")
        strict_file = Path("output/ga_pool_strict_both.csv")
        if not terminal_file.exists() or not strict_file.exists():
            self.skipTest("output/ pool files not present -- run prepare_dataset.py first")
        reserved: set = set()
        terminal_keys = set(load_pool_from_file(terminal_file, reserved))
        strict_keys = set(load_pool_from_file(strict_file, reserved))
        self.assertTrue(strict_keys.issubset(terminal_keys),
                        f"strict_both keys not a subset of terminal_complete: "
                        f"{strict_keys - terminal_keys}")

    def test_completeness_overlap_partitions_the_union(self):
        """both + only_coding + only_terminal must partition coding_complete |
        terminal_complete -- every key in the union appears in exactly one of the
        three sets. Uses fixture isolates so it runs without the output/ snapshot."""
        from types import SimpleNamespace as NS
        coding_only = NS(isolate_key="c1", has_complete_triplet=True)
        terminal_only = NS(isolate_key="t1", has_complete_triplet=False)
        both_ruler = NS(isolate_key="b1", has_complete_triplet=True)
        reference = NS(isolates=[coding_only, terminal_only, both_ruler])
        reserved: set = set()

        with mock.patch.object(
                rc, "wild_classified_pool",
                return_value=[coding_only, terminal_only, both_ruler]), \
             mock.patch.object(
                rc, "usable_pool", return_value=[terminal_only, both_ruler]):
            overlap = completeness_overlap(reference, reserved)

        coding_complete = {"c1", "b1"}
        terminal_complete = {"t1", "b1"}
        union = coding_complete | terminal_complete
        self.assertEqual(overlap["both"] + overlap["only_coding"] + overlap["only_terminal"],
                          len(union))
        self.assertEqual(overlap, {"both": 1, "only_coding": 1, "only_terminal": 1})

    def test_completeness_overlap_ignores_isolates_outside_shared_universe(self):
        """Regression for I4: an isolate with has_complete_triplet=True that is NOT in
        rc.wild_classified_pool's universe (e.g. non-wild, unclassified, or reserved)
        must not inflate only_coding -- it was never eligible for either ruler's pool,
        so it isn't a real disagreement between the two rulers."""
        from types import SimpleNamespace as NS
        eligible_and_coding = NS(isolate_key="e1", has_complete_triplet=True)
        ineligible_but_coding = NS(isolate_key="x1", has_complete_triplet=True)
        reference = NS(isolates=[eligible_and_coding, ineligible_but_coding])
        reserved: set = set()

        with mock.patch.object(
                rc, "wild_classified_pool", return_value=[eligible_and_coding]), \
             mock.patch.object(rc, "usable_pool", return_value=[]):
            overlap = completeness_overlap(reference, reserved)

        self.assertEqual(overlap, {"both": 0, "only_coding": 1, "only_terminal": 0})


class ConfirmWinnersTest(unittest.TestCase):
    """§7e Finding 2: confirm_winners used to hardcode mafft_bin=None, iqtree_bin=None,
    so gate4_check -> run_mafft -> subprocess.run([None, ...]) raised TypeError for every
    winner, while the function still returned the same "attempted" string it returns on
    a real, working run -- a success message over a value that failed (the §8 defect
    shape). These tests pin the real behavior: binaries are actually resolved, and total
    failure to resolve or score is reported as failure, not as "attempted"."""

    def _oracle_and_winners(self):
        from types import SimpleNamespace as NS
        reference = NS(isolates=[NS(isolate_key=k, has_complete_triplet=True,
                                    is_lab_derived=False, is_vaccine_strain=False)
                                 for k in ("a", "b", "c")])
        oracle = Oracle(reference, precomputed={}, args=object(),
                        candidate_factory=lambda cid, sources: {"id": cid, "src": sources})
        winners = [{"triplet": ["a", "b", "c"], "fitness": 0.5}]
        return oracle, winners

    def test_not_requested_defers_without_touching_binaries(self):
        oracle, winners = self._oracle_and_winners()
        with mock.patch.object(ph, "resolve_binary") as spy:
            result = confirm_winners(oracle, winners, Path("/tmp/out"), with_phylogeny=False)
        spy.assert_not_called()
        self.assertIn("deferred", result)

    def test_missing_binaries_reports_failure_not_attempted(self):
        oracle, winners = self._oracle_and_winners()
        with mock.patch.object(ph, "resolve_binary", side_effect=SystemExit("no mafft")):
            result = confirm_winners(oracle, winners, Path("/tmp/out"), with_phylogeny=True)
        self.assertIn("failed", result)
        self.assertNotIn("attempted", result)

    def test_resolved_binaries_are_threaded_into_score_candidate_args(self):
        oracle, winners = self._oracle_and_winners()
        captured = {}

        def fake_resolve(explicit, name):
            return f"/resolved/{name}"

        def fake_score_candidate(candidate, reference, precomputed, args):
            captured["mafft_bin"] = args.mafft_bin
            captured["iqtree_bin"] = args.iqtree_bin
            return {"gate4_phylo_status": "ok", "gate4_reassortment_call": "no_conflict"}

        with mock.patch.object(ph, "resolve_binary", side_effect=fake_resolve), \
             mock.patch.object(rc, "score_candidate", side_effect=fake_score_candidate):
            result = confirm_winners(oracle, winners, Path("/tmp/out"), with_phylogeny=True)

        self.assertEqual(captured["mafft_bin"], "/resolved/mafft")
        self.assertEqual(captured["iqtree_bin"], "/resolved/iqtree3")
        self.assertIn("attempted", result)
        self.assertEqual(winners[0]["gate4_phylo_status"], "ok")

    def test_every_winner_erroring_reports_failure_not_attempted(self):
        oracle, winners = self._oracle_and_winners()
        with mock.patch.object(ph, "resolve_binary", side_effect=lambda e, n: f"/bin/{n}"), \
             mock.patch.object(rc, "score_candidate", side_effect=RuntimeError("boom")):
            result = confirm_winners(oracle, winners, Path("/tmp/out"), with_phylogeny=True)
        self.assertIn("failed", result)
        self.assertNotIn("attempted", result)

    def test_checkpoint_path_writes_one_line_per_winner(self):
        """A multi-hour real gate4 run must not lose completed work to a crash --
        checkpoint_path makes confirm_winners persist each winner's result as it
        finishes, not only in the in-memory list returned at the end."""
        oracle, winners = self._oracle_and_winners()
        winners.append({"triplet": ["a", "b", "c"], "fitness": 0.4})
        with tempfile.TemporaryDirectory() as tmp:
            checkpoint = Path(tmp) / "checkpoint.jsonl"
            with mock.patch.object(ph, "resolve_binary", side_effect=lambda e, n: f"/bin/{n}"), \
                 mock.patch.object(rc, "score_candidate",
                                   return_value={"gate4_phylo_status": "insufficient_comparators",
                                                "gate4_reassortment_call": "not_analyzed"}):
                confirm_winners(oracle, winners, Path("/tmp/out"), with_phylogeny=True,
                                checkpoint_path=checkpoint)
            lines = checkpoint.read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(lines), 2)
        recorded = [json.loads(line) for line in lines]
        self.assertEqual([r["gate4_phylo_status"] for r in recorded],
                         ["insufficient_comparators", "insufficient_comparators"])

    def test_checkpoint_survives_a_hard_crash_mid_run(self):
        """The actual failure mode this exists to prevent, stated directly: if the
        process dies partway (here simulated with KeyboardInterrupt, which the
        per-winner `except Exception` deliberately does not swallow), everything
        completed before the crash must already be on disk."""
        oracle, winners = self._oracle_and_winners()
        winners.append({"triplet": ["a", "b", "c"], "fitness": 0.4})
        with tempfile.TemporaryDirectory() as tmp:
            checkpoint = Path(tmp) / "checkpoint.jsonl"
            with mock.patch.object(ph, "resolve_binary", side_effect=lambda e, n: f"/bin/{n}"), \
                 mock.patch.object(rc, "score_candidate",
                                   side_effect=[{"gate4_phylo_status": "ok",
                                                "gate4_reassortment_call": "no_conflict"},
                                               KeyboardInterrupt()]):
                with self.assertRaises(KeyboardInterrupt):
                    confirm_winners(oracle, winners, Path("/tmp/out"), with_phylogeny=True,
                                    checkpoint_path=checkpoint)
            lines = checkpoint.read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(lines), 1)
        self.assertEqual(json.loads(lines[0])["gate4_phylo_status"], "ok")

    def test_no_checkpoint_path_behaves_exactly_as_before(self):
        oracle, winners = self._oracle_and_winners()
        with mock.patch.object(ph, "resolve_binary", side_effect=lambda e, n: f"/bin/{n}"), \
             mock.patch.object(rc, "score_candidate",
                               return_value={"gate4_phylo_status": "ok",
                                            "gate4_reassortment_call": "no_conflict"}):
            result = confirm_winners(oracle, winners, Path("/tmp/out"), with_phylogeny=True)
        self.assertIn("attempted", result)

    def test_each_winner_gets_a_distinct_candidate_id(self):
        """Finding 9 (2026-09-05 review): every winner reused candidate_id="ga-confirm",
        so with --with-phylogeny all winners wrote gate4 artifacts to the same directory
        -- only the last survived and none could be attributed. Each winner must get its
        own id so per-winner gate4 trees don't collide."""
        oracle, winners = self._oracle_and_winners()
        winners.append({"triplet": ["a", "b", "c"], "fitness": 0.4})
        winners.append({"triplet": ["a", "b", "c"], "fitness": 0.3})
        captured_ids = []

        def fake_candidate_factory(cid, sources):
            captured_ids.append(cid)
            return {"id": cid, "src": sources}

        oracle.candidate_factory = fake_candidate_factory
        with mock.patch.object(ph, "resolve_binary", side_effect=lambda e, n: f"/bin/{n}"), \
             mock.patch.object(rc, "score_candidate",
                               return_value={"gate4_phylo_status": "ok",
                                            "gate4_reassortment_call": "no_conflict"}):
            confirm_winners(oracle, winners, Path("/tmp/out"), with_phylogeny=True)
        self.assertEqual(captured_ids, ["ga-confirm-0", "ga-confirm-1", "ga-confirm-2"])
        self.assertEqual(len(set(captured_ids)), len(captured_ids))


class DumpScoredIsUnconditionalTest(unittest.TestCase):
    """The all-scored dump is unconditional -- there is NO flag governing it (2026-09-26).

    History: it was action="store_true" (opt-in), and README's documented Step 3b command
    did not pass it, so the canonical run wrote only its 10 winners per pool and discarded
    the 1366/1758 triplets it had scored -- the denominator every "the GA found X" claim
    rests on, unrecoverable except by re-running the search. It was briefly made
    default-on with a --no-dump-scored opt-out; the opt-out was then removed outright so
    no invocation can produce a run whose search space went unrecorded.

    These tests pin the ABSENCE of the knob. A future convenience flag would silently
    reintroduce exactly the gap this removal closed.
    """

    def test_no_dump_scored_knob_on_the_parser(self):
        parsed = build_parser().parse_args([])
        self.assertFalse(hasattr(parsed, "dump_scored"),
                         "dump_scored must not be a parsed option -- the dump is unconditional")

    def test_old_opt_in_flag_is_rejected(self):
        with self.assertRaises(SystemExit):
            build_parser().parse_args(["--dump-scored"])

    def test_opt_out_flag_is_rejected(self):
        with self.assertRaises(SystemExit):
            build_parser().parse_args(["--no-dump-scored"])

    def test_dump_call_site_is_not_guarded_by_a_condition(self):
        """Guards the call site itself, not just the parser: re-adding `if args.<flag>:`
        around dump_scored_table would restore the old behaviour while leaving every
        parser-level test above still passing."""
        import inspect
        src = inspect.getsource(_run_one_pool)
        call_line = next(l for l in src.splitlines() if "dump_scored_table(" in l)
        self.assertEqual(call_line.strip().find("if "), -1,
                         f"dump_scored_table call must be unconditional, got: {call_line.strip()}")
        indent = len(call_line) - len(call_line.lstrip())
        self.assertEqual(indent, 4, "call must sit at function-body level, not inside a block")


class DumpScoredTableTest(unittest.TestCase):
    """Every unique scored triplet (not just winners) is written alongside the existing
    seed_<N>.json output. Scoring, selection, and existing outputs must stay unchanged --
    this only adds a read of oracle.cache/oracle.rows after the search is done."""

    def test_dump_file_exists_parses_and_row_count_matches_unique_triplets_scored(self):
        oracle = Oracle.__new__(Oracle)
        oracle.cache = {
            ("a", "b", "c"): 0.5,
            ("d", "e", "f"): 0.0,
            ("g", "h", "i"): 0.2,
        }
        oracle.rows = {
            ("a", "b", "c"): {"verdict": "scored", "rpi_generation": 0.5},
            ("d", "e", "f"): {"verdict": "scoring_error", "scoring_error": "boom"},
            # ("g", "h", "i") deliberately has no rows entry -- must not crash.
        }
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            jsonl_path, csv_path, n = dump_scored_table(oracle, out_dir, seed=99)

            self.assertEqual(n, 3)
            self.assertTrue(jsonl_path.exists())
            self.assertTrue(csv_path.exists())

            lines = jsonl_path.read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(lines), 3)
            parsed = [json.loads(line) for line in lines]
            self.assertEqual(
                [(r["triplet_s"], r["triplet_m"], r["triplet_l"]) for r in parsed],
                [("a", "b", "c"), ("d", "e", "f"), ("g", "h", "i")])
            self.assertEqual(parsed[0]["verdict"], "scored")
            self.assertEqual(parsed[1]["scoring_error"], "boom")
            self.assertEqual(parsed[2]["fitness"], 0.2)

            with csv_path.open(encoding="utf-8", newline="") as f:
                csv_rows = list(csv.DictReader(f))
            self.assertEqual(len(csv_rows), 3)
            self.assertEqual(csv_rows[0]["fitness"], "0.5")


class ResolveIsolateSegmentsAccessionBaseTest(unittest.TestCase):
    """2026-09-30 /code-review finding: an isolate's s/m/l_accession field can carry a
    different NCBI version suffix than the FASTA record genomic_index indexed it under
    -- every other genomic_index lookup in the codebase (rpi_core.py's
    genomic_sequences_by_segment and three other call sites) falls back to the
    version-stripped base accession for exactly this reason. _resolve_isolate_segments
    must do the same, or an isolate that resolves everywhere else in the pipeline
    silently counts as 'missing' here."""

    def test_falls_back_to_version_stripped_accession(self):
        from unittest.mock import Mock
        isolate = Mock(s_accession="NC_077669.2", m_accession="", l_accession="")
        # load_reference_data dual-keys genomic_index on a FASTA record's own accession
        # (.1, the version actually deposited) AND its base -- never on the isolate's
        # accession field's version (.2) directly. Only the base key is shared ground.
        record = Mock(sequence="ACGT")
        genomic_index = {"NC_077669.1": record, "NC_077669": record}
        reference = Mock(genomic_index=genomic_index)
        resolved = _resolve_isolate_segments(isolate, reference)
        self.assertEqual(resolved, [("S", "NC_077669.2", "ACGT")])


class WriteWinnerSequencesTest(unittest.TestCase):
    """Winners carry isolate_keys, not sequence data (see dump_scored_table); this is
    the one place actual FASTA sequences join the GA's output, resolved from the
    reference snapshot's isolates + genomic_index by the same two-hop lookup
    phylogeny.py uses to build sequences.fasta."""

    def _reference(self):
        from unittest.mock import Mock
        isolate_a = Mock(isolate_key="iso-a", species_name="Orthohantavirus puumalaense",
                         s_accession="ACC_A_S", m_accession="ACC_A_M", l_accession="ACC_A_L")
        isolate_b = Mock(isolate_key="iso-b", species_name="Orthohantavirus tulaense",
                         s_accession="ACC_B_S", m_accession="ACC_B_M", l_accession="ACC_B_L")
        genomic_index = {
            "ACC_A_S": Mock(sequence="AAAA"), "ACC_A_M": Mock(sequence="MMMM"),
            "ACC_A_L": Mock(sequence="LLLL"), "ACC_B_S": Mock(sequence="AAAA2"),
            # ACC_B_M and ACC_B_L deliberately absent -- the missing-sequence case.
        }
        return Mock(isolates=[isolate_a, isolate_b], genomic_index=genomic_index)

    def test_writes_one_fasta_record_per_segment_per_winner(self):
        winners = [{"triplet": ["iso-a", "iso-a", "iso-a"], "fitness": 0.5}]
        with tempfile.TemporaryDirectory() as tmp:
            fasta_path, missing = write_winner_sequences(
                self._reference(), winners, Path(tmp), seed=99)
            self.assertEqual(missing, 0)
            text = fasta_path.read_text(encoding="utf-8")
            self.assertEqual(text.count(">"), 3)
            self.assertIn(">winner00_S|iso-a|ACC_A_S|Orthohantavirus puumalaense\nAAAA\n", text)
            self.assertIn(">winner00_M|iso-a|ACC_A_M|Orthohantavirus puumalaense\nMMMM\n", text)
            self.assertIn(">winner00_L|iso-a|ACC_A_L|Orthohantavirus puumalaense\nLLLL\n", text)

    def test_missing_genomic_sequence_is_counted_not_raised(self):
        winners = [{"triplet": ["iso-b", "iso-b", "iso-b"], "fitness": 0.3}]
        with tempfile.TemporaryDirectory() as tmp:
            fasta_path, missing = write_winner_sequences(
                self._reference(), winners, Path(tmp), seed=99)
            # M and L accessions aren't in genomic_index; only S resolves.
            self.assertEqual(missing, 2)
            text = fasta_path.read_text(encoding="utf-8")
            self.assertEqual(text.count(">"), 1)
            self.assertIn(">winner00_S|iso-b|ACC_B_S|Orthohantavirus tulaense\nAAAA2\n", text)

    def test_multiple_winners_get_distinct_indices_in_headers(self):
        winners = [{"triplet": ["iso-a", "iso-a", "iso-a"], "fitness": 0.5},
                  {"triplet": ["iso-a", "iso-a", "iso-a"], "fitness": 0.4}]
        with tempfile.TemporaryDirectory() as tmp:
            fasta_path, _missing = write_winner_sequences(
                self._reference(), winners, Path(tmp), seed=99)
            text = fasta_path.read_text(encoding="utf-8")
            self.assertIn("winner00_S", text)
            self.assertIn("winner01_S", text)


class WritePoolSequencesTest(unittest.TestCase):
    """Covers every isolate the GA search could draw a triplet from, not just the top
    --top winners -- see write_pool_sequences' docstring for why that's cheap (the pool
    is a few hundred isolates, not the thousands of triplet rows built from it)."""

    def _reference(self):
        from unittest.mock import Mock
        isolate_a = Mock(isolate_key="iso-a", species_name="Orthohantavirus puumalaense",
                         s_accession="ACC_A_S", m_accession="ACC_A_M", l_accession="ACC_A_L")
        isolate_b = Mock(isolate_key="iso-b", species_name="Orthohantavirus tulaense",
                         s_accession="ACC_B_S", m_accession="ACC_B_M", l_accession="ACC_B_L")
        genomic_index = {
            "ACC_A_S": Mock(sequence="AAAA"), "ACC_A_M": Mock(sequence="MMMM"),
            "ACC_A_L": Mock(sequence="LLLL"), "ACC_B_S": Mock(sequence="AAAA2"),
            # ACC_B_M and ACC_B_L deliberately absent -- the missing-sequence case.
        }
        return Mock(isolates=[isolate_a, isolate_b], genomic_index=genomic_index)

    def test_writes_one_record_per_isolate_per_resolved_segment(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta_path, missing = write_pool_sequences(
                self._reference(), ["iso-a", "iso-b"], Path(tmp), seed=99)
            # iso-a: S, M, L all resolve. iso-b: only S resolves (M, L missing).
            self.assertEqual(missing, 2)
            text = fasta_path.read_text(encoding="utf-8")
            self.assertEqual(text.count(">"), 4)
            self.assertIn(">iso-a_S|ACC_A_S|Orthohantavirus puumalaense\nAAAA\n", text)
            self.assertIn(">iso-a_M|ACC_A_M|Orthohantavirus puumalaense\nMMMM\n", text)
            self.assertIn(">iso-a_L|ACC_A_L|Orthohantavirus puumalaense\nLLLL\n", text)
            self.assertIn(">iso-b_S|ACC_B_S|Orthohantavirus tulaense\nAAAA2\n", text)

    def test_isolate_appears_once_regardless_of_how_many_triplets_used_it(self):
        """Unlike winner_sequences (one record per winner appearance), pool_sequences
        writes one record per isolate per segment -- the pool list itself has no
        duplicates by construction (load_pool_from_file sorts unique isolate_keys), so
        this just confirms there's no separate per-appearance multiplication anywhere."""
        with tempfile.TemporaryDirectory() as tmp:
            fasta_path, _missing = write_pool_sequences(
                self._reference(), ["iso-a"], Path(tmp), seed=99)
            text = fasta_path.read_text(encoding="utf-8")
            self.assertEqual(text.count("iso-a_S"), 1)

    def test_isolate_key_not_in_reference_counts_all_three_segments_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta_path, missing = write_pool_sequences(
                self._reference(), ["iso-ghost"], Path(tmp), seed=99)
            self.assertEqual(missing, 3)
            self.assertEqual(fasta_path.read_text(encoding="utf-8"), "")


class TotalScoringFailureGuardTest(unittest.TestCase):
    """§7e Finding 10: a systemic scoring failure (bad --reference-input, a factory
    KeyError on every candidate) must abort instead of silently reporting a
    successful run with best fitness 0.0 -- indistinguishable from a real flat
    landscape otherwise."""

    def _oracle(self, cache, errors):
        oracle = Oracle.__new__(Oracle)
        oracle.cache = cache
        oracle.errors = errors
        return oracle

    def test_raises_when_every_scored_triplet_errored(self):
        oracle = self._oracle(
            cache={("a", "b", "c"): 0.0, ("d", "e", "f"): 0.0},
            errors=[(("a", "b", "c"), "KeyError"), (("d", "e", "f"), "KeyError")])
        with self.assertRaises(SystemExit):
            _check_no_total_scoring_failure(oracle)

    def test_does_not_raise_with_partial_errors(self):
        oracle = self._oracle(
            cache={("a", "b", "c"): 0.5, ("d", "e", "f"): 0.0},
            errors=[(("d", "e", "f"), "KeyError")])
        _check_no_total_scoring_failure(oracle)  # must not raise

    def test_does_not_raise_with_no_errors(self):
        oracle = self._oracle(cache={("a", "b", "c"): 0.5}, errors=[])
        _check_no_total_scoring_failure(oracle)  # must not raise

    def test_does_not_raise_when_cache_empty(self):
        oracle = self._oracle(cache={}, errors=[])
        _check_no_total_scoring_failure(oracle)  # must not raise


class OracleTest(unittest.TestCase):
    def _oracle(self):
        from types import SimpleNamespace as NS
        reference = NS(isolates=[NS(isolate_key=k, has_complete_triplet=True,
                                    is_lab_derived=False, is_vaccine_strain=False)
                                 for k in ("a", "b", "c")])
        oracle = Oracle(reference, precomputed={}, args=object(),
                        candidate_factory=lambda cid, sources: {"id": cid, "src": sources})
        return oracle

    def test_memoization_scores_each_triplet_once(self):
        oracle = self._oracle()
        row = {"verdict": "noconflict", "rpi_generation": 0.42}
        with mock.patch.object(rc, "score_candidate", return_value=row) as spy:
            oracle.score(("a", "b", "c"))
            oracle.score(("a", "b", "c"))
        self.assertEqual(spy.call_count, 1)
        self.assertEqual(oracle.cache[("a", "b", "c")], 0.42)

    def test_junk_verdict_maps_to_zero(self):
        oracle = self._oracle()
        row = {"verdict": "junk_sequence", "rpi_generation": None}
        with mock.patch.object(rc, "score_candidate", return_value=row):
            self.assertEqual(oracle.score(("a", "b", "c")), 0.0)

    def test_rpi_e_unknown_components_reaches_oracle_rows(self):
        """rpi_core.py's rpi_e_unknown_components diagnostic (which winner's Factor E
        hit a missing-metadata sentinel) must survive Oracle.score()'s field whitelist,
        or it never reaches all_scored_seed_<N>.csv, seed_<N>.json's winners, or
        figures.py's ga_winners_parents.csv -- the only places a reader could see it."""
        oracle = self._oracle()
        row = {"verdict": "no_conflict", "rpi_generation": 0.47, "rpi_e_unknown_components": "host,geo"}
        with mock.patch.object(rc, "score_candidate", return_value=row):
            oracle.score(("a", "b", "c"))
        self.assertEqual(oracle.rows[("a", "b", "c")]["rpi_e_unknown_components"], "host,geo")

    def test_factor_statuses_and_neutral_summary_reach_oracle_rows(self):
        """A neutral 1.0 from a missing-data fallback must be distinguishable from a
        measured value in the GA outputs, for every factor. The statuses rpi_core.py
        already computes (M_pack, M_pol) and K's provisional marker must survive the
        field whitelist, and rpi_neutral_factors summarises which factors are not
        measured. Purely descriptive: the returned fitness must be unchanged."""
        oracle = self._oracle()
        row = {"verdict": "no_conflict", "rpi_generation": 0.47, "rpi_m_pack": 1.0,
               "rpi_m_pack_status": "no_parental_sequence", "rpi_m_pol": 1.0,
               "rpi_m_pol_status": "no_signal_beyond_tag", "rpi_k_stratum": "close",
               "rpi_k_stratum_provisional": True, "rpi_e_unknown_components": "host,geo"}
        with mock.patch.object(rc, "score_candidate", return_value=row):
            value = oracle.score(("a", "b", "c"))
        self.assertEqual(value, 0.47)
        stored = oracle.rows[("a", "b", "c")]
        self.assertEqual(stored["rpi_m_pack_status"], "no_parental_sequence")
        self.assertEqual(stored["rpi_m_pol_status"], "no_signal_beyond_tag")
        self.assertIs(stored["rpi_k_stratum_provisional"], True)
        self.assertEqual(stored["rpi_neutral_factors"],
                         "E:host+geo;M_pack:no_parental_sequence;M_pol:no_signal_beyond_tag")

    def test_neutral_factors_is_empty_when_every_factor_was_measured(self):
        self.assertEqual(neutral_factors({
            "rpi_e_unknown_components": "", "rpi_m_pack_status": "scored",
            "rpi_m_pol_status": "scored"}), "")

    def test_neutral_factors_does_not_count_a_scored_one_as_neutral(self):
        """A measured M_pack of exactly 1.0 (saturated, status 'scored') is a real value."""
        self.assertEqual(neutral_factors({
            "rpi_m_pack": 1.0, "rpi_m_pack_status": "scored",
            "rpi_m_pol_status": "no_signal_beyond_tag"}), "M_pol:no_signal_beyond_tag")

    def test_score_candidate_exception_scores_zero_and_is_recorded_not_raised(self):
        """An unattended multi-hour search must survive one bad triplet -- see errors
        docstring on Oracle.__init__. This is the behavior that makes that possible."""
        oracle = self._oracle()
        with mock.patch.object(rc, "score_candidate", side_effect=ValueError("boom")):
            value = oracle.score(("a", "b", "c"))
        self.assertEqual(value, 0.0)
        self.assertEqual(oracle.cache[("a", "b", "c")], 0.0)
        self.assertEqual(len(oracle.errors), 1)
        self.assertEqual(oracle.errors[0][0], ("a", "b", "c"))
        self.assertIn("boom", oracle.errors[0][1])
        self.assertEqual(oracle.rows[("a", "b", "c")]["verdict"], "scoring_error")

    def test_score_candidate_exception_is_memoized_like_any_other_score(self):
        oracle = self._oracle()
        with mock.patch.object(rc, "score_candidate", side_effect=ValueError("boom")) as spy:
            oracle.score(("a", "b", "c"))
            oracle.score(("a", "b", "c"))
        self.assertEqual(spy.call_count, 1)
        self.assertEqual(len(oracle.errors), 1)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        sys.argv = [sys.argv[0]]
        unittest.main()
    else:
        raise SystemExit(main())
