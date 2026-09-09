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


def build_pool(reference: Any, reserved_keys: set) -> List[str]:
    """Use the canonical usable_pool from rpi_core rather than a reimplementation.

    usable_pool's docstring explicitly forbids a second implementation: it filters on
    terminally_complete_isolates (all three segments full-length AND both panhandle tags)
    to protect panhandle-adjacent computation like M_pack. The old build_pool only
    checked has_complete_triplet and omitted terminal-completeness, contaminating 51% of
    what it searched and missing 65% of the real pool.
    """
    pool_records = rc.usable_pool(reference, reserved_keys)
    pool = sorted([isolate_key_of(i) for i in pool_records])
    if not pool:
        raise SystemExit("search pool is empty after filters -- nothing to sample from")
    return pool


def default_candidate_factory(reference: Any) -> Callable[[str, Dict[str, Any]], Any]:
    factory = getattr(rc, "_candidate_from_isolates", None)
    if factory is None:
        raise SystemExit(
            "rpi_core._candidate_from_isolates not found; grep rpi_core.py for "
            "'def candidate' and wire the correct factory name here")
    return lambda candidate_id, sources: factory(candidate_id, reference, sources)


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
            "rpi_e_isolate_basis", "rpi_m_pack", "rpi_m_pol", "rpi_k_stratum",
            "gate1_pattern", "gate1_discordant_segment", "gate1_implicated_species")}
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


def main(argv: Optional[Sequence[str]] = None) -> int:
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
    parser.add_argument("--dump-scored", action="store_true",
                        help="write every unique triplet the search scored (not just "
                             "the top --top winners) to "
                             "all_scored_seed_<seed>.{jsonl,csv} in --output")
    args = parser.parse_args(argv)

    reference = rc.prep.load_reference_data(Path(args.reference_input))
    precomputed = rc._build_precomputed(reference)
    fast_args = rc._score_args(args.output)  # positional: keyword name is internal
    oracle = Oracle(reference, precomputed, fast_args)
    reserved = load_reserved_isolate_keys()
    pool = build_pool(reference, reserved)

    result = run_search(pool, oracle.score, args.seed, args.population,
                        args.generations, args.mutation_rate, args.tournament_size,
                        args.elitism, args.plateau, args.diversity)
    _check_no_total_scoring_failure(oracle)

    if args.dump_scored:
        jsonl_path, csv_path, n_scored = dump_scored_table(oracle, args.output, args.seed)
        print(f"wrote {n_scored} scored triplet(s) to {jsonl_path} and {csv_path}")

    seen, winners = set(), []
    for triplet, value in result["final_population"]:
        if triplet in seen:
            continue
        seen.add(triplet)
        winners.append({"triplet": list(triplet), "fitness": value,
                        **oracle.rows.get(triplet, {})})
        if len(winners) >= args.top:
            break
    confirmation = confirm_winners(oracle, winners, Path(args.output), args.with_phylogeny)

    payload = {"config": vars(args), "calibration_reserved_excluded": len(reserved),
               "pool_size": len(pool), "generations_run": result["generations_run"],
               "stopped_on_plateau": result["stopped_on_plateau"],
               "log": result["log"], "winners": winners, "confirmation": confirmation,
               "scoring_errors": [{"triplet": list(t), "error": e} for t, e in oracle.errors]}
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"seed_{args.seed}.json"
    out_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    if result["log"]:
        first, last = result["log"][0], result["log"][-1]
        print(f"generations: {result['generations_run']} "
              f"(plateau stop: {result['stopped_on_plateau']})")
        print(f"best fitness: first {first['best']:.4f} -> final {last['best']:.4f}")
        print(f"winners written: {len(winners)}; confirmation: {confirmation}")
    if oracle.errors:
        print(f"scoring errors (scored 0.0, run continued): {len(oracle.errors)}",
              file=sys.stderr)
    print(f"wrote {out_path}")
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

    def test_pool_matches_usable_pool(self):
        """build_pool now calls the canonical rc.usable_pool, so test that it produces
        the same result. This guards against accidental divergence."""
        from pathlib import Path
        if not Path("output/ortho_reference_isolates.json").exists():
            self.skipTest("output/ reference snapshot not present -- run prepare_dataset.py first")
        reference = rc.prep.load_reference_data(Path("output"))
        reserved = rc.reserved_calibration_isolate_keys()
        # Both should produce the same pool
        alg_pool = set(build_pool(reference, reserved))
        canonical = set(i.isolate_key for i in rc.usable_pool(reference, reserved))
        self.assertEqual(alg_pool, canonical)


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


class DumpScoredTableTest(unittest.TestCase):
    """--dump-scored writes every unique scored triplet (not just winners) alongside
    the existing seed_<N>.json output. Scoring, selection, and existing outputs must
    stay unchanged -- this only adds a read of oracle.cache/oracle.rows after the
    search is done."""

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
