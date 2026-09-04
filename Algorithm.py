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
import json
import random
import statistics
import unittest
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple
from unittest import mock

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
    pool = []
    for isolate in reference.isolates:
        if not _field(isolate, "has_complete_triplet", "hascompletetriplet"):
            continue
        if _field(isolate, "is_lab_derived", "islabderived"):
            continue
        if _field(isolate, "is_vaccine_strain", "isvaccinestrain"):
            continue
        key = isolate_key_of(isolate)
        if key in reserved_keys:
            continue
        pool.append(key)
    if not pool:
        raise SystemExit("search pool is empty after filters -- nothing to sample from")
    return sorted(pool)


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

    def score(self, triplet: Triplet) -> float:
        if triplet in self.cache:
            return self.cache[triplet]
        sources = {"S": self.by_key[triplet[0]], "M": self.by_key[triplet[1]],
                   "L": self.by_key[triplet[2]]}
        candidate = self.candidate_factory("ga", sources)
        row = rc.score_candidate(candidate, self.reference, self.precomputed, self.args)
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
                    with_phylogeny: bool) -> str:
    if not with_phylogeny:
        return ("deferred: full-gate confirmation not requested "
                "(re-run with --with-phylogeny on a machine with mafft/iqtree)")
    args = argparse.Namespace(output=str(output_dir), mafft_bin=None, iqtree_bin=None,
                              phylo_min_taxa=4, bootstrap=1000,
                              min_bootstrap_support=95.0, phylo_seed=42, phylo_threads=1)
    for winner in winners:
        triplet = tuple(winner["triplet"])
        sources = {"S": oracle.by_key[triplet[0]], "M": oracle.by_key[triplet[1]],
                   "L": oracle.by_key[triplet[2]]}
        try:
            row = rc.score_candidate(oracle.candidate_factory("ga-confirm", sources),
                                     oracle.reference, oracle.precomputed, args)
            winner["gate4_phylo_status"] = row.get("gate4_phylo_status")
            winner["gate4_reassortment_call"] = row.get("gate4_reassortment_call")
        except Exception as exc:  # binaries missing etc. -- recorded, never hidden
            winner["gate4_error"] = str(exc)
    return "attempted: per-winner gate4 fields recorded (errors included, if any)"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference-input", default="all_hanta")
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
    args = parser.parse_args(argv)

    reference = rc.load_reference_data(Path(args.reference_input))
    precomputed = rc._build_precomputed(reference)
    fast_args = rc._score_args(args.output)  # positional: keyword name is internal
    oracle = Oracle(reference, precomputed, fast_args)
    reserved = load_reserved_isolate_keys()
    pool = build_pool(reference, reserved)

    result = run_search(pool, oracle.score, args.seed, args.population,
                        args.generations, args.mutation_rate, args.tournament_size,
                        args.elitism, args.plateau, args.diversity)

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
               "log": result["log"], "winners": winners, "confirmation": confirmation}
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

    def test_pool_excludes_reserved_lab_and_incomplete(self):
        from types import SimpleNamespace as NS
        reference = NS(isolates=[
            NS(isolate_key="keep1", has_complete_triplet=True, is_lab_derived=False,
               is_vaccine_strain=False),
            NS(isolate_key="drop_reserved", has_complete_triplet=True,
               is_lab_derived=False, is_vaccine_strain=False),
            NS(isolate_key="drop_lab", has_complete_triplet=True, is_lab_derived=True,
               is_vaccine_strain=False),
            NS(isolate_key="drop_incomplete", has_complete_triplet=False,
               is_lab_derived=False, is_vaccine_strain=False)])
        pool = build_pool(reference, {"drop_reserved"})
        self.assertEqual(pool, ["keep1"])


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


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        sys.argv = [sys.argv[0]]
        unittest.main()
    else:
        raise SystemExit(main())
