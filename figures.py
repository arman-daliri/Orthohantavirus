"""Before/after phylogenetic-tree figures over Stage 2's COMMITTED outputs only.

No new scoring, no new gate calls, no new sequences. Every tree drawn here is either
read directly from a `tree.treefile` already committed under `output/phylogenetics*/`
(Stage 2's own MAFFT+IQ-TREE output) or is that same file pruned to a local neighborhood
for legibility -- never rebuilt, realigned, or reinferred.

The one piece of real computation this module does is reconstructing the synthetic tip
label -> isolate_key mapping (`I0000`, `I0001`, ... -> `isolate_key`) that
`phylogeny._analyze_group` assigned when it built the committed trees, since that mapping
is not itself persisted to disk. This is done by calling phylogeny.py's own
`collect_phylo_candidates` / `group_candidates` (deterministic, sorted by isolate_key --
see `group_candidates`'s own docstring) with the same inputs and defaults Stage 2's
committed run used. This re-derives a label ASSIGNMENT, not a tree -- no mafft/iqtree
subprocess runs anywhere in this module.

BEFORE panel: the 8 `supported_reassortment_candidate` isolates from
`output/ortho_reassortment_table.csv` (prepare_dataset.py's own genome-wide phylogenetic
surveillance). For each, render its S/M/L segment trees' LOCAL NEIGHBORHOOD (the
smallest ancestor clade with at least `--min-tips` taxa, since several species trees
here run into the hundreds of taxa and a full tree would be illegible), highlighting
that isolate's tip and connecting it across the three panels.

AFTER panel: the GA search's top-fitness winners (see `--ga-search-json`), each a
triplet of existing reference isolates assigned to S/M/L slots. Renders each slot
donor's own segment tree local neighborhood, highlighted in a color distinct per slot,
plus a caption stating this is an in-silico combination of existing isolates -- no new
sequence is created. A companion CSV records candidate -> (S, M, L) parent isolate,
species, country, host.

Run:
    python figures.py before --output output/figures/before
    python figures.py after  --ga-search-json output/ga_search_corrected_pool_long/seed_20260904.json \\
                              --output output/figures/after
    python -m unittest figures
"""

import argparse
import copy
import csv
import json
import logging
import unittest
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import prepare_dataset as prep
import phylogeny as ph

SEGMENTS = ("S", "M", "L")
SLOT_COLOR = {"S": "#1b6ca8", "M": "#2e8b57", "L": "#c1440e"}
DEFAULT_MIN_TIPS = 8


# ============================================================================
# Ingest + label-map reconstruction (re-derives phylogeny.py's own deterministic
# grouping; never rebuilds a tree)
# ============================================================================

def load_ingest(
    input_dir: str = ".", cluster_threshold: float = 0.02, min_species_size: int = 8,
) -> Tuple[List["prep.IsolateRecord"], Dict[str, "prep.SegmentDivergenceProfile"],
          Dict[str, List[Tuple[str, str]]], "prep.InputPaths"]:
    """Mirrors phylogeny.run_phylogeny's ingest chain through cds_index -- the same
    ~15-line duplication the three-stage split spec already accepts (§3.2), here for a
    fourth reason: figures need the same isolate population Stage 2 scored, not a
    filtered output table.
    """
    paths = prep.resolve_inputs(Path(input_dir))
    data_records = list(prep.iter_jsonl(paths.data_report))
    annotation_index = prep.build_annotation_index(prep.iter_jsonl(paths.annotation_report))
    biosample_index = prep.build_biosample_index(prep.iter_jsonl(paths.biosample_report))
    genomic_index = prep.parse_fasta(paths.genomic_fasta)
    protein_index = prep.parse_fasta(paths.protein_fasta)

    accessions = []
    for record in data_records:
        transformed, _qc_row = prep.transform_one_record(
            record, annotation_index, biosample_index, genomic_index, protein_index)
        if transformed is not None:
            accessions.append(transformed)

    isolates = prep.build_isolate_records(accessions)
    prep.assign_kmer_clusters(isolates, cluster_threshold)
    divergence = prep.compute_segment_divergence(isolates, min_species_size=min_species_size)
    divergence_by_key = {d.isolate_key: d for d in divergence}
    cds_index = ph.parse_cds_fasta(paths.cds_fasta)
    return isolates, divergence_by_key, cds_index, paths


def group_for_scope(
    isolates: Sequence["prep.IsolateRecord"], divergence_by_key: Dict[str, Any],
    cds_index: Dict[str, Any], scope: str, phylo_min_taxa: int = 4,
) -> Dict[str, List["ph.PhyloCandidate"]]:
    """by_group exactly as phylogeny._run_phylogenetics_scope built it: same candidate
    collection, same min-taxa filter, members sorted by isolate_key (group_candidates'
    own ordering guarantee -- this is what makes label reconstruction deterministic)."""
    candidates, _ = ph.collect_phylo_candidates(isolates, divergence_by_key, cds_index, scope)
    by_group, _ = ph.group_candidates(candidates, phylo_min_taxa, scope)
    return by_group


def label_map_for_members(members: Sequence["ph.PhyloCandidate"]) -> Dict[str, "ph.PhyloCandidate"]:
    """The exact label assignment phylogeny._analyze_group makes: I0000, I0001, ... in
    the members' given order. Members must already be in the same order _analyze_group
    received them (group_candidates' sorted-by-isolate_key output) for this to line up
    with a committed tree.treefile."""
    labels = [f"I{i:04d}" for i in range(len(members))]
    return dict(zip(labels, members))


def find_group_for_isolate(
    by_group: Dict[str, List["ph.PhyloCandidate"]], isolate_key: str,
) -> Tuple[Optional[str], Optional[List["ph.PhyloCandidate"]]]:
    for group_key, members in by_group.items():
        if any(m.isolate.isolate_key == isolate_key for m in members):
            return group_key, members
    return None, None


def load_tree_relabeled(group_dir_root: Path, group_key: str, segment: str,
                        label_map: Dict[str, "ph.PhyloCandidate"]) -> "Bio.Phylo.BaseTree.Tree":
    from Bio import Phylo
    path = group_dir_root / prep.normalize_attr_name(group_key) / segment / "tree.treefile"
    if not path.exists():
        raise FileNotFoundError(f"committed tree not found: {path}")
    tree = Phylo.read(str(path), "newick")
    for tip in tree.get_terminals():
        candidate = label_map.get(tip.name)
        if candidate is not None:
            tip.name = candidate.isolate.isolate_key
    return tree


# ============================================================================
# Local-neighborhood pruning (legibility only -- never touches the underlying tree file)
# ============================================================================

def local_neighborhood(tree: "Bio.Phylo.BaseTree.Tree", target_name: str,
                       min_tips: int = DEFAULT_MIN_TIPS) -> "Bio.Phylo.BaseTree.Clade":
    """Smallest ancestor clade of `target_name` with at least `min_tips` terminals, or
    the whole tree if it never reaches that size. A deep copy so callers can mutate
    (ladderize) the returned clade without touching the cached full tree."""
    path = tree.get_path(target_name)
    if not path:
        return copy.deepcopy(tree.root)
    for ancestor in reversed(path[:-1]):
        if ancestor.count_terminals() >= min_tips:
            return copy.deepcopy(ancestor)
    return copy.deepcopy(tree.root)


# ============================================================================
# Minimal cladogram renderer (figure-drawing code -- not phylogenetic inference; the
# instruction to never hand-roll a tree refers to alignment/tree-BUILDING, which this
# module never does)
# ============================================================================

def _layout(clade: "Bio.Phylo.BaseTree.Clade") -> Tuple[Dict[int, float], List["Bio.Phylo.BaseTree.Clade"]]:
    clade.ladderize()
    tips = clade.get_terminals()
    y: Dict[int, float] = {id(t): i for i, t in enumerate(tips)}

    def assign(c):
        if c.is_terminal():
            return y[id(c)]
        ys = [assign(child) for child in c.clades]
        y[id(c)] = sum(ys) / len(ys)
        return y[id(c)]

    assign(clade)
    return y, tips


def _x_positions(clade: "Bio.Phylo.BaseTree.Clade") -> Dict[int, float]:
    x: Dict[int, float] = {}

    def assign(c, parent_x):
        branch_length = c.branch_length if c.branch_length is not None else 0.0
        cx = parent_x + branch_length
        x[id(c)] = cx
        for child in c.clades:
            assign(child, cx)

    assign(clade, 0.0)
    return x


def draw_cladogram(ax, clade: "Bio.Phylo.BaseTree.Clade", highlight_name: Optional[str] = None,
                   highlight_color: str = "crimson") -> Dict[str, Tuple[float, float]]:
    """Draws a rectangular cladogram of `clade` onto `ax`. Returns {tip_name: (x, y)} in
    axes data coordinates, so callers can draw cross-panel connecting lines."""
    y, tips = _layout(clade)
    x = _x_positions(clade)

    def draw(c):
        if not c.is_terminal():
            child_ys = [y[id(child)] for child in c.clades]
            ax.plot([x[id(c)], x[id(c)]], [min(child_ys), max(child_ys)], color="#555555", lw=0.8)
            for child in c.clades:
                ax.plot([x[id(c)], x[id(child)]], [y[id(child)], y[id(child)]], color="#555555", lw=0.8)
                draw(child)

    draw(clade)

    tip_pos: Dict[str, Tuple[float, float]] = {}
    for tip in tips:
        cx, cy = x[id(tip)], y[id(tip)]
        tip_pos[tip.name] = (cx, cy)
        is_hl = tip.name == highlight_name
        ax.plot([cx], [cy], marker="o", markersize=6 if is_hl else 2,
                color=highlight_color if is_hl else "#9aa5b1", zorder=5)
        if is_hl:
            ax.text(cx, cy + 0.4, tip.name, fontsize=7, fontweight="bold",
                    color=highlight_color, va="bottom", ha="left")

    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.set_xlabel("substitutions/site (local neighborhood)")
    return tip_pos


# ============================================================================
# BEFORE panel
# ============================================================================

def supported_reassortment_rows(reassortment_table: Path) -> List[Dict[str, str]]:
    with reassortment_table.open(encoding="utf-8") as handle:
        rows = [r for r in csv.DictReader(handle)
                if r.get("reassortment_call") == "supported_reassortment_candidate"]
    # One row per isolate: species-scope preferred over phylogroup-scope when both
    # exist for the same isolate (single-species tree is the simpler, more directly
    # interpretable figure; phylogroup is the fallback for isolates species-scope
    # structurally excludes, e.g. unclassified taxa).
    by_isolate: Dict[str, Dict[str, str]] = {}
    for row in rows:
        key = row["isolate_key"]
        existing = by_isolate.get(key)
        if existing is None or (existing.get("analysis_scope") == "phylogroup"
                                and row.get("analysis_scope") == "species"):
            by_isolate[key] = row
    return sorted(by_isolate.values(), key=lambda r: r["isolate_key"])


def render_before_figures(output_dir: Path, input_dir: str = ".",
                          reassortment_table: Optional[Path] = None,
                          min_tips: int = DEFAULT_MIN_TIPS) -> List[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import ConnectionPatch

    isolates, divergence_by_key, cds_index, paths = load_ingest(input_dir)
    by_group = {
        "species": group_for_scope(isolates, divergence_by_key, cds_index, "species"),
        "phylogroup": group_for_scope(isolates, divergence_by_key, cds_index, "phylogroup"),
    }
    group_dir_root = {
        "species": paths.output_dir / "phylogenetics",
        "phylogroup": paths.output_dir / "phylogenetics_phylogroup",
    }

    table_path = reassortment_table or (paths.output_dir / "ortho_reassortment_table.csv")
    rows = supported_reassortment_rows(table_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for row in rows:
        isolate_key = row["isolate_key"]
        scope = row["analysis_scope"]
        discordant = row.get("discordant_segment") or ""
        fig, axes = plt.subplots(1, 3, figsize=(13, 5))
        tip_pos: Dict[str, Tuple[float, float]] = {}
        for ax, segment in zip(axes, SEGMENTS):
            group_key, members = find_group_for_isolate(by_group[scope], isolate_key)
            if group_key is None:
                ax.set_axis_off()
                ax.text(0.5, 0.5, f"{segment}: not in a\n{scope}-scope tree",
                        ha="center", va="center", wrap=True, transform=ax.transAxes)
                continue
            label_map = label_map_for_members(members)
            tree = load_tree_relabeled(group_dir_root[scope], group_key, segment, label_map)
            clade = local_neighborhood(tree, isolate_key, min_tips=min_tips)
            color = "#d4380d" if segment == discordant else "#1b6ca8"
            positions = draw_cladogram(ax, clade, highlight_name=isolate_key, highlight_color=color)
            title = segment + (" -- discordant" if segment == discordant else "")
            ax.set_title(title, color=color if segment == discordant else "black",
                        fontweight="bold" if segment == discordant else "normal")
            if isolate_key in positions:
                tip_pos[segment] = positions[isolate_key]

        for seg_a, seg_b in zip(SEGMENTS, SEGMENTS[1:]):
            if seg_a in tip_pos and seg_b in tip_pos:
                idx_a, idx_b = SEGMENTS.index(seg_a), SEGMENTS.index(seg_b)
                con = ConnectionPatch(
                    xyA=tip_pos[seg_a], coordsA=axes[idx_a].transData,
                    xyB=tip_pos[seg_b], coordsB=axes[idx_b].transData,
                    color="#d4380d" if discordant in (seg_a, seg_b) else "#9aa5b1",
                    lw=1.2, linestyle="--", alpha=0.8)
                fig.add_artist(con)

        fig.suptitle(f"{isolate_key}  ({row['species_name']}) -- IQ-TREE-supported "
                    f"topology conflict, discordant segment: {discordant or 'n/a'}",
                    fontsize=10)
        fig.text(0.5, 0.01,
                 "Real, unmodified S/M/L sequences from a single deposited isolate; "
                 "trees pruned to a local neighborhood for legibility (full tree in "
                 f"output/phylogenetics{'_phylogroup' if scope == 'phylogroup' else ''}/).",
                 ha="center", fontsize=7, color="#555555")
        fig.tight_layout(rect=(0, 0.03, 1, 1))
        out_path = output_dir / f"{prep.normalize_attr_name(isolate_key)}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        written.append(out_path)
    return written


# ============================================================================
# AFTER panel (GA winners)
# ============================================================================

def parent_metadata(isolate_key: str, by_key: Dict[str, "prep.IsolateRecord"]) -> Dict[str, str]:
    isolate = by_key.get(isolate_key)
    if isolate is None:
        return {"isolate_key": isolate_key, "species": "", "country": "", "host": ""}
    return {"isolate_key": isolate_key, "species": isolate.species_name,
            "country": isolate.country, "host": isolate.host_name}


def render_after_figures(ga_search_json: Path, output_dir: Path, input_dir: str = ".",
                         min_tips: int = DEFAULT_MIN_TIPS) -> Tuple[List[Path], Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    payload = json.loads(ga_search_json.read_text(encoding="utf-8"))
    winners = payload["winners"]

    isolates, divergence_by_key, cds_index, paths = load_ingest(input_dir)
    by_key = {i.isolate_key: i for i in isolates}
    by_group = {
        "species": group_for_scope(isolates, divergence_by_key, cds_index, "species"),
        "phylogroup": group_for_scope(isolates, divergence_by_key, cds_index, "phylogroup"),
    }
    group_dir_root = {
        "species": paths.output_dir / "phylogenetics",
        "phylogroup": paths.output_dir / "phylogenetics_phylogroup",
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    table_rows: List[Dict[str, str]] = []

    for idx, winner in enumerate(winners):
        triplet = winner["triplet"]
        candidate_id = f"ga-winner-{idx:02d}"
        fig, axes = plt.subplots(1, 3, figsize=(13, 5))
        for ax, segment, donor_key in zip(axes, SEGMENTS, triplet):
            group_key = members = scope_used = None
            for scope in ("species", "phylogroup"):
                group_key, members = find_group_for_isolate(by_group[scope], donor_key)
                if group_key is not None:
                    scope_used = scope
                    break
            if group_key is None:
                ax.set_axis_off()
                ax.text(0.5, 0.5, f"{segment} donor:\nnot in a committed tree",
                        ha="center", va="center", wrap=True, transform=ax.transAxes)
                continue
            label_map = label_map_for_members(members)
            tree = load_tree_relabeled(group_dir_root[scope_used], group_key, segment, label_map)
            clade = local_neighborhood(tree, donor_key, min_tips=min_tips)
            draw_cladogram(ax, clade, highlight_name=donor_key, highlight_color=SLOT_COLOR[segment])
            ax.set_title(f"{segment} donor", color=SLOT_COLOR[segment], fontweight="bold")

        fig.suptitle(f"{candidate_id} -- fitness {winner.get('fitness'):.4f}", fontsize=10)
        fig.text(0.5, 0.01,
                 "In-silico combination of existing, independently-deposited isolates; "
                 "no new sequence is created or synthesized.",
                 ha="center", fontsize=7, color="#555555", fontweight="bold")
        fig.tight_layout(rect=(0, 0.03, 1, 1))
        out_path = output_dir / f"{candidate_id}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        written.append(out_path)

        row = {"candidate_id": candidate_id, "fitness": winner.get("fitness")}
        for segment, donor_key in zip(SEGMENTS, triplet):
            meta = parent_metadata(donor_key, by_key)
            for field in ("isolate_key", "species", "country", "host"):
                row[f"{segment.lower()}_parent_{field}"] = meta[field]
        table_rows.append(row)

    table_path = output_dir / "ga_winners_parents.csv"
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table_rows[0].keys()) if table_rows else [])
        writer.writeheader()
        writer.writerows(table_rows)

    return written, table_path


# ============================================================================
# CLI
# ============================================================================

def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    before = sub.add_parser("before")
    before.add_argument("--input", default=".")
    before.add_argument("--output", default="output/figures/before")
    before.add_argument("--min-tips", type=int, default=DEFAULT_MIN_TIPS)

    after = sub.add_parser("after")
    after.add_argument("--input", default=".")
    after.add_argument("--ga-search-json", required=True)
    after.add_argument("--output", default="output/figures/after")
    after.add_argument("--min-tips", type=int, default=DEFAULT_MIN_TIPS)

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if args.command == "before":
        written = render_before_figures(Path(args.output), args.input, min_tips=args.min_tips)
        print(f"wrote {len(written)} figure(s) to {args.output}")
    else:
        written, table_path = render_after_figures(
            Path(args.ga_search_json), Path(args.output), args.input, min_tips=args.min_tips)
        print(f"wrote {len(written)} figure(s) and {table_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# ============================================================================
# Tests
# ============================================================================

class LabelMapReconstructionTest(unittest.TestCase):
    """The one piece of real computation this module does: reconstructing which
    synthetic tip label a committed tree assigned to which isolate. Wrong here means
    every highlighted tip in every figure is silently the wrong isolate."""

    def _candidate(self, isolate_key):
        from types import SimpleNamespace as NS
        isolate = NS(isolate_key=isolate_key)
        return ph.PhyloCandidate(isolate=isolate, tags=["complete_triplet"], sequences={})

    def test_labels_assigned_in_member_order(self):
        members = [self._candidate(k) for k in ("iso_a", "iso_b", "iso_c")]
        label_map = label_map_for_members(members)
        self.assertEqual(label_map["I0000"].isolate.isolate_key, "iso_a")
        self.assertEqual(label_map["I0001"].isolate.isolate_key, "iso_b")
        self.assertEqual(label_map["I0002"].isolate.isolate_key, "iso_c")

    def test_matches_group_candidates_own_sort_order(self):
        """Guard the guard: group_candidates's docstring guarantees members sorted by
        isolate_key ("so re-runs assign identical tip labels"). label_map_for_members
        must consume that same order -- pin it against group_candidates itself, using
        real IsolateRecord-shaped candidates (it reads .species_name internally)."""
        from types import SimpleNamespace as NS
        unsorted = ["iso_c", "iso_a", "iso_b"]
        candidates = [ph.PhyloCandidate(isolate=NS(isolate_key=k, species_name="sp"),
                                        tags=["complete_triplet"], sequences={})
                     for k in unsorted]
        accepted, _ = ph.group_candidates(candidates, min_taxa=1, scope="species")
        label_map = label_map_for_members(accepted["sp"])
        self.assertEqual([label_map[f"I{i:04d}"].isolate.isolate_key for i in range(3)],
                         ["iso_a", "iso_b", "iso_c"])

    def test_find_group_for_isolate_locates_the_right_group(self):
        members_x = [self._candidate("a"), self._candidate("b")]
        members_y = [self._candidate("c")]
        by_group = {"group_x": members_x, "group_y": members_y}
        key, members = find_group_for_isolate(by_group, "c")
        self.assertEqual(key, "group_y")
        self.assertIs(members, members_y)

    def test_find_group_for_isolate_absent_returns_none(self):
        key, members = find_group_for_isolate({"g": [self._candidate("a")]}, "nope")
        self.assertIsNone(key)
        self.assertIsNone(members)


class LocalNeighborhoodTest(unittest.TestCase):
    def _tree(self):
        from io import StringIO
        from Bio import Phylo
        newick = "((((t1:0.1,t2:0.1):0.1,t3:0.1):0.1,(t4:0.1,t5:0.1):0.1):0.1,(t6:0.1,t7:0.1):0.1);"
        return Phylo.read(StringIO(newick), "newick")

    def test_returns_ancestor_with_at_least_min_tips(self):
        tree = self._tree()
        clade = local_neighborhood(tree, "t1", min_tips=3)
        self.assertGreaterEqual(clade.count_terminals(), 3)
        names = {t.name for t in clade.get_terminals()}
        self.assertIn("t1", names)

    def test_min_tips_larger_than_tree_returns_whole_tree(self):
        tree = self._tree()
        clade = local_neighborhood(tree, "t1", min_tips=100)
        self.assertEqual(clade.count_terminals(), 7)

    def test_does_not_mutate_the_source_tree(self):
        tree = self._tree()
        before = tree.count_terminals()
        local_neighborhood(tree, "t1", min_tips=3).ladderize()
        self.assertEqual(tree.count_terminals(), before)
        self.assertIsNotNone(tree.get_path("t1"))


class DrawCladogramTest(unittest.TestCase):
    def test_returns_position_for_every_tip(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from io import StringIO
        from Bio import Phylo
        newick = "((t1:0.1,t2:0.2):0.1,t3:0.3);"
        tree = Phylo.read(StringIO(newick), "newick")
        fig, ax = plt.subplots()
        try:
            positions = draw_cladogram(ax, tree.root, highlight_name="t2")
            self.assertEqual(set(positions), {"t1", "t2", "t3"})
            self.assertGreater(positions["t2"][0], 0)
        finally:
            plt.close(fig)


class RowSelectionTest(unittest.TestCase):
    def test_species_scope_preferred_over_phylogroup_for_same_isolate(self):
        import tempfile
        rows = [
            {"isolate_key": "iso1", "reassortment_call": "supported_reassortment_candidate",
             "analysis_scope": "phylogroup", "species_name": "sp", "discordant_segment": "L"},
            {"isolate_key": "iso1", "reassortment_call": "supported_reassortment_candidate",
             "analysis_scope": "species", "species_name": "sp", "discordant_segment": "L"},
            {"isolate_key": "iso2", "reassortment_call": "no_conflict",
             "analysis_scope": "species", "species_name": "sp", "discordant_segment": ""},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ortho_reassortment_table.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
                writer.writeheader()
                writer.writerows(rows)
            result = supported_reassortment_rows(path)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["analysis_scope"], "species")


class RenderBeforeFiguresIntegrationTest(unittest.TestCase):
    """Smoke test against the real committed all_hanta/ package and output/ -- confirms
    the whole chain (ingest -> grouping -> label reconstruction -> tree load -> prune ->
    draw) actually produces a non-empty figure for a real supported_reassortment_candidate
    isolate, not just that each piece works in isolation."""

    def test_renders_at_least_one_real_figure(self):
        if not Path("all_hanta/ncbi_dataset/data").exists():
            self.skipTest("all_hanta/ package not present in this environment")
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            written = render_before_figures(Path(tmp), input_dir=".", min_tips=6)
            self.assertGreaterEqual(len(written), 1)
            for path in written:
                self.assertTrue(path.exists())
                self.assertGreater(path.stat().st_size, 1000)
