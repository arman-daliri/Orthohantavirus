"""Before/after phylogenetic-tree figures over Stage 2's COMMITTED outputs only.

No new scoring, no new gate calls, no new sequences. Every tree drawn here is either
read directly from a `tree.treefile` already committed under `output/phylogenetics*/`
(Stage 2's own MAFFT+IQ-TREE output) or is that same file pruned to a local neighborhood
for legibility -- never rebuilt, realigned, or reinferred.

This module never reconstructs which tree group a segment belongs to: tree paths are
read straight off `ortho_reassortment_table.csv`'s `s_tree_path`/`m_tree_path`/
`l_tree_path` columns, and each committed tree's synthetic tip names (`I0000`, `I0001`,
...) are relabeled back to `isolate_key` using the `labels.json` file `phylogeny.py`
persists alongside each tree group's S/M/L segment directories. No mafft/iqtree
subprocess runs anywhere in this module, and no grouping logic is re-derived.

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
from typing import Dict, List, Optional, Sequence, Tuple

import prepare_dataset as prep

SEGMENTS = ("S", "M", "L")
SLOT_COLOR = {"S": "#1b6ca8", "M": "#2e8b57", "L": "#c1440e"}
DEFAULT_MIN_TIPS = 8


# ============================================================================
# Ingest + tip relabeling (reads phylogeny.py's persisted labels.json; never
# reconstructs group membership, never rebuilds a tree)
# ============================================================================

def load_ingest(output_dir: str = "output") -> List["prep.IsolateRecord"]:
    """Isolate metadata only (species/country/host, for the AFTER panel's parent-metadata
    table) -- loaded from Stage 1's reference snapshot, never all_hanta/ncbi_dataset/
    directly. Grouping/divergence/cds_index are not needed here: phylogeny.py now
    persists tip-label -> isolate_key maps as labels.json per tree group, and tree paths
    are read straight off ortho_reassortment_table.csv, so this module never reconstructs
    group membership itself (see docs/superpowers/specs/2026-09-09-stage1-reference-
    snapshot-design.md's amendment).
    """
    return prep.load_reference_data(Path(output_dir)).isolates


def load_tree_relabeled(tree_path: Path) -> "Bio.Phylo.BaseTree.Tree":
    """Loads a committed tree.treefile and relabels its synthetic tip names (I0000, ...)
    back to isolate_key, using the labels.json phylogeny.py wrote alongside the group's
    S/M/L segment directories (tree_path.parent is the segment dir, tree_path.parent.parent
    is the group dir labels.json lives in).
    """
    from Bio import Phylo
    if not tree_path.exists():
        raise FileNotFoundError(f"committed tree not found: {tree_path}")
    labels_path = tree_path.parent.parent / "labels.json"
    if not labels_path.exists():
        raise FileNotFoundError(f"labels.json not found next to {tree_path} "
                                f"(expected at {labels_path})")
    label_to_isolate_key = json.loads(labels_path.read_text(encoding="utf-8"))
    tree = Phylo.read(str(tree_path), "newick")
    for tip in tree.get_terminals():
        isolate_key = label_to_isolate_key.get(tip.name)
        if isolate_key is not None:
            tip.name = isolate_key
    return tree


# ============================================================================
# Local-neighborhood pruning (legibility only -- never touches the underlying tree file)
# ============================================================================

def local_neighborhood(tree: "Bio.Phylo.BaseTree.Tree", target_name: str,
                       min_tips: int = DEFAULT_MIN_TIPS) -> "Bio.Phylo.BaseTree.Clade":
    """Smallest ancestor clade of `target_name` with at least `min_tips` terminals, or
    the whole tree if it never reaches that size. A deep copy so callers can mutate
    (ladderize) the returned clade without touching the cached full tree.

    Raises ValueError if `target_name` is not a tip of `tree` -- this means the tree
    and the caller's isolate_key have drifted (stale CSV row, regenerated tree, or a
    tip missing from labels.json), and silently falling back to the full tree would
    render a wrong, unhighlighted figure instead of failing loudly.
    """
    path = tree.get_path(target_name)
    if not path:
        raise ValueError(f"{target_name!r} is not a tip of the given tree -- "
                         f"tree and isolate_key/labels.json have likely drifted")
    for ancestor in reversed(path[:-1]):
        if ancestor.count_terminals() >= min_tips:
            return copy.deepcopy(ancestor)
    return copy.deepcopy(tree.root)


# ============================================================================
# Minimal cladogram renderer (figure-drawing code -- not phylogenetic inference; the
# instruction to never hand-roll a tree refers to alignment/tree-BUILDING, which this
# module never does)
# ============================================================================

def _layout(clade: "Bio.Phylo.BaseTree.Clade",
            ladderize: bool = True) -> Tuple[Dict[int, float], List["Bio.Phylo.BaseTree.Clade"]]:
    if ladderize:
        clade.ladderize()
    tips = clade.get_terminals()
    y: Dict[int, float] = {id(t): float(i) for i, t in enumerate(tips)}

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


def tree_paths_by_isolate_from_table(reassortment_table: Path) -> Dict[str, Dict[str, str]]:
    """isolate_key -> {segment: tree_path}, for every analyzed isolate regardless of its
    reassortment_call -- a GA winner need not be a supported_reassortment_candidate to
    have a committed tree. Species-scope preferred over phylogroup-scope per isolate,
    same precedence supported_reassortment_rows uses, and only rows with at least one
    non-empty tree path are kept (phylo_status='species_below_min_taxa'/'failed' rows
    have none).
    """
    with reassortment_table.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    by_isolate: Dict[str, Dict[str, str]] = {}
    scope_used: Dict[str, str] = {}
    for row in rows:
        paths = {seg: row.get(f"{seg.lower()}_tree_path") or "" for seg in SEGMENTS}
        if not any(paths.values()):
            continue
        key = row["isolate_key"]
        scope = row["analysis_scope"]
        if key not in by_isolate or (scope_used[key] == "phylogroup" and scope == "species"):
            by_isolate[key] = paths
            scope_used[key] = scope
    return by_isolate


def render_before_figures(output_dir: Path, data_dir: str = "output",
                          reassortment_table: Optional[Path] = None,
                          min_tips: int = DEFAULT_MIN_TIPS) -> List[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import ConnectionPatch

    stage_dir = Path(data_dir)
    table_path = reassortment_table or (stage_dir / "ortho_reassortment_table.csv")
    rows = supported_reassortment_rows(table_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    seen_names: Dict[str, int] = {}
    for row in rows:
        isolate_key = row["isolate_key"]
        discordant = row.get("discordant_segment") or ""
        fig, axes = plt.subplots(1, 3, figsize=(13, 5))
        tip_pos: Dict[str, Tuple[float, float]] = {}
        for ax, segment in zip(axes, SEGMENTS):
            tree_path_str = row.get(f"{segment.lower()}_tree_path") or ""
            if not tree_path_str:
                ax.set_axis_off()
                ax.text(0.5, 0.5, f"{segment}: not in a committed tree",
                        ha="center", va="center", wrap=True, transform=ax.transAxes)
                continue
            tree = load_tree_relabeled(stage_dir / tree_path_str)
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
        row_scope = row["analysis_scope"]
        fig.text(0.5, 0.01,
                 "Real, unmodified S/M/L sequences from a single deposited isolate; "
                 "trees pruned to a local neighborhood for legibility (full tree in "
                 f"output/phylogenetics{'_phylogroup' if row_scope == 'phylogroup' else ''}/).",
                 ha="center", fontsize=7, color="#555555")
        fig.tight_layout(rect=(0, 0.03, 1, 1))
        base_name = prep.normalize_attr_name(isolate_key)
        if base_name in seen_names:
            seen_names[base_name] += 1
            file_name = f"{base_name}_{seen_names[base_name]}"
        else:
            seen_names[base_name] = 0
            file_name = base_name
        out_path = output_dir / f"{file_name}.png"
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


def render_after_figures(ga_search_json: Path, output_dir: Path, data_dir: str = "output",
                         min_tips: int = DEFAULT_MIN_TIPS) -> Tuple[List[Path], Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    payload = json.loads(ga_search_json.read_text(encoding="utf-8"))
    winners = payload["winners"]

    stage_dir = Path(data_dir)
    isolates = load_ingest(data_dir)
    by_key = {i.isolate_key: i for i in isolates}
    tree_paths_by_isolate = tree_paths_by_isolate_from_table(stage_dir / "ortho_reassortment_table.csv")

    output_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    table_rows: List[Dict[str, str]] = []

    for idx, winner in enumerate(winners):
        triplet = winner["triplet"]
        candidate_id = f"ga-winner-{idx:02d}"
        fig, axes = plt.subplots(1, 3, figsize=(13, 5))
        for ax, segment, donor_key in zip(axes, SEGMENTS, triplet):
            tree_path_str = tree_paths_by_isolate.get(donor_key, {}).get(segment, "")
            if not tree_path_str:
                ax.set_axis_off()
                ax.text(0.5, 0.5, f"{segment} donor:\nnot in a committed tree",
                        ha="center", va="center", wrap=True, transform=ax.transAxes)
                continue
            tree = load_tree_relabeled(stage_dir / tree_path_str)
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

        # rpi_e_unknown_components (2026-09-26): comma-joined subset of {host,geo,time}
        # whose Factor E term hit a "missing metadata" neutral sentinel rather than a
        # measured value, for THIS winner's actual parent pair. Surfaced here, not just
        # in the raw all_scored dump, because this is the file a reader opens to see
        # who a winner's donors are -- a winner whose fitness rests partly on unscored
        # metadata should say so at the point where that fitness is being read, not only
        # in a 1000+-row CSV three steps upstream.
        row = {"candidate_id": candidate_id, "fitness": winner.get("fitness"),
              "rpi_e_unknown_components": winner.get("rpi_e_unknown_components") or "",
              "rpi_neutral_factors": winner.get("rpi_neutral_factors") or ""}
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
# TREES panel: every group's S/M/L trees as pairwise tanglegrams. Display only -- no
# reassortment calls, no highlighting, no distances. Reads the committed tree files.
# ============================================================================

SEGMENT_PAIRS = (("S", "M"), ("M", "L"), ("S", "L"))
TREE_WIDTH_IN = 5.0
LINK_WIDTH_IN = 2.0
ROW_HEIGHT_IN = 0.11
LABEL_FONT_PT = 6.0
TREE_COLOR = "#444444"
LINK_COLOR = "#7a8794"


def short_tip_label(isolate_key: str) -> str:
    """Drops the leading 'tier|species|' fields of an isolate_key; display only."""
    return isolate_key.split("|", 2)[-1]


def order_by_reference(clade: "Bio.Phylo.BaseTree.Clade", ref_y: Dict[str, float]) -> None:
    """Reorders children (in place) so each clade sits near the mean position its tips
    have in a reference tree -- fewer crossing lines in the tanglegram. Changes drawing
    order only; the tree's topology and branch lengths are untouched."""
    neutral = (len(ref_y) - 1) / 2.0
    mean_y: Dict[int, float] = {}

    def walk(c) -> Tuple[float, int]:
        if c.is_terminal():
            if c.name in ref_y:
                mean_y[id(c)] = ref_y[c.name]
                return ref_y[c.name], 1
            mean_y[id(c)] = neutral
            return 0.0, 0
        total, count = 0.0, 0
        for child in c.clades:
            child_total, child_count = walk(child)
            total += child_total
            count += child_count
        mean_y[id(c)] = total / count if count else neutral
        return total, count

    walk(clade)
    for c in clade.find_clades(order="level"):
        if not c.is_terminal():
            c.clades.sort(key=lambda child: mean_y[id(child)])


def _scale_bar_length(max_depth: float) -> float:
    import math
    target = max_depth / 4.0
    if target <= 0:
        return 0.0
    exponent = math.floor(math.log10(target))
    for mult in (5, 2, 1):
        if mult * 10 ** exponent <= target:
            return mult * 10 ** exponent


def _draw_tree_side(ax, clade, y, x_offset: float, mirrored: bool):
    """Draws one tree into `ax` (axes in inches). Returns ({tip_name: (tip_x, y)},
    max_depth). Tips of a tree run ragged by branch length; mirrored=True flips it so
    its tips face left."""
    from matplotlib.collections import LineCollection
    x = _x_positions(clade)
    max_depth = max(x.values()) if x else 0.0
    scale = TREE_WIDTH_IN / max_depth if max_depth > 0 else 0.0

    def fx(depth: float) -> float:
        offset = depth * scale
        return x_offset + (TREE_WIDTH_IN - offset if mirrored else offset)

    segments = []

    def walk(c):
        if c.is_terminal():
            return
        child_ys = [y[id(child)] for child in c.clades]
        segments.append(((fx(x[id(c)]), min(child_ys)), (fx(x[id(c)]), max(child_ys))))
        for child in c.clades:
            segments.append(((fx(x[id(c)]), y[id(child)]), (fx(x[id(child)]), y[id(child)])))
            walk(child)

    walk(clade)
    ax.add_collection(LineCollection(segments, colors=TREE_COLOR, linewidths=0.5))
    tip_pos = {t.name: (fx(x[id(t)]), y[id(t)]) for t in clade.get_terminals()}
    return tip_pos, max_depth


def draw_tanglegram(left_tree, right_tree, left_name: str, right_name: str, title: str):
    """One figure: left tree (tips facing right) and right tree (mirrored), the tip
    labels between them, and a thin line joining each isolate present in both trees.
    The only things drawn are the trees as written by IQ-TREE and those joins."""
    import copy as _copy
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    left = _copy.deepcopy(left_tree.root)
    right = _copy.deepcopy(right_tree.root)
    left_y, left_tips = _layout(left)
    left_pos_by_name = {t.name: left_y[id(t)] for t in left_tips}
    order_by_reference(right, left_pos_by_name)
    right_y, right_tips = _layout(right, ladderize=False)

    n_rows = max(len(left_tips), len(right_tips))
    label_chars = max((len(short_tip_label(t.name or "")) for t in left_tips + right_tips), default=0)
    label_in = label_chars * LABEL_FONT_PT * 0.55 / 72.0
    fig_w = 2 * TREE_WIDTH_IN + 2 * label_in + LINK_WIDTH_IN + 0.6
    fig_h = max(4.0, n_rows * ROW_HEIGHT_IN + 2.2)

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, fig_w)
    ax.set_ylim(n_rows + 3.2, -2.5)
    ax.axis("off")

    left_x0 = 0.3
    right_x0 = fig_w - 0.3 - TREE_WIDTH_IN
    left_pos, left_depth = _draw_tree_side(ax, left, left_y, left_x0, mirrored=False)
    right_pos, right_depth = _draw_tree_side(ax, right, right_y, right_x0, mirrored=True)

    left_label_x = left_x0 + TREE_WIDTH_IN + 0.05
    right_label_x = right_x0 - 0.05
    leaders, links = [], []
    for name, (tx, ty) in left_pos.items():
        leaders.append(((tx, ty), (left_label_x - 0.02, ty)))
        ax.text(left_label_x, ty, short_tip_label(name), fontsize=LABEL_FONT_PT,
                va="center", ha="left")
    for name, (tx, ty) in right_pos.items():
        leaders.append(((tx, ty), (right_label_x + 0.02, ty)))
        ax.text(right_label_x, ty, short_tip_label(name), fontsize=LABEL_FONT_PT,
                va="center", ha="right")
    shared = [n for n in left_pos if n in right_pos]
    for name in shared:
        links.append(((left_label_x + label_in + 0.05, left_pos[name][1]),
                      (right_label_x - label_in - 0.05, right_pos[name][1])))
    ax.add_collection(LineCollection(leaders, colors="#cccccc", linewidths=0.3,
                                     linestyles="dotted"))
    ax.add_collection(LineCollection(links, colors=LINK_COLOR, linewidths=0.4, alpha=0.7))

    bar_y = n_rows + 0.8
    for x0, mirrored, depth in ((left_x0, False, left_depth), (right_x0, True, right_depth)):
        bar = _scale_bar_length(depth)
        if bar > 0 and depth > 0:
            width = bar / depth * TREE_WIDTH_IN
            start = x0 + TREE_WIDTH_IN - width if mirrored else x0
            ax.plot([start, start + width], [bar_y, bar_y], color="black", lw=1)
            ax.text(start + width / 2, bar_y + 0.5, f"{bar:g} subst./site", ha="center",
                    va="top", fontsize=7)
    ax.text(left_x0, -1.6, f"{left_name} segment", fontsize=10, fontweight="bold", va="center")
    ax.text(fig_w - 0.3, -1.6, f"{right_name} segment", fontsize=10, fontweight="bold",
            va="center", ha="right")
    ax.text(fig_w / 2, -1.6, title, fontsize=9, ha="center", va="center")
    ax.text(fig_w / 2, -0.7,
            f"{len(left_tips)} tips | {len(right_tips)} tips | {len(shared)} shared (joined by lines)",
            fontsize=7, ha="center", va="center", color="#555555")
    ax.text(fig_w / 2, n_rows + 2.7,
            "IQ-TREE output, unrooted (no outgroup), drawn from the file's basal node. Each tree is "
            "scaled to the same drawn width; see scale bars. Right-hand tree order rearranged to "
            "reduce line crossings.",
            fontsize=6, ha="center", va="center", color="#555555")
    return fig, len(shared)


def discover_tree_groups(data_dir: str = "output") -> List[Tuple[str, Path]]:
    """[(scope_label, group_dir)] for every group under output/phylogenetics*/ that has a
    labels.json and at least one segment tree on disk."""
    found = []
    for scope in ("phylogenetics", "phylogenetics_phylogroup"):
        scope_dir = Path(data_dir) / scope
        if not scope_dir.is_dir():
            continue
        for group_dir in sorted(p for p in scope_dir.iterdir() if p.is_dir()):
            if (group_dir / "labels.json").exists() and any(
                    (group_dir / seg / "tree.treefile").exists() for seg in SEGMENTS):
                found.append((scope, group_dir))
    return found


def render_tree_figures(output_dir: Path, data_dir: str = "output", fmt: str = "pdf",
                        only_group: Optional[str] = None) -> List[Path]:
    """One multi-page PDF per group (pages S-M, M-L, S-L), or one PNG per pair with
    fmt='png'. A pair is skipped, with a logged note, if either segment's tree is missing."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    if fmt not in ("pdf", "png"):
        raise ValueError(f"fmt must be 'pdf' or 'png', got {fmt!r}")
    written: List[Path] = []
    for scope, group_dir in discover_tree_groups(data_dir):
        if only_group and group_dir.name != only_group:
            continue
        trees = {}
        for seg in SEGMENTS:
            tree_path = group_dir / seg / "tree.treefile"
            if tree_path.exists():
                trees[seg] = load_tree_relabeled(tree_path)
        out_dir = Path(output_dir) / scope
        out_dir.mkdir(parents=True, exist_ok=True)
        pdf_path = out_dir / f"{group_dir.name}.pdf"
        pdf = None
        try:
            for left_seg, right_seg in SEGMENT_PAIRS:
                if left_seg not in trees or right_seg not in trees:
                    logging.info("skip %s %s-%s: segment tree missing", group_dir.name, left_seg, right_seg)
                    continue
                fig, _shared = draw_tanglegram(trees[left_seg], trees[right_seg], left_seg, right_seg,
                                               f"{group_dir.name} ({scope})")
                if fmt == "pdf":
                    if pdf is None:
                        pdf = PdfPages(pdf_path)
                    pdf.savefig(fig)
                else:
                    path = out_dir / f"{group_dir.name}_{left_seg}-{right_seg}.png"
                    fig.savefig(path, dpi=100)
                    written.append(path)
                plt.close(fig)
        finally:
            if pdf is not None:
                pdf.close()
                if pdf_path.exists() and pdf_path.stat().st_size > 0:
                    written.append(pdf_path)
    return written


# ============================================================================
# PHYLOGROUP OVERVIEW — species-condensed cladogram per host family
# ============================================================================

_PHYLOGROUP_LABELS = {
    "arvicolinae_borne": "Arvicolinae-borne (voles, lemmings)",
    "murinae_borne": "Murinae-borne (Old World rats & mice)",
    "sigmodontinae_borne": "Sigmodontinae-borne (New World rats & mice)",
    "soricomorpha_borne": "Soricomorpha-borne (shrews & moles)",
}


def _species_short(leaf_name: str) -> str:
    """'t2|Orthohantavirus puumalaense|foo' -> 'puumalaense'; falls back to full name."""
    parts = leaf_name.split("|")
    if len(parts) >= 2:
        words = parts[1].split()
        return words[-1] if words else leaf_name
    return leaf_name


def _prune_to_species_representatives(tree: "Bio.Phylo.BaseTree.Tree") -> "Bio.Phylo.BaseTree.Tree":
    """Return a deep copy of tree with all but one representative per species pruned.
    The kept representative is the alphabetically first isolate_key per species."""
    import copy
    from collections import defaultdict
    t = copy.deepcopy(tree)
    by_species: Dict[str, List[str]] = defaultdict(list)
    for tip in t.get_terminals():
        by_species[_species_short(tip.name)].append(tip.name)
    to_prune = []
    for names in by_species.values():
        to_prune.extend(sorted(names)[1:])  # keep sorted-first; remove the rest
    for name in to_prune:
        t.prune(name)
    # Rename surviving tips to the short species name for readability
    for tip in t.get_terminals():
        tip.name = _species_short(tip.name)
    return t


def render_phylogroup_overview(output_dir: Path, data_dir: str = "output",
                               segment: str = "S") -> List[Path]:
    """One figure per phylogroup showing the species-condensed cladogram (one
    representative per species, labelled by short species name), plus a combined
    2×2 panel figure. Reads committed trees under output/phylogenetics_phylogroup/;
    raises FileNotFoundError if the directory is absent."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    phylogroup_dir = Path(data_dir) / "phylogenetics_phylogroup"
    if not phylogroup_dir.is_dir():
        raise FileNotFoundError(
            f"phylogroup output directory not found: {phylogroup_dir}\n"
            f"Run phylogeny.py (Stage 2) first.")

    groups = sorted(p for p in phylogroup_dir.iterdir()
                    if p.is_dir()
                    and (p / "labels.json").exists()
                    and (p / segment / "tree.treefile").exists())
    if not groups:
        raise FileNotFoundError(
            f"no {segment}-segment trees found under {phylogroup_dir}")

    written: List[Path] = []
    panel_items: List[Tuple[str, "Bio.Phylo.BaseTree.Tree", int]] = []

    for g in groups:
        raw = load_tree_relabeled(g / segment / "tree.treefile")
        n_isolates = len(raw.get_terminals())
        condensed = _prune_to_species_representatives(raw)
        panel_items.append((g.name, condensed, n_isolates))

    # Combined 2×2 figure
    n_panels = len(panel_items)
    ncols = 2
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 5 * nrows))
    import numpy as _np
    axes_flat = list(axes.flat) if isinstance(axes, _np.ndarray) else [axes]

    for idx, (g_name, tree, n_iso) in enumerate(panel_items):
        ax = axes_flat[idx]
        tip_pos = draw_cladogram(ax, tree.root)
        n_species = len(tip_pos)
        host_label = _PHYLOGROUP_LABELS.get(g_name, g_name.replace("_", " "))
        ax.set_title(f"{host_label}\n"
                     f"{n_species} species · {n_iso} isolates · {segment}-segment IQ-TREE",
                     fontsize=8, pad=6)
        ax.set_xlabel("substitutions/site", fontsize=7)
        # Label each tip
        x_positions = _x_positions(tree.root)
        y_positions, tips = _layout(tree.root)
        x_max = max(x_positions.values()) if x_positions else 1.0
        for tip in tips:
            cx = x_positions[id(tip)]
            cy = y_positions[id(tip)]
            ax.text(cx + x_max * 0.01, cy, tip.name,
                    fontsize=6, va="center", ha="left")
        # Extra right margin for labels
        ax.set_xlim(right=x_max * 1.55)

    # Hide unused panels (if panel count is odd)
    for idx in range(n_panels, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    fig.suptitle(
        f"Orthohantavirus species-condensed phylogenies by host family "
        f"({segment} segment, one representative per species)",
        fontsize=10, y=1.01)
    fig.tight_layout()
    combined_path = output_dir / f"phylogroup_overview_{segment}.png"
    fig.savefig(combined_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    written.append(combined_path)

    return written


# ============================================================================
# CLI
# ============================================================================

def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    before = sub.add_parser("before")
    before.add_argument("--input", default="output",
                        help="Stage 1/2 output directory (reference snapshot + "
                             "ortho_reassortment_table.csv + phylogenetics*/)")
    before.add_argument("--output", default="output/figures/before")
    before.add_argument("--min-tips", type=int, default=DEFAULT_MIN_TIPS)

    after = sub.add_parser("after")
    after.add_argument("--input", default="output",
                        help="Stage 1/2 output directory (reference snapshot + "
                             "ortho_reassortment_table.csv + phylogenetics*/)")
    after.add_argument("--ga-search-json", required=True)
    after.add_argument("--output", default="output/figures/after")
    after.add_argument("--min-tips", type=int, default=DEFAULT_MIN_TIPS)

    trees = sub.add_parser("trees", help="every group's S/M/L trees as tanglegram pairs "
                                         "(display only, no interpretation)")
    trees.add_argument("--input", default="output",
                       help="Stage 2 output directory containing phylogenetics*/")
    trees.add_argument("--output", default="output/figures/trees")
    trees.add_argument("--format", choices=("pdf", "png"), default="pdf")
    trees.add_argument("--group", default=None,
                       help="only this group directory name, e.g. orthohantavirus_tulaense")

    overview = sub.add_parser("overview",
                               help="species-condensed phylogroup overview: one representative "
                                    "per species per host family, coloured cladogram")
    overview.add_argument("--input", default="output",
                          help="Stage 2 output directory containing phylogenetics_phylogroup/")
    overview.add_argument("--output", default="output/figures/overview")
    overview.add_argument("--segment", choices=("S", "M", "L"), default="S",
                          help="which segment's tree to draw (default: S)")

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if args.command == "overview":
        written = render_phylogroup_overview(Path(args.output), args.input, args.segment)
        print(f"wrote {len(written)} file(s) to {args.output}")
        return 0
    if args.command == "trees":
        written = render_tree_figures(Path(args.output), args.input, args.format, args.group)
        print(f"wrote {len(written)} file(s) to {args.output}")
        return 0
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

class LoadTreeRelabeledTest(unittest.TestCase):
    """load_tree_relabeled must relabel a committed tree's synthetic tip names using
    labels.json persisted one directory up from the segment dir (tree_path.parent.parent).
    Wrong here means every highlighted tip in every figure is silently the wrong isolate."""

    def test_relabels_tips_from_sibling_labels_json(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            group_dir = Path(tmp) / "some_group"
            segment_dir = group_dir / "S"
            segment_dir.mkdir(parents=True)
            (segment_dir / "tree.treefile").write_text(
                "((I0000:0.1,I0001:0.1):0.1,I0002:0.1);", encoding="utf-8")
            (group_dir / "labels.json").write_text(
                json.dumps({"I0000": "iso_a", "I0001": "iso_b", "I0002": "iso_c"}),
                encoding="utf-8")
            tree = load_tree_relabeled(segment_dir / "tree.treefile")
            names = {t.name for t in tree.get_terminals()}
            self.assertEqual(names, {"iso_a", "iso_b", "iso_c"})

    def test_missing_tree_raises_file_not_found(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                load_tree_relabeled(Path(tmp) / "group" / "S" / "tree.treefile")

    def test_missing_labels_json_raises_file_not_found(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            segment_dir = Path(tmp) / "some_group" / "S"
            segment_dir.mkdir(parents=True)
            (segment_dir / "tree.treefile").write_text("(I0000:0.1,I0001:0.1);", encoding="utf-8")
            with self.assertRaises(FileNotFoundError):
                load_tree_relabeled(segment_dir / "tree.treefile")


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

    def test_missing_tip_raises_instead_of_returning_whole_tree(self):
        tree = self._tree()
        with self.assertRaises(ValueError):
            local_neighborhood(tree, "not_a_real_tip", min_tips=3)


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


class TreePathsByIsolateFromTableTest(unittest.TestCase):
    def test_species_scope_preferred_and_rows_without_any_tree_path_dropped(self):
        import tempfile
        rows = [
            {"isolate_key": "iso1", "analysis_scope": "phylogroup",
             "s_tree_path": "phylogenetics_phylogroup/g/S/tree.treefile",
             "m_tree_path": "", "l_tree_path": ""},
            {"isolate_key": "iso1", "analysis_scope": "species",
             "s_tree_path": "phylogenetics/g/S/tree.treefile",
             "m_tree_path": "phylogenetics/g/M/tree.treefile", "l_tree_path": ""},
            {"isolate_key": "iso2", "analysis_scope": "species",
             "s_tree_path": "", "m_tree_path": "", "l_tree_path": ""},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ortho_reassortment_table.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
                writer.writeheader()
                writer.writerows(rows)
            result = tree_paths_by_isolate_from_table(path)
        self.assertEqual(set(result), {"iso1"})
        self.assertEqual(result["iso1"]["S"], "phylogenetics/g/S/tree.treefile")
        self.assertEqual(result["iso1"]["M"], "phylogenetics/g/M/tree.treefile")


class RenderBeforeFiguresFilenameCollisionTest(unittest.TestCase):
    """Two isolate_keys that normalize to the same filename (case/punctuation only)
    must not silently overwrite one another -- each must get a distinct file, and
    `written` must accurately list exactly one path per input row."""

    def test_normalized_name_collision_gets_disambiguated(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            group_dir = tmp_path / "phylogenetics" / "g"
            seg_dir = group_dir / "S"
            seg_dir.mkdir(parents=True)
            newick = ("((((t1:0.1,t2:0.1):0.1,t3:0.1):0.1,(t4:0.1,t5:0.1):0.1):0.1,"
                      "(t6:0.1,t7:0.1):0.1);")
            (seg_dir / "tree.treefile").write_text(newick, encoding="utf-8")
            (group_dir / "labels.json").write_text(
                json.dumps({"t1": "ISO-1", "t2": "iso_1"}), encoding="utf-8")

            table_path = tmp_path / "ortho_reassortment_table.csv"
            rows = [
                {"isolate_key": "ISO-1", "reassortment_call": "supported_reassortment_candidate",
                 "analysis_scope": "species", "species_name": "sp", "discordant_segment": "",
                 "s_tree_path": "phylogenetics/g/S/tree.treefile", "m_tree_path": "", "l_tree_path": ""},
                {"isolate_key": "iso_1", "reassortment_call": "supported_reassortment_candidate",
                 "analysis_scope": "species", "species_name": "sp", "discordant_segment": "",
                 "s_tree_path": "phylogenetics/g/S/tree.treefile", "m_tree_path": "", "l_tree_path": ""},
            ]
            with table_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
                writer.writeheader()
                writer.writerows(rows)

            out_dir = tmp_path / "figs"
            written = render_before_figures(out_dir, data_dir=str(tmp_path),
                                            reassortment_table=table_path, min_tips=3)

            self.assertEqual(len(written), 2)
            self.assertEqual(len(set(written)), 2)
            for path in written:
                self.assertTrue(path.exists())


class RenderAfterFiguresUnknownComponentsColumnTest(unittest.TestCase):
    """ga_winners_parents.csv must carry rpi_e_unknown_components (2026-09-26) so a
    reader of THIS file -- the one that names a winner's actual donors -- can see
    whether that winner's fitness rests partly on a missing-metadata sentinel rather
    than a measured value, without having to cross-reference the all_scored dump.
    Isolated from real output/ data (empty load_ingest, no tree paths) so this checks
    the CSV-writing logic only, not tree rendering."""

    def test_column_present_and_carries_the_winner_flag(self):
        import tempfile
        from unittest import mock

        winners = [
            {"triplet": ["s1", "m1", "l1"], "fitness": 0.4716981132,
             "rpi_e_unknown_components": "host,geo",
             "rpi_neutral_factors": "E:host+geo;M_pol:no_signal_beyond_tag"},
            {"triplet": ["s2", "m2", "l2"], "fitness": 0.01,
             "rpi_e_unknown_components": ""},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ga_json = tmp_path / "seed.json"
            ga_json.write_text(json.dumps({"winners": winners}), encoding="utf-8")
            out_dir = tmp_path / "after"
            with mock.patch("figures.load_ingest", return_value=[]), \
                 mock.patch("figures.tree_paths_by_isolate_from_table", return_value={}):
                _written, table_path = render_after_figures(ga_json, out_dir, data_dir=str(tmp_path))
            with table_path.open(newline="", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh))
            self.assertIn("rpi_e_unknown_components", rows[0])
            self.assertEqual(rows[0]["rpi_e_unknown_components"], "host,geo")
            self.assertEqual(rows[1]["rpi_e_unknown_components"], "")
            self.assertEqual(rows[0]["rpi_neutral_factors"], "E:host+geo;M_pol:no_signal_beyond_tag")
            self.assertEqual(rows[1]["rpi_neutral_factors"], "")


class RenderBeforeFiguresIntegrationTest(unittest.TestCase):
    """Smoke test against the real committed reference snapshot and phylogenetics output
    under output/ -- confirms the whole chain (read CSV tree paths -> relabel from
    labels.json -> prune -> draw) actually produces a non-empty figure for a real
    supported_reassortment_candidate isolate, not just that each piece works in isolation."""

    def test_renders_at_least_one_real_figure(self):
        if not Path("output/phylo_rpi_ga_isolates.json").exists():
            self.skipTest("output/ reference snapshot not present -- run prepare_dataset.py first")
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            written = render_before_figures(Path(tmp), data_dir="output", min_tips=6)
            self.assertGreaterEqual(len(written), 1)
            for path in written:
                self.assertTrue(path.exists())
                self.assertGreater(path.stat().st_size, 1000)


class TreeTanglegramTest(unittest.TestCase):
    """The trees panel is display-only: it must draw every tip of both trees, join only
    isolates present in both, and never alter topology or branch lengths."""

    def _write_group(self, root: Path, newicks: Dict[str, str]) -> Path:
        group = root / "phylogenetics" / "grp"
        for seg, newick in newicks.items():
            (group / seg).mkdir(parents=True)
            (group / seg / "tree.treefile").write_text(newick, encoding="utf-8")
        (group / "labels.json").write_text(json.dumps(
            {"I0": "t2|sp|a", "I1": "t2|sp|b", "I2": "t2|sp|c", "I3": "t2|sp|d"}), encoding="utf-8")
        return group

    def test_order_by_reference_keeps_topology_and_tips(self):
        import io
        from Bio import Phylo
        tree = Phylo.read(io.StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
        before = sorted(t.name for t in tree.get_terminals())
        order_by_reference(tree.root, {"c": 0.0, "d": 1.0, "a": 2.0, "b": 3.0})
        self.assertEqual([t.name for t in tree.get_terminals()], ["c", "d", "a", "b"])
        self.assertEqual(sorted(t.name for t in tree.get_terminals()), before)
        self.assertAlmostEqual(tree.total_branch_length(), 6.0)

    def test_short_tip_label_drops_tier_and_species(self):
        self.assertEqual(short_tip_label("t1|Orthohantavirus x|SAMN1|name"), "SAMN1|name")
        self.assertEqual(short_tip_label("t2|Orthohantavirus x|name"), "name")

    def test_pdf_has_three_pages_and_png_has_three_files(self):
        import re
        import tempfile
        newicks = {"S": "((I0:0.1,I1:0.1):0.1,(I2:0.1,I3:0.1):0.1);",
                   "M": "((I0:0.2,I2:0.1):0.1,(I1:0.1,I3:0.3):0.1);",
                   "L": "((I0:0.1,I3:0.1):0.1,(I1:0.1,I2:0.1):0.1);"}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write_group(root, newicks)
            pdfs = render_tree_figures(root / "out", data_dir=str(root), fmt="pdf")
            self.assertEqual(len(pdfs), 1)
            self.assertGreater(pdfs[0].stat().st_size, 1000)
            self.assertEqual(len(re.findall(rb"/Type\s*/Page(?!s)", pdfs[0].read_bytes())), 3)
            pngs = render_tree_figures(root / "png", data_dir=str(root), fmt="png")
            self.assertEqual(len(pngs), 3)

    def test_pair_with_missing_segment_tree_is_skipped_not_raised(self):
        import tempfile
        newicks = {"S": "((I0:0.1,I1:0.1):0.1,(I2:0.1,I3:0.1):0.1);",
                   "M": "((I0:0.2,I2:0.1):0.1,(I1:0.1,I3:0.3):0.1);"}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write_group(root, newicks)
            pngs = render_tree_figures(root / "png", data_dir=str(root), fmt="png")
            self.assertEqual([p.name for p in pngs], ["grp_S-M.png"])

    def test_group_with_a_single_segment_tree_writes_nothing_and_does_not_raise_for_pdf(self):
        """/code-review finding: PdfPages was opened before any page existed, so a group
        with only one segment tree (no drawable pair) crashed the whole run on unlink."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write_group(root, {"S": "((I0:0.1,I1:0.1):0.1,(I2:0.1,I3:0.1):0.1);"})
            written = render_tree_figures(root / "out", data_dir=str(root), fmt="pdf")
            self.assertEqual(written, [])
            self.assertFalse((root / "out" / "phylogenetics" / "grp.pdf").exists())

    def test_only_isolates_in_both_trees_are_joined(self):
        import io
        from Bio import Phylo
        left = Phylo.read(io.StringIO("((a:1,b:1):1,c:1);"), "newick")
        right = Phylo.read(io.StringIO("((a:1,b:1):1,d:1);"), "newick")
        fig, shared = draw_tanglegram(left, right, "S", "M", "t")
        self.assertEqual(shared, 2)


class PhylogroupOverviewTest(unittest.TestCase):
    def _write_phylogroup(self, root: Path, group_name: str,
                          newick: str, labels: Dict[str, str]) -> None:
        g = root / "phylogenetics_phylogroup" / group_name
        seg = g / "S"
        seg.mkdir(parents=True)
        (seg / "tree.treefile").write_text(newick, encoding="utf-8")
        (g / "labels.json").write_text(
            json.dumps(labels), encoding="utf-8")

    def test_species_short_parses_isolate_key(self):
        self.assertEqual(
            _species_short("t2|Orthohantavirus puumalaense|foo bar"),
            "puumalaense")

    def test_species_short_falls_back_for_bare_name(self):
        self.assertEqual(_species_short("bare"), "bare")

    def test_prune_keeps_one_per_species(self):
        from Bio.Phylo.BaseTree import Tree, Clade
        # Build tree directly: two puumalaense tips, one hantanense, one tulaense
        def tip(name): return Clade(branch_length=0.1, name=name)
        root = Clade(clades=[
            Clade(branch_length=0.1, clades=[
                tip("t2|Orthohantavirus puumalaense|iso1"),
                tip("t2|Orthohantavirus puumalaense|iso2"),
            ]),
            Clade(branch_length=0.1, clades=[
                tip("t2|Orthohantavirus hantanense|iso3"),
                tip("t2|Orthohantavirus tulaense|iso4"),
            ]),
        ])
        tree = Tree(root=root)
        condensed = _prune_to_species_representatives(tree)
        tip_names = {t.name for t in condensed.get_terminals()}
        self.assertEqual(tip_names, {"puumalaense", "hantanense", "tulaense"})

    def test_render_phylogroup_overview_writes_combined_png(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write_phylogroup(
                root, "murinae_borne",
                "((I0:0.1,I1:0.1):0.1,(I2:0.1,I3:0.1):0.1);",
                {"I0": "t2|Orthohantavirus hantanense|a",
                 "I1": "t2|Orthohantavirus hantanense|b",
                 "I2": "t2|Orthohantavirus seoulense|c",
                 "I3": "t2|Orthohantavirus dobravaense|d"})
            out = root / "out"
            written = render_phylogroup_overview(out, data_dir=str(root), segment="S")
            self.assertEqual(len(written), 1)
            self.assertTrue(written[0].exists())
            self.assertIn("overview", written[0].name)
