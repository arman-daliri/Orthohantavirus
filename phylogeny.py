#!/usr/bin/env python3
"""Orthohantavirus reassortment phylogenetics stage (MAFFT + IQ-TREE).

Stage 2 of the three-stage pipeline documented in README.md's "Pipeline stages" section.
Analyzes the reference isolate set ingested from an NCBI Datasets Orthohantavirus package
for real, tree-based reassortment evidence -- the rigorous replacement for the k-mer
divergence diagnostics all_hanta.py (Stage 1) still emits as informational-only columns.

Independent of all_hanta.py's output files: re-derives isolates, the genomic/protein FASTA
index, and the CDS index from --input itself, the same way all_hanta.py's own run_pipeline
did before this code was split out (docs/superpowers/specs/2026-08-13-three-stage-pipeline-
split-design.md has the split rationale and the exact symbol inventory this file was
extracted from). Reuses all_hanta.py's ingest primitives via `import all_hanta as ah`.
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, DefaultDict, Dict, List, Optional, Sequence, Tuple

import numpy as np

import all_hanta as ah


# Host-order fallback for isolates ah.PHYLOGROUP can't place -- species not yet formally
# named or classified (species_status != 'classified'), which have no curated reservoir
# to look up. Only usable where a WHOLE taxonomic order maps to one phylogroup: every
# Eulipotyphla genus in ah.GENUS_PHYLOGROUP is soricomorpha_borne, so an unclassified
# isolate observed in any shrew or mole is unambiguous even without a species-level
# identification. Chiroptera is included for the same reason, even though no curated
# ah.SPECIES_RISK reservoir is bat-borne -- this fallback is the only route into
# chiroptera_borne. Rodentia is deliberately absent: it spans Murinae, Arvicolinae and
# Sigmodontinae (three different phylogroups) all at once, so the order alone can't
# pick one -- see _phylogroup_for_isolate, below.
HOST_ORDER_PHYLOGROUP_FALLBACK = {
    "Eulipotyphla": "soricomorpha_borne",
    "Soricomorpha": "soricomorpha_borne",
    "Insectivora": "soricomorpha_borne",
    "Chiroptera": "chiroptera_borne",
}


@dataclass
class PhyloCandidate:
    """One isolate selected for rigorous (MAFFT + IQ-TREE) reassortment analysis.

    ``tags`` records every reason the isolate was selected (it can be more than one of
    human_affecting_or_likely / divergence_candidate / complete_triplet) so the output
    table can be filtered by which goal it serves without re-deriving the criteria.
    ``sequences`` holds the CDS nucleotide sequence actually used per segment -- never
    the full genomic sequence and never the amino acid translation, see
    select_primary_cds.
    """
    isolate: "ah.IsolateRecord"
    tags: List[str]
    sequences: Dict[str, str]


def _segment_rank_metrics(
    position: Dict[str, Any], span: int, index: int, nn: Dict[str, int],
) -> Dict[str, float]:
    """Three-segment generalization of _sl_m_rank_metrics, used only by the tree-distance
    path (tree_divergence_for_species). _sl_m_rank_metrics itself is untouched and still
    backs compute_segment_divergence's k-mer diagnostics.

    offset[X] -- how poorly X's own nearest-neighbour pick (nn[X]) ranks in the OTHER two
    segments' orderings, taking the more favourable (min) of the two. High offset[X] means
    X disagrees with both other segments about who is close to this isolate: "X is the odd
    one out". This is _sl_m_rank_metrics's m_rank_offset, generalized from X=M to any
    segment -- for X="M" it is the same value by construction (min over the other two,
    which are S and L either way).

    separation[Y,Z] -- how far apart Y and Z's own nearest-neighbour picks land in each
    other's ordering (the worse of the two directions). High separation[Y,Z] means Y and Z
    don't even agree with each other about who is close to this isolate. This is
    _sl_m_rank_metrics's sl_rank_separation, generalized from the fixed (S, L) pair to any
    pair -- for (Y, Z)=(S, L) it is the same value by construction.
    """
    segments = ("S", "M", "L")
    # rank[(other, seg)]: seg's own nearest-neighbour pick, looked up in other's ordering.
    rank = {(other, seg): float(position[other][index, nn[seg]]) / span
            for other in segments for seg in segments if other != seg}
    offset = {seg: min(rank[(other, seg)] for other in segments if other != seg)
              for seg in segments}
    separation = {
        "sm": max(rank[("M", "S")], rank[("S", "M")]),
        "sl": max(rank[("L", "S")], rank[("S", "L")]),
        "ml": max(rank[("L", "M")], rank[("M", "L")]),
    }
    return {
        "tree_s_rank_offset": offset["S"],
        "tree_m_rank_offset": offset["M"],
        "tree_l_rank_offset": offset["L"],
        "tree_sm_rank_separation": separation["sm"],
        "tree_sl_rank_separation": separation["sl"],
        "tree_ml_rank_separation": separation["ml"],
    }

# --------------------------------------------------------------------------------------
# Rigorous phylogenetic reassortment analysis (tree-based MAFFT + IQ-TREE stage)
# --------------------------------------------------------------------------------------
#
# all_hanta.py's k-mer diagnostics (assign_kmer_clusters, compute_segment_divergence) are
# alignment-free and deliberately NOT a reassortment call -- see compute_segment_divergence's docstring
# for the validation failure that got that classifier withdrawn. This section is the
# real thing it pointed to: per-segment CDS-nucleotide alignments (MAFFT) and
# maximum-likelihood trees with ultrafast bootstrap support (IQ-TREE), compared for
# topological incongruence between S, M and L.
#
# This is the phylogenetics stage (Stage 2 of the three-stage pipeline) and requires
# mafft and iqtree on PATH (or --mafft-bin/--iqtree-bin), and Biopython + dendropy
# (imported lazily below so a default run never needs them). It only ever touches the
# three named candidate groups -- human-affecting-or-likely, divergence candidates,
# complete triplets -- and every candidate isolate is reported exactly once PER
# REQUESTED SCOPE (--phylo-scope) in the output table, tagged by why it was selected
# and, if analysis didn't happen, why not. Nothing is silently dropped. Species-only
# trees (scope='species') can only ever show reassortment WITHIN a species, because a
# segment acquired from a different species has no possible donor in a single-species
# tree; phylogroup trees (scope='phylogroup') group candidates by reservoir-host
# phylogroup (ah.PHYLOGROUP) instead, spanning every species that shares one, so
# interspecies exchange has somewhere to show up (Laenen et al. 2019,
# 10.3390/v11090788). --phylo-scope=both (the default) runs both, so a candidate can
# appear up to twice, once per scope, distinguished by the analysis_scope column -- the
# old "reported exactly once, full stop" invariant now applies per scope.
#
# Setup:
#   conda env create -f environment.yml
#   conda activate hanta_phylo   # or pass --mafft-bin / --iqtree-bin with full paths

HUMAN_AFFECTING_CLASSES = {"pathogenic_hcps", "pathogenic_hfrs", "rare_or_opportunistic", "suspected"}

# CDS sequences shorter than this cannot be meaningfully aligned or placed in a tree --
# a floor against corrupt/truncated annotations, not a completeness filter.
MIN_CDS_LENGTH = 100

# Rank-separation thresholds for classify_reassortment_call. Not M-specific: any segment
# can be the odd one out (Park et al. 2025, PLOS NTD 10.1371/journal.pntd.0012859, report
# the same pattern in HTNV's L segment). PAIR_SEPARATION_LOW asks "do the two segments
# NOT under suspicion still agree with each other's neighbourhood" (should be LOW);
# SEGMENT_OFFSET_HIGH asks "is the suspect segment's answer different from the other
# two's" (should be HIGH for a genuine reassortant). Neither threshold is fitted to this
# data -- they mark the qualitative pattern the literature describes (M-segment exchange
# with S/L cosegregation, per PMC6153690, is the most-reported instance of it, not the
# only possible one), and every candidate above them still has to clear the bootstrap
# gate below.
SEGMENT_OFFSET_HIGH = 0.5
PAIR_SEPARATION_LOW = 0.35

# Which pairwise separation column gates a given segment's candidacy as the odd one out:
# if M is the suspect, S and L are the segments that must still agree with each other.
OTHER_PAIR_SEPARATION_COLUMN = {
    "S": "tree_ml_rank_separation",
    "M": "tree_sl_rank_separation",
    "L": "tree_sm_rank_separation",
}


def is_human_affecting_or_likely(pathogenicity_class: str, host_role: str) -> bool:
    """Curated pathogenic/suspected species, OR a species not yet curated but observed
    directly in a human host in this package -- the 'likely to affect humans' half of
    the brief that a static curated table alone cannot express.
    """
    if pathogenicity_class in HUMAN_AFFECTING_CLASSES:
        return True
    return pathogenicity_class == "unknown" and host_role == "spillover_host_human"


def parse_cds_fasta(path: Optional[Path]) -> Dict[str, List[Tuple[str, str]]]:
    """Index cds.fna by NUCLEOTIDE accession -- not protein accession.

    Unlike protein.faa (keyed by protein accessionVersion), cds.fna's defline carries the
    parent nucleotide accession plus a coordinate suffix (e.g. 'PX492530.1:1-430'),
    because a CDS is a coordinate range on its parent nucleotide record. A single
    nucleotide accession can carry more than one CDS -- 134 in this package, mostly
    ambisense S segments encoding both N and NSs on overlapping reading frames, plus a
    few M segments annotated with separate G1/G2 CDS alongside a combined precursor
    entry. parse_fasta's setdefault-first-wins index would silently keep whichever CDS
    happened to be read first and drop the other, so this keeps every distinct
    (product, sequence) pair per accession and leaves selection to select_primary_cds.
    """
    index: DefaultDict[str, List[Tuple[str, str]]] = defaultdict(list)
    seen: set = set()
    for desc, seq in ah._iter_fasta_records(path):
        if not desc:
            continue
        tokens = desc.split(None, 1)
        acc = tokens[0].split(":")[0]
        base = ah.accession_base(acc)
        product = tokens[1].split("[")[0].strip() if len(tokens) > 1 else ""
        seq = seq.upper()
        if not seq or not base:
            continue
        dedupe_key = (base, product, seq)
        if dedupe_key in seen:
            continue
        seen.add(dedupe_key)
        index[base].append((product, seq))
    return dict(index)


def select_primary_cds(records: Optional[Sequence[Tuple[str, str]]], segment: str) -> str:
    """Pick the CDS nucleotide sequence matching the segment's primary product.

    Mirrors select_primary_protein's role-matching (ah.classify_protein_role), longest-wins
    among matches -- for M segments annotated with separate G1/G2 CDS *and* a combined
    G1+G2 precursor entry, this correctly prefers the full precursor over either half.
    """
    if not records:
        return ""
    wanted = ah.SEGMENT_PRIMARY_PRODUCT.get(segment)
    if not wanted:
        return ""
    candidates = [seq for product, seq in records if ah.classify_protein_role(product) == wanted]
    return max(candidates, key=len) if candidates else ""


def collect_phylo_candidates(
    isolates: Sequence[ah.IsolateRecord],
    divergence_by_key: Dict[str, ah.SegmentDivergenceProfile],
    cds_index: Dict[str, List[Tuple[str, str]]],
    scope: str,
) -> Tuple[List[PhyloCandidate], Counter]:
    """Extract the three named candidate groups and require full S+M+L CDS for each.

    The three groups -- human-affecting-or-likely, divergence candidates (already
    flagged 'requires_phylogenetic_confirmation' by the k-mer diagnostic), complete
    triplets -- are UNIONED, not intersected: each is independently a reason this
    isolate matters for the reassortment/pathogenicity goal. But a topology-incongruence
    test needs every isolate placed in all three segment trees, so 'CDS available for S,
    M and L' is a hard requirement on top of the union, not an optional nicety.

    Lab-derived and vaccine-strain isolates are excluded here specifically because a
    passage/attenuation history is a tree-topology risk (branch-length distortion from
    cell-culture or lab-host adaptation can bias a reassortment call). This is a
    different reason than build_ml_tables now uses -- that function keeps is_lab_derived
    isolates (flagged, not filtered) because the ML task doesn't depend on branch
    lengths, and blanket-excluding them there discarded real pathogenicity signal (see
    its docstring). Only is_vaccine_strain is excluded from both stages.

    ``scope`` governs how unclassified isolates (species_status != 'classified') are
    treated -- this function is called once per requested --phylo-scope, not once total,
    specifically so the two scopes can disagree here:

    scope='species': excluded outright, same as before scope existed. 'unclassified
    Orthohantavirus' and genus-only records are a catch-all bucket, not one coherent
    taxon -- pooling them under species_name would build a tree across genuinely
    different, unrelated lineages that just happen to share a placeholder label.

    scope='phylogroup': admitted if _phylogroup_for_isolate can place them by observed
    host order (see HOST_ORDER_PHYLOGROUP_FALLBACK). These are exactly the novel and
    divergent lineages -- unnamed species, genus-only records -- most relevant to
    emergence, and species-scope's blanket exclusion meant they never entered a tree at
    all, on EITHER scope, even though pooling by phylogroup rather than by species name
    doesn't have the "unrelated lineages sharing a label" problem. An isolate that still
    can't be placed (host order not in HOST_ORDER_PHYLOGROUP_FALLBACK, or absent) is
    excluded with its own reason, same as every other exclusion here.
    """
    exclusions: Counter = Counter()
    candidates: List[PhyloCandidate] = []
    for isolate in isolates:
        if isolate.is_lab_derived or isolate.is_vaccine_strain:
            exclusions["lab_or_vaccine_excluded"] += 1
            continue
        if isolate.species_status != "classified":
            if scope == "species":
                exclusions["species_not_classified"] += 1
                continue
            _, basis = _phylogroup_for_isolate(isolate)
            if basis == "unassigned":
                exclusions["species_not_classified_and_no_resolvable_phylogroup"] += 1
                continue
            # else: basis == 'observed_host_order' -- admitted, falls through below
            # exactly like a classified isolate.

        risk = ah.species_risk_profile(isolate.species_name)
        tags: List[str] = []
        if is_human_affecting_or_likely(risk["pathogenicity_class"], isolate.host_role):
            tags.append("human_affecting_or_likely")
        profile = divergence_by_key.get(isolate.isolate_key)
        if profile is not None and profile.reassortment_status == "requires_phylogenetic_confirmation":
            tags.append("divergence_candidate")
        if isolate.has_complete_triplet:
            tags.append("complete_triplet")
        if not tags:
            continue

        sequences: Dict[str, str] = {}
        for segment, accession in (("S", isolate.s_accession), ("M", isolate.m_accession),
                                   ("L", isolate.l_accession)):
            records = cds_index.get(ah.accession_base(accession)) if accession else None
            seq = select_primary_cds(records, segment)
            if len(seq) >= MIN_CDS_LENGTH:
                sequences[segment] = seq

        if len(sequences) < 3:
            exclusions["missing_cds_for_one_or_more_segments"] += 1
            continue

        candidates.append(PhyloCandidate(isolate=isolate, tags=tags, sequences=sequences))

    return candidates, exclusions


def _phylogroup_for_isolate(isolate: ah.IsolateRecord) -> Tuple[str, str]:
    """Phylogroup and assignment basis for one isolate at scope='phylogroup'.

    Classified species use the curated ah.PHYLOGROUP lookup (basis='curated_reservoir') --
    every ah.SPECIES_RISK entry resolves here, so this branch never falls through to
    'unassigned' for a classified isolate.

    Unclassified isolates (species_status != 'classified') have no curated reservoir to
    look up, so they fall back to their own OBSERVED host order
    (basis='observed_host_order'). That fallback only resolves an order that maps
    entirely to one phylogroup in this dataset -- true for Eulipotyphla (shrews/moles,
    all soricomorpha_borne) and Chiroptera (bats, all chiroptera_borne), see
    HOST_ORDER_PHYLOGROUP_FALLBACK -- but NOT for Rodentia, which spans Murinae,
    Arvicolinae and Sigmodontinae all at once. A Rodentia-hosted (or host-unknown)
    unclassified isolate therefore has no fallback and stays 'unassigned': host order
    alone cannot tell which rodent subfamily it belongs to, and guessing one would
    misrepresent where the isolate actually sits, which is worse than excluding it.
    """
    curated = ah.PHYLOGROUP.get(isolate.species_name)
    if curated is not None:
        return curated, "curated_reservoir"
    fallback = HOST_ORDER_PHYLOGROUP_FALLBACK.get(isolate.host_order)
    if fallback is not None:
        return fallback, "observed_host_order"
    return "unassigned", "unassigned"


def _phylogroup_for_candidate(candidate: PhyloCandidate) -> Tuple[str, str]:
    return _phylogroup_for_isolate(candidate.isolate)


def _assignment_basis(candidate: PhyloCandidate, scope: str) -> str:
    """phylogroup_assignment_basis for one candidate's output row -- empty at
    scope='species', where the concept doesn't apply.
    """
    return _phylogroup_for_candidate(candidate)[1] if scope == "phylogroup" else ""


def _candidate_group_key(candidate: PhyloCandidate, scope: str) -> str:
    if scope == "species":
        return candidate.isolate.species_name
    if scope == "phylogroup":
        # Every candidate reaching this point already has a resolvable phylogroup --
        # collect_phylo_candidates excluded anything that would resolve to 'unassigned'
        # before it became a candidate at all -- so the [0] here never actually returns
        # 'unassigned' in practice.
        return _phylogroup_for_candidate(candidate)[0]
    raise ValueError(f"unknown phylogenetics scope: {scope!r}")


def group_candidates(
    candidates: Sequence[PhyloCandidate], min_taxa: int, scope: str,
) -> Tuple[Dict[str, List[PhyloCandidate]], Counter]:
    """Partition candidates into tree groups at the given scope. Only groups with enough
    candidates for a meaningful tree get analyzed.

    scope='species': one group per canonical species. Can only ever show reassortment
    WITHIN a species -- a segment acquired from a different species has no possible donor
    in a single-species tree.

    scope='phylogroup': one group per reservoir-host phylogroup (ah.PHYLOGROUP), spanning
    every species that shares one. Makes interspecies reassortment visible, which
    species-scope structurally cannot (Laenen et al. 2019, 10.3390/v11090788). These
    groups are larger and taxonomically heterogeneous, so unlike species trees they need
    an explicit outgroup to be interpretable -- see --phylo-outgroup in
    run_phylogenetics_stage.

    An unrooted tree needs at least 4 taxa to have any internal structure at all; below
    min_taxa there is nothing for MAFFT/IQ-TREE to resolve. Sorted by isolate_key so
    re-runs assign identical tip labels and produce byte-identical alignments/trees given
    the same inputs.
    """
    grouped: DefaultDict[str, List[PhyloCandidate]] = defaultdict(list)
    for candidate in candidates:
        grouped[_candidate_group_key(candidate, scope)].append(candidate)

    exclusions: Counter = Counter()
    accepted: Dict[str, List[PhyloCandidate]] = {}
    for key, members in grouped.items():
        if len(members) < min_taxa:
            # Species-scope message text is unchanged from the pre-scope version (no
            # "species:" prefix) so a species-only run's log output stays byte-identical
            # to before this function existed (this stage's own log, ortho_phylogeny.log,
            # not all_hanta.py's ortho_pipeline.log -- they're separate files).
            label = key if scope == "species" else f"phylogroup {key}"
            exclusions[f"{label}: only {len(members)} candidate(s), need {min_taxa}"] += 1
            continue
        accepted[key] = sorted(members, key=lambda c: c.isolate.isolate_key)
    return accepted, exclusions


def group_candidates_by_species(
    candidates: Sequence[PhyloCandidate], min_taxa: int,
) -> Tuple[Dict[str, List[PhyloCandidate]], Counter]:
    """Thin species-scope wrapper, kept so call sites that only ever want species-level
    grouping don't need to know the scope parameter exists.
    """
    return group_candidates(candidates, min_taxa, scope="species")


def resolve_binary(explicit: Optional[str], name: str) -> str:
    """Resolve a required external binary to an absolute path, or fail with setup help."""
    candidate = explicit or shutil.which(name)
    resolved = shutil.which(candidate) if candidate else None
    if not resolved:
        raise SystemExit(
            f"Phylogenetics stage requires '{name}' on PATH (or pass --{name}-bin=/path/to/{name}).\n"
            "Install dependencies with:\n"
            "  conda env create -f environment.yml\n"
            f"then either 'conda activate hanta_phylo' before running this script, or pass "
            f"--{name}-bin pointing at that environment's bin/{name}."
        )
    return resolved


def write_segment_fasta(path: Path, sequences: Dict[str, str]) -> None:
    """Write label -> nucleotide sequence. Labels are already synthetic (I0000, I0001,
    ...), never isolate_key directly -- isolate_key contains '|' and other punctuation
    that is legal in FASTA headers but risks corruption in Newick tip names downstream,
    so a synthetic label sidesteps the question entirely rather than trying to escape it.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for label, seq in sequences.items():
            handle.write(f">{label}\n")
            for i in range(0, len(seq), 70):
                handle.write(seq[i:i + 70] + "\n")


def run_mafft(mafft_bin: str, input_fasta: Path, output_fasta: Path, threads: int) -> None:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8") as out:
        result = subprocess.run(
            [mafft_bin, "--auto", "--thread", str(threads), str(input_fasta)],
            stdout=out, stderr=subprocess.PIPE, text=True, check=False,
        )
    if result.returncode != 0 or not output_fasta.exists() or output_fasta.stat().st_size == 0:
        raise RuntimeError(f"mafft failed on {input_fasta}: {(result.stderr or '')[-2000:]}")


def run_iqtree(
    iqtree_bin: str, alignment_fasta: Path, prefix: Path, threads: int, bootstrap: int, seed: int,
    outgroup: Optional[str] = None,
) -> Path:
    cmd = [
        # --seqtype DNA is explicit, not left to IQ-TREE's auto-detection: these are CDS
        # nucleotide sequences by construction (select_primary_cds), and heavily-N'd
        # partial CDS can otherwise make the ACGT-vs-amino-acid heuristic misfire.
        iqtree_bin, "-s", str(alignment_fasta), "--seqtype", "DNA", "-m", "MFP", "-B", str(bootstrap),
        "-T", str(threads), "--seed", str(seed), "--prefix", str(prefix), "--quiet", "--redo",
    ]
    if outgroup:
        cmd += ["-o", outgroup]
    prefix.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    treefile = Path(f"{prefix}.treefile")
    if result.returncode != 0 or not treefile.exists():
        raise RuntimeError(
            f"iqtree failed on {alignment_fasta}: {(result.stderr or result.stdout or '')[-2000:]}")
    return treefile


def _patristic_distance_matrix(tree: Any, labels: Sequence[str]) -> Any:
    """All-pairs patristic (branch-length, i.e. real evolutionary) distance."""
    clades = {clade.name: clade for clade in tree.get_terminals()}
    n = len(labels)
    dist = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            d = tree.distance(clades[labels[i]], clades[labels[j]])
            dist[i, j] = dist[j, i] = d
    return dist


def tree_divergence_for_species(labels: Sequence[str], trees: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Per-isolate nearest-neighbour and rank-separation diagnostics from REAL segment
    trees, replacing compute_segment_divergence's k-mer/cosine version with the same
    statistical shape (same span/rank normalization) but on maximum-likelihood
    phylogenetic distance instead of alignment-free composition distance -- the actual
    fix for the failure mode documented there.
    """
    n = len(labels)
    span = max(1, n - 1)
    distance: Dict[str, Any] = {}
    order: Dict[str, Any] = {}
    for segment, tree in trees.items():
        dist = _patristic_distance_matrix(tree, labels)
        np.fill_diagonal(dist, np.inf)
        order[segment] = np.argsort(dist, axis=1)
        distance[segment] = dist

    position = ah._rank_positions(order)

    rows_out: List[Dict[str, Any]] = []
    for i, label in enumerate(labels):
        s_nn, m_nn, l_nn = order["S"][i, 0], order["M"][i, 0], order["L"][i, 0]
        metrics = _segment_rank_metrics(position, span, i, {"S": s_nn, "M": m_nn, "L": l_nn})
        rows_out.append({
            "label": label,
            "tree_s_nearest": labels[s_nn], "tree_m_nearest": labels[m_nn], "tree_l_nearest": labels[l_nn],
            "tree_s_distance": float(distance["S"][i, s_nn]),
            "tree_m_distance": float(distance["M"][i, m_nn]),
            "tree_l_distance": float(distance["L"][i, l_nn]),
            **{key: round(value, 4) for key, value in metrics.items()},
        })
    return rows_out


def _nearest_supported_descendant(node: Any, target_name: str) -> Optional[float]:
    """Walk from `node` down toward the tip named `target_name`, returning the first
    confidence value found (None if the walk reaches a tip without ever finding one).
    Helper for bootstrap_support_for_pair's root fallback -- see there for why this walk
    is needed at all.
    """
    while True:
        children = [c for c in node.clades if any(t.name == target_name for t in c.get_terminals())]
        if not children:
            return None
        node = children[0]
        if node.confidence is not None:
            return float(node.confidence)
        if node.is_terminal():
            return None


def bootstrap_support_for_pair(tree: Any, label_a: str, label_b: str) -> Tuple[float, str]:
    """Ultrafast-bootstrap support at the smallest clade containing both tips -- how well
    the discordant segment's own tree supports grouping this isolate with the neighbour
    it disagrees with the other two segments about. IQ-TREE writes support as the
    internal-node name, which Biopython's Newick parser reads into Clade.confidence for
    non-terminal clades.

    Returns (support, source). source is 'mrca' in the overwhelming majority of calls:
    every internal node in an IQ-TREE treefile carries a UFBoot value EXCEPT the root --
    an unrooted tree's root is not a real bipartition, so IQ-TREE writes nothing there,
    and Biopython reads that absence as confidence=None. The previous version of this
    function returned 0.0 whenever the pair's MRCA happened to BE that root, which reads
    as "definitively no support" and is indistinguishable from a genuine 0% UFBoot value
    -- a real reassortant whose discordant-segment neighbour sits at the tree's deepest
    split got silently zeroed out and demoted to weak_signal for a reason that has
    nothing to do with how well-supported the grouping actually is.

    When the MRCA is unsupported, this walks from it toward each tip in turn (there is no
    smaller clade containing BOTH tips once the MRCA is unsupported -- by definition
    that only happens at the root -- so "toward the tips" necessarily means down into the
    two separate subtrees that lead to label_a and to label_b) and takes the MORE
    CONSERVATIVE (lower) of whichever supported ancestors it finds on each side --
    source='nearest_supported_ancestor'. Consistent with the bootstrap gate's whole
    purpose: don't let a lucky pick of the better-supported side overstate the evidence.
    If neither walk ever finds a support value (only possible in a tree too small to have
    any internal structure between the root and its tips), returns
    (0.0, 'unsupported_root') -- now an honest "not measured", not a disguised "measured
    and it was zero".
    """
    clade_a = next(tree.find_clades(name=label_a))
    clade_b = next(tree.find_clades(name=label_b))
    mrca = tree.common_ancestor([clade_a, clade_b])
    if mrca.confidence is not None:
        return float(mrca.confidence), "mrca"

    candidates = [s for s in (_nearest_supported_descendant(mrca, label_a),
                              _nearest_supported_descendant(mrca, label_b)) if s is not None]
    if candidates:
        return min(candidates), "nearest_supported_ancestor"
    return 0.0, "unsupported_root"


def classify_reassortment_call(
    row: Dict[str, Any],
    trees: Dict[str, Any],
    min_bootstrap_support: float,
) -> Tuple[str, float, str, str]:
    """Identify the discordant segment (if any), then require the OTHER two to still
    cosegregate, then gate on bootstrap support in the discordant segment's own tree.

    Generalizes the old M-only check: that version could only ever find the PMC6153690
    M-exchange/S-L-cosegregation pattern, so an isolate whose S and M trees agree but
    whose L tree disagrees -- the pattern Park et al. 2025 (PLOS NTD,
    10.1371/journal.pntd.0012859) report for three candidate HTNV reassortants --
    returned 'no_conflict' unconditionally. The M-permissive prior from PMC6153690 is
    preserved downstream as matches_canonical_m_exchange, not as a filter here.

    The discordant segment is whichever has the highest rank offset (ties broken toward
    S, then M, then L, by dict iteration order -- deterministic, and ties are already rare
    since offsets are continuous). Everything else stays 'no_conflict' or 'weak_signal' --
    the whole point of doing this with real trees is to not repeat the withdrawn
    classifier's mistake of reporting an unsupported pattern as a finding.

    Returns (call, support, discordant_segment, bootstrap_support_source).
    discordant_segment is "" for no_conflict, and set (even if the call is weak_signal)
    whenever a rank pattern was found -- the pattern existed structurally regardless of
    whether bootstrap confirmed it. bootstrap_support_source is "" for no_conflict (the
    bootstrap gate is never reached), otherwise see bootstrap_support_for_pair.
    """
    offsets = {"S": row["tree_s_rank_offset"], "M": row["tree_m_rank_offset"], "L": row["tree_l_rank_offset"]}
    discordant = max(offsets, key=offsets.get)
    other_pair_separation = row[OTHER_PAIR_SEPARATION_COLUMN[discordant]]
    if offsets[discordant] < SEGMENT_OFFSET_HIGH or other_pair_separation > PAIR_SEPARATION_LOW:
        return "no_conflict", 0.0, "", ""
    tree = trees[discordant]
    nearest = row[f"tree_{discordant.lower()}_nearest"]
    support, support_source = bootstrap_support_for_pair(tree, row["label"], nearest)
    if support >= min_bootstrap_support:
        return "supported_reassortment_candidate", support, discordant, support_source
    return "weak_signal_low_bootstrap_support", support, discordant, support_source


def _phylo_row(candidate: PhyloCandidate, species: str, n_taxa_in_tree: int,
              status: str, error: str, scope: str, phylogroup_assignment_basis: str = "") -> Dict[str, Any]:
    """Default (not-yet-analyzed / failed) row shape. Every field a successful analysis
    fills in is declared here too, so ortho_reassortment_table.csv has one stable
    schema regardless of status -- required for a clean, joinable table.

    ``species`` is always the isolate's own true canonical species -- never the group
    label a phylogroup-scope tree was built from, since a phylogroup tree spans several
    species. Which phylogroup that species belongs to is recoverable from ah.PHYLOGROUP.

    ``phylogroup_assignment_basis`` is only meaningful at scope='phylogroup' (empty
    string at scope='species', where the concept doesn't apply): 'curated_reservoir' for
    a classified species looked up in ah.PHYLOGROUP, 'observed_host_order' for an
    unclassified isolate placed via HOST_ORDER_PHYLOGROUP_FALLBACK instead -- see
    _phylogroup_for_isolate.
    """
    return {
        "isolate_key": candidate.isolate.isolate_key,
        "species_name": species,
        "candidate_tags": ";".join(candidate.tags),
        "analysis_scope": scope,
        "phylogroup_assignment_basis": phylogroup_assignment_basis,
        # Name predates scope='phylogroup' and is now imprecise for those rows (the tree
        # spans several species, not one) -- kept as-is rather than renamed, so a
        # species-only run's column set stays identical to before this task. Read as
        # "n_taxa_in_this_row's_tree" regardless of scope.
        "n_taxa_in_species_tree": n_taxa_in_tree,
        "phylo_status": status,
        "phylo_error": error,
        "tree_s_nearest_isolate": "", "tree_m_nearest_isolate": "", "tree_l_nearest_isolate": "",
        "tree_s_distance": 0.0, "tree_m_distance": 0.0, "tree_l_distance": 0.0,
        "tree_s_rank_offset": 0.0, "tree_m_rank_offset": 0.0, "tree_l_rank_offset": 0.0,
        "tree_sm_rank_separation": 0.0, "tree_sl_rank_separation": 0.0, "tree_ml_rank_separation": 0.0,
        "discordant_segment": "",
        "matches_canonical_m_exchange": 0,
        # None, not 0.0: a genuine 0% bootstrap support is a real (weak) measurement, and
        # must stay distinguishable from "no discordant clade exists to measure at all"
        # (not_analyzed/failed rows, or an analyzed row that called no_conflict). Same
        # reasoning m_clade_bootstrap_support inherits for "M wasn't the discordant segment".
        "m_clade_bootstrap_support": None,
        "discordant_clade_bootstrap_support": None,
        # "" for not-yet-analyzed rows and for no_conflict (bootstrap gate never reached);
        # 'mrca' / 'nearest_supported_ancestor' / 'unsupported_root' otherwise -- see
        # bootstrap_support_for_pair.
        "bootstrap_support_source": "",
        "reassortment_call": "not_analyzed",
        # species_* here means "constant across every isolate row sharing this group's
        # tree" -- under scope='phylogroup' that's a phylogroup, not literally one
        # species, same imprecision n_taxa_in_species_tree already carries and for the
        # same reason (not renamed further to avoid a second unrequested rename).
        "species_rf_s_vs_m": None, "species_rf_m_vs_l": None, "species_rf_s_vs_l": None,
        "species_rf_s_vs_m_normalized": None, "species_rf_m_vs_l_normalized": None,
        "species_rf_s_vs_l_normalized": None,
        "s_alignment_path": "", "m_alignment_path": "", "l_alignment_path": "",
        "s_tree_path": "", "m_tree_path": "", "l_tree_path": "",
    }


def _resolve_outgroup_label(label_map: Dict[str, "PhyloCandidate"], outgroup: str) -> Optional[str]:
    """Map --phylo-outgroup (an accession or isolate_key) to the synthetic tip label
    (I0000, ...) used in THIS group's alignment/tree, or None if that isolate isn't one
    of this group's members -- one global outgroup will not be present in every tree.
    """
    base = ah.accession_base(outgroup)
    for label, candidate in label_map.items():
        isolate = candidate.isolate
        if outgroup == isolate.isolate_key:
            return label
        if base and base in {ah.accession_base(a) for a in
                             (isolate.s_accession, isolate.m_accession, isolate.l_accession) if a}:
            return label
    return None


def normalized_rf(raw_rf: int, n_taxa: int) -> float:
    """raw_rf as a fraction of the maximum possible Robinson-Foulds distance for n_taxa,
    clamped to at most 1.0.

    2*(n-3) is the maximum RF distance between two STRICTLY BIFURCATING unrooted trees
    on n taxa. IQ-TREE can emit a multifurcating tree on ties/polytomies -- a real,
    non-error output on some data, not a crash -- and the raw RF distance between a
    multifurcating tree and a fully-resolved one can then exceed this bound. Clamped
    here rather than left to silently read as "more disagreement than two trees can
    possibly have"; the raw, un-normalized species_rf_*_vs_* column still carries the
    true value regardless.
    """
    max_rf = 2 * (n_taxa - 3) if n_taxa > 3 else 1
    return min(1.0, raw_rf / max_rf) if max_rf else 0.0


def _analyze_group(
    group_key: str, scope: str, members: Sequence[PhyloCandidate], group_dir_root: Path, output_dir: Path,
    args: argparse.Namespace, mafft_bin: str, iqtree_bin: str,
) -> List[Dict[str, Any]]:
    """Align + build trees for one group's (species or phylogroup) three segments, then
    score every member.

    Every step -- MAFFT/IQ-TREE, tree parsing, patristic-distance/Robinson-Foulds
    computation, and per-isolate classification -- runs inside one try/except: a failure
    anywhere fails only this group (isolates reported with status='failed' and the
    reason) and never the run, the same per-unit error isolation transform_one_record
    uses for accessions. An earlier version of this function wrapped only the MAFFT/
    IQ-TREE calls, so a failure in the tree-parsing or distance code that followed them
    (e.g. IQ-TREE dropping a duplicate tip, a malformed treefile) escaped uncaught.

    ``group_key`` is the species name at scope='species' or the phylogroup name at
    scope='phylogroup' -- it names the tree, not any one member's own species (members of
    a phylogroup tree can each have a different species_name; see _phylo_row).
    """
    # Re-imported here rather than threaded through as parameters: run_phylogenetics_stage
    # already proved these are importable (its own ImportError check, with the install
    # hint, runs before any work starts), and Python caches modules in sys.modules, so
    # this costs nothing per group.
    from Bio import Phylo
    import dendropy
    from dendropy.calculate import treecompare

    group_dir = group_dir_root / ah.normalize_attr_name(group_key)
    labels = [f"I{i:04d}" for i in range(len(members))]
    label_map: Dict[str, PhyloCandidate] = dict(zip(labels, members))

    outgroup_label = _resolve_outgroup_label(label_map, args.phylo_outgroup) if args.phylo_outgroup else None
    if args.phylo_outgroup and outgroup_label is None:
        logging.info("Phylogenetics (%s): --phylo-outgroup %s is not a member of group %s, "
                    "its tree(s) are left unrooted", scope, args.phylo_outgroup, group_key)

    try:
        trees: Dict[str, Any] = {}
        tree_paths: Dict[str, Path] = {}
        for segment in ("S", "M", "L"):
            seg_dir = group_dir / segment
            fasta_path = seg_dir / "sequences.fasta"
            aligned_path = seg_dir / "aligned.fasta"
            prefix = seg_dir / "tree"
            sequences = {label: candidate.sequences[segment] for label, candidate in label_map.items()}
            write_segment_fasta(fasta_path, sequences)
            run_mafft(mafft_bin, fasta_path, aligned_path, args.phylo_threads)
            treefile = run_iqtree(iqtree_bin, aligned_path, prefix, args.phylo_threads,
                                  args.bootstrap, args.phylo_seed, outgroup=outgroup_label)
            trees[segment] = Phylo.read(str(treefile), "newick")
            tree_paths[segment] = treefile

        divergence_by_label = {row["label"]: row for row in tree_divergence_for_species(labels, trees)}

        # Robinson-Foulds between each pair of segment trees on the shared taxon set
        # (identical by construction) -- a per-group overall measure of how much the three
        # trees disagree, independent of any single isolate's placement.
        taxon_namespace = dendropy.TaxonNamespace()
        dendro_trees = {
            segment: dendropy.Tree.get(path=str(tree_paths[segment]), schema="newick",
                                       taxon_namespace=taxon_namespace)
            for segment in ("S", "M", "L")
        }
        for dtree in dendro_trees.values():
            dtree.encode_bipartitions()
        rf: Dict[str, Any] = {}
        for a, b in (("s", "m"), ("m", "l"), ("s", "l")):
            raw = treecompare.symmetric_difference(dendro_trees[a.upper()], dendro_trees[b.upper()])
            rf[f"species_rf_{a}_vs_{b}"] = raw
            rf[f"species_rf_{a}_vs_{b}_normalized"] = normalized_rf(raw, len(labels))

        rows: List[Dict[str, Any]] = []
        for label, candidate in label_map.items():
            row = divergence_by_label[label]
            call, support, discordant, support_source = classify_reassortment_call(
                row, trees, args.min_bootstrap_support)
            out_row = _phylo_row(candidate, candidate.isolate.species_name, len(members), "analyzed", "", scope,
                                _assignment_basis(candidate, scope))
            out_row.update({
                "tree_s_nearest_isolate": label_map[row["tree_s_nearest"]].isolate.isolate_key,
                "tree_m_nearest_isolate": label_map[row["tree_m_nearest"]].isolate.isolate_key,
                "tree_l_nearest_isolate": label_map[row["tree_l_nearest"]].isolate.isolate_key,
                "tree_s_distance": row["tree_s_distance"],
                "tree_m_distance": row["tree_m_distance"],
                "tree_l_distance": row["tree_l_distance"],
                "tree_s_rank_offset": row["tree_s_rank_offset"],
                "tree_m_rank_offset": row["tree_m_rank_offset"],
                "tree_l_rank_offset": row["tree_l_rank_offset"],
                "tree_sm_rank_separation": row["tree_sm_rank_separation"],
                "tree_sl_rank_separation": row["tree_sl_rank_separation"],
                "tree_ml_rank_separation": row["tree_ml_rank_separation"],
                "discordant_segment": discordant,
                "matches_canonical_m_exchange": int(discordant == "M"),
                # None (not 0.0) when there's no discordant clade to measure at all
                # (no_conflict, discordant == ""), so a genuine 0% bootstrap support stays
                # distinguishable from "not applicable". m_clade_bootstrap_support is kept
                # for existing consumers, populated only when M is actually the discordant
                # segment -- see discordant_clade_bootstrap_support for the generalized
                # value regardless of which segment it is.
                "m_clade_bootstrap_support": support if discordant == "M" else None,
                "discordant_clade_bootstrap_support": support if discordant else None,
                "bootstrap_support_source": support_source,
                "reassortment_call": call,
                **rf,
                "s_alignment_path": str((group_dir / "S" / "aligned.fasta").relative_to(output_dir)),
                "m_alignment_path": str((group_dir / "M" / "aligned.fasta").relative_to(output_dir)),
                "l_alignment_path": str((group_dir / "L" / "aligned.fasta").relative_to(output_dir)),
                "s_tree_path": str(tree_paths["S"].relative_to(output_dir)),
                "m_tree_path": str(tree_paths["M"].relative_to(output_dir)),
                "l_tree_path": str(tree_paths["L"].relative_to(output_dir)),
            })
            rows.append(out_row)
        return rows
    except Exception as exc:  # noqa: BLE001 - one group's failure must not kill the run
        logging.exception("Phylogenetics failed for %s (%s scope)", group_key, scope)
        error = f"{type(exc).__name__}:{exc}"
        return [_phylo_row(c, c.isolate.species_name, len(members), "failed", error, scope,
                          _assignment_basis(c, scope)) for c in members]


def _run_phylogenetics_scope(
    scope: str, candidates: Sequence[PhyloCandidate], output_dir: Path, args: argparse.Namespace,
    mafft_bin: str, iqtree_bin: str,
) -> List[Dict[str, Any]]:
    """One scope's worth of grouping + tree-building + classification. Factored out of
    run_phylogenetics_stage so --phylo-scope=both is just this called twice; the
    scope='species' call alone is exactly what run_phylogenetics_stage did before this
    function existed.
    """
    by_group, group_drop_reasons = group_candidates(candidates, args.phylo_min_taxa, scope)
    for reason in group_drop_reasons:
        logging.info("Phylogenetics: excluded group below --phylo-min-taxa: %s", reason)
    logging.info("Phylogenetics (%s): analyzing %s group(s): %s", scope, len(by_group), sorted(by_group))

    grouped_keys = {c.isolate.isolate_key for members in by_group.values() for c in members}
    candidates_per_group: DefaultDict[str, int] = defaultdict(int)
    for c in candidates:
        candidates_per_group[_candidate_group_key(c, scope)] += 1
    rows: List[Dict[str, Any]] = [
        _phylo_row(c, c.isolate.species_name, candidates_per_group[_candidate_group_key(c, scope)],
                  "species_below_min_taxa", "", scope, _assignment_basis(c, scope))
        for c in candidates if c.isolate.isolate_key not in grouped_keys
    ]

    # Phylogroup output lives under a sibling directory, never inside phylogenetics/, so a
    # species-scope-only run's directory tree (and therefore the alignment/tree path
    # columns it writes) stays exactly what it was before scope existed.
    group_dir_root = output_dir / "phylogenetics" if scope == "species" else output_dir / "phylogenetics_phylogroup"
    for group_key, members in sorted(by_group.items()):
        logging.info("Phylogenetics (%s): %s (%s taxa) -- aligning + building trees",
                     scope, group_key, len(members))
        rows.extend(_analyze_group(
            group_key, scope, members, group_dir_root, output_dir, args, mafft_bin, iqtree_bin))
    return rows


def run_phylogenetics_stage(
    isolates: Sequence[ah.IsolateRecord],
    divergence_by_key: Dict[str, ah.SegmentDivergenceProfile],
    cds_index: Dict[str, List[Tuple[str, str]]],
    output_dir: Path,
    args: argparse.Namespace,
) -> List[Dict[str, Any]]:
    """Entry point for the phylogenetics (tree-based) analysis stage.

    Runs once per scope named by --phylo-scope ('both', the default, runs species then
    phylogroup). Each candidate is reported once per scope it was analyzed under, so at
    --phylo-scope=both a candidate can appear up to twice in the returned rows, tagged by
    analysis_scope.
    """
    try:
        import Bio.Phylo  # noqa: F401 - fail fast with the install hint below, before any work starts
        import dendropy  # noqa: F401
        import dendropy.calculate.treecompare  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "Phylogenetics stage requires biopython and dendropy. Install with:\n"
            "  conda env create -f environment.yml\n"
            f"(missing: {exc.name})"
        ) from exc

    mafft_bin = resolve_binary(args.mafft_bin, "mafft")
    iqtree_bin = resolve_binary(args.iqtree_bin, "iqtree")

    scopes = ("species", "phylogroup") if args.phylo_scope == "both" else (args.phylo_scope,)
    if "phylogroup" in scopes and not args.phylo_outgroup:
        logging.warning(
            "Phylogroup trees are unrooted (no --phylo-outgroup given): they can show "
            "THAT segments from different species cluster together, but not which "
            "species was donor and which was recipient.")

    rows: List[Dict[str, Any]] = []
    for scope in scopes:
        # Candidate collection is scope-dependent (Task 3): unclassified isolates are
        # excluded outright at scope='species' but admitted at scope='phylogroup' if
        # observed host order places them -- see collect_phylo_candidates. So this runs
        # once per scope, not once total.
        candidates, drop_reasons = collect_phylo_candidates(isolates, divergence_by_key, cds_index, scope)
        logging.info("Phylogenetics (%s): %s candidate isolates (%s excluded before grouping: %s)",
                     scope, len(candidates), sum(drop_reasons.values()), dict(drop_reasons))
        rows.extend(_run_phylogenetics_scope(scope, candidates, output_dir, args, mafft_bin, iqtree_bin))
    return rows


# Ranks reassortment_call outcomes so phylo_call_by_isolate can pick the most informative
# one when an isolate has more than one row (up to one per --phylo-scope).
_REASSORTMENT_CALL_RANK = {
    "supported_reassortment_candidate": 2,
    "weak_signal_low_bootstrap_support": 1,
    "no_conflict": 0,
}


def phylo_call_by_isolate(phylo_rows: Sequence[Dict[str, Any]]) -> Dict[str, str]:
    """Reduce ortho_reassortment_table.csv's rows to one reassortment_call per isolate.

    Public utility, not called anywhere in this pipeline's own scripts: before the
    three-stage split, this fed all_hanta.py's surveillance_priority(phylo_call=...)
    directly, in the same process. Post-split, ortho_risk_table.csv (Stage 1) is
    unconditionally phylo-blind by design (see README's Pipeline stages section and
    docs/superpowers/specs/2026-08-13-three-stage-pipeline-split-design.md section 5),
    so nothing in Stage 1 or Stage 2 calls this anymore. It's kept, public, and tested
    for exactly the manual isolate_key join the README documents -- call it yourself on
    ortho_reassortment_table.csv's rows if you want this same scope-reduction before
    joining onto ortho_risk_table.csv or your own analysis.

    At --phylo-scope=both an isolate can have up to two rows, one per scope, and they can
    disagree WITHOUT contradicting each other: species-only and phylogroup-only trees
    answer different questions (within-species vs across-species topology), so one
    scope's 'no_conflict' never invalidates the other's 'supported_reassortment_candidate'
    -- it just means that scope didn't find (or couldn't look for) the same signal. So
    this takes the single most informative call across an isolate's rows, not the last
    one written or an arbitrary pick. Rows with phylo_status != 'analyzed' (not_analyzed,
    species_below_min_taxa, failed) carry reassortment_call='not_analyzed', which isn't in
    _REASSORTMENT_CALL_RANK and is simply never picked over anything analyzed.
    """
    best: Dict[str, str] = {}
    best_rank: Dict[str, int] = {}
    for row in phylo_rows:
        key = row["isolate_key"]
        call = row.get("reassortment_call", "")
        rank = _REASSORTMENT_CALL_RANK.get(call, -1)
        if rank > best_rank.get(key, -1):
            best[key] = call
            best_rank[key] = rank
    return best


def run_phylogeny(paths: "ah.InputPaths", args: argparse.Namespace) -> None:
    """Mirrors all_hanta.py's run_pipeline() ingest chain exactly through
    compute_segment_divergence (needed for divergence_by_key, which
    collect_phylo_candidates uses to tag divergence_candidate isolates) -- see spec
    section 3.2 for why this ~duplication is a deliberate choice, not an oversight.
    """
    logging.info("Verifying package integrity against NCBI md5sum.txt")
    checksums = ah.verify_inputs(paths)
    ah.report_checksums(checksums, allow_mismatch=args.allow_checksum_mismatch)

    logging.info("Loading reports")
    data_records = list(ah.iter_jsonl(paths.data_report))
    annotation_index = ah.build_annotation_index(ah.iter_jsonl(paths.annotation_report))
    biosample_index = ah.build_biosample_index(ah.iter_jsonl(paths.biosample_report))

    logging.info("Parsing FASTA")
    genomic_index = ah.parse_fasta(paths.genomic_fasta)
    protein_index = ah.parse_fasta(paths.protein_fasta)

    logging.info("Transforming accession records")
    accessions = []
    for record in data_records:
        transformed, _qc_row = ah.transform_one_record(
            record, annotation_index, biosample_index, genomic_index, protein_index)
        if transformed is not None:
            accessions.append(transformed)

    logging.info("Reconstructing isolates")
    isolates = ah.build_isolate_records(accessions)
    logging.info("  %s isolates, %s complete triplets",
                 len(isolates), sum(1 for i in isolates if i.has_complete_triplet))

    ah.assign_kmer_clusters(isolates, args.cluster_threshold)
    divergence = ah.compute_segment_divergence(isolates, min_species_size=args.min_species_size)
    divergence_by_key = {d.isolate_key: d for d in divergence}

    logging.info("Loading CDS")
    cds_index = parse_cds_fasta(paths.cds_fasta)
    logging.info("  cds: %s accessions with at least one CDS", len(cds_index))

    logging.info("Running phylogenetics stage")
    phylo_rows = run_phylogenetics_stage(isolates, divergence_by_key, cds_index, paths.output_dir, args)
    ah.safe_write_csv(paths.output_dir / "ortho_reassortment_table.csv", phylo_rows)
    logging.info("Done. %s candidate rows written to ortho_reassortment_table.csv", len(phylo_rows))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Orthohantavirus phylogenetics stage: MAFFT + IQ-TREE reassortment analysis "
                    "over the reference isolate set (Stage 2 of the three-stage pipeline)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", default=".",
                        help="Package root containing all_hanta/ncbi_dataset/data")
    parser.add_argument("--output", help="Output directory; defaults to <input>/output")
    parser.add_argument("--allow-checksum-mismatch", action="store_true",
                        help="Continue even if inputs do not match NCBI's md5sum.txt")
    parser.add_argument("--cluster-threshold", type=float, default=0.02,
                        help="Cosine k-mer distance for single-linkage clustering within a species")
    parser.add_argument("--min-species-size", type=int, default=8,
                        help="Minimum assessable isolates in a species before divergence "
                             "diagnostics (and divergence_candidate tagging) are computed")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    parser.add_argument("--mafft-bin", default=None, help="Path to mafft (default: PATH lookup)")
    parser.add_argument("--iqtree-bin", default=None, help="Path to iqtree (default: PATH lookup)")
    parser.add_argument("--phylo-min-taxa", type=int, default=4,
                        help="Minimum candidate isolates in a species before building segment trees")
    parser.add_argument("--bootstrap", type=int, default=1000,
                        help="IQ-TREE ultrafast bootstrap replicates (>=1000)")
    parser.add_argument("--min-bootstrap-support", type=float, default=95.0,
                        help="Minimum ultrafast bootstrap support to call a topology conflict supported")
    parser.add_argument("--phylo-seed", type=int, default=42,
                        help="Fixed IQ-TREE random seed, for reproducible trees")
    parser.add_argument("--phylo-threads", type=int, default=max(1, (os.cpu_count() or 4) - 1),
                        help="Threads for mafft and iqtree")
    parser.add_argument("--phylo-scope", choices=("species", "phylogroup", "both"), default="both",
                        help="species: within-species trees only. phylogroup: trees spanning "
                             "every species in a reservoir-host phylogroup (see ah.PHYLOGROUP), "
                             "for interspecies signal a species-only tree structurally cannot "
                             "show. both: run each independently.")
    parser.add_argument("--phylo-outgroup", default=None,
                        help="Accession or isolate_key to root trees on, via IQ-TREE's -o.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    paths = ah.resolve_inputs(Path(args.input), output=Path(args.output) if args.output else None)
    ah.configure_logging(paths.output_dir / "ortho_phylogeny.log", verbose=args.verbose)
    logging.info("Orthohantavirus phylogenetics stage (Stage 2)")
    logging.info("Package: %s", paths.data_dir)
    logging.info("Output : %s", paths.output_dir)
    run_phylogeny(paths, args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
