# detector_gates.py
"""Gates 0-4 and the verdict rubric for the reassortment candidate detector.

See docs/superpowers/specs/2026-08-12-reassortment-detector-design.md for the full
framework this implements.
"""

from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

import all_hanta as ah
import reassortment_reference_data as rd

_STOP_CODONS = {"TAA", "TAG", "TGA"}
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "CAT": "H",
    "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R",
    "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass
class CandidateInput:
    candidate_id: str
    s_genomic: str
    m_genomic: str
    l_genomic: str


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def _translate(seq: str) -> str:
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        protein.append(_CODON_TABLE.get(codon, "X"))
    return "".join(protein)


def _longest_orf_span(genomic_seq: str) -> str:
    """Longest ATG...stop nucleotide span (stop codon included) across all 6 reading
    frames. Shared by longest_orf_protein (Gate 0, translates it) and
    candidate_cds_from_orf (Gate 4, uses the nucleotide span directly) so the ORF-search
    loop exists exactly once. Returns '' if no in-frame stop codon terminates any
    ATG-initiated ORF.
    """
    genomic_seq = genomic_seq.upper()
    best = ""
    for strand_seq in (genomic_seq, _reverse_complement(genomic_seq)):
        for frame in range(3):
            i = frame
            while i < len(strand_seq) - 2:
                if strand_seq[i:i + 3] == "ATG":
                    j = i + 3
                    while j < len(strand_seq) - 2:
                        codon = strand_seq[j:j + 3]
                        if codon in _STOP_CODONS:
                            span = strand_seq[i:j + 3]
                            if len(span) > len(best):
                                best = span
                            break
                        j += 3
                i += 3
    return best


def longest_orf_protein(genomic_seq: str) -> str:
    """Longest ATG...stop ORF across all 6 reading frames, translated to protein (start
    and stop codons both excluded from the returned protein -- the initiator Met is not
    counted). Returns '' if no in-frame stop codon terminates any ATG-initiated ORF --
    Gate 0 treats that as a structural failure, not a crash.
    """
    span = _longest_orf_span(genomic_seq)
    return _translate(span[3:-3]) if span else ""


def empirical_protein_length_range(protein_lengths: List[int]) -> Tuple[int, int]:
    """(min, max) observed protein length. Wide-open (0, 10**9) when no reference data
    is available, so an empty reference set never spuriously fails Gate 0.
    """
    if not protein_lengths:
        return (0, 10**9)
    return (min(protein_lengths), max(protein_lengths))


def empirical_panhandle_consensus(genomic_sequences: List[str], n: int = 11) -> Tuple[str, str]:
    """Majority-vote consensus of the first and last n nucleotides across reference
    genomic sequences for one segment. Sequences shorter than n are skipped.
    """
    five_columns: List[Counter] = [Counter() for _ in range(n)]
    three_columns: List[Counter] = [Counter() for _ in range(n)]
    for seq in genomic_sequences:
        if len(seq) < n:
            continue
        for idx, base in enumerate(seq[:n]):
            five_columns[idx][base] += 1
        for idx, base in enumerate(seq[-n:]):
            three_columns[idx][base] += 1
    five = "".join(col.most_common(1)[0][0] if col else "N" for col in five_columns)
    three = "".join(col.most_common(1)[0][0] if col else "N" for col in three_columns)
    return five, three


@dataclass
class Gate0Result:
    passed_structure: bool
    passed_identity_floor: bool
    failures: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    orf_proteins: Dict[str, str] = field(default_factory=dict)
    best_kmer_identity: float = 0.0


def _kmer_identity_to_reference(protein: str, reference_profiles: List[Tuple[str, str]]) -> float:
    """1 - minimum cosine distance from `protein` to any reference isolate's protein for
    this segment. Computes only the candidate's own row (not a full NxN pairwise
    matrix), since only one comparison target (the candidate) is ever needed here.
    """
    if not protein or not reference_profiles:
        return 0.0
    candidate_profile = ah.kmer_profile(protein)
    if not candidate_profile:
        return 0.0
    best_identity = 0.0
    for _species, ref_protein in reference_profiles:
        ref_profile = ah.kmer_profile(ref_protein)
        matrix, _vocab = ah.profiles_to_matrix([candidate_profile, ref_profile])
        distance = ah.cosine_distance_matrix(matrix)[0, 1]
        identity = 1.0 - float(distance)
        if identity > best_identity:
            best_identity = identity
    return best_identity


def gate0_check(
    candidate: CandidateInput,
    protein_len_ranges: Dict[str, Tuple[int, int]],
    panhandle_consensus: Dict[str, Tuple[str, str]],
    reference_protein_profiles: Dict[str, List[Tuple[str, str]]],
    identity_floor: float = 0.35,
) -> Gate0Result:
    failures: List[str] = []
    warnings: List[str] = []
    orf_proteins: Dict[str, str] = {}
    genomic_by_segment = {"S": candidate.s_genomic, "M": candidate.m_genomic, "L": candidate.l_genomic}

    for segment, genomic in genomic_by_segment.items():
        low, high = ah.SEGMENT_LENGTH_RANGE[segment]
        if not (low <= len(genomic) <= high):
            failures.append(f"{segment}: genomic length {len(genomic)} outside expected {low}-{high}")

        protein = longest_orf_protein(genomic)
        orf_proteins[segment] = protein
        plen_low, plen_high = protein_len_ranges.get(segment, (0, 10**9))
        # Allow a wide tolerance band (empirical range +/- 15%) since a real but
        # divergent/novel candidate can legitimately sit outside the exact observed
        # reference range without being junk.
        tolerance = max(1, int(0.15 * (plen_high - plen_low + 1)))
        if not protein or not (plen_low - tolerance <= len(protein) <= plen_high + tolerance):
            failures.append(
                f"{segment}: no plausible ORF found (got {len(protein)} aa, "
                f"expected roughly {plen_low}-{plen_high})")

        consensus5, consensus3 = panhandle_consensus.get(segment, ("", ""))
        if consensus5 and consensus3 and len(genomic) >= len(consensus5):
            observed5 = genomic[:len(consensus5)].upper()
            observed3 = genomic[-len(consensus3):].upper()
            # Consensus region only (positions 2-12, i.e. index 1 onward): archived
            # "complete" sequences can carry 1-2 artificially added terminal bases
            # (Liu et al. 2024) -- checking the very first/last base as a hard
            # requirement would false-positive on real, correctly-submitted isolates.
            core5_match = observed5[1:] == consensus5[1:]
            core3_match = observed3[:-1] == consensus3[:-1]
            if not (core5_match and core3_match):
                failures.append(f"{segment}: terminal panhandle motif not recognizable")
            elif observed5[0] != consensus5[0] or observed3[-1] != consensus3[-1]:
                warnings.append(f"{segment}: terminal base mismatch at position 1 (or last), "
                                f"consistent with GenBank submission padding, not a failure")

    passed_structure = not failures

    identities = [
        _kmer_identity_to_reference(orf_proteins[segment], reference_protein_profiles.get(segment, []))
        for segment in ("S", "M", "L")
    ]
    best_kmer_identity = max(identities) if identities else 0.0
    passed_identity_floor = best_kmer_identity >= identity_floor

    return Gate0Result(
        passed_structure=passed_structure,
        passed_identity_floor=passed_identity_floor,
        failures=failures,
        warnings=warnings,
        orf_proteins=orf_proteins,
        best_kmer_identity=best_kmer_identity,
    )


@dataclass
class SpeciesHit:
    species: str
    distance: float


def nearest_species_by_kmer(
    candidate_protein: str, reference_profiles: List[Tuple[str, str]], top_n: int = 3,
) -> List[SpeciesHit]:
    """Top-n distinct species by minimum per-isolate cosine distance to
    candidate_protein. Deliberately per-isolate (not pre-aggregated per species) so a
    single close outlier isolate is still found, matching how compute_segment_divergence
    already works in all_hanta.py.

    Builds one candidate-vs-all matrix (not the full NxN pairwise matrix
    cosine_distance_matrix would build if handed all reference isolates together) --
    the candidate is always row 0.
    """
    if not candidate_protein or not reference_profiles:
        return []
    profiles = [ah.kmer_profile(candidate_protein)] + [
        ah.kmer_profile(protein) for _species, protein in reference_profiles]
    matrix, _vocab = ah.profiles_to_matrix(profiles)
    distances = ah.cosine_distance_matrix(matrix)[0, 1:]

    best_per_species: Dict[str, float] = {}
    for (species, _protein), distance in zip(reference_profiles, distances.tolist()):
        if species not in best_per_species or distance < best_per_species[species]:
            best_per_species[species] = distance

    ranked = sorted(best_per_species.items(), key=lambda item: item[1])
    return [SpeciesHit(species=species, distance=distance) for species, distance in ranked[:top_n]]


@dataclass
class Gate1Result:
    top_hits: Dict[str, List[SpeciesHit]] = field(default_factory=dict)
    implicated_species: List[str] = field(default_factory=list)
    discordant_segment: str = ""
    pattern: str = "concordant"


def gate1_check(
    orf_proteins: Dict[str, str],
    reference_profiles: Dict[str, List[Tuple[str, str]]],
    undetectable_floor: float = 0.02,
    too_divergent_ceiling: float = 0.9,
) -> Gate1Result:
    top_hits = {
        segment: nearest_species_by_kmer(orf_proteins.get(segment, ""), reference_profiles.get(segment, []))
        for segment in ("S", "M", "L")
    }
    best_species = {
        segment: hits[0].species if hits else "" for segment, hits in top_hits.items()
    }
    best_distance = {
        segment: hits[0].distance if hits else 1.0 for segment, hits in top_hits.items()
    }

    if len({s for s in best_species.values() if s}) <= 1:
        implicated = [s for s in {best_species["S"]} if s]
        return Gate1Result(top_hits=top_hits, implicated_species=implicated,
                            discordant_segment="", pattern="concordant")

    # Discordant segment = whichever disagrees with the other two (S/L majority vs M,
    # or S vs M/L, etc.) -- with only 3 segments, the odd one out is whichever species
    # doesn't match at least one other segment's top hit.
    counts = Counter(best_species.values())
    majority_species, majority_count = counts.most_common(1)[0]
    discordant_segment = next(
        (segment for segment, species in best_species.items() if species != majority_species),
        "M")
    pattern = "m_only" if discordant_segment == "M" else "s_or_l_involving"

    # too_divergent: the discordant segment doesn't convincingly match ANYTHING in the
    # reference set (best_distance is itself close to 1.0) -- no plausible phylogenetic
    # bridge, more consistent with an artifact than a real donor. This uses the
    # candidate's own best-match distance, which is the right quantity for "no good
    # match exists at all".
    #
    # undetectable_by_design: a *different* question -- can the top-2 candidate parent
    # species for this segment even be told apart? A tiny GAP between the best and
    # second-best hit's distances means the segment is genuinely ambiguous between two
    # near-identical parents (Reassortment_Report.md's "reassortment between
    # near-identical parents is invisible to any tree-based method" case). Deliberately
    # NOT based on best_distance itself: a small best_distance means a confident,
    # high-quality match -- the opposite of ambiguous -- so gating on that quantity
    # would mislabel every clean match as "undetectable".
    discordant_distance = best_distance[discordant_segment]
    if discordant_distance >= too_divergent_ceiling:
        pattern = "too_divergent"
    else:
        discordant_hits = top_hits[discordant_segment]
        if len(discordant_hits) >= 2 and (
                discordant_hits[1].distance - discordant_hits[0].distance) <= undetectable_floor:
            pattern = "undetectable_by_design"

    implicated_species = list({s for s in best_species.values() if s})
    return Gate1Result(top_hits=top_hits, implicated_species=implicated_species,
                        discordant_segment=discordant_segment, pattern=pattern)


@dataclass
class Gate2TrackAResult:
    anchor_species: str = ""
    anchor_basis: str = ""
    waasa_intact: Optional[bool] = None
    tmis_conserved: Optional[bool] = None
    ph_switch_conserved: Optional[bool] = None
    structural_check_not_available: bool = True


def gate2_track_a(candidate_m_protein: str, implicated_species: List[str]) -> Gate2TrackAResult:
    anchor_species = next((s for s in implicated_species if s in rd.STRUCTURAL_ANCHORS), None)
    if anchor_species is None or not candidate_m_protein:
        return Gate2TrackAResult(structural_check_not_available=True)

    anchor = rd.STRUCTURAL_ANCHORS[anchor_species]
    ref_seq = anchor["sequence"]

    tmis_ok = True
    tmis_checked_any = False
    for _label, position in anchor["tmis"].items():
        if position is None:
            continue
        tmis_checked_any = True
        residue, aligned = rd.find_aligned_residue(candidate_m_protein, ref_seq, position)
        ref_residue = ref_seq[position - 1]
        if not aligned or residue != ref_residue:
            tmis_ok = False

    ph_ok = True
    for _label, position in anchor["ph_switch"].items():
        residue, aligned = rd.find_aligned_residue(candidate_m_protein, ref_seq, position)
        ref_residue = ref_seq[position - 1]
        if not aligned or residue != ref_residue:
            ph_ok = False

    return Gate2TrackAResult(
        anchor_species=anchor_species,
        anchor_basis=anchor["basis"],
        waasa_intact="WAASA" in candidate_m_protein,
        tmis_conserved=tmis_ok if tmis_checked_any else None,
        ph_switch_conserved=ph_ok,
        structural_check_not_available=False,
    )
