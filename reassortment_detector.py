# reassortment_detector.py
"""Reassortment candidate detector: reference data, reference loading, Gates 0-4 and the
verdict rubric, and the CLI entry point -- all in one module.

Scores a brand-new candidate hantavirus genome (S+M+L genomic FASTA) against the
five-gate framework in docs/superpowers/specs/2026-08-12-reassortment-detector-design.md.
Imports all_hanta.py's ingest and phylogenetics machinery; makes no changes to it.

This module was originally split into reassortment_detector.py (CLI),
detector_gates.py (gates 0-4 and the verdict rubric), detector_reference.py (reference
dataset loading), and reassortment_reference_data.py (curated structural/literature
reference data) per docs/superpowers/specs/2026-08-13-three-stage-pipeline-split-design.md.
Reconsolidated back into this single file per
docs/superpowers/notes/2026-08-14-reassortment-detector-reconsolidation.md.
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Tuple

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

import all_hanta as ah
import phylogeny as ph


# ---------------------------------------------------------------------------
# Structural/literature reference data for Gate 2 Track A.
#
# Every sequence and residue position here was fetched and verified during the design
# session -- see docs/superpowers/specs/2026-08-12-reassortment-detector-design.md sections
# 7 and 15 for the primary-source citations and the empirical WAASA-anchored verification
# that found PUUV's own accession does NOT match Tscherne et al. 2025's literal integers.
# ---------------------------------------------------------------------------

# Fetched via NCBI eutils efetch (db=protein, rettype=fasta) during design, 2026-08-12.
PUUV_GPC_SEQUENCE = (
    "MGKSSPVCLYLILQGLLLFDTVNARNLNELKMECPHTIGLGQGLVVGSVELPPVPIQQIESLKLESSCNF"
    "DLHTSTAGQQSFTKWTWETKGDLAENTQASSTSFQTKSSEVNLRGLCLIPTLVVETAARMRKTIACYDLS"
    "CNQTVCQPTVYLMGPIQTCLTTKSCLLGLGDQRIQVNYERTYCVSGQLVEGVCFNPIHTMALSQPSHTYD"
    "IVTIMVRCFLVIKKVTSGDSMKIEKNFETLVQKTGCTANGFQGYYICLIGSSSEPLYVPTLDDYRSAEVL"
    "SRMAFAPHGEDHDIEKNAVSALRIAGKVTGKAPSTESSDTVQGIAFSGSPLYTSTGVLTAKDDPVYVWAP"
    "GIIMEGNHSVCEKKTLPLTWTGFIPLPGEIEKTTQCTVFCTLAGPGADCEAYSETGIFNISSPTCLINRV"
    "QRFRGAEQQIKFVCQRVDMDITVYCNGVKKVILTKTLVIGQCIYTFTSIFSLIPGIAHSLAVELCVPGLH"
    "GWATVLLLLTFCSGWVLIPTITMILLKILIAFAHLCSKYNTDSKFKILVEKVKKEYQKTMGSMVCEVCQY"
    "ECETAKELESHRKSCSIGSCPYCLNPSEATPSALQAHFKVCKLTSRFQENLKKSLTMYEPMQGCYRTLSL"
    "FRYRSRFFVGLVWCMLLVLELIVWAASAETQNLNAGWTDTAHGSGIIPMKADLELDFSLPSSASYTYRRQ"
    "LQNPANEQEKIPFHLQISKQVIHAEIQHLGHWMDATFNLKTAFHCYGSCEKYAYPWQTAGCFVEKDYEYE"
    "TGWGCNPPDCPGVGTGCTACGVYLDKLKSVGKVFKIVSLRYTRKVCIQLGTEQTCKTVDSNDCLITTSVK"
    "VCLIGTISKFQPSDTLLFLGPLQQGGLIFKQWCTTTCQFGDPGDIMSTPTGMKCPELNGSFRKKCAFATT"
    "PVCQFDGNTISGYKRMVATKDSFQSFNVTEPHISTSALEWIDPDSSLRDHINVIVSRDLSFQDLSETPCQ"
    "VDLTTSAIDGAWGSGVGFNLVCTVSLTECSAFLTSIKACDAAMCYGSTTTNLVRGQNTIHVVGKGGHSGS"
    "KFMCCHDTKCSSTGLVAAAPHLDRVTGFNQADSDKIFDDGAPECGMSCWFKKSGEWVLGVLNGNWMVVAV"
    "LIALLILSIFLFALCCPRRPSYKKEHKP"
)  # QOC68549.1, PUUV Sotkamo GPC precursor, 1148 aa

ANDV_GPC_SEQUENCE = (
    "MEGWYLVVLGVCYTLTLAMPKTIYELKMECPHTVGLGQGYIIGSTELGLISIEAASDIKLESSCNFDLHT"
    "TSMAQKSFTQVEWRKKSDTTDTTNAASTTFEAQTKTVNLRGTCILAPELYDTLKKVKKTVLCYDLTCNQT"
    "HCQPTVYLIAPVLTCMSIRSCMASVFTSRIQVIYEKTHCVTGQLIEGQCFNPAHTLTLSQPAHTYDTVTL"
    "PISCFFTPKKSEQLKVIKTFEGILTKTGCTENALQGYYVCFLGSHSEPLIVPSLEDIRSAEVVSRMLVHP"
    "RGEDHDAIQNSQSHLRIVGPITAKVPSTSSTDTLKGTAFAGVPMYSSLSTLVRNADPEFVFSPGIVPESN"
    "HSTCDKKTVPITWTGYLPISGEMEKVTGCTVFCTLAGPGASCEAYSENGIFNISSPTCLVNKVQRFRGSE"
    "QKINFICQRVDQDVVVYCNGQKKVILTKTLVIGQCIYTFTSLFSLMPDVAHSLAVELCVPGLHGWATVML"
    "LSTFCFGWVLIPAVTLIILKCLRVLTFSCSHYTNESKFKFILEKVKIEYQKTMGSMVCDVCHHECETAKE"
    "LESHRQSCINGQCPYCMTITEATESALQAHYSICKLTGRFQEALKKSLKKPEVKKGCYRTLGVFRYKSRC"
    "YVGLVWCLLLTCEIVIWAASAETPLMESGWSDTAHGVGEIPMKTDLELDFSLPSSSSYSYRRKLTNPANK"
    "EESIPFHFQMEKQVIHAEIQPLGHWMDATFNIKTAFHCYGACQKYSYPWQTSKCFFEKDYQYETGWGCNP"
    "GDCPGVGTGCTACGVYLDKLKSVGKAYKIISLKYTRKVCIQLGTEQTCKHIDANDCLVTPSVKVCIVGTV"
    "SKLQPSDTLLFLGPLEQGGIILKQWCTTSCAFGDPGDIMSTPSGMRCPEHTGSFRKICGFATTPVCEYQG"
    "NTISGYKRMMATKDSFQSFNLTEPHITTNKLEWIDPDGNTRDHVNLVLNRDVSFQDLSDNPCKVDLHTQA"
    "IEGAWGSGVGFTLTCTVGLTECPSFMTSIKACDLAMCYGSTVTNLARGSNTVKVVGKGGHSGSSFKCCHD"
    "TDCSSEGLLASAPHLERVTGFNQIDSDKVYDDGAPPCTFKCWFTKSGEWLLGILNGNWIVVVVLVVILIL"
    "SIIMFSVLCPRRGHKKTV"
)  # NP_604472.1, ANDV Chile-9717869 GPC precursor, 1138 aa -- the numerical anchor

HTNV_GPC_SEQUENCE = (
    "MGIWKWLVMASLVWPVLTLRNVYDMKIECPHTVSFGENSVIGYVELPPVPLADTAQMVPESSCNMDNHQS"
    "LNTITKYTQVSWRGKADQSQSSQNSFETVSTEVDLKGTCVLKHKMVEESYRSRKSVTCYDLSCNSTYCKP"
    "TLYMIVPIHACNMMKSCLIALGPYRVQVVYERSYCMTGVLIEGKCFVPDQSVVSIIKHGIFDIASVHIVC"
    "FFVAVKGNTYKIFEQVKKSFESTCNDTENKVQGYYICIVGGNSAPIYVPTLDDFRSMEAFTGIFRSPHGE"
    "DHDLAGEEIASYSIVGPANAKVPHSASSDTLSLIAYSGIPSYSSLSILTSSTEAKHVFSPGLFPKLNHTN"
    "CDKSAIPLIWTGMIDLPGYYEAVHPCTVFCVLSGPGASCEAFSEGGIFNITSPMCLVSKQNRFRLTEQQV"
    "NFVCQRVDMDIVVYCNGQRKVILTKTLVIGQCIYTITSLFSLLPGVAHSIAVELCVPGFHGWATAALLVT"
    "FCFGWVLIPAITFIILTVLKFIANIFHTSNQENRLKSVLRKIKEEFEKTKGSMVCDVCKYECETYKELKA"
    "HGVSCPQSQCPYCFTHCEPTEAAFQAHYKVCQVTHRFRDDLKKTVTPQNFTPGCYRTLNLFRYKSRCYIF"
    "TMWIFLLVLESILWAASASETPLTPVWNDNAHGVGSVPMHTDLELDFSLTSSSKYTYRRKLTNPLEEAQS"
    "IDLHIEIEEQTIGVDVHALGHWFDGRLNLKTSFHCYGACTKYEYPWHTAKCHYERDYQYETSWGCNPSDC"
    "PGVGTGCTACGLYLDQLKPVGSAYKIITIRYSRRVCVQFGEENLCKIIDMNDCFVSRHVKVCIIGTVSKF"
    "SQGDTLLFFGPLEGGGLIFKHWCTSTCQFGDPGDIMSPRDKGFLCPEFPGSFRKKCNFATTPICEYDGNM"
    "VSGYKKVMATIDSFQSFNTSTMHFTDERIEWKDPDGMLRDHINILVTKDIDFDNLGENPCKIGLQTSSIE"
    "GAWGSGVGFTLTCLVSLTECPTFLTSIKACDKAICYGAESVTLTRGQNTVKVSGKGGHSGSTFRCCHGED"
    "CSQIGLHAAAPHLDKVNGISEIENSKVYDDGAPQCGIKCWFVKSGEWISGIFSGNWIVLIVLCVFLLFSL"
    "VLLSILCPVRKHKKS"
)  # P08668.1, HTNV strain 76-118 GPC precursor, 1135 aa

# Positions verified empirically during design by anchoring on each sequence's own
# WAASA motif and checking whether the ANDV-relative offset reproduces the expected
# residue identity at each position -- see spec section 15 for the full table. ANDV
# matches at Tscherne et al. 2025's literal cited integers (745/757/759/766/900);
# PUUV only matches after a -7 offset; HTNV matches 4 of 5 after a +3 offset (F900's
# equivalent is a local exception, recorded as None rather than a guessed number).
STRUCTURAL_ANCHORS: Dict[str, Dict[str, Any]] = {
    "Orthohantavirus andesense": {
        "accession": "NP_604472.1",
        "sequence": ANDV_GPC_SEQUENCE,
        "basis": "solved",
        "waasa_position": 647,
        "tmis": {"Y745": 745, "W766": 766, "F900": 900},
        "ph_switch": {"E757": 757, "D759": 759},
    },
    "Orthohantavirus puumalaense": {
        "accession": "QOC68549.1",
        "sequence": PUUV_GPC_SEQUENCE,
        "basis": "solved",
        "waasa_position": 654,
        "tmis": {"Y745": 752, "W766": 773, "F900": 907},
        "ph_switch": {"E757": 764, "D759": 766},
    },
    "Orthohantavirus hantanense": {
        "accession": "P08668.1",
        "sequence": HTNV_GPC_SEQUENCE,
        "basis": "homology_inferred_verified",
        "waasa_position": 644,
        "tmis": {"Y745": 742, "W766": 763, "F900": None},
        "ph_switch": {"E757": 754, "D759": 756},
    },
}

# Empty by design -- see spec section 7 Track B implementation note. No citable reference
# accession was found for Mydocs/Orthohantavirus_01.md's L179 NP / V71/G700/D967 GP markers.
# Shape when populated: {species: {"segment": "S"|"M", "accession": str, "sequence": str,
# "markers": {label: position_1based}}}.
LINEAGE_MARKER_TABLE: Dict[str, Dict[str, Any]] = {}

# By genus clade (New World vs. Old World), per Jangra et al. 2018 (Nature,
# doi:10.1038/s41586-018-0702-1) establishing PCDH1 dependence specifically for New
# World hantaviruses -- see spec section 7 Track C for the correction from
# Mydocs/Orthohantavirus_01.md's oversimplified beta-3-integrin framing.
NEW_WORLD_SPECIES: FrozenSet[str] = frozenset({
    "Orthohantavirus andesense",
    "Orthohantavirus sinnombreense",
    "Orthohantavirus maporalense",
    "Orthohantavirus montanoense",
    "Orthohantavirus carrizalense",
    "Orthohantavirus delgaditoense",
})

_ALIGNER = PairwiseAligner()
_ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
_ALIGNER.mode = "global"
_ALIGNER.open_gap_score = -10
_ALIGNER.extend_gap_score = -0.5


@lru_cache(maxsize=None)
def _cached_aligned_pair(reference_sequence: str, candidate_protein: str) -> Tuple[str, str]:
    """Memoized global alignment of one (reference_sequence, candidate_protein) pair,
    returning just the two aligned strings. Callers checking several residue positions
    for the same pair (gate2_track_a checks up to 5 positions per species, gate2_track_b
    checks one per marker) previously re-ran this full alignment once per position --
    the alignment itself doesn't depend on which position is being looked up.
    """
    alignment = _ALIGNER.align(reference_sequence, candidate_protein)[0]
    return str(alignment[0]), str(alignment[1])


def find_aligned_residue(
    candidate_protein: str, reference_sequence: str, reference_position_1based: int,
) -> Tuple[str, bool]:
    """Map reference_position_1based (1-based, into reference_sequence) to the aligned
    residue in candidate_protein, via a global alignment of the two sequences.

    Returns (residue, aligned) -- residue is '' and aligned is False if the reference
    position falls on a gap in the candidate, or either sequence is too short/empty to
    align meaningfully.
    """
    if not candidate_protein or not reference_sequence:
        return "", False
    if not (1 <= reference_position_1based <= len(reference_sequence)):
        return "", False

    ref_aligned, cand_aligned = _cached_aligned_pair(reference_sequence, candidate_protein)

    ref_index = 0  # 0-based index into reference_sequence, i.e. reference_position_1based - 1
    for ref_char, cand_char in zip(ref_aligned, cand_aligned):
        if ref_char != "-":
            if ref_index == reference_position_1based - 1:
                if cand_char == "-":
                    return "", False
                return cand_char, True
            ref_index += 1
    return "", False


# ---------------------------------------------------------------------------
# Reference dataset loading.
#
# load_reference_data replicates the ingest portion of all_hanta.py's run_pipeline (see
# that function, lines ~3960-4001) without any of the table-building, ML-matrix, or
# output-writing steps this detector doesn't need.
# ---------------------------------------------------------------------------

@dataclass
class ReferenceData:
    isolates: List["ah.IsolateRecord"]
    cds_index: Dict[str, List[Tuple[str, str]]]
    genomic_index: Dict[str, "ah.FastaRecord"]


def load_reference_data(root: Path) -> ReferenceData:
    paths = ah.resolve_inputs(root)

    data_records = list(ah.iter_jsonl(paths.data_report))
    annotation_index = ah.build_annotation_index(ah.iter_jsonl(paths.annotation_report))
    biosample_index = ah.build_biosample_index(ah.iter_jsonl(paths.biosample_report))
    genomic_index = ah.parse_fasta(paths.genomic_fasta)
    protein_index = ah.parse_fasta(paths.protein_fasta)

    accessions = []
    for record in data_records:
        transformed, _qc_row = ah.transform_one_record(
            record, annotation_index, biosample_index, genomic_index, protein_index)
        if transformed is not None:
            accessions.append(transformed)

    isolates = ah.build_isolate_records(accessions)
    cds_index = ph.parse_cds_fasta(paths.cds_fasta)

    return ReferenceData(isolates=isolates, cds_index=cds_index, genomic_index=genomic_index)


def protein_lengths_by_segment(isolates: List["ah.IsolateRecord"]) -> Dict[str, List[int]]:
    lengths: Dict[str, List[int]] = {"S": [], "M": [], "L": []}
    for isolate in isolates:
        for segment, attribute in (("S", "s_protein_aa"), ("M", "m_protein_aa"), ("L", "l_protein_aa")):
            protein = getattr(isolate, attribute)
            if protein:
                lengths[segment].append(len(protein))
    return lengths


def genomic_sequences_by_segment(
    isolates: List["ah.IsolateRecord"], genomic_index: Dict[str, "ah.FastaRecord"],
) -> Dict[str, List[str]]:
    sequences: Dict[str, List[str]] = {"S": [], "M": [], "L": []}
    for isolate in isolates:
        for segment, accession in (("S", isolate.s_accession), ("M", isolate.m_accession),
                                    ("L", isolate.l_accession)):
            if not accession:
                continue
            record = genomic_index.get(accession) or genomic_index.get(ah.accession_base(accession))
            if record is not None and record.sequence:
                sequences[segment].append(record.sequence)
    return sequences


def protein_profiles_by_segment(
    isolates: List["ah.IsolateRecord"],
) -> Dict[str, List[Tuple[str, str]]]:
    profiles: Dict[str, List[Tuple[str, str]]] = {"S": [], "M": [], "L": []}
    for isolate in isolates:
        for segment, attribute in (("S", "s_protein_aa"), ("M", "m_protein_aa"), ("L", "l_protein_aa")):
            protein = getattr(isolate, attribute)
            if protein:
                profiles[segment].append((isolate.species_name, protein))
    return profiles


# ---------------------------------------------------------------------------
# Gates 0-4 and the verdict rubric.
#
# See docs/superpowers/specs/2026-08-12-reassortment-detector-design.md for the full
# framework this implements.
# ---------------------------------------------------------------------------

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


@lru_cache(maxsize=None)
def _cached_reference_kmer_profile(protein: str) -> Dict[str, int]:
    """Memoized ah.kmer_profile for REFERENCE proteins only -- the reference set is
    identical across every candidate scored in one run (unlike the candidate's own
    protein, which changes every call and gains nothing from caching). Safe to cache a
    shared dict here: profiles_to_matrix only reads its inputs, never mutates them.
    """
    return ah.kmer_profile(protein)


def _kmer_identity_to_reference(protein: str, reference_profiles: List[Tuple[str, str]]) -> float:
    """1 - minimum cosine distance from `protein` to any reference isolate's protein for
    this segment. Vectorized like nearest_species_by_kmer: one candidate-vs-all matrix,
    one batch of dot products, instead of a fresh 2-row profiles_to_matrix +
    cosine_distance_matrix call per reference isolate. The per-isolate loop this
    replaced paid that full pairwise-matrix cost once per reference isolate (thousands
    per segment), measured at several minutes per candidate against the real reference
    set -- this pays it once, using the same _cached_reference_kmer_profile memoization
    the loop already relied on.
    """
    if not protein or not reference_profiles:
        return 0.0
    candidate_profile = ah.kmer_profile(protein)
    if not candidate_profile:
        return 0.0
    profiles = [candidate_profile] + [
        _cached_reference_kmer_profile(ref_protein) for _species, ref_protein in reference_profiles]
    matrix, _vocab = ah.profiles_to_matrix(profiles)
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    unit = matrix / norms
    similarity = np.clip(unit[0] @ unit[1:].T, -1.0, 1.0)
    distances = 1.0 - similarity
    return float(max(0.0, 1.0 - distances.min()))


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
            def _core_match(candidate_seq: str) -> Tuple[bool, bool]:
                """(core_match, exact_match) for one orientation of candidate_seq
                against consensus5/consensus3. Consensus region only (positions 2-12,
                i.e. index 1 onward): archived "complete" sequences can carry 1-2
                artificially added terminal bases (Liu et al. 2024) -- checking the
                very first/last base as a hard requirement would false-positive on
                real, correctly-submitted isolates.
                """
                observed5 = candidate_seq[:len(consensus5)].upper()
                observed3 = candidate_seq[-len(consensus3):].upper()
                core = observed5[1:] == consensus5[1:] and observed3[:-1] == consensus3[:-1]
                exact = observed5[0] == consensus5[0] and observed3[-1] == consensus3[-1]
                return core, exact

            core_match, exact_match = _core_match(genomic)
            if not core_match:
                # Real deposits aren't always oriented the same way -- e.g. this
                # dataset's own puumalaense isolate OR915064.1 (vindeln) starts
                # "GGGAGTCTACTACTA", whose reverse complement is
                # "TAGTAGTAGACTCCC", a near-exact panhandle match. The ORF search
                # already checks both strands (_longest_orf_span); this check should
                # too, rather than penalizing a real isolate for strand orientation.
                core_match, exact_match = _core_match(_reverse_complement(genomic))

            if core_match and not exact_match:
                warnings.append(f"{segment}: terminal base mismatch at position 1 (or last), "
                                f"consistent with GenBank submission padding, not a failure")
            elif not core_match:
                # A downgrade from Gate 0's original design: measured directly against
                # this project's own reference dataset, only ~22% of genomic S records
                # start with the recognizable panhandle motif even after checking both
                # orientations -- many real, non-junk NCBI submissions simply don't
                # include the untranslated terminal region. Treating this as a hard
                # failure would call real reference isolates "junk_sequence". ORF
                # length and overall genomic length remain the load-bearing structural
                # checks; a missing panhandle is now informational only.
                warnings.append(f"{segment}: terminal panhandle motif not recognizable "
                                f"in either orientation")

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

    Computes only the candidate-vs-all row of distances directly (O(N) dot products),
    not the full NxN pairwise matrix ah.cosine_distance_matrix builds -- that function
    always computes unit @ unit.T internally regardless of which rows the caller
    actually reads afterward, so calling it here and slicing row 0 would still pay the
    full O(N^2) cost up front. Confirmed by direct measurement against the real
    reference set (thousands of profiles per segment): the O(N^2) call took several
    seconds per invocation where this takes a fraction of that.
    """
    if not candidate_protein or not reference_profiles:
        return []
    profiles = [ah.kmer_profile(candidate_protein)] + [
        ah.kmer_profile(protein) for _species, protein in reference_profiles]
    matrix, _vocab = ah.profiles_to_matrix(profiles)
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    unit = matrix / norms
    similarity = np.clip(unit[0] @ unit[1:].T, -1.0, 1.0)
    distances = 1.0 - similarity

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
    anchor_species = next((s for s in implicated_species if s in STRUCTURAL_ANCHORS), None)
    if anchor_species is None or not candidate_m_protein:
        return Gate2TrackAResult(structural_check_not_available=True)

    anchor = STRUCTURAL_ANCHORS[anchor_species]
    ref_seq = anchor["sequence"]

    tmis_ok = True
    tmis_checked_any = False
    for _label, position in anchor["tmis"].items():
        if position is None:
            continue
        tmis_checked_any = True
        residue, aligned = find_aligned_residue(candidate_m_protein, ref_seq, position)
        ref_residue = ref_seq[position - 1]
        if not aligned or residue != ref_residue:
            tmis_ok = False

    ph_ok = True
    for _label, position in anchor["ph_switch"].items():
        residue, aligned = find_aligned_residue(candidate_m_protein, ref_seq, position)
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


@dataclass
class Gate2TrackBResult:
    lineage_marker_match: Optional[bool] = None
    matched_species: str = ""
    marker_set_used: str = ""


def gate2_track_b(candidate_proteins: Dict[str, str], implicated_species: List[str]) -> Gate2TrackBResult:
    for species in implicated_species:
        entry = LINEAGE_MARKER_TABLE.get(species)
        if entry is None:
            continue
        candidate_protein = candidate_proteins.get(entry["segment"], "")
        if not candidate_protein:
            continue
        all_match = True
        for label, position in entry["markers"].items():
            residue, aligned = find_aligned_residue(candidate_protein, entry["sequence"], position)
            ref_residue = entry["sequence"][position - 1]
            if not aligned or residue != ref_residue:
                all_match = False
        return Gate2TrackBResult(
            lineage_marker_match=all_match, matched_species=species,
            marker_set_used=f"{species}:{entry['accession']}")
    return Gate2TrackBResult(lineage_marker_match=None, matched_species="", marker_set_used="")


@dataclass
class Gate2TrackCResult:
    sl_mismatch_high_risk: bool = False
    receptor_class_note: str = "unknown"


def gate2_track_c(gate1: Gate1Result) -> Gate2TrackCResult:
    sl_mismatch = gate1.discordant_segment in ("S", "L")

    receptor_note = "unknown"
    if gate1.discordant_segment == "M" and len(gate1.implicated_species) >= 2:
        clades = {
            species in NEW_WORLD_SPECIES for species in gate1.implicated_species
        }
        if clades == {True, False}:
            # Ambiguous about which side donor/recipient are without more info -- record
            # the crossing happened, not a specific direction claim.
            receptor_note = "new_world_to_old_world"
        else:
            receptor_note = "same_clade"

    return Gate2TrackCResult(sl_mismatch_high_risk=sl_mismatch, receptor_class_note=receptor_note)


@dataclass
class Gate3Result:
    status: str = "ecology_unknown"
    detail: str = ""


def gate3_check(implicated_species: List[str], isolates: List["ah.IsolateRecord"]) -> Gate3Result:
    if len(implicated_species) < 2:
        return Gate3Result(status="ecology_unknown", detail="single implicated species, no overlap to test")

    by_species: Dict[str, List["ah.IsolateRecord"]] = defaultdict(list)
    for isolate in isolates:
        if isolate.is_lab_derived or isolate.is_vaccine_strain:
            continue
        if isolate.species_name in implicated_species:
            by_species[isolate.species_name].append(isolate)

    missing = [s for s in implicated_species if not by_species.get(s)]
    if missing:
        return Gate3Result(status="ecology_unknown",
                           detail=f"no wild reference isolate metadata for: {', '.join(missing)}")

    host_orders = {s: {i.host_order for i in members if i.host_order} for s, members in by_species.items()}
    countries = {s: {i.country for i in members if i.country} for s, members in by_species.items()}

    species_list = implicated_species[:2]
    a, b = species_list[0], species_list[1]
    host_overlap = bool(host_orders[a] & host_orders[b])
    country_overlap = bool(countries[a] & countries[b])

    if host_overlap and country_overlap:
        return Gate3Result(status="ecological_overlap_plausible",
                           detail=f"shared host order and country between {a} and {b}")
    if not host_orders[a] or not host_orders[b] or not countries[a] or not countries[b]:
        return Gate3Result(status="ecology_unknown", detail="insufficient metadata to judge overlap")
    return Gate3Result(status="no_ecological_route",
                       detail=f"no shared host order/country found between {a} and {b}")


def candidate_cds_from_orf(genomic_seq: str) -> str:
    """The ORF's own nucleotide span (ATG..stop, stop included) -- what Gate 4 aligns,
    since _analyze_group works on nucleotide CDS, not translated protein. Reuses
    _longest_orf_span (Task 2) rather than re-implementing the same 6-frame ORF search --
    that function already returns exactly this span; longest_orf_protein is the one that
    additionally translates it.
    """
    return _longest_orf_span(genomic_seq)


def build_candidate_isolate_record(candidate_id: str, species_label: str) -> "ah.IsolateRecord":
    """Minimal synthetic IsolateRecord for Gate 4's PhyloCandidate -- only the fields
    _analyze_group, group_candidates and _phylo_row actually read are given meaningful
    values (confirmed by reading all_hanta.py directly during design: isolate_key,
    species_name, host_order, is_lab_derived, is_vaccine_strain, species_status,
    has_complete_triplet). Everything else is a valid placeholder of the right type,
    following the same pattern test_all_hanta.py's _make_isolate already uses.
    """
    return ah.IsolateRecord(
        isolate_key=f"CANDIDATE:{candidate_id}", species_name=species_label, species_taxid="",
        species_status="classified", biosample="", isolate=candidate_id, strain="",
        host_name="", host_taxid="", host_order="", host_family="", host_role="",
        collection_date="", year="", location="", country="", region="",
        grouping_tier="", grouping_confidence=0.0,
        s_accession="", m_accession="", l_accession="",
        s_completeness="", m_completeness="", l_completeness="",
        s_nt_length=0, m_nt_length=0, l_nt_length=0,
        s_protein_aa="", m_protein_aa="", l_protein_aa="", concatenated_protein_aa="",
        segments_present=["S", "M", "L"], coding_complete_segments=["S", "M", "L"],
        triplet_status="coding_complete_triplet", has_complete_triplet=True,
        is_lab_derived=False, is_vaccine_strain=False, has_nss_orf=False,
        metadata_conflict=False, taxonomic_incongruence=False, cross_species_partners=[],
        kmer_cluster_id="", s_protein_accessions=[], m_protein_accessions=[],
        l_protein_accessions=[], all_accessions=[], warnings=[],
    )


def build_comparator_members(
    implicated_species: List[str], reference: "ReferenceData", phylo_min_taxa: int,
) -> Tuple[List["ph.PhyloCandidate"], str]:
    """Union of reference isolates (with full S+M+L CDS) from every implicated species --
    not a single-species tree, so an M-only exchange candidate has both the S/L-matching
    and M-matching species' isolates available to be visible against.
    """
    members: List["ph.PhyloCandidate"] = []
    for isolate in reference.isolates:
        if isolate.species_name not in implicated_species:
            continue
        if isolate.is_lab_derived or isolate.is_vaccine_strain:
            continue
        sequences: Dict[str, str] = {}
        for segment, accession in (("S", isolate.s_accession), ("M", isolate.m_accession),
                                   ("L", isolate.l_accession)):
            records = reference.cds_index.get(ah.accession_base(accession)) if accession else None
            seq = ph.select_primary_cds(records, segment)
            if len(seq) >= ph.MIN_CDS_LENGTH:
                sequences[segment] = seq
        if len(sequences) == 3:
            members.append(ph.PhyloCandidate(isolate=isolate, tags=["comparator"], sequences=sequences))

    if len(members) < phylo_min_taxa:
        return [], "insufficient_comparators"
    return sorted(members, key=lambda m: m.isolate.isolate_key), ""


@dataclass
class Gate4Result:
    phylo_status: str = "not_analyzed"
    reassortment_call: str = "not_analyzed"
    discordant_segment: str = ""
    bootstrap_support: Optional[float] = None
    detail: Dict[str, Any] = field(default_factory=dict)


def gate4_check(
    candidate: CandidateInput,
    candidate_cds: Dict[str, str],
    implicated_species: List[str],
    reference: "ReferenceData",
    output_dir: Path,
    mafft_bin: str,
    iqtree_bin: str,
    phylo_min_taxa: int = 4,
    bootstrap: int = 1000,
    min_bootstrap_support: float = 95.0,
    phylo_seed: int = 42,
    phylo_threads: int = 1,
) -> Gate4Result:
    comparator_members, status = build_comparator_members(implicated_species, reference, phylo_min_taxa)
    if status:
        return Gate4Result(phylo_status=status, reassortment_call="not_analyzed")

    if any(not candidate_cds.get(segment) for segment in ("S", "M", "L")):
        return Gate4Result(phylo_status="failed", reassortment_call="not_analyzed",
                           detail={"error": "candidate is missing CDS for one or more segments"})

    species_label = implicated_species[0] if implicated_species else "unclassified"
    candidate_isolate = build_candidate_isolate_record(candidate.candidate_id, species_label)
    candidate_member = ph.PhyloCandidate(
        isolate=candidate_isolate, tags=["candidate"], sequences=dict(candidate_cds))

    members = [candidate_member] + comparator_members
    args = argparse.Namespace(
        phylo_outgroup=None, bootstrap=bootstrap, min_bootstrap_support=min_bootstrap_support,
        phylo_seed=phylo_seed, phylo_threads=phylo_threads,
    )
    # build_comparator_members unions isolates from every implicated species, not one --
    # so this is a "species" tree only when exactly one species is implicated. With two
    # or more, it's the same multi-species situation phylogeny.py's own docstrings call
    # "phylogroup" scope (group_key names a group, not any one member's real species;
    # see _analyze_group's docstring and group_candidates's). Mislabeling a multi-species
    # tree "species" would mean a downstream reader of analysis_scope=='species' wrongly
    # assumes every member shares one species, when the "+"-joined group_key already
    # signals otherwise.
    scope = "species" if len(implicated_species) == 1 else "phylogroup"
    group_dir_root = output_dir / "candidate_phylogenetics" / candidate.candidate_id
    rows = ph._analyze_group(
        "+".join(sorted(implicated_species)), scope, members, group_dir_root, output_dir,
        args, mafft_bin, iqtree_bin,
    )
    candidate_row = next(r for r in rows if r["isolate_key"] == candidate_isolate.isolate_key)

    return Gate4Result(
        phylo_status=candidate_row["phylo_status"],
        reassortment_call=candidate_row["reassortment_call"],
        discordant_segment=candidate_row["discordant_segment"],
        bootstrap_support=candidate_row["discordant_clade_bootstrap_support"],
        detail=candidate_row,
    )


def combine_verdict(
    gate0: Gate0Result,
    gate1: Gate1Result,
    gate2a: Gate2TrackAResult,
    gate2b: Gate2TrackBResult,
    gate2c: Gate2TrackCResult,
    gate3: Gate3Result,
    gate4: Gate4Result,
) -> str:
    """Combines every gate's finding into one verdict, per the rubric in
    docs/superpowers/specs/2026-08-12-reassortment-detector-design.md section 10.

    No min_bootstrap_support parameter here: gate4.reassortment_call already reflects
    the threshold (all_hanta.py's classify_reassortment_call applies it when gate4_check
    calls _analyze_group), so re-checking it here would be a second, redundant copy of
    the same decision rather than an independent gate.
    """
    if not gate0.passed_structure:
        return "junk_sequence"
    if not gate0.passed_identity_floor:
        return "divergent_novel_candidate"
    if gate1.pattern == "undetectable_by_design":
        return "undetectable_by_design"
    # gate1.pattern == "too_divergent" is deliberately NOT gated here. It has no verdict
    # of its own in the spec's rubric (section 10, a closed 10-verdict list) -- unlike
    # undetectable_by_design, which does. This matches the same principle the spec uses
    # to retire dead_end_non_viable: a gate's finding can ride as an annotation column
    # (gate1_pattern, already in the output row) without holding hard-verdict veto
    # power. A too_divergent candidate still gets scored through Gates 3/4 normally;
    # filter on the gate1_pattern column yourself if you want to exclude these.
    if gate4.phylo_status == "insufficient_comparators":
        return "insufficient_comparators"
    if gate3.status == "no_ecological_route":
        return "no_ecological_route"
    if gate4.reassortment_call == "no_conflict" or gate4.phylo_status != "analyzed":
        return "no_conflict"
    if gate4.reassortment_call == "weak_signal_low_bootstrap_support":
        return "insufficient_evidence"
    if gate4.reassortment_call == "supported_reassortment_candidate":
        structural_red_flag = (
            gate2a.tmis_conserved is False
            or gate2a.ph_switch_conserved is False
            or gate2a.waasa_intact is False
            or gate2b.lineage_marker_match is False
            or gate2c.sl_mismatch_high_risk
        )
        if structural_red_flag:
            return "supported_reassortment_candidate_with_structural_warning"
        return "supported_reassortment_candidate"
    return "no_conflict"


# ---------------------------------------------------------------------------
# CLI entry point.
# ---------------------------------------------------------------------------

def read_candidate(candidate_id: str, s_path: Path, m_path: Path, l_path: Path) -> CandidateInput:
    def read_one(path: Path) -> str:
        record = next(iter(ah.parse_fasta(Path(path)).values()), None)
        return record.sequence if record else ""

    return CandidateInput(
        candidate_id=candidate_id,
        s_genomic=read_one(s_path), m_genomic=read_one(m_path), l_genomic=read_one(l_path),
    )


def read_manifest(path: Path) -> List[CandidateInput]:
    candidates = []
    base_dir = Path(path).parent
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            candidates.append(read_candidate(
                row["candidate_id"],
                base_dir / row["s_fasta"], base_dir / row["m_fasta"], base_dir / row["l_fasta"],
            ))
    return candidates


def _build_precomputed(reference: ReferenceData) -> Dict[str, Any]:
    protein_lengths = protein_lengths_by_segment(reference.isolates)
    genomic_sequences = genomic_sequences_by_segment(reference.isolates, reference.genomic_index)
    return {
        "protein_len_ranges": {
            segment: empirical_protein_length_range(lengths)
            for segment, lengths in protein_lengths.items()
        },
        "panhandle_consensus": {
            segment: empirical_panhandle_consensus(sequences)
            for segment, sequences in genomic_sequences.items()
        },
        "protein_profiles": protein_profiles_by_segment(reference.isolates),
    }


def score_candidate(
    candidate: CandidateInput, reference: ReferenceData, precomputed: Dict[str, Any],
    args: argparse.Namespace,
) -> Dict[str, Any]:
    gate0 = gate0_check(
        candidate, precomputed["protein_len_ranges"], precomputed["panhandle_consensus"],
        precomputed["protein_profiles"])

    if not gate0.passed_structure:
        return {
            "candidate_id": candidate.candidate_id, "verdict": "junk_sequence",
            "gate0_failures": ";".join(gate0.failures), "gate0_warnings": ";".join(gate0.warnings),
        }

    gate1 = gate1_check(gate0.orf_proteins, precomputed["protein_profiles"])
    gate2a = gate2_track_a(gate0.orf_proteins.get("M", ""), gate1.implicated_species)
    gate2b = gate2_track_b(gate0.orf_proteins, gate1.implicated_species)
    gate2c = gate2_track_c(gate1)
    gate3 = gate3_check(gate1.implicated_species, reference.isolates)

    candidate_cds = {
        "S": candidate_cds_from_orf(candidate.s_genomic),
        "M": candidate_cds_from_orf(candidate.m_genomic),
        "L": candidate_cds_from_orf(candidate.l_genomic),
    }
    gate4 = gate4_check(
        candidate, candidate_cds, gate1.implicated_species, reference, Path(args.output),
        args.mafft_bin, args.iqtree_bin, phylo_min_taxa=args.phylo_min_taxa,
        bootstrap=args.bootstrap, min_bootstrap_support=args.min_bootstrap_support,
        phylo_seed=args.phylo_seed, phylo_threads=args.phylo_threads,
    )

    verdict = combine_verdict(gate0, gate1, gate2a, gate2b, gate2c, gate3, gate4)

    return {
        "candidate_id": candidate.candidate_id,
        "verdict": verdict,
        "gate0_warnings": ";".join(gate0.warnings),
        "gate0_best_kmer_identity": round(gate0.best_kmer_identity, 4),
        "gate1_pattern": gate1.pattern,
        "gate1_discordant_segment": gate1.discordant_segment,
        "gate1_implicated_species": ";".join(gate1.implicated_species),
        "gate2_waasa_intact": gate2a.waasa_intact,
        "gate2_tmis_conserved": gate2a.tmis_conserved,
        "gate2_ph_switch_conserved": gate2a.ph_switch_conserved,
        "gate2_structural_check_not_available": gate2a.structural_check_not_available,
        "gate2_lineage_marker_match": gate2b.lineage_marker_match,
        "gate2_matched_species": gate2b.matched_species,
        "gate2_sl_mismatch_high_risk": gate2c.sl_mismatch_high_risk,
        "gate2_receptor_class_note": gate2c.receptor_class_note,
        "gate3_status": gate3.status,
        "gate3_detail": gate3.detail,
        "gate4_phylo_status": gate4.phylo_status,
        "gate4_reassortment_call": gate4.reassortment_call,
        "gate4_discordant_segment": gate4.discordant_segment,
        "gate4_bootstrap_support": gate4.bootstrap_support,
    }


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-input", required=True,
                        help="Package root containing all_hanta/ncbi_dataset/data")
    parser.add_argument("--output", required=True, help="Output directory")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--candidate", help="candidate_id:s_fasta,m_fasta,l_fasta")
    group.add_argument("--candidates-manifest", help="CSV: candidate_id,s_fasta,m_fasta,l_fasta")
    parser.add_argument("--mafft-bin", default=None)
    parser.add_argument("--iqtree-bin", default=None)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--min-bootstrap-support", type=float, default=95.0)
    parser.add_argument("--phylo-min-taxa", type=int, default=4)
    parser.add_argument("--phylo-seed", type=int, default=42)
    parser.add_argument("--phylo-threads", type=int, default=1)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    ah.configure_logging(Path(args.output) / "reassortment_detector.log", verbose=args.verbose)

    if args.candidates_manifest:
        candidates = read_manifest(Path(args.candidates_manifest))
    else:
        candidate_id, rest = args.candidate.split(":", 1)
        s_fasta, m_fasta, l_fasta = rest.split(",")
        candidates = [read_candidate(candidate_id, Path(s_fasta), Path(m_fasta), Path(l_fasta))]

    reference = load_reference_data(Path(args.reference_input))
    precomputed = _build_precomputed(reference)

    mafft_bin = ph.resolve_binary(args.mafft_bin, "mafft")
    iqtree_bin = ph.resolve_binary(args.iqtree_bin, "iqtree")
    args.mafft_bin, args.iqtree_bin = mafft_bin, iqtree_bin

    rows = [score_candidate(c, reference, precomputed, args) for c in candidates]

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "candidate_reassortment_table.csv"
    fieldnames = sorted({key for row in rows for key in row})
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(f"{row['candidate_id']}: {row['verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())


# ============================================================================
# Tests
# ============================================================================
#
# Merged from test_detector_gates.py, test_detector_reference.py,
# test_reassortment_detector.py, and test_reassortment_reference_data.py.
# Run with: python3 -m unittest reassortment_detector -v

import random
import tempfile
import unittest


class LongestOrfProteinTest(unittest.TestCase):
    def test_finds_orf_on_forward_frame(self):
        # ATG ... stop, 30 nt = 9 codons incl. stop -> 8 aa protein
        seq = "ATG" + "GCT" * 7 + "TAA"
        self.assertEqual(longest_orf_protein(seq), "A" * 7)

    def test_picks_longest_among_multiple_candidate_orfs(self):
        short = "ATG" + "GCT" * 2 + "TAA"
        long_ = "ATG" + "GCT" * 10 + "TAA"
        seq = short + "CCCCC" + long_
        self.assertEqual(longest_orf_protein(seq), "A" * 10)

    def test_no_stop_codon_returns_empty(self):
        self.assertEqual(longest_orf_protein("ATG" + "GCT" * 5), "")

    def test_empty_sequence_returns_empty(self):
        self.assertEqual(longest_orf_protein(""), "")


class EmpiricalDerivationTest(unittest.TestCase):
    def test_protein_length_range_is_min_max(self):
        self.assertEqual(empirical_protein_length_range([430, 428, 435, 431]), (428, 435))

    def test_protein_length_range_empty_list_is_wide_open(self):
        self.assertEqual(empirical_protein_length_range([]), (0, 10**9))

    def test_panhandle_consensus_is_majority_vote_per_position(self):
        seqs = ["TAGTAGTAGAC" + "N" * 20, "TAGTAGTAGAC" + "N" * 20, "TAGTAGTAGAT" + "N" * 20]
        five_prime, three_prime = empirical_panhandle_consensus(seqs, n=11)
        self.assertEqual(five_prime, "TAGTAGTAGAC")  # majority at the last position


class Gate0CheckTest(unittest.TestCase):
    def setUp(self):
        self.protein_len_ranges = {"S": (420, 440), "M": (1100, 1160), "L": (2100, 2200)}
        self.panhandle = {
            "S": ("TAGTAGTAGAC", "GTCTACTACTA"),
            "M": ("TAGTAGTAGAC", "GTCTACTACTA"),
            "L": ("TAGTAGTAGAC", "GTCTACTACTA"),
        }
        # One reference protein profile per segment, used only for the identity floor.
        self.reference_profiles = {
            "S": [("Orthohantavirus puumalaense", "M" * 430)],
            "M": [("Orthohantavirus puumalaense", "M" * 1130)],
            "L": [("Orthohantavirus puumalaense", "M" * 2150)],
        }

    def _valid_orf_genomic(self, aa_len, panhandle5, panhandle3):
        body = "ATG" + "GCT" * aa_len + "TAA"
        return panhandle5 + body + panhandle3

    def test_correct_length_and_terminals_and_high_identity_passes_structure(self):
        candidate = CandidateInput(
            candidate_id="C1",
            s_genomic=self._valid_orf_genomic(428, "TAGTAGTAGAC", "GTCTACTACTA"),
            m_genomic=self._valid_orf_genomic(1128, "TAGTAGTAGAC", "GTCTACTACTA"),
            l_genomic=self._valid_orf_genomic(2148, "TAGTAGTAGAC", "GTCTACTACTA"),
        )
        result = gate0_check(
            candidate, self.protein_len_ranges, self.panhandle, self.reference_profiles)
        self.assertTrue(result.passed_structure)
        self.assertEqual(result.failures, [])

    def test_wrong_orf_length_fails_structure_with_named_reason(self):
        candidate = CandidateInput(
            candidate_id="C2",
            s_genomic=self._valid_orf_genomic(5, "TAGTAGTAGAC", "GTCTACTACTA"),  # way too short
            m_genomic=self._valid_orf_genomic(1128, "TAGTAGTAGAC", "GTCTACTACTA"),
            l_genomic=self._valid_orf_genomic(2148, "TAGTAGTAGAC", "GTCTACTACTA"),
        )
        result = gate0_check(
            candidate, self.protein_len_ranges, self.panhandle, self.reference_profiles)
        self.assertFalse(result.passed_structure)
        self.assertTrue(any("S" in f and "orf" in f.lower() for f in result.failures))

    def test_terminal_mismatch_at_literal_first_base_is_warning_not_failure(self):
        # Positions 2-12 (the checked consensus window) still match; only the very
        # first base differs -- per Liu et al. 2024's note about GenBank-added terminal
        # bases, this must be a warning, not a Gate 0 failure.
        mutated = "C" + "AGTAGTAGAC"
        candidate = CandidateInput(
            candidate_id="C3",
            s_genomic=self._valid_orf_genomic(428, mutated, "GTCTACTACTA"),
            m_genomic=self._valid_orf_genomic(1128, "TAGTAGTAGAC", "GTCTACTACTA"),
            l_genomic=self._valid_orf_genomic(2148, "TAGTAGTAGAC", "GTCTACTACTA"),
        )
        result = gate0_check(
            candidate, self.protein_len_ranges, self.panhandle, self.reference_profiles)
        self.assertTrue(result.passed_structure)
        self.assertTrue(any("terminal" in w.lower() for w in result.warnings))

    def test_low_identity_to_reference_fails_floor_but_not_structure(self):
        candidate = CandidateInput(
            candidate_id="C4",
            s_genomic=self._valid_orf_genomic(428, "TAGTAGTAGAC", "GTCTACTACTA"),
            m_genomic=self._valid_orf_genomic(1128, "TAGTAGTAGAC", "GTCTACTACTA"),
            l_genomic=self._valid_orf_genomic(2148, "TAGTAGTAGAC", "GTCTACTACTA"),
        )
        # Reference profiles totally unlike the candidate's all-alanine ORFs.
        far_reference = {
            "S": [("X", "WYWYWYWYWY" * 43)],
            "M": [("X", "WYWYWYWYWY" * 113)],
            "L": [("X", "WYWYWYWYWY" * 215)],
        }
        result = gate0_check(candidate, self.protein_len_ranges, self.panhandle, far_reference)
        self.assertTrue(result.passed_structure)
        self.assertFalse(result.passed_identity_floor)


class NearestSpeciesByKmerTest(unittest.TestCase):
    def test_returns_top_n_distinct_species_sorted_by_distance(self):
        reference = [
            ("SpA", "MMMMMMMMMM"),
            ("SpA", "MMMMMMMMMA"),  # near-duplicate of SpA, must not count as a 2nd species
            ("SpB", "WYWYWYWYWY"),
            ("SpC", "MMMMMMMMMW"),
        ]
        hits = nearest_species_by_kmer("MMMMMMMMMM", reference, top_n=3)
        self.assertEqual(len(hits), 3)
        self.assertEqual(hits[0].species, "SpA")
        self.assertLessEqual(hits[0].distance, hits[1].distance)
        self.assertLessEqual(hits[1].distance, hits[2].distance)

    def test_empty_reference_returns_empty(self):
        self.assertEqual(nearest_species_by_kmer("MMM", [], top_n=3), [])


class Gate1CheckTest(unittest.TestCase):
    def setUp(self):
        self.reference = {
            "S": [("SpA", "M" * 50), ("SpB", "W" * 50)],
            "M": [("SpA", "M" * 50), ("SpB", "W" * 50)],
            "L": [("SpA", "M" * 50), ("SpB", "W" * 50)],
        }

    def test_all_segments_concordant(self):
        orf_proteins = {"S": "M" * 50, "M": "M" * 50, "L": "M" * 50}
        result = gate1_check(orf_proteins, self.reference)
        self.assertEqual(result.pattern, "concordant")
        self.assertEqual(result.discordant_segment, "")
        self.assertEqual(result.implicated_species, ["SpA"])

    def test_m_only_discordant_flags_m_pattern(self):
        orf_proteins = {"S": "M" * 50, "M": "W" * 50, "L": "M" * 50}
        result = gate1_check(orf_proteins, self.reference)
        self.assertEqual(result.discordant_segment, "M")
        self.assertEqual(result.pattern, "m_only")
        self.assertEqual(set(result.implicated_species), {"SpA", "SpB"})

    def test_s_discordant_flags_s_or_l_involving_pattern(self):
        orf_proteins = {"S": "W" * 50, "M": "M" * 50, "L": "M" * 50}
        result = gate1_check(orf_proteins, self.reference)
        self.assertEqual(result.discordant_segment, "S")
        self.assertEqual(result.pattern, "s_or_l_involving")


class Gate2TrackATest(unittest.TestCase):
    def test_candidate_identical_to_andv_anchor_conserves_everything(self):
        andv_seq = STRUCTURAL_ANCHORS["Orthohantavirus andesense"]["sequence"]
        result = gate2_track_a(andv_seq, ["Orthohantavirus andesense"])
        self.assertEqual(result.anchor_species, "Orthohantavirus andesense")
        self.assertEqual(result.anchor_basis, "solved")
        self.assertTrue(result.waasa_intact)
        self.assertTrue(result.tmis_conserved)
        self.assertTrue(result.ph_switch_conserved)
        self.assertFalse(result.structural_check_not_available)

    def test_candidate_identical_to_puuv_anchor_also_conserves_everything(self):
        puuv_seq = STRUCTURAL_ANCHORS["Orthohantavirus puumalaense"]["sequence"]
        result = gate2_track_a(puuv_seq, ["Orthohantavirus puumalaense"])
        self.assertTrue(result.waasa_intact)
        self.assertTrue(result.tmis_conserved)
        self.assertTrue(result.ph_switch_conserved)

    def test_disrupted_tmis_residue_is_detected(self):
        andv_seq = STRUCTURAL_ANCHORS["Orthohantavirus andesense"]["sequence"]
        mutated = list(andv_seq)
        mutated[766 - 1] = "A"  # W766 -> A, breaking the TMIS hydrophobic residue
        result = gate2_track_a("".join(mutated), ["Orthohantavirus andesense"])
        self.assertFalse(result.tmis_conserved)

    def test_unimplicated_species_without_a_table_entry_is_not_available(self):
        result = gate2_track_a("M" * 1100, ["Orthohantavirus seoulense"])
        self.assertTrue(result.structural_check_not_available)
        self.assertIsNone(result.waasa_intact)


class Gate2TrackBTest(unittest.TestCase):
    def test_no_marker_set_for_implicated_species_is_not_informative(self):
        result = gate2_track_b({"S": "M" * 400, "M": "M" * 1100}, ["Orthohantavirus hantanense"])
        self.assertIsNone(result.lineage_marker_match)
        self.assertEqual(result.marker_set_used, "")

    def test_empty_marker_table_never_raises(self):
        # LINEAGE_MARKER_TABLE ships empty by design (spec section 7) -- this must not
        # crash regardless of which species Gate 1 implicates.
        for species in ["Orthohantavirus andesense", "Orthohantavirus puumalaense", "Unlisted sp."]:
            result = gate2_track_b({"S": "M" * 400, "M": "M" * 1100}, [species])
            self.assertIsNone(result.lineage_marker_match)


class Gate2TrackCTest(unittest.TestCase):
    def test_no_flag_when_gate1_concordant(self):
        gate1 = Gate1Result(implicated_species=["Orthohantavirus andesense"],
                            discordant_segment="", pattern="concordant")
        result = gate2_track_c(gate1)
        self.assertFalse(result.sl_mismatch_high_risk)

    def test_s_discordant_flags_sl_mismatch_high_risk(self):
        gate1 = Gate1Result(implicated_species=["Orthohantavirus andesense", "Orthohantavirus puumalaense"],
                            discordant_segment="S", pattern="s_or_l_involving")
        result = gate2_track_c(gate1)
        self.assertTrue(result.sl_mismatch_high_risk)

    def test_m_only_never_flags_sl_mismatch(self):
        gate1 = Gate1Result(implicated_species=["Orthohantavirus andesense", "Orthohantavirus hantanense"],
                            discordant_segment="M", pattern="m_only")
        result = gate2_track_c(gate1)
        self.assertFalse(result.sl_mismatch_high_risk)

    def test_m_only_new_world_to_old_world_crossing_is_noted(self):
        gate1 = Gate1Result(implicated_species=["Orthohantavirus andesense", "Orthohantavirus hantanense"],
                            discordant_segment="M", pattern="m_only")
        result = gate2_track_c(gate1)
        self.assertIn(result.receptor_class_note, {"new_world_to_old_world", "old_world_to_new_world"})

    def test_m_only_within_same_clade_is_same_clade(self):
        gate1 = Gate1Result(implicated_species=["Orthohantavirus hantanense", "Orthohantavirus puumalaense"],
                            discordant_segment="M", pattern="m_only")
        result = gate2_track_c(gate1)
        self.assertEqual(result.receptor_class_note, "same_clade")


def _make_ecology_isolate(species, country, host_order, year, lab=False, vaccine=False):
    return ah.IsolateRecord(
        isolate_key=f"{species}|{country}|{year}", species_name=species, species_taxid="999",
        species_status="classified", biosample="", isolate="", strain="", host_name="Test host",
        host_taxid="1", host_order=host_order, host_family="", host_role="reservoir_rodent",
        collection_date=year, year=year, location="", country=country, region="",
        grouping_tier="tier_2_species_isolate", grouping_confidence=1.0,
        s_accession="", m_accession="", l_accession="",
        s_completeness="", m_completeness="", l_completeness="",
        s_nt_length=0, m_nt_length=0, l_nt_length=0,
        s_protein_aa="", m_protein_aa="", l_protein_aa="", concatenated_protein_aa="",
        segments_present=[], coding_complete_segments=[], triplet_status="", has_complete_triplet=False,
        is_lab_derived=lab, is_vaccine_strain=vaccine, has_nss_orf=False,
        metadata_conflict=False, taxonomic_incongruence=False, cross_species_partners=[],
        kmer_cluster_id="", s_protein_accessions=[], m_protein_accessions=[],
        l_protein_accessions=[], all_accessions=[], warnings=[],
    )


class Gate3CheckTest(unittest.TestCase):
    def test_single_implicated_species_needs_no_overlap_check(self):
        isolates = [_make_ecology_isolate("SpA", "Finland", "Rodentia", "2010")]
        result = gate3_check(["SpA"], isolates)
        self.assertEqual(result.status, "ecology_unknown")

    def test_sympatric_host_order_and_country_is_plausible_overlap(self):
        isolates = [
            _make_ecology_isolate("SpA", "Finland", "Rodentia", "2010"),
            _make_ecology_isolate("SpB", "Finland", "Rodentia", "2011"),
        ]
        result = gate3_check(["SpA", "SpB"], isolates)
        self.assertEqual(result.status, "ecological_overlap_plausible")

    def test_different_continents_and_host_orders_is_no_route(self):
        isolates = [
            _make_ecology_isolate("SpA", "Finland", "Rodentia", "2010"),
            _make_ecology_isolate("SpC", "Guinea", "Eulipotyphla", "1990"),
        ]
        result = gate3_check(["SpA", "SpC"], isolates)
        self.assertEqual(result.status, "no_ecological_route")

    def test_missing_metadata_for_a_parent_is_unknown_not_no_route(self):
        isolates = [_make_ecology_isolate("SpA", "Finland", "Rodentia", "2010")]
        # SpB has no reference isolates at all -- metadata absent, not evidence of
        # non-overlap.
        result = gate3_check(["SpA", "SpB"], isolates)
        self.assertEqual(result.status, "ecology_unknown")

    def test_lab_and_vaccine_isolates_excluded_from_consideration(self):
        isolates = [
            _make_ecology_isolate("SpA", "Finland", "Rodentia", "2010", lab=True),
            _make_ecology_isolate("SpB", "Finland", "Rodentia", "2011"),
        ]
        result = gate3_check(["SpA", "SpB"], isolates)
        self.assertEqual(result.status, "ecology_unknown")  # SpA has no WILD isolates


def _make_isolate(key, species, s_protein, m_protein, l_protein):
    """Shared by both the Gate 0-4 tests (below) and the reference-loading tests --
    identical fixture, kept as one copy after the merge that used to duplicate it across
    test_detector_gates.py and test_detector_reference.py.
    """
    return ah.IsolateRecord(
        isolate_key=key, species_name=species, species_taxid="999", species_status="classified",
        biosample="", isolate=key, strain="", host_name="Test host", host_taxid="1",
        host_order="Rodentia", host_family="", host_role="reservoir_rodent",
        collection_date="", year="2020", location="", country="Testland", region="",
        grouping_tier="tier_2_species_isolate", grouping_confidence=1.0,
        s_accession="SACC", m_accession="MACC", l_accession="LACC",
        s_completeness="complete", m_completeness="complete", l_completeness="complete",
        s_nt_length=len(s_protein) * 3, m_nt_length=len(m_protein) * 3, l_nt_length=len(l_protein) * 3,
        s_protein_aa=s_protein, m_protein_aa=m_protein, l_protein_aa=l_protein,
        concatenated_protein_aa=s_protein + m_protein + l_protein,
        segments_present=["S", "M", "L"], coding_complete_segments=["S", "M", "L"],
        triplet_status="coding_complete_triplet", has_complete_triplet=True,
        is_lab_derived=False, is_vaccine_strain=False, has_nss_orf=False,
        metadata_conflict=False, taxonomic_incongruence=False, cross_species_partners=[],
        kmer_cluster_id="", s_protein_accessions=[], m_protein_accessions=[],
        l_protein_accessions=[], all_accessions=["SACC", "MACC", "LACC"], warnings=[],
    )


def _make_reference_isolate_with_cds(key, species, seed):
    rng = random.Random(seed)
    bases = "ACGT"

    def rand_seq(n):
        return "".join(rng.choice(bases) for _ in range(n))

    isolate = _make_isolate(key, species, "M" * 60, "M" * 60, "M" * 60)
    isolate.s_accession, isolate.m_accession, isolate.l_accession = f"{key}S", f"{key}M", f"{key}L"
    isolate.all_accessions = [isolate.s_accession, isolate.m_accession, isolate.l_accession]
    cds = {"S": rand_seq(150), "M": rand_seq(150), "L": rand_seq(150)}
    return isolate, cds


class CandidateCdsFromOrfTest(unittest.TestCase):
    def test_returns_the_orfs_own_nucleotide_span(self):
        genomic = "TTTTT" + "ATG" + "GCT" * 5 + "TAA" + "GGGGG"
        cds = candidate_cds_from_orf(genomic)
        self.assertTrue(cds.startswith("ATG"))
        self.assertTrue(cds.endswith("TAA"))
        self.assertEqual(len(cds), 3 + 15 + 3)


class BuildComparatorMembersTest(unittest.TestCase):
    def test_builds_members_from_reference_isolates_with_cds(self):
        isolates, cds_index = [], {}
        for i in range(5):
            isolate, cds = _make_reference_isolate_with_cds(f"SpA_{i}", "SpA", seed=i)
            isolates.append(isolate)
            for segment, accession in (("S", isolate.s_accession), ("M", isolate.m_accession),
                                       ("L", isolate.l_accession)):
                cds_index[accession] = [("nucleocapsid" if segment == "S" else
                                         "glycoprotein" if segment == "M" else "polymerase",
                                         cds[segment])]
        reference = ReferenceData(isolates=isolates, cds_index=cds_index, genomic_index={})
        members, status = build_comparator_members(["SpA"], reference, phylo_min_taxa=4)
        self.assertEqual(status, "")
        self.assertEqual(len(members), 5)
        self.assertTrue(all(set(m.sequences) == {"S", "M", "L"} for m in members))

    def test_too_few_comparators_reports_insufficient(self):
        isolates, cds_index = [], {}
        for i in range(2):
            isolate, cds = _make_reference_isolate_with_cds(f"SpA_{i}", "SpA", seed=i)
            isolates.append(isolate)
            for segment, accession in (("S", isolate.s_accession), ("M", isolate.m_accession),
                                       ("L", isolate.l_accession)):
                cds_index[accession] = [("nucleocapsid" if segment == "S" else
                                         "glycoprotein" if segment == "M" else "polymerase",
                                         cds[segment])]
        reference = ReferenceData(isolates=isolates, cds_index=cds_index, genomic_index={})
        members, status = build_comparator_members(["SpA"], reference, phylo_min_taxa=4)
        self.assertEqual(status, "insufficient_comparators")
        self.assertEqual(members, [])


class Gate4CheckEndToEndTest(unittest.TestCase):
    """Mirrors test_all_hanta.py's PhylogeneticsScopeEndToEndTest fake-binary pattern --
    mafft/iqtree are not installed in this environment (confirmed via `which` during
    design), so this uses the same structurally-valid-but-meaningless fake scripts.
    """

    _FAKE_MAFFT = '''#!/usr/bin/env python3
import sys
path = sys.argv[-1]
seqs, label, buf = [], None, []
with open(path) as f:
    for line in f:
        line = line.rstrip("\\n")
        if line.startswith(">"):
            if label is not None:
                seqs.append((label, "".join(buf)))
            label, buf = line[1:], []
        else:
            buf.append(line)
    if label is not None:
        seqs.append((label, "".join(buf)))
maxlen = max(len(s) for _, s in seqs)
for label, s in seqs:
    print(f">{label}")
    print(s.ljust(maxlen, "-"))
'''

    _FAKE_IQTREE = '''#!/usr/bin/env python3
import random, sys
argv = sys.argv[1:]
def get_opt(flag):
    return argv[argv.index(flag) + 1] if flag in argv else None
alignment, prefix, seed = get_opt("-s"), get_opt("--prefix"), get_opt("--seed")
random.seed(int(seed) if seed else 0)
labels = []
with open(alignment) as f:
    for line in f:
        if line.startswith(">"):
            labels.append(line[1:].strip())
nodes = list(labels)
random.shuffle(nodes)
cur = f"{nodes[0]}:{random.random():.4f}"
for label in nodes[1:-1]:
    cur = f"({cur},{label}:{random.random():.4f}){random.randint(50, 100)}:{random.random():.4f}"
cur = f"({cur},{nodes[-1]}:{random.random():.4f})"
with open(f"{prefix}.treefile", "w") as f:
    f.write(cur + ";\\n")
'''

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        tmp_path = Path(self.tmpdir.name)
        self.mafft_bin = tmp_path / "fake_mafft.py"
        self.iqtree_bin = tmp_path / "fake_iqtree.py"
        self.mafft_bin.write_text(self._FAKE_MAFFT)
        self.iqtree_bin.write_text(self._FAKE_IQTREE)
        self.mafft_bin.chmod(0o755)
        self.iqtree_bin.chmod(0o755)
        self.output_dir = tmp_path / "output"
        self.output_dir.mkdir()

        isolates, cds_index = [], {}
        for i in range(5):
            isolate, cds = _make_reference_isolate_with_cds(f"SpA_{i}", "SpA", seed=i)
            isolates.append(isolate)
            for segment, accession in (("S", isolate.s_accession), ("M", isolate.m_accession),
                                       ("L", isolate.l_accession)):
                cds_index[accession] = [("nucleocapsid" if segment == "S" else
                                         "glycoprotein" if segment == "M" else "polymerase",
                                         cds[segment])]
        self.reference = ReferenceData(isolates=isolates, cds_index=cds_index, genomic_index={})

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_candidate_gets_analyzed_and_placed_in_the_tree(self):
        candidate = CandidateInput(candidate_id="CAND1", s_genomic="", m_genomic="", l_genomic="")
        candidate_cds = {"S": "ACGT" * 40, "M": "TGCA" * 40, "L": "AATT" * 40}
        result = gate4_check(
            candidate, candidate_cds, ["SpA"], self.reference, self.output_dir,
            str(self.mafft_bin), str(self.iqtree_bin), phylo_min_taxa=4, bootstrap=10,
            min_bootstrap_support=95.0, phylo_seed=42, phylo_threads=1,
        )
        self.assertEqual(result.phylo_status, "analyzed")
        self.assertIn(result.reassortment_call, {"no_conflict", "supported_reassortment_candidate",
                                                  "weak_signal_low_bootstrap_support"})

    def test_missing_cds_for_one_segment_fails_without_running_mafft(self):
        # candidate_cds always has all three keys (reassortment_detector.py builds it as
        # a dict literal with S/M/L always present) -- the guard must check for an empty
        # VALUE, not key membership, or a genuinely missing CDS silently reaches mafft
        # with an empty-sequence FASTA record instead of failing with this diagnostic.
        candidate = CandidateInput(candidate_id="CAND1", s_genomic="", m_genomic="", l_genomic="")
        candidate_cds = {"S": "ACGT" * 40, "M": "", "L": "AATT" * 40}
        result = gate4_check(
            candidate, candidate_cds, ["SpA"], self.reference, self.output_dir,
            str(self.mafft_bin), str(self.iqtree_bin), phylo_min_taxa=4, bootstrap=10,
            min_bootstrap_support=95.0, phylo_seed=42, phylo_threads=1,
        )
        self.assertEqual(result.phylo_status, "failed")
        self.assertEqual(result.reassortment_call, "not_analyzed")
        self.assertIn("missing CDS", result.detail.get("error", ""))

    def test_insufficient_comparators_short_circuits_without_running_mafft(self):
        candidate = CandidateInput(candidate_id="CAND1", s_genomic="", m_genomic="", l_genomic="")
        candidate_cds = {"S": "ACGT" * 40, "M": "TGCA" * 40, "L": "AATT" * 40}
        empty_reference = ReferenceData(isolates=[], cds_index={}, genomic_index={})
        result = gate4_check(
            candidate, candidate_cds, ["SpUnknown"], empty_reference, self.output_dir,
            str(self.mafft_bin), str(self.iqtree_bin), phylo_min_taxa=4,
        )
        self.assertEqual(result.phylo_status, "insufficient_comparators")


class CombineVerdictTest(unittest.TestCase):
    def setUp(self):
        self.gate0_pass = Gate0Result(passed_structure=True, passed_identity_floor=True)
        self.gate2a_clean = Gate2TrackAResult(structural_check_not_available=True)
        self.gate2b_clean = Gate2TrackBResult()
        self.gate2c_clean = Gate2TrackCResult()

    def test_gate0_structural_failure_is_junk_sequence(self):
        gate0 = Gate0Result(passed_structure=False, passed_identity_floor=True, failures=["S: bad"])
        verdict = combine_verdict(
            gate0, Gate1Result(), self.gate2a_clean, self.gate2b_clean, self.gate2c_clean,
            Gate3Result(), Gate4Result())
        self.assertEqual(verdict, "junk_sequence")

    def test_gate0_identity_floor_failure_is_divergent_novel(self):
        gate0 = Gate0Result(passed_structure=True, passed_identity_floor=False)
        verdict = combine_verdict(
            gate0, Gate1Result(), self.gate2a_clean, self.gate2b_clean, self.gate2c_clean,
            Gate3Result(), Gate4Result())
        self.assertEqual(verdict, "divergent_novel_candidate")

    def test_undetectable_by_design_pattern_wins_regardless_of_gate4(self):
        gate1 = Gate1Result(pattern="undetectable_by_design", discordant_segment="M")
        gate4 = Gate4Result(phylo_status="analyzed", reassortment_call="no_conflict")
        verdict = combine_verdict(
            self.gate0_pass, gate1, self.gate2a_clean, self.gate2b_clean, self.gate2c_clean,
            Gate3Result(), gate4)
        self.assertEqual(verdict, "undetectable_by_design")

    def test_insufficient_comparators_from_gate4(self):
        gate4 = Gate4Result(phylo_status="insufficient_comparators")
        verdict = combine_verdict(
            self.gate0_pass, Gate1Result(), self.gate2a_clean, self.gate2b_clean,
            self.gate2c_clean, Gate3Result(), gate4)
        self.assertEqual(verdict, "insufficient_comparators")

    def test_no_ecological_route_from_gate3(self):
        gate4 = Gate4Result(phylo_status="analyzed", reassortment_call="no_conflict")
        verdict = combine_verdict(
            self.gate0_pass, Gate1Result(pattern="m_only", discordant_segment="M"),
            self.gate2a_clean, self.gate2b_clean, self.gate2c_clean,
            Gate3Result(status="no_ecological_route"), gate4)
        self.assertEqual(verdict, "no_ecological_route")

    def test_ecology_unknown_does_not_block_a_supported_call(self):
        gate4 = Gate4Result(phylo_status="analyzed", reassortment_call="supported_reassortment_candidate",
                            bootstrap_support=98.0)
        verdict = combine_verdict(
            self.gate0_pass, Gate1Result(pattern="m_only", discordant_segment="M"),
            self.gate2a_clean, self.gate2b_clean, self.gate2c_clean,
            Gate3Result(status="ecology_unknown"), gate4)
        self.assertEqual(verdict, "supported_reassortment_candidate")

    def test_weak_bootstrap_is_insufficient_evidence(self):
        gate4 = Gate4Result(phylo_status="analyzed", reassortment_call="weak_signal_low_bootstrap_support",
                            bootstrap_support=60.0)
        verdict = combine_verdict(
            self.gate0_pass, Gate1Result(pattern="m_only", discordant_segment="M"),
            self.gate2a_clean, self.gate2b_clean, self.gate2c_clean,
            Gate3Result(status="ecological_overlap_plausible"), gate4)
        self.assertEqual(verdict, "insufficient_evidence")

    def test_supported_with_structural_warning_when_track_a_flags_a_red_flag(self):
        gate4 = Gate4Result(phylo_status="analyzed", reassortment_call="supported_reassortment_candidate",
                            bootstrap_support=99.0)
        gate2a = Gate2TrackAResult(structural_check_not_available=False, tmis_conserved=False)
        verdict = combine_verdict(
            self.gate0_pass, Gate1Result(pattern="m_only", discordant_segment="M"), gate2a,
            self.gate2b_clean, self.gate2c_clean,
            Gate3Result(status="ecological_overlap_plausible"), gate4)
        self.assertEqual(verdict, "supported_reassortment_candidate_with_structural_warning")

    def test_no_conflict_when_gate4_finds_nothing(self):
        gate4 = Gate4Result(phylo_status="analyzed", reassortment_call="no_conflict")
        verdict = combine_verdict(
            self.gate0_pass, Gate1Result(pattern="concordant"), self.gate2a_clean,
            self.gate2b_clean, self.gate2c_clean, Gate3Result(), gate4)
        self.assertEqual(verdict, "no_conflict")


class ProteinLengthsBySegmentTest(unittest.TestCase):
    def test_collects_lengths_per_segment_skipping_empty(self):
        isolates = [
            _make_isolate("I1", "SpA", "M" * 430, "M" * 1130, ""),
            _make_isolate("I2", "SpA", "M" * 428, "M" * 1135, "M" * 2150),
        ]
        lengths = protein_lengths_by_segment(isolates)
        self.assertEqual(sorted(lengths["S"]), [428, 430])
        self.assertEqual(sorted(lengths["M"]), [1130, 1135])
        self.assertEqual(lengths["L"], [2150])


class ProteinProfilesBySegmentTest(unittest.TestCase):
    def test_pairs_species_with_protein_per_segment(self):
        isolates = [_make_isolate("I1", "SpA", "MMM", "GGG", "LLL")]
        profiles = protein_profiles_by_segment(isolates)
        self.assertEqual(profiles["S"], [("SpA", "MMM")])
        self.assertEqual(profiles["M"], [("SpA", "GGG")])
        self.assertEqual(profiles["L"], [("SpA", "LLL")])


class LoadReferenceDataIntegrationTest(unittest.TestCase):
    """Exercises load_reference_data against the real all_hanta/ NCBI package already
    present in this repo -- skipped if it isn't there (e.g. a fresh checkout before
    Step 0 of the README has been run). This is plumbing around already-tested
    all_hanta.py functions, so a synthetic fixture package would just re-test
    all_hanta.py itself; using the real package matches how the rest of this repo
    validates the pipeline (see README's own reproducibility instructions).
    """

    def test_loads_real_package_when_present(self):
        root = Path(__file__).parent / "all_hanta"
        if not (root / "ncbi_dataset" / "data" / "data_report.jsonl").exists():
            self.skipTest("all_hanta/ NCBI package not present in this checkout")
        data = load_reference_data(Path(__file__).parent)
        self.assertGreater(len(data.isolates), 1000)
        self.assertGreater(len(data.cds_index), 1000)
        self.assertTrue(any(i.has_complete_triplet for i in data.isolates))


def _write_fasta(path: Path, header: str, sequence: str) -> None:
    path.write_text(f">{header}\n{sequence}\n")


class ReadCandidateTest(unittest.TestCase):
    def test_reads_three_fasta_files_into_candidate_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            s_path, m_path, l_path = tmp_path / "s.fa", tmp_path / "m.fa", tmp_path / "l.fa"
            _write_fasta(s_path, "s_header", "ACGT")
            _write_fasta(m_path, "m_header", "TGCA")
            _write_fasta(l_path, "l_header", "AATT")
            candidate = read_candidate("CAND1", s_path, m_path, l_path)
            self.assertEqual(candidate.candidate_id, "CAND1")
            self.assertEqual(candidate.s_genomic, "ACGT")
            self.assertEqual(candidate.m_genomic, "TGCA")
            self.assertEqual(candidate.l_genomic, "AATT")


class ReadManifestTest(unittest.TestCase):
    def test_reads_multiple_candidates_from_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            for cid in ("A", "B"):
                for seg, seq in (("s", "ACGT"), ("m", "TGCA"), ("l", "AATT")):
                    _write_fasta(tmp_path / f"{cid}_{seg}.fa", f"{cid}_{seg}", seq)
            manifest = tmp_path / "manifest.csv"
            with manifest.open("w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["candidate_id", "s_fasta", "m_fasta", "l_fasta"])
                for cid in ("A", "B"):
                    writer.writerow([cid, f"{cid}_s.fa", f"{cid}_m.fa", f"{cid}_l.fa"])
            candidates = read_manifest(manifest)
            self.assertEqual([c.candidate_id for c in candidates], ["A", "B"])
            self.assertEqual(candidates[0].s_genomic, "ACGT")


class ParseArgsTest(unittest.TestCase):
    def test_single_candidate_form(self):
        args = parse_args([
            "--reference-input", "all_hanta", "--output", "output",
            "--candidate", "CAND01:s.fa,m.fa,l.fa",
        ])
        self.assertEqual(args.candidate, "CAND01:s.fa,m.fa,l.fa")
        self.assertEqual(args.phylo_seed, 42)

    def test_batch_manifest_form(self):
        args = parse_args([
            "--reference-input", "all_hanta", "--output", "output",
            "--candidates-manifest", "manifest.csv",
        ])
        self.assertEqual(args.candidates_manifest, "manifest.csv")

    def test_requires_one_of_candidate_or_manifest(self):
        with self.assertRaises(SystemExit):
            parse_args(["--reference-input", "all_hanta", "--output", "output"])


class EndToEndRealDataTest(unittest.TestCase):
    """Builds two synthetic candidates from REAL held-out reference sequences (pulled
    from the actual all_hanta/ package, not fabricated) and scores them with fake
    mafft/iqtree binaries (same rationale as Gate4CheckEndToEndTest above -- neither
    binary is installed in this environment). Skipped if all_hanta/ isn't present in
    this checkout.

    Candidate 1 ("no-conflict control"): all three segments' CDS come from one real
    isolate -- Gate 1/4 should find no discordance.
    Candidate 2 ("artificial M-swap"): S and L CDS come from one isolate of species A,
    M CDS comes from a different isolate of species B -- Gate 1 should flag M as
    discordant, and Gate 4 should have both species' isolates in its comparator set.
    """

    def setUp(self):
        root = Path(__file__).parent / "all_hanta"
        if not (root / "ncbi_dataset" / "data" / "data_report.jsonl").exists():
            self.skipTest("all_hanta/ NCBI package not present in this checkout")
        self.reference = load_reference_data(root)

        self.tmpdir = tempfile.TemporaryDirectory()
        tmp_path = Path(self.tmpdir.name)
        self.mafft_bin = tmp_path / "fake_mafft.py"
        self.iqtree_bin = tmp_path / "fake_iqtree.py"
        # Reuses the exact fake-binary scripts from Gate4CheckEndToEndTest above -- see
        # that class for the source. Both classes live in this same module now, so no
        # cross-module import is needed to reach it.
        self.mafft_bin.write_text(Gate4CheckEndToEndTest._FAKE_MAFFT)
        self.iqtree_bin.write_text(Gate4CheckEndToEndTest._FAKE_IQTREE)
        self.mafft_bin.chmod(0o755)
        self.iqtree_bin.chmod(0o755)
        self.output_dir = tmp_path / "output"
        self.output_dir.mkdir()

    def tearDown(self):
        self.tmpdir.cleanup()

    def _args(self):
        return argparse.Namespace(
            output=str(self.output_dir), mafft_bin=str(self.mafft_bin), iqtree_bin=str(self.iqtree_bin),
            bootstrap=10, min_bootstrap_support=95.0, phylo_min_taxa=4, phylo_seed=42, phylo_threads=1,
        )

    def _find_complete_triplet_isolates(self, min_count):
        by_species = {}
        for isolate in self.reference.isolates:
            if not isolate.has_complete_triplet:
                continue
            by_species.setdefault(isolate.species_name, []).append(isolate)
        return next((v for v in by_species.values() if len(v) >= min_count), None)

    def _genomic(self, accession):
        record = self.reference.genomic_index.get(accession) or self.reference.genomic_index.get(
            accession.split(".")[0])
        return record.sequence if record else ""

    def test_no_conflict_control_from_one_real_isolate(self):
        members = self._find_complete_triplet_isolates(min_count=5)
        if members is None:
            self.skipTest("no species in this dataset has >=5 complete-triplet isolates")
        source_isolate = members[0]
        candidate = CandidateInput(
            candidate_id="NOCONFLICT",
            s_genomic=self._genomic(source_isolate.s_accession),
            m_genomic=self._genomic(source_isolate.m_accession),
            l_genomic=self._genomic(source_isolate.l_accession),
        )
        precomputed = _build_precomputed(self.reference)
        result = score_candidate(candidate, self.reference, precomputed, self._args())
        self.assertIn(result["verdict"], {"no_conflict", "insufficient_evidence",
                                          "supported_reassortment_candidate",
                                          "supported_reassortment_candidate_with_structural_warning"})
        self.assertEqual(result["gate1_discordant_segment"], "")

    def test_artificial_m_swap_is_flagged_discordant_by_gate1(self):
        by_species = {}
        for isolate in self.reference.isolates:
            if isolate.has_complete_triplet:
                by_species.setdefault(isolate.species_name, []).append(isolate)
        species_with_enough = [sp for sp, members in by_species.items() if len(members) >= 5]
        if len(species_with_enough) < 2:
            self.skipTest("need at least 2 species with >=5 complete-triplet isolates")
        species_a, species_b = species_with_enough[0], species_with_enough[1]
        backbone = by_species[species_a][0]
        m_donor = by_species[species_b][0]
        candidate = CandidateInput(
            candidate_id="MSWAP",
            s_genomic=self._genomic(backbone.s_accession),
            m_genomic=self._genomic(m_donor.m_accession),
            l_genomic=self._genomic(backbone.l_accession),
        )
        precomputed = _build_precomputed(self.reference)
        result = score_candidate(candidate, self.reference, precomputed, self._args())
        self.assertEqual(result["gate1_discordant_segment"], "M")
        self.assertEqual(result["gate1_pattern"], "m_only")


class StructuralAnchorDataTest(unittest.TestCase):
    def test_three_species_present_with_expected_basis(self):
        self.assertEqual(STRUCTURAL_ANCHORS["Orthohantavirus andesense"]["basis"], "solved")
        self.assertEqual(STRUCTURAL_ANCHORS["Orthohantavirus puumalaense"]["basis"], "solved")
        self.assertEqual(
            STRUCTURAL_ANCHORS["Orthohantavirus hantanense"]["basis"],
            "homology_inferred_verified",
        )

    def test_andv_tmis_and_ph_switch_at_cited_literal_positions(self):
        seq = STRUCTURAL_ANCHORS["Orthohantavirus andesense"]["sequence"]
        # 1-based positions from Guardado-Calvo et al. (PMC7572791), matching the raw
        # sequence directly -- ANDV is the numerical anchor, see spec section 15.
        self.assertEqual(seq[744], "Y")  # Y745
        self.assertEqual(seq[756], "E")  # E757
        self.assertEqual(seq[758], "D")  # D759
        self.assertEqual(seq[765], "W")  # W766
        self.assertEqual(seq[899], "F")  # F900

    def test_puuv_tmis_positions_are_offset_not_literal(self):
        anchor = STRUCTURAL_ANCHORS["Orthohantavirus puumalaense"]
        seq = anchor["sequence"]
        # PUUV's own QOC68549.1 does NOT match Tscherne et al. 2025's literal integers
        # (confirmed empirically during design, see spec section 15) -- the table stores
        # the corrected, offset positions instead.
        self.assertEqual(seq[anchor["tmis"]["Y745"] - 1], "Y")
        self.assertEqual(seq[anchor["tmis"]["W766"] - 1], "W")
        self.assertEqual(seq[anchor["tmis"]["F900"] - 1], "F")
        self.assertEqual(seq[anchor["ph_switch"]["E757"] - 1], "E")
        self.assertEqual(seq[anchor["ph_switch"]["D759"] - 1], "D")

    def test_htnv_f900_equivalent_is_none_known_exception(self):
        anchor = STRUCTURAL_ANCHORS["Orthohantavirus hantanense"]
        self.assertIsNone(anchor["tmis"]["F900"])
        self.assertIsNotNone(anchor["tmis"]["Y745"])
        self.assertIsNotNone(anchor["tmis"]["W766"])

    def test_waasa_motif_present_at_recorded_position_in_all_three(self):
        for species, anchor in STRUCTURAL_ANCHORS.items():
            pos = anchor["waasa_position"]
            self.assertEqual(anchor["sequence"][pos - 1:pos - 1 + 5], "WAASA", species)

    def test_lineage_marker_table_ships_empty_by_design(self):
        self.assertEqual(LINEAGE_MARKER_TABLE, {})

    def test_new_world_species_membership(self):
        self.assertIn("Orthohantavirus andesense", NEW_WORLD_SPECIES)
        self.assertNotIn("Orthohantavirus hantanense", NEW_WORLD_SPECIES)
        self.assertNotIn("Orthohantavirus puumalaense", NEW_WORLD_SPECIES)


class FindAlignedResidueTest(unittest.TestCase):
    def test_identical_sequence_maps_to_same_position(self):
        seq = STRUCTURAL_ANCHORS["Orthohantavirus andesense"]["sequence"]
        residue, aligned = find_aligned_residue(seq, seq, 745)
        self.assertEqual(residue, "Y")
        self.assertTrue(aligned)

    def test_cross_species_alignment_recovers_the_conserved_residue(self):
        andv = STRUCTURAL_ANCHORS["Orthohantavirus andesense"]["sequence"]
        puuv_anchor = STRUCTURAL_ANCHORS["Orthohantavirus puumalaense"]
        residue, aligned = find_aligned_residue(
            puuv_anchor["sequence"], andv, 745)
        self.assertEqual(residue, "Y")
        self.assertTrue(aligned)

    def test_short_unrelated_sequence_reports_not_aligned_or_wrong_residue(self):
        residue, aligned = find_aligned_residue("MAAAAAAAAAA", "M" * 50, 745)
        self.assertFalse(aligned)
