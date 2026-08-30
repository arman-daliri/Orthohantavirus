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
import copy
import csv
import itertools
import json
import math
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Tuple

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
# Reassortment Plausibility Index (RPI) -- see
# docs/superpowers/specs/2026-08-24-reassortment-plausibility-index-design.md
#
# This is a RELATIVE index, never a probability. Anchored on primary experimental
# measurements, not on any published model's fitted coefficients -- the pipeline has
# only 8 positive labels, which cannot support fitting anything.
# ---------------------------------------------------------------------------

RPI_WEIGHTS: Dict[str, Any] = {
    # Jeffreys prior, Dirichlet(1/2, 1/2, 1/2). NOT arbitrary: Rizvanov's S-odd = 0 and
    # L-odd = 0 are UNOBSERVED events in n = 10, not measured impossibilities. Swept in
    # the sensitivity analysis (spec section 12).
    "jeffreys_alpha": 0.5,

    # Raw in vitro constellation counts -- which segment was the odd one out -- by
    # divergence stratum. Spec section 7.
    #   close   = Rodriguez 1998, SNV NMR11 x CC107, MOI 1.0-1.5, n = 25 (12 M, 12 S, 1 L)
    #   distant = Rizvanov 2004, SNV x ANDV, MOI 3, n = 10 (all M-odd)
    #           + Rodriguez 1998, SNV x BCCV, MOI 1.0-1.5, n = 1 (M-odd)
    #
    # CONFOUNDED AND UNTESTABLE WITH CURRENT LITERATURE: divergence is confounded with
    # pair identity AND with MOI. No distant-pair MOI-1 experiment exists. Leave-one-
    # pair-out (Task 8) and the sensitivity sweep (Task 9) probe this; they cannot
    # resolve it.
    "constellation_counts": {
        "close": {"M": 12, "S": 12, "L": 1},
        "distant": {"M": 11, "S": 0, "L": 0},
    },

    # Pairwise nucleotide identity below which a pair is treated as "distant".
    # PROVISIONAL, literature-anchored. All three anchors below are TIER 1 PRIMARY --
    # the numeric values were read off and verified against the source PDFs during
    # design, and are labelled here so they are not the only unlabelled figures in this
    # module (every other citation in this file carries its verification state):
    #   ~0.85 identity / ~15% divergence -- Rodriguez et al. 1998, Virology 242:99-106,
    #       Table 7 (SNV NMR11 x CC107, the close pair). TIER 1 PRIMARY, value verified
    #       against the source PDF.
    #   0.548 -- Rizvanov et al. 2004, Virology 327:225-232 (SNV x ANDV).
    #       TIER 1 PRIMARY, value verified against the source PDF.
    #   0.707 -- Rizvanov et al. 2004, Virology 327:225-232 (SNV x BCCV).
    #       TIER 1 PRIMARY, value verified against the source PDF.
    # 0.80 separates them with margin on both sides.
    # Spec section 7 calls for replacing this with a mixture-model fit to the project's
    # own PED distribution -- a "DEmARC-style analysis", NOT the DEmARC tool. Swept in
    # the sensitivity analysis.
    "stratum_identity_boundary": 0.80,

    # STRUCTURAL ASSUMPTION -- NOT A MEASUREMENT. These decay weights cannot be
    # estimated from 8 positives. The ladder's ORDERING is supported (Rodriguez 1998
    # measured 8.5% within host species vs 0.6% across host species, both MOI 1); the
    # INTERMEDIATE VALUES and the functional form are assumed. Both the values and the
    # shape of the decay are swept in the sensitivity analysis (spec section 5, 12).
    "host_decay": {
        "same_species": 1.0,
        "same_genus": 0.5,
        "same_family": 0.25,
        "same_order": 0.10,
        "different_order": 0.01,
        # Neither isolate has any host taxonomy at all (host_name, host_family, and
        # host_order all empty). Unknown, not separated -- must not collapse to
        # 'different_order'. See host_rank_proximity's top-of-function guard.
        "unknown_host_metadata": 1.0,
    },
    # Floors keeping every factor strictly positive. See the no-literal-zeros rule.
    "eco_floor": 0.01,
    # Generation -> retention shift, normalised to M = 1.0.
    #   Generation (Rodriguez 1998, in vitro, plaque-level, no amplification): M:S:L = 12:12:1
    #   Retention  (Henderson 1995 via Rodriguez 1998, natural; TIER 3 SECONDHAND): 14:2:1
    # WEAKLY DETERMINED: the natural denominators are 2 S-odd and 1 L-odd. L carries no
    # signal at n = 1 and is set equal to S rather than read off one observation.
    # Swept in the sensitivity analysis.
    "retention_by_odd_segment": {"M": 1.0, "S": 0.15, "L": 0.15},
    # Gate 2 structural integrity (WAASA, TM/IS, pH switch) failing does not veto -- it
    # discounts. Provisional; not fitted.
    "retention_structural_penalty": 0.5,

    # M_pack rescaling (D1, decided WS3.1). `boltzmann_compatibility`'s raw exp(-ddG/RT)
    # form spans a physically-scaled range (RPI_RT_KCAL_PER_MOL ~= 0.6163) incommensurate
    # with the other RPI factors -- see docs/superpowers/specs/2026-08-29-rpi-fitness-
    # function-design.md, "Option A". `m_pack_tau` REPLACES RT at the call site (Task
    # 10 / WS3.3), so M_pack becomes exp(-ddG / m_pack_tau) rather than a physical
    # Boltzmann weight -- the index is explicitly RELATIVE, per compute_rpi's docstring.
    #
    # tau = 3.914, DERIVED (not chosen) from data/rpi_calibration.json: the robust
    # p5-p95 span of DeltaDeltaG over 756 cognate/non-cognate folds (42 isolates, 14
    # species, ViennaRNA Python bindings 2.7.2), divided by ln(100) -- 100 matching
    # host_decay, the widest existing RPI factor. Reproduce with
    # `python rpi_calibration.py`; re-derive across seeds/windows with
    # `python rpi_calibration.py --stability`.
    #
    # THE D1 BRANCH: bimodal-with-a-clean-gap DeltaDeltaG would have selected a
    # logistic (Option C) instead. It measured NOT bimodal -- BIC prefers one
    # component (delta +15.1 over two), component means only 1.24 SD apart against a
    # 2.0 requirement, one density mode. So Option A stands. Stable across seeds
    # (tau = 3.90/4.28/3.90) and per-species sampling caps (3.93-4.39) at this window.
    #
    # KNOWN WEAK, NOT A BUG: DeltaDeltaG separates cognate from non-cognate only
    # weakly (Cliff's delta = +0.222, "small"; 38% of the calibration panel has
    # DeltaDeltaG <= 0, where M_pack saturates at exactly 1.0 and contributes no
    # signal). The genus's terminal 14 nt are invariant across species AND across
    # S/M/L, which structurally caps what a terminal duplex can distinguish --
    # M_pack alone is not expected to separate the WS4 b2 (cross-species sympatric)
    # group; M_pol carries more of that weight. See data/rpi_calibration.json ->
    # "discrimination" for the full measurement. This is evidence to weigh, not a
    # defect in the rescaling to fix by retuning tau.
    "m_pack_tau": 3.914078624686338,
    # The panhandle window (D2) tau was measured under. Swept 25/75 in WS4's
    # sensitivity analysis (Task 15) -- NOT re-derived here, so a caller changing
    # RPI_PANHANDLE_WINDOW_NT without re-running rpi_calibration.py would silently
    # score against a tau measured at a different window than the one folding.
    "m_pack_window_nt": 50,
}


def _weights(weights: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """The weight table to score against: an explicitly-passed one, else the global.

    Exists so the sensitivity sweep can score a perturbed table WITHOUT reassigning the
    module-level RPI_WEIGHTS. That reassignment was process-global: anything scoring
    concurrently -- the WS-API scorer the downstream team calls in a loop -- would have
    silently read perturbed weights and returned numbers attributable to no declared
    weights version. Threading the table through is what makes "never adjust tau or any
    weight inside the caller's loop" an enforced property rather than a documented
    request.
    """
    return RPI_WEIGHTS if weights is None else weights


def smoothed_prior(counts: Dict[str, float], alpha: float) -> Dict[str, float]:
    """Jeffreys-smoothed distribution over "S", "M", "L" from raw odd-segment counts.

    Smoothing is applied AFTER any count arithmetic the caller has done (e.g. the
    leave-one-pair-out subtraction in retrospective_reproduction), so no category can
    reach zero even when the subtraction empties it -- a literal zero inside a
    multiplicative model is a hidden hard gate, which is precisely what spec section 7
    forbids.
    """
    smoothed = {segment: max(0.0, counts.get(segment, 0)) + alpha for segment in ("S", "M", "L")}
    total = sum(smoothed.values())
    return {segment: value / total for segment, value in smoothed.items()}


def constellation_prior(stratum: str, alpha: Optional[float] = None,
                        weights: Optional[Dict[str, Any]] = None) -> Dict[str, float]:
    """Jeffreys-smoothed prior over which segment is the odd one out.

    Returns a distribution over "S", "M", "L" summing to 1.0, with every category
    strictly positive. Raises KeyError for an unknown stratum rather than silently
    defaulting -- a wrong stratum must be loud.
    """
    table = _weights(weights)
    alpha = table["jeffreys_alpha"] if alpha is None else alpha
    return smoothed_prior(table["constellation_counts"][stratum], alpha)


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
    # Segment -> the REAL isolate that contributed that segment, when a caller already
    # knows it (e.g. rpi_discrimination.py's group a/b1/b2 builders, which pick specific
    # isolates on purpose). None (the default) for ordinary production candidates, whose
    # true parent isolate is exactly what reassortment detection exists to resolve.
    # Consumed only by resolve_e_parent_isolates -- see rpi_e_isolate_basis.
    known_segment_isolates: Optional[Dict[str, "ah.IsolateRecord"]] = None


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def _translate(seq: str) -> str:
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        protein.append(_CODON_TABLE.get(codon, "X"))
    return "".join(protein)


def _longest_orf_span_with_terminus(genomic_seq: str) -> Tuple[str, bool]:
    """(longest ATG-initiated nucleotide span, whether a stop codon terminated it).

    Searches all 6 reading frames. A span that reaches the end of the record WITHOUT
    hitting an in-frame stop is a valid candidate, and is returned truncated to whole
    codons with stop_terminated=False.

    WHY OPEN-ENDED SPANS COUNT. Requiring a terminal stop silently discarded real
    full-length ORFs, because GenBank records are routinely trimmed to (or just inside)
    the CDS. Measured case: Tatenale virus Upton Heath M, NC_055637.1 -- a RefSeq
    record this project's own upstream annotates as coding-complete with an 1118 aa
    glycoprotein. The record is 3355 nt; that CDS needs 3360 nt with start and stop, so
    the deposit stops 5 nt short of its own stop codon. The real ORF was therefore
    thrown away and a spurious 62 aa ORF elsewhere in the record won instead --
    whereupon gate1 could match nothing and reported the segment "too_divergent" on an
    ordinary, unmodified isolate. Translating the first codons of that record gives
    MEQRSTVCLCLFLTGLLLWNPSESRN, an unmistakable hantavirus GPC signal peptide.

    This does not soften junk detection. An open-ended span is only long when a frame
    genuinely runs stop-free, which is the signature of coding sequence: in random
    sequence a stop appears every ~21 codons, so the probability of 400 stop-free
    codons is about 4e-9. Stop-terminated spans still win ties, so nothing that has a
    real terminus is reinterpreted.
    """
    genomic_seq = genomic_seq.upper()
    best, best_stop_terminated = "", False
    for strand_seq in (genomic_seq, _reverse_complement(genomic_seq)):
        for frame in range(3):
            i = frame
            while i < len(strand_seq) - 2:
                if strand_seq[i:i + 3] == "ATG":
                    j = i + 3
                    stop_found = False
                    while j < len(strand_seq) - 2:
                        if strand_seq[j:j + 3] in _STOP_CODONS:
                            stop_found = True
                            break
                        j += 3
                    if stop_found:
                        span = strand_seq[i:j + 3]
                    else:
                        # Run to the last COMPLETE codon; a trailing partial codon is
                        # not translatable and must not inflate the comparison.
                        span = strand_seq[i:i + 3 * ((len(strand_seq) - i) // 3)]
                    if len(span) > len(best) or (
                            len(span) == len(best) and stop_found and not best_stop_terminated):
                        best, best_stop_terminated = span, stop_found
                i += 3
    return best, best_stop_terminated


def _longest_orf_span(genomic_seq: str) -> str:
    """Longest ATG-initiated nucleotide span (stop codon included when there is one)
    across all 6 reading frames. Shared by longest_orf_protein (Gate 0, translates it)
    and candidate_cds_from_orf (Gate 4, uses the nucleotide span directly) so the
    ORF-search loop exists exactly once. Returns '' when no ATG-initiated span exists.
    """
    return _longest_orf_span_with_terminus(genomic_seq)[0]


def longest_orf_protein(genomic_seq: str) -> str:
    """Longest ATG-initiated ORF across all 6 reading frames, translated to protein.

    The initiator Met and the stop codon are both excluded from the returned protein.
    An ORF running off the end of a CDS-trimmed record is translated to its last
    complete codon -- see _longest_orf_span_with_terminus for why that case is real and
    what it cost when it was not handled. Returns '' when no ATG-initiated span exists;
    Gate 0 treats that as a structural failure, not a crash.
    """
    span, stop_terminated = _longest_orf_span_with_terminus(genomic_seq)
    if not span:
        return ""
    return _translate(span[3:-3] if stop_terminated else span[3:])


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
        # Same tolerance principle as the protein-length check below, for the same
        # reason, and it was missing here: SEGMENT_LENGTH_RANGE is ICTV's band for
        # DECIDING WHICH SEGMENT a record is (all_hanta.py:86), not a validity
        # criterion, so enforcing it exactly turned "unusual length" into "junk".
        #
        # Measured against this project's own reference set, that rejected real NCBI
        # RefSeq complete genomes at the margins: Luxi virus LX309 S (NC_038530.1) is
        # 2107 nt against a 1200-2100 band -- junk_sequence by SEVEN nucleotides -- and
        # Tatenale virus Upton Heath M (NC_055637.1) is 3355 nt against 3400-3800.
        # junk_sequence is terminal (score_candidate returns immediately), so a
        # divergent novel hantavirus was discarded unscored for being marginally off a
        # band derived from better-sampled relatives. That inverts the tool's purpose.
        #
        # The band is still enforced, just with room at the edges. Genuinely truncated
        # records are unaffected: a 1014 nt M is 70% below the band and fails either
        # way. all_hanta.py's own SEGMENT_LENGTH_ABSURD_MAX stays a hard ceiling --
        # beyond it a record "is almost certainly a multi-segment concatenation or a
        # mis-label" (all_hanta.py:90), which is junk in the sense this gate means.
        tolerance = max(1, int(0.15 * (high - low + 1)))
        absurd_max = ah.SEGMENT_LENGTH_ABSURD_MAX.get(segment, high + tolerance)
        if len(genomic) > absurd_max:
            failures.append(
                f"{segment}: genomic length {len(genomic)} exceeds the plausible maximum "
                f"{absurd_max} -- likely a concatenation or mis-labelled record")
        elif not (low - tolerance <= len(genomic) <= high + tolerance):
            failures.append(
                f"{segment}: genomic length {len(genomic)} outside expected {low}-{high} "
                f"(+/- {tolerance} nt tolerance)")
        elif not (low <= len(genomic) <= high):
            warnings.append(
                f"{segment}: genomic length {len(genomic)} outside the ICTV band "
                f"{low}-{high} but within tolerance -- unusual, not junk")

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

    # sorted(), not list(): this was list({...}) over a SET, and Python randomises string
    # hashing per process, so the order varied between runs of identical input. That was
    # not cosmetic -- gate3_check reads implicated_species[:2], so WHICH two species gate
    # 3 examined was chosen by the hash seed. Measured on three real species, 2 of the 6
    # orderings returned "ecological_overlap_plausible" where the rest returned
    # "no_ecological_route", and combine_verdict returns that status directly.
    #
    # Alphabetical order carries no meaning and is not claimed to; it is simply stable.
    # See the KNOWN LIMITATION on gate3_check for the 2-of-N question this does not fix.
    implicated_species = sorted({s for s in best_species.values() if s})
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

    # KNOWN LIMITATION -- gate 3 judges only the FIRST TWO of N implicated species.
    #
    # With three implicated species the other two pairings are never examined, and they
    # can disagree. Measured on Orthohantavirus hantanense / sinnombreense / luxiense:
    #     hantanense    + sinnombreense -> no_ecological_route
    #     hantanense    + luxiense      -> ecological_overlap_plausible
    #     sinnombreense + luxiense      -> no_ecological_route
    # so the reported status depends on which pair happens to sort first. Until
    # 2026-08-29 that choice was made by the hash seed and varied run to run; gate1_check
    # now sorts, which makes the choice STABLE but does not make it PRINCIPLED.
    #
    # DELIBERATELY NOT FIXED HERE. The right answer is the worst-pair-over-all-pairs
    # treatment that ecological_opportunity already implements for Factor E, and moving
    # gate 3 onto it changes a verdict-bearing gate -- out of scope for the RPI work.
    # Deferred to the Factor E rebuild. Until then: a "no_ecological_route" verdict on a
    # 3-species candidate means "the alphabetically first two do not overlap", not "no
    # pair overlaps". Check gate1_implicated_species before quoting such a row.
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


def divergence_stratum(identity: float, boundary: Optional[float] = None,
                       weights: Optional[Dict[str, Any]] = None) -> str:
    """Classify a parental pair as 'close' or 'distant' by nucleotide identity.

    `identity` is a fraction in [0, 1], NOT a percentage. The stratum selects which
    constellation prior applies -- see constellation_prior.
    """
    if not 0.0 <= identity <= 1.0:
        raise ValueError(f"identity must be a fraction in [0, 1], got {identity!r}")
    boundary = (_weights(weights)["stratum_identity_boundary"]
                if boundary is None else boundary)
    return "close" if identity >= boundary else "distant"


def host_genus_of(isolate: "ah.IsolateRecord") -> str:
    """Genus from the host binomial. IsolateRecord carries no host_genus field; host_name
    is a binomial in the curated data ('Myodes glareolus'), so the first token is the
    genus. Returns '' when host_name is absent, which callers treat as unknown rather
    than as a mismatch.
    """
    return isolate.host_name.split()[0] if isolate.host_name else ""


def host_rank_proximity(a: "ah.IsolateRecord", b: "ah.IsolateRecord") -> str:
    """Finest shared host taxonomic rank between two isolates' reservoirs.

    Missing metadata never yields 'different_order' by default: an empty field means
    unknown, and unknown must not masquerade as ecological separation.
    """
    def has_any_host_field(isolate: "ah.IsolateRecord") -> bool:
        return bool(isolate.host_name or isolate.host_family or isolate.host_order)

    if not has_any_host_field(a) or not has_any_host_field(b):
        # Either side entirely lacking host taxonomy (e.g. resolve_host_taxonomy
        # failed and no biosample host was available) is unknown, not separated --
        # mirrors ecological_opportunity's own geo term, which stays neutral when
        # either side's country set is empty rather than assuming the worst.
        return "unknown_host_metadata"
    if a.host_name and a.host_name == b.host_name:
        return "same_species"
    genus_a, genus_b = host_genus_of(a), host_genus_of(b)
    if genus_a and genus_a == genus_b:
        return "same_genus"
    if a.host_family and a.host_family == b.host_family:
        return "same_family"
    if a.host_order and a.host_order == b.host_order:
        return "same_order"
    return "different_order"


def _year_span(isolates: Sequence["ah.IsolateRecord"]) -> Optional[Tuple[int, int]]:
    years = [int(i.year) for i in isolates if str(i.year).isdigit()]
    return (min(years), max(years)) if years else None


def ecological_opportunity(
    implicated_species: List[str], isolates: List["ah.IsolateRecord"],
    parent_isolates: Dict[str, "ah.IsolateRecord"], isolate_basis: str,
    weights: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Factor E: relative opportunity for the implicated lineages to co-infect one cell.

    E = h_host * g_geo * t_time, over the WORST pair among all implicated species --
    deliberately unlike gate3_check, which reads implicated_species[:2] and has no
    temporal term. Disagreement between the two is expected and informative.

    FIX (2026-08-30, docs/superpowers/specs/2026-08-30-rpi-e-geo-isolate-level-fix-
    design.md): host/geo/time are computed from the ONE isolate that stands for each
    implicated species IN THIS CANDIDATE (`parent_isolates`), never from that species'
    full country/host/year union across the whole reference. The old species-wide union
    let widely-distributed or border-spanning species (Puumala, Seoul) register "geo
    opportunity" against every other species that shares ANY of their many countries,
    even when the two specific parent isolates being scored were sampled continents
    apart -- see the design doc and notes_diagnosis.md for the confirmed mechanism.

    `parent_isolates` and `isolate_basis` are REQUIRED, not defaulted (same reasoning
    `m_pack_status`/`m_pol_status` already use): a caller must state which single
    isolate stands for each species and where it came from
    (see resolve_e_parent_isolates) rather than let this function silently guess one.
    `isolate_basis` is echoed back verbatim as `rpi_e_isolate_basis` -- every row must
    declare which basis its E actually used, never leave it implicit.

    Returns 1.0 (neutral) when fewer than two species are implicated: there is no
    co-infection requirement to test, so E must not penalise.
    """
    table = _weights(weights)
    floor = table["eco_floor"]
    if len(implicated_species) < 2:
        return {"rpi_e_host": 1.0, "rpi_e_geo": 1.0, "rpi_e_time": 1.0, "rpi_e_total": 1.0,
                "rpi_e_isolate_basis": isolate_basis}

    all_species_present = {isolate.species_name for isolate in isolates}
    wild_species_present = {
        isolate.species_name for isolate in isolates
        if not isolate.is_lab_derived and not isolate.is_vaccine_strain}

    worst = {"rpi_e_host": 1.0, "rpi_e_geo": 1.0, "rpi_e_time": 1.0}
    for index, species_a in enumerate(implicated_species):
        for species_b in implicated_species[index + 1:]:
            isolate_a = parent_isolates.get(species_a)
            isolate_b = parent_isolates.get(species_b)
            if isolate_a is None or isolate_b is None:
                # Same two absence cases the old species-wide union distinguished,
                # checked against the UNFILTERED isolates list so they aren't conflated
                # -- now evaluated per species that failed to resolve, since a resolved
                # isolate on one side and an absent species on the other is still
                # informative:
                #  (a) the species IS present in the reference, but every one of its
                #      isolates is lab- or vaccine-derived (so resolution correctly
                #      found none): there is genuinely no wild population, so no
                #      natural ecological route -- the host term takes the floor rather
                #      than defaulting to full opportunity.
                #  (b) the species does not appear in the reference AT ALL (novel, or
                #      simply never sampled): that is unknown, not separated. Scoring
                #      it as impossible would contradict gate3_check's own
                #      "ecology_unknown" treatment of this exact situation and would
                #      punish a candidate for lack of data rather than lack of
                #      plausibility -- stay neutral and leave `worst` untouched.
                for species, isolate in ((species_a, isolate_a), (species_b, isolate_b)):
                    if isolate is not None:
                        continue
                    if species not in all_species_present:
                        continue
                    if species not in wild_species_present:
                        worst["rpi_e_host"] = min(worst["rpi_e_host"], floor)
                continue

            host = table["host_decay"][host_rank_proximity(isolate_a, isolate_b)]
            country_a, country_b = isolate_a.country, isolate_b.country
            # Missing-means-neutral takes precedence: an empty country on either side
            # is unknown, not separated, so it must be checked before the floor branch
            # could otherwise apply.
            if not country_a or not country_b:
                geo = 1.0
            elif country_a == country_b:
                geo = 1.0
            else:
                geo = floor

            span_a, span_b = _year_span([isolate_a]), _year_span([isolate_b])
            if span_a and span_b:
                # +1: inclusive-year spans that merely touch (e.g. 2001 vs 2002, a
                # one-year gap between point samples) count as overlapping -- point
                # samples a year apart are effectively contemporaneous for a
                # multi-year reservoir infection cycle, not evidence of separation.
                overlap = min(span_a[1], span_b[1]) - max(span_a[0], span_b[0]) + 1
                time = 1.0 if overlap >= 0 else max(floor, 1.0 / (1.0 + abs(overlap)))
            else:
                time = 1.0

            worst["rpi_e_host"] = min(worst["rpi_e_host"], host)
            worst["rpi_e_geo"] = min(worst["rpi_e_geo"], geo)
            worst["rpi_e_time"] = min(worst["rpi_e_time"], time)

    worst["rpi_e_total"] = worst["rpi_e_host"] * worst["rpi_e_geo"] * worst["rpi_e_time"]
    worst["rpi_e_isolate_basis"] = isolate_basis
    return worst


def nearest_reference_isolate(
    candidate_protein: str, segment: str, isolates: List["ah.IsolateRecord"],
) -> Optional["ah.IsolateRecord"]:
    """The isolate (already filtered to one species) whose own `segment` protein is
    nearest `candidate_protein` by k-mer cosine distance -- the nearest-reference-isolate
    proxy `resolve_e_parent_isolates` falls back to for Factor E when a candidate's true
    parent isolate is unknown (design doc section 1a). Mirrors `nearest_species_by_kmer`'s
    method, at isolate instead of species granularity.

    Returns None if `candidate_protein` is empty or no candidate isolate has a usable
    protein for `segment` -- the caller must treat that as unresolved (`insufficient_
    metadata`), never guess an isolate.
    """
    attr = {"S": "s_protein_aa", "M": "m_protein_aa", "L": "l_protein_aa"}[segment]
    candidates = [(isolate, getattr(isolate, attr)) for isolate in isolates]
    candidates = [(isolate, protein) for isolate, protein in candidates if protein]
    if not candidate_protein or not candidates:
        return None
    profiles = [ah.kmer_profile(candidate_protein)] + [
        ah.kmer_profile(protein) for _isolate, protein in candidates]
    matrix, _vocab = ah.profiles_to_matrix(profiles)
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    unit = matrix / norms
    similarity = np.clip(unit[0] @ unit[1:].T, -1.0, 1.0)
    distances = 1.0 - similarity
    best_index = int(np.argmin(distances))
    return candidates[best_index][0]


def resolve_e_parent_isolates(
    implicated_species: List[str], isolates: List["ah.IsolateRecord"],
    gate1: "Gate1Result", orf_proteins: Dict[str, str],
    known_isolates: Optional[Dict[str, "ah.IsolateRecord"]] = None,
) -> Tuple[Dict[str, "ah.IsolateRecord"], str]:
    """(species -> the one isolate standing for it, basis) for `ecological_opportunity`
    (design doc section 1a).

    `known_isolates` (segment -> the REAL isolate that contributed that segment) is
    supplied only by a caller that already knows the truth -- e.g. rpi_discrimination.py's
    group a/b1/b2 builders, which construct candidates from specific isolates on
    purpose. When every implicated species is covered by it, basis is "actual_parents".

    Otherwise (ordinary production candidates, and WS4 group c's real positives, whose
    true OTHER parent is exactly what reassortment detection exists to resolve) this
    falls back to the k-mer-nearest wild reference isolate of each implicated species,
    on the segment gate1 top-matched that species against -- the same "nearest
    well-defined ground truth, since the real one is unknown" reasoning
    species_consensus_panhandle_window (D3) already uses for M_pack/M_pol, applied per
    isolate instead of per species-consensus. basis is "nearest_reference" when every
    species resolves this way (or a mix of actual+fallback), "insufficient_metadata"
    when at least one species resolves neither way.
    """
    if known_isolates:
        by_species = {isolate.species_name: isolate for isolate in known_isolates.values()}
        if all(species in by_species for species in implicated_species):
            return ({species: by_species[species] for species in implicated_species},
                    "actual_parents")

    resolved: Dict[str, "ah.IsolateRecord"] = {}
    basis = "nearest_reference"
    for species in implicated_species:
        segment = next(
            (seg for seg in ("S", "M", "L")
             if gate1.top_hits.get(seg) and gate1.top_hits[seg][0].species == species), "")
        pool = _wild_isolates_of_species(isolates, species)
        found = (nearest_reference_isolate(orf_proteins.get(segment, ""), segment, pool)
                 if segment else None)
        if found is None:
            basis = "insufficient_metadata"
            continue
        resolved[species] = found
    return resolved, basis


# Gas constant x RNAcofold's default temperature (37 C = 310.15 K), in kcal/mol so the
# units match RNAcofold's reported dG. Both the temperature and the energy-parameter set
# are pinned with the ViennaRNA version in environment.yml -- dG is parameter-set
# dependent. Spec section 6.1.
RPI_RT_KCAL_PER_MOL = 0.0019872 * 310.15  # ~= 0.6163


def boltzmann_compatibility(
    dg_non_cognate: float, dg_cognate: float, rt: Optional[float] = None,
) -> float:
    """Relative Boltzmann weight of the non-cognate duplex against the cognate one.

    ddG = dG_non_cognate - dG_cognate. A WEAKER non-cognate duplex has a LESS NEGATIVE
    dG, so ddG > 0 and the result is < 1. REVERSING THIS SUBTRACTION INVERTS THE ENTIRE
    INDEX -- see test_sign_convention_is_not_flipped.

    Clamped to 1.0 above: this factor is a penalty for non-cognate competition, never a
    bonus. Strictly positive below: a zero would be a hidden hard gate.

    Methods sentence: the exp(-ddG/RT) form is the relative Boltzmann weight of the
    non-cognate versus cognate duplex -- the equilibrium expression of the competition
    argument of McDonald et al. (2016), in which low-affinity non-cognate RNA-RNA
    interactions are outcompeted by higher-affinity cognate ones.
    """
    rt = RPI_RT_KCAL_PER_MOL if rt is None else rt
    ddg = dg_non_cognate - dg_cognate
    # A non-finite ddG means RNAcofold returned something unusable. It must NOT be
    # scored, and the two silent answers available here are opposites -- which is what
    # makes the silence dangerous rather than merely imprecise. max(-700.0, nan) returns
    # -700.0, so nan scored as MAXIMAL penalty (9.86e-305); before the underflow clamp,
    # min(0.0, nan) returned 0.0 and the same input scored as NO penalty (1.0). Raising
    # is the only honest option. The call site in score_candidate catches this and
    # records rpi_m_pack_status = "undefined_energy" rather than aborting the run.
    if not math.isfinite(ddg):
        raise ValueError(
            f"non-finite ddG ({ddg!r}) from dg_non_cognate={dg_non_cognate!r}, "
            f"dg_cognate={dg_cognate!r} -- RNAcofold returned an unusable energy"
        )
    # Clamp the EXPONENT, not just the result: math.exp overflows above ~709, so a
    # far-stronger non-cognate duplex (ddg very negative, exponent very positive) would
    # raise OverflowError before an outer min(1.0, ...) ever got a chance to clamp it.
    # min(0.0, ...) here means math.exp's argument is never positive, so its result is
    # never above 1.0 -- the outer min(1.0, ...) below is now redundant but kept as a
    # cheap, harmless belt-and-suspenders guard against floating-point edge cases.
    #
    # The exponent is clamped BELOW as well, and that half is not symmetry for its own
    # sake. math.exp does not raise on underflow the way it raises on overflow -- it
    # returns 0.0 silently, so a large POSITIVE ddg produced exactly the hidden hard
    # gate the docstring above promises never to produce: boltzmann_compatibility(500.0,
    # 0.0) was 0.0.
    # -700 keeps exp() in normal range (its smallest normal result is ~1e-308) while
    # still being far below any physically meaningful ddg.
    #
    # SCOPE NOTE (D2, 2026-08-29): this comment previously justified the clamp by saying
    # cofolding took WHOLE-SEGMENT duplexes whose MFEs run to hundreds of kcal/mol.
    # That is no longer true and was never runnable -- cofolding is O(n^3), so an L x M
    # duplex (~10 kb) will not complete. The scored window is the 50-nt panhandle
    # (RPI_PANHANDLE_WINDOW_NT), where
    # observed ddg runs ~0-22 kcal/mol and the clamp is never reached in practice. The
    # clamp is retained anyway: it is cheap, and it is what keeps the no-literal-zeros
    # guarantee true for adversarial input rather than merely true for expected input.
    exponent = min(0.0, max(-700.0, -ddg / rt))
    return min(1.0, math.exp(exponent))


# Panhandle cofold window, in nt, per design doc D2. Swept 25 / 50 / 75 in the
# sensitivity analysis. 50 comfortably contains the invariant 14-nt tag plus the
# lineage-variable region beyond it (RPI_N_FOOTPRINT_NT runs to 23), folds in ~7 ms, and
# keeps the fold tractable -- whole-segment duplexes are O(n^3) and do not complete.
RPI_PANHANDLE_WINDOW_NT = 50


def panhandle_window(genomic_seq: str, end: str, n: int = RPI_PANHANDLE_WINDOW_NT) -> str:
    """The terminal `n` nt of one end of a segment, for panhandle cofolding (D2).

    `end` is "5" or "3". The hantavirus panhandle is formed by the complementary 5' and
    3' termini, so a duplex is scored as one segment's 5' window against another's 3'
    window.

    WHY A WINDOW AND NOT THE WHOLE SEGMENT: RNAcofold is O(n^3) in time and O(n^2) in
    memory. S is ~1.8 kb, M ~3.6 kb and L ~6.5 kb, so an L x M duplex is a ~10 kb fold --
    on the order of 1e12 operations. It does not complete. It is also the wrong quantity:
    packaging compatibility is a terminal interaction, not a distributed one.
    """
    if end not in ("5", "3"):
        raise ValueError(f"end must be '5' or '3', got {end!r}")
    sequence = genomic_seq.strip().upper()
    return sequence[:n] if end == "5" else sequence[-n:]


def cofold_dg(seq_a: str, seq_b: str) -> float:
    """Minimum free energy of the seq_a/seq_b duplex, in kcal/mol, via ViennaRNA.

    USES THE VIENNARNA PYTHON BINDINGS, NOT A SUBPROCESS (design doc D2b). This replaced
    an `RNAcofold` CLI call, and the replacement is deliberate on three counts:

      1. The pip `viennarna` package ships the BINDINGS ONLY. There is no RNAcofold
         executable in the interpreter tree and none on PATH, so the subprocess form
         could not run here at all -- it raised FileNotFoundError.
      2. One pinned version governs the energy-parameter set. dG is parameter-set
         dependent, so a separately-installed CLI could drift from the bindings and
         silently change every downstream score. There is exactly one pin to check:
         viennarna==2.7.2, asserted by EnvironmentReconciliationTest.
      3. No process spawn. The WS-API scorer calls this in a loop over thousands of
         triplets; a spawn per fold would dominate the measured ~7 ms fold cost.

    THERE IS DELIBERATELY NO CLI FALLBACK. A half-working alternative path could diverge
    from this one without announcing it, which is the failure the single-mechanism rule
    exists to prevent.

    Raises rather than returning a sentinel on unusable input: a silently-defaulted dG
    would propagate a fabricated energy into every downstream score.
    """
    import RNA
    for name, seq in (("seq_a", seq_a), ("seq_b", seq_b)):
        if not seq:
            raise ValueError(f"{name} is empty -- cannot cofold")
    _structure, energy = RNA.fold_compound(f"{seq_a}&{seq_b}").mfe_dimer()
    energy = float(energy)
    if not math.isfinite(energy):
        raise RuntimeError(
            f"ViennaRNA returned a non-finite dG ({energy!r}) for a "
            f"{len(seq_a)}nt x {len(seq_b)}nt duplex"
        )
    return energy


def viennarna_provenance() -> Dict[str, str]:
    """Mechanism and version, for the calibration artifact to record (D2b).

    The entire tau derivation depends on WHICH implementation produced the energies, so
    the artifact records the mechanism alongside the version rather than the version
    alone -- a bare "2.7.2" would not distinguish these energies from CLI-produced ones.
    """
    import RNA
    return {"mechanism": "ViennaRNA Python bindings", "version": str(RNA.__version__)}


# The genus-conserved GENOMIC 5'/3' termini that form the intramolecular panhandle
# (D2, WS3.1). NOT the same thing as HANTAVIRUS_TERMINAL_TAG_NT below -- that is a
# CDS-relative invariant LENGTH used by polymerase_compatibility's N-recognition
# scoring; these are the literal genomic RNA sequence at each end. A deposit missing
# them is missing its actual termini: GenBank submissions are routinely trimmed, and
# `SEGMENT_LENGTH_RANGE` does not catch it -- a segment 30 nt short at each end is still
# comfortably inside the ICTV band. Measured on a real deposit (rpi_calibration.py, WS3.1):
# a truncated segment's "terminal window" folded to -10.0 kcal/mol against -42.2 for a
# terminally complete one -- arbitrary internal sequence, not a panhandle.
HANTAVIRUS_GENOMIC_5P_TAG = "TAGTAGTA"
HANTAVIRUS_GENOMIC_3P_TAG = "TACTACTA"  # == _reverse_complement(HANTAVIRUS_GENOMIC_5P_TAG)


def has_genomic_panhandle_termini(sequence: str) -> bool:
    """Does this genomic deposit carry both true termini, or is it trimmed?

    48% of the wild full-length triplet pool in this dataset carries NEITHER tag
    (measured in rpi_calibration.py, WS3.1). Folding a trimmed deposit's "terminal
    window" scores noise, not a panhandle -- this predicate is what
    `species_consensus_panhandle_window` filters on before building a consensus, for
    the same reason the WS3.1 calibration panel filters on it before deriving tau: a
    consensus built from mostly-trimmed sequences would be a consensus of arbitrary
    internal sequence at the columns nearest each end.
    """
    sequence = sequence.strip().upper()
    return (sequence.startswith(HANTAVIRUS_GENOMIC_5P_TAG)
            and sequence.endswith(HANTAVIRUS_GENOMIC_3P_TAG))


# The terminal 14 nt of S, M and L are highly conserved across the genus AND shared
# between the three segments (J Virol 1995;69(12):8132-8136, PMC189768 -- content
# verified, author list NOT verified; resolve before any bibliography). M_pol must
# therefore NEVER be scored on this window: it would saturate at ~1.0 for every pair.
HANTAVIRUS_TERMINAL_TAG_NT = 14

# N-recognition footprint lengths. Mir MA & Panganiban AT, "Hantavirus N Protein Exhibits
# Genus-Specific Recognition of the Viral RNA Panhandle", J Virol 2006,
# doi:10.1128/JVI.00820-06 -- existence and abstract verified, FULL TEXT NOT READ.
# SNV N requires the first 23 nt; most hantavirus N proteins need only the terminal 9.
# This asymmetry VARIES BETWEEN LINEAGES, which is why it carries signal where the
# invariant tag does not.
RPI_N_FOOTPRINT_NT = {"Orthohantavirus sinnombreense": 23}
RPI_N_FOOTPRINT_DEFAULT_NT = 9

# Below this variance across sampled pairs, M_pol is declared uninformative and fixed.
RPI_M_POL_MIN_VARIANCE = 1e-4

# CDS-relative window length M_pol compares, on both the candidate and the reference
# consensus side (D3, WS3.4). Comfortably contains HANTAVIRUS_TERMINAL_TAG_NT (14) plus
# every entry in RPI_N_FOOTPRINT_NT (23 at most), so _beyond_tag's slice never runs off
# the end of either window. A NAMED constant rather than an inline 64 at both call
# sites -- score_candidate's candidate-side window and species_consensus_cds_window's
# reference-side window MUST agree on this number, or the two windows being compared
# would not even be the same length.
RPI_MPOL_CDS_WINDOW_NT = 64


def n_recognition_footprint(species_name: str) -> int:
    """Length of the N-recognition window for a species, in nucleotides."""
    return RPI_N_FOOTPRINT_NT.get(species_name, RPI_N_FOOTPRINT_DEFAULT_NT)


def _beyond_tag(sequence: str, footprint: int) -> str:
    """The scoring window: the part of the N-recognition footprint lying BEYOND the tag.

    Mir & Panganiban's footprints are measured FROM THE TERMINUS -- "the first 23 nt"
    for SNV N, "the terminal 9" for most other hantavirus N proteins. So the footprint
    spans sequence[0:footprint], and the part of it that is not the genus-invariant
    14-nt tag is sequence[14:footprint].

    An earlier version returned sequence[14:14+footprint], which for the 9-nt default
    scored nt 15-23 -- entirely OUTSIDE a footprint that ends at nt 9. That window had
    no source behind it.

    For the 9-nt default this slice is EMPTY, and correctly so: the whole default
    footprint lies inside the invariant tag, so there is nothing left to measure. Only
    the SNV 23-vs-9 asymmetry -- the part that actually varies between lineages -- can
    produce a scoring window at all. See polymerase_compatibility's empty-window branch.
    """
    return sequence[HANTAVIRUS_TERMINAL_TAG_NT:footprint]


def polymerase_compatibility(
    candidate_termini: Dict[str, str], parent_termini: Dict[str, str],
    footprint: int = RPI_N_FOOTPRINT_DEFAULT_NT,
) -> Tuple[float, str]:
    """Factor M_pol: protein-RNA recognition compatibility, scored BEYOND the conserved tag.

    Returns (value, status). Status is one of:
      "scored"               -- a non-empty window beyond the tag was compared;
      "no_signal_beyond_tag" -- the footprint lies entirely inside the invariant 14-nt
                                tag (the 9-nt default) or no sequence was supplied, so
                                there is nothing to measure. Value is 1.0, neutral.
    m_pol_variance_check may later downgrade the factor to "fixed_constant" across a run.

    Reporting "no_signal_beyond_tag" rather than a scored 1.0 is the honest half of the
    _beyond_tag fix: a neutral value that came from an EMPTY comparison must not be
    reported in the same status as a neutral value that came from a perfect match.

    FRAMING REQUIREMENT (spec 6.2): hantavirus L-polymerase terminal recognition is a
    STRUCTURALLY MOTIVATED HYPOTHESIS, by analogy with the rotavirus 3'-terminal
    seven-nucleotide polymerase recognition described by McDonald et al. (2016). It is
    NOT an experimentally demonstrated hantavirus mechanism. Methods text must say so.
    """
    matches = total = 0
    for segment in ("S", "M", "L"):
        window_a = _beyond_tag(candidate_termini.get(segment, ""), footprint)
        window_b = _beyond_tag(parent_termini.get(segment, ""), footprint)
        for base_a, base_b in zip(window_a, window_b):
            total += 1
            matches += int(base_a == base_b)
    if total == 0:
        return 1.0, "no_signal_beyond_tag"
    identity = matches / total
    return max(RPI_WEIGHTS["eco_floor"], identity), "scored"


def m_pol_variance_check(
    observed: Sequence[Tuple[float, Any]], min_variance: float = RPI_M_POL_MIN_VARIANCE,
) -> Tuple[float, str]:
    """Mandatory gate from spec 6.2.

    If M_pol shows negligible variance across a representative sample of pairs, it is
    measuring nothing and must be FIXED AT 1.0 AND DECLARED SO. A constant factor that
    announces itself is honest; a constant dressed as a measurement is not.
    """
    values = [value for value, _ in observed]
    if not values or float(np.var(values)) < min_variance:
        return 1.0, "fixed_constant"
    return float(np.mean(values)), "scored"


def _wild_isolates_of_species(
        isolates: List["ah.IsolateRecord"], species: str) -> List["ah.IsolateRecord"]:
    """Real (not lab-derived, not vaccine-strain) isolates of one species (D3).

    Shared by every "implicated parent/donor consensus" builder -- M_pack's genomic
    termini and M_pol's CDS windows both start from this same wild-isolate pool, they
    just do different things to each sequence afterward.
    """
    return [i for i in isolates if i.species_name == species
            and not i.is_lab_derived and not i.is_vaccine_strain]


def _wild_species_genomic_sequences(
        reference: ReferenceData, species: str, segment: str) -> List[str]:
    """Terminally complete wild genomic sequences for one species' segment (D3, M_pack).

    Filters on `has_genomic_panhandle_termini` for the same reason WS3.1's calibration
    panel does: a trimmed deposit's terminal window is arbitrary internal sequence, and
    folding it would corrupt the consensus this feeds into M_pack.
    """
    isolates = _wild_isolates_of_species(reference.isolates, species)
    sequences = genomic_sequences_by_segment(isolates, reference.genomic_index)[segment]
    return [seq for seq in sequences if has_genomic_panhandle_termini(seq)]


def _wild_species_cds_windows(
        reference: ReferenceData, species: str, segment: str, window: int) -> List[str]:
    """CDS-derived terminal windows for one species' segment (D3, M_pol) -- the ORF's
    own nucleotide span (`candidate_cds_from_orf`), same as `terminal_windows` is built
    from on the candidate side in score_candidate.

    NOT filtered on `has_genomic_panhandle_termini`: that predicate is about the
    GENOMIC 5'/3' UTR termini M_pack folds, which is a different region from the CDS/ORF
    this scores. A deposit trimmed of its UTR can still carry a fully intact ORF, so
    applying that filter here would exclude good CDS candidates for the wrong reason.
    Isolates whose segment has no findable ORF at all are dropped (`candidate_cds_from_orf`
    returns "" in that case).
    """
    isolates = _wild_isolates_of_species(reference.isolates, species)
    sequences = genomic_sequences_by_segment(isolates, reference.genomic_index)[segment]
    windows = [candidate_cds_from_orf(seq)[:window] for seq in sequences]
    return [w for w in windows if w]


def species_consensus_panhandle_window(
        reference: ReferenceData, species: str, segment: str, end: str,
        n: int = RPI_PANHANDLE_WINDOW_NT) -> str:
    """Majority-vote consensus terminal window for one species' segment (D3) -- the
    "implicated parent/donor consensus" M_pack is scored against.

    A SPECIES-level substitute for the isolate-level cognate pairing
    rpi_calibration.py measures tau from: at scoring time the candidate's true parent
    isolate is exactly the unknown the whole exercise is trying to resolve, so the
    nearest well-defined ground truth is the species' own consensus, not one isolate's
    sequence picked arbitrarily.

    Returns "" if the species has no usable wild reference sequence for this segment --
    the caller MUST treat that as `no_parental_sequence`, never fold an empty window.
    """
    sequences = _wild_species_genomic_sequences(reference, species, segment)
    if not sequences:
        return ""
    five, three = empirical_panhandle_consensus(sequences, n=n)
    return five if end == "5" else three


def species_consensus_cds_window(
        reference: ReferenceData, species: str, segment: str,
        window: int = RPI_MPOL_CDS_WINDOW_NT) -> str:
    """Majority-vote consensus CDS-start window for one species' segment (D3, M_pol) --
    the "implicated parent/donor consensus" score_m_pol compares the candidate against.

    Same species-level-substitute reasoning as species_consensus_panhandle_window, but
    over ORF-anchored windows (`_wild_species_cds_windows`) rather than genomic-terminus
    ones: correctly annotated ORFs across isolates of one species start at the same
    relative position (barring rare indels), so majority-voting per column is just as
    meaningful here as it is over genomic termini.

    Returns "" if the species has no usable wild reference sequence for this segment --
    the caller MUST treat that as `no_parental_sequence`, never compare against an
    empty window.
    """
    windows = _wild_species_cds_windows(reference, species, segment, window)
    if not windows:
        return ""
    five, _three = empirical_panhandle_consensus(windows, n=window)
    return five


def rpi_species_assignment(
        gate1: "Gate1Result", odd_segment: str) -> Tuple[str, str]:
    """(donor_species, backbone_species) for the candidate's odd segment (D3).

    SHARED by score_m_pack and score_m_pol: both need the same answer to "which
    species is this candidate's odd segment implicated against, and which species is
    its backbone" -- there is exactly one such assignment per candidate, not one per
    factor, so both call the same function rather than risk the two factors silently
    disagreeing about who the backbone is.

    donor_species is gate1's own top hit for `odd_segment` -- the species whose
    reference consensus stands in for "the segment this candidate's odd segment
    resembles". backbone_species is the OTHER implicated species, when exactly one
    exists.

    Two well-defined outcomes and one degenerate one:
      - Exactly one OTHER implicated species -> (donor, that species). The ordinary
        single-discordant-segment case.
      - Zero other implicated species (gate1 pattern "concordant": every segment,
        including `odd_segment`, best-matches the SAME species) -> (donor, donor).
        Comparing that species' own consensus against itself is the CORRECT neutral
        case here, not the self-comparison bug D3 forbids -- see score_m_pack's
        docstring for why these are different things.
      - Two or more OTHER implicated species (every segment disagrees, no single
        backbone) -> ("", "") -- gate1's majority-vote heuristic does not resolve a
        single backbone here, and forcing a choice would be arbitrary. Caller must
        treat this as `no_parental_sequence`.
    """
    hits = gate1.top_hits.get(odd_segment) or []
    donor = hits[0].species if hits else ""
    if not donor:
        return "", ""
    others = sorted({s for s in gate1.implicated_species if s != donor})
    if not others:
        return donor, donor
    if len(others) == 1:
        return donor, others[0]
    return donor, ""


def score_m_pack(
        candidate: CandidateInput, reference: ReferenceData, gate1: "Gate1Result",
        odd_segment: str, weights: Optional[Dict[str, Any]] = None,
) -> Tuple[float, str]:
    """Factor M_pack (D3): does the candidate's actual cross-segment panhandle duplex
    fold as tightly as a genuine single-species (backbone) pairing would?

    Returns (value, status). Status is one of:
      "scored"               -- at least one candidate-vs-backbone-consensus duplex pair
                                was folded successfully.
      "no_parental_sequence" -- the donor/backbone species could not be determined
                                unambiguously (rpi_species_assignment), or neither has
                                a usable wild reference sequence for the relevant
                                segment (species_consensus_panhandle_window).
      "undefined_energy"     -- ViennaRNA (or boltzmann_compatibility, N4) could not
                                produce a usable ddG for ANY pair -- caught here and
                                logged via the status, per the global "log and continue"
                                constraint. compute_rpi floors m_pack at `eco_floor`
                                regardless, so this never propagates a fabricated value.
      "fixed_constant"       -- reserved for a future m_pack_variance_check analogous to
                                m_pol_variance_check; not produced by this function. The
                                WS3.1 calibration measured real, if small, signal
                                (Cliff's delta = +0.222) -- not the negligible variance
                                m_pol_variance_check exists to catch -- so fixing M_pack
                                at a constant is not warranted by that evidence.

    NEVER SELF-COMPARISON (D3). This is the exact defect class M_pol currently has
    (`terminal_windows` passed as BOTH arguments to `polymerase_compatibility`, so it is
    identically 1.0 for every candidate) -- and it is avoided here because the two sides
    of every fold come from genuinely different sources: the "actual" side is the
    CANDIDATE's own genomic sequence; the "reference" side is the BACKBONE SPECIES' OWN
    CONSENSUS, sourced from `reference.genomic_index`. The rpi_species_assignment
    donor==backbone case (gate1 pattern "concordant") looks superficially similar but is
    NOT the same bug: both sides still come from different sources (this candidate's
    actual sequence vs that species' independently-built reference consensus), and a
    candidate that genuinely matches one species on every segment SHOULD score near 1.0
    -- that is the correct answer, not an artifact of comparing something to itself.

    Averages ddG over every inter-segment duplex touching `odd_segment` -- both
    orientations against each of the two backbone segments, up to 4 observations -- the
    same exhaustive treatment rpi_calibration.py's SEGMENT_DUPLEXES uses. There is no
    single canonical pairing among three segments (the real intramolecular panhandle is
    per-segment; the cross-segment duplex tested here is D2's modelling proxy), so
    averaging over all of them avoids an arbitrary choice of just one.
    """
    table = _weights(weights)
    donor_species, backbone_species = rpi_species_assignment(gate1, odd_segment)
    if not donor_species or not backbone_species:
        return 1.0, "no_parental_sequence"

    candidate_genomic = {
        "S": candidate.s_genomic, "M": candidate.m_genomic, "L": candidate.l_genomic}
    backbone_segments = tuple(s for s in ("S", "M", "L") if s != odd_segment)

    # Local cache: within one call, the odd segment's consensus/candidate windows are
    # reused across both backbone segments, and each backbone segment's own windows are
    # otherwise computed twice (once per orientation). Scoped to this call only -- the
    # module-level memoisation a production WS-API loop would want is WS-API's D9, not
    # this task's.
    consensus_cache: Dict[Tuple[str, str, str], str] = {}
    candidate_window_cache: Dict[Tuple[str, str], str] = {}

    def consensus(species: str, segment: str, end: str) -> str:
        key = (species, segment, end)
        if key not in consensus_cache:
            consensus_cache[key] = species_consensus_panhandle_window(
                reference, species, segment, end, n=int(table["m_pack_window_nt"]))
        return consensus_cache[key]

    def candidate_window(segment: str, end: str) -> str:
        key = (segment, end)
        if key not in candidate_window_cache:
            candidate_window_cache[key] = panhandle_window(
                candidate_genomic[segment], end, n=int(table["m_pack_window_nt"]))
        return candidate_window_cache[key]

    ddgs: List[float] = []
    undefined_energy = False
    for backbone_segment in backbone_segments:
        # Two orientations: odd-segment-5' against backbone-3', and backbone-5' against
        # odd-segment-3'. Both are real candidates for a panhandle duplex (D2), and
        # neither is more canonical than the other.
        for five_segment, three_segment in (
                (odd_segment, backbone_segment), (backbone_segment, odd_segment)):
            cognate_five = consensus(backbone_species, five_segment, "5")
            cognate_three = consensus(backbone_species, three_segment, "3")
            if not cognate_five or not cognate_three:
                continue
            try:
                dg_candidate = cofold_dg(
                    candidate_window(five_segment, "5"), candidate_window(three_segment, "3"))
                dg_cognate = cofold_dg(cognate_five, cognate_three)
                ddgs.append(dg_candidate - dg_cognate)
            except (ValueError, RuntimeError):
                undefined_energy = True
                continue

    if not ddgs:
        return (1.0, "undefined_energy") if undefined_energy else (1.0, "no_parental_sequence")

    mean_ddg = float(np.mean(ddgs))
    try:
        value = boltzmann_compatibility(mean_ddg, 0.0, rt=table["m_pack_tau"])
    except ValueError:
        return 1.0, "undefined_energy"
    return value, "scored"


def score_m_pol(
        candidate_cds: Dict[str, str], reference: ReferenceData, gate1: "Gate1Result",
        odd_segment: str, window: int = RPI_MPOL_CDS_WINDOW_NT,
) -> Tuple[float, str]:
    """Factor M_pol (D3, WS3.4): does the candidate's actual N-recognition footprint,
    beyond the genus-invariant tag, match what the BACKBONE species' own reference
    consensus looks like?

    Returns (value, status). Status is whatever `polymerase_compatibility` itself
    returns ("scored" or "no_signal_beyond_tag" -- see its docstring), UNLESS the
    backbone species or its reference sequence cannot be determined at all, in which
    case this returns "no_parental_sequence" -- the same vocabulary score_m_pack uses
    for the same reason: this function never got as far as attempting a comparison.

    "fixed_constant" is DELIBERATELY NOT produced here. Spec 6.2's mandatory variance
    check (`m_pol_variance_check`) is a RUN-LEVEL check over "a representative sample of
    lineage pairs" -- it decides once, for a whole deployment, whether M_pol is
    informative enough to admit to the product; it is not a per-candidate computation,
    and score_candidate scores exactly one candidate. Wiring that run-level check in is
    future work (spec 6.2 also notes "no real sample of lineage pairs to run over in
    production" existed until this task's real parental termini did), not something
    Task 11 asked for -- Task 11's own checklist is "stop self-comparison" and "stop
    forcing status", both of which this function does.

    NEVER SELF-COMPARISON (D3), for the exact reason this factor existed as a
    self-comparison bug before this task: `candidate_termini` is built from the
    CANDIDATE's own CDS; `parent_termini` is built from the BACKBONE SPECIES' OWN
    REFERENCE CONSENSUS (`species_consensus_cds_window`), sourced from
    `reference.genomic_index` -- two genuinely different sources.

    ONE backbone species covers all three segments here, same reasoning as
    score_m_pack: the question is "does the candidate look like a genuine backbone
    genome," so the yardstick for every segment -- including the truly-donated one --
    is what the backbone's OWN native segment would look like, not the donor's. Using
    the donor's own consensus for the donated segment would make M_pol insensitive to
    reassortment by construction (a genuinely-donated segment would always match its
    own true species closely).

    `footprint` is sourced from `n_recognition_footprint(backbone_species)` -- the
    species-specific N-recognition length (Mir & Panganiban's SNV-23-vs-genus-9
    asymmetry) that was previously never wired to any real species at any production
    call site, because no real backbone species existed to look it up for.
    """
    _donor, backbone_species = rpi_species_assignment(gate1, odd_segment)
    if not backbone_species:
        return 1.0, "no_parental_sequence"

    candidate_termini = {seg: candidate_cds.get(seg, "")[:window] for seg in ("S", "M", "L")}
    parent_termini = {
        seg: species_consensus_cds_window(reference, backbone_species, seg, window=window)
        for seg in ("S", "M", "L")
    }
    if not any(parent_termini.values()):
        return 1.0, "no_parental_sequence"

    footprint = n_recognition_footprint(backbone_species)
    return polymerase_compatibility(candidate_termini, parent_termini, footprint=footprint)


def retention_factor(odd_segment: str, structural_intact: bool,
                     weights: Optional[Dict[str, Any]] = None) -> float:
    """Factor F: scales rpi_generation into rpi_retention.

    SCOPE (spec 0.7, 8): models IMMEDIATE post-reassortment fitness only. McDonald's
    Fig. 4 mode (c) -- compensatory mutations restoring interaction interfaces over time
    -- is OUT OF SCOPE and is a named limitation, not a term here. A reassortant scoring
    low may still emerge after adaptation.
    """
    table = _weights(weights)
    base = table["retention_by_odd_segment"].get(odd_segment, 0.15)
    if not structural_intact:
        base *= table["retention_structural_penalty"]
    return max(table["eco_floor"], min(1.0, base))


# Tag identifying which revision of RPI_WEIGHTS produced a given rpi_* row -- not itself
# a weight, so it lives at module scope rather than inside RPI_WEIGHTS. Bump this any
# time RPI_WEIGHTS changes so historical output rows remain attributable.
RPI_WEIGHTS_VERSION = "2026-08-29.1"

# Emitted as the rpi_k_stratum_provisional column on every row.
#
# TRUE because the `identity` compute_rpi is currently handed by the production call
# site is gate0.best_kmer_identity -- PROTEIN k-mer novelty against the nearest
# reference, maximised over segments -- and NOT the parental NUCLEOTIDE divergence that
# spec section 7's close/distant stratification is defined on. The two are different
# quantities on different scales (best_kmer_identity runs ~0.35-0.46 in this dataset;
# the boundary's calibration points are 0.548 / 0.707 / 0.85 nucleotide identity), so
# rpi_k_stratum must NOT be read as a real divergence stratification.
#
# Computing real parental identity is a documented follow-up, not a rename: it needs
# parental sequences sourced per implicated species from the reference index. Flip this
# constant to False in the same change that lands it, and bump RPI_WEIGHTS_VERSION.
RPI_K_STRATUM_PROVISIONAL = True


def _rpi_decimals(floor_value: float) -> int:
    """Decimal places at which `floor_value` still survives round() as nonzero.

    A column whose smallest legitimate value is 1e-6 cannot be reported at 4 dp: the
    floor exists precisely so no factor is ever exactly zero, and rounding coarser than
    the floor throws that guarantee away at the last step, in the CSV a reader actually
    sees. Deriving the precision from the floor keeps the two in step automatically if
    eco_floor is ever retuned, instead of leaving a hand-picked literal to drift.

    Never returns fewer than 4, so the human-readable columns stay readable.
    """
    return max(4, int(math.ceil(-math.log10(floor_value))))


def compute_rpi(
    implicated_species: List[str], isolates: List["ah.IsolateRecord"],
    parent_isolates: Dict[str, "ah.IsolateRecord"], isolate_basis: str,
    odd_segment: str, identity: float, m_pack: float, m_pack_status: str, m_pol: float,
    m_pol_status: str, structural_intact: bool,
    weights: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Assemble the Reassortment Plausibility Index.

        rpi_generation = E * M_pack * M_pol * K
        rpi_retention  = rpi_generation * F

    MULTIPLICATIVE, not additive -- a deliberate departure from surveillance_priority's
    house style. Lowen (2018) Fig. 1 makes each stage a NECESSARY gate; McDonald's
    cognate-competition result is inherently a ratio; and an additive score would let
    molecular elegance compensate for ecological impossibility.

    This index is RELATIVE. It is not a probability, and must never be reported as one.

    `m_pack_status` is REQUIRED, not defaulted, on the same reasoning `m_pol_status`
    already is: a caller must state where its m_pack value came from (see
    score_m_pack's four statuses) rather than let a constant land in the report dressed
    as a measurement. `parent_isolates`/`isolate_basis` (design doc
    2026-08-30-rpi-e-geo-isolate-level-fix-design.md) are REQUIRED for the same reason --
    see resolve_e_parent_isolates and ecological_opportunity's own docstring.
    """
    # Every lookup below goes through `table`, and every helper is handed it
    # explicitly. Nothing on this path reads the RPI_WEIGHTS global directly, which is
    # what lets the sensitivity sweep score a perturbed table without reassigning it.
    table = _weights(weights)
    floor = table["eco_floor"]
    eco = ecological_opportunity(
        implicated_species, isolates, parent_isolates, isolate_basis, weights=table)
    stratum = divergence_stratum(identity, weights=table)
    k_value = constellation_prior(stratum, weights=table)[odd_segment]
    f_value = retention_factor(odd_segment, structural_intact, weights=table)

    generation = eco["rpi_e_total"] * max(m_pack, floor) * max(m_pol, floor) * k_value
    # floor ** 4: the product of four factors each individually floored at `floor` can
    # still underflow toward zero in float terms as the factors compound. Re-floor the
    # PRODUCT at floor**4 (not floor) so the no-literal-zeros guarantee holds for the
    # combined quantity too, not just for each input factor in isolation.
    generation = max(floor ** 4, min(1.0, generation))

    row = dict(eco)
    row.update({
        "rpi_m_pack": round(m_pack, 4),
        "rpi_m_pack_status": m_pack_status,
        "rpi_m_pol": round(m_pol, 4),
        "rpi_m_pol_status": m_pol_status,
        "rpi_k_stratum": stratum,
        # Explicit honesty marker -- see RPI_K_STRATUM_PROVISIONAL. Without it a reader
        # would take rpi_k_stratum for a real parental-divergence stratification.
        "rpi_k_stratum_provisional": RPI_K_STRATUM_PROVISIONAL,
        "rpi_k_value": round(k_value, 4),
        "rpi_f": round(f_value, 4),
        "rpi_smoothing": table["jeffreys_alpha"],
        "rpi_weights_version": RPI_WEIGHTS_VERSION,
        # Precision per column is set by THAT COLUMN'S OWN FLOOR, not picked by eye.
        # Each is a product of floored factors, so each has a different smallest
        # legitimate value, and rounding any of them coarser than its floor collapses a
        # small-but-nonzero score to a literal 0.0 -- the hidden hard gate spec 7
        # forbids. See _rpi_decimals.
        #
        #   rpi_generation = E * M_pack * M_pol * K, clamped at floor**4 = 1e-8
        #   rpi_retention  = rpi_generation * F, and F is itself floored at eco_floor,
        #                    so the product bottoms out one factor lower, at
        #                    floor**5 = 1e-10. At 8 dp that rounded to 0.0: generation
        #                    at its clamp with F = 0.075 (odd segment S/L,
        #                    structural_intact False) gives 7.5e-10.
        #
        # BOTH use the deeper of the two precisions, and must. F <= 1, so
        # rpi_retention <= rpi_generation is an invariant of the model
        # (test_retention_never_exceeds_generation). Rounding is only order-preserving
        # at EQUAL precision -- rounding retention finer than generation let a retention
        # of 1.57233e-05 be reported above a generation of 1.572e-05, inverting the pair
        # by pure formatting.
        "rpi_generation": round(generation, _rpi_decimals(floor ** 5)),
        "rpi_retention": round(generation * f_value, _rpi_decimals(floor ** 5)),
    })
    # rpi_e_host/geo/time are single floored terms (min eco_floor = 0.01), but
    # rpi_e_total is their PRODUCT and bottoms out at floor**3 = 1e-6. Rounding it to 4
    # dp alongside its own components reported a literal 0.0 for entirely ordinary
    # candidates -- two implicated species in different host orders, different
    # countries, a few years apart gives 3.33e-05, which read as 0.0 in the CSV.
    # EcologicalOpportunityTest.test_no_factor_is_ever_exactly_zero did not catch it
    # because it asserts on ecological_opportunity's RAW output, not on this rounded row.
    #
    # THE WHOLE FAMILY SHARES THE DEEPEST PRECISION, for the same reason rpi_generation
    # and rpi_retention do. rpi_e_total = h * g * t with every term <= 1, so
    # rpi_e_total <= min(h, g, t) is an invariant, and round() only preserves order at
    # EQUAL precision. Reporting the total at 6 dp beside its components at 4 dp put the
    # PRODUCT above its own factor: h = g = 1.0 with t = 1/3 gave rpi_e_total = 0.333333
    # against rpi_e_time = 0.3333. Arithmetically impossible, and it would not survive a
    # reader recomputing h * g * t from the reported columns.
    e_decimals = _rpi_decimals(floor ** 3)
    for key in ("rpi_e_host", "rpi_e_geo", "rpi_e_time", "rpi_e_total"):
        row[key] = round(row[key], e_decimals)
    return row


# ---------------------------------------------------------------------------
# Retrospective in vitro reproduction harness -- the PRIMARY validation
# evidence for the RPI (spec 11.1), ranking above the sensitivity sweep in
# Task 9 because the pipeline carries only 8 positive labels, far too few to
# validate anything statistically.
#
# Ecology is held at 1.0 throughout: every experiment below FORCED
# co-infection by co-inoculating cells, so ecological opportunity was never a
# variable in them. That is exactly why this data can calibrate the
# molecular factors (constellation_prior) but cannot calibrate Factor E --
# Factor E is never computed here.
# ---------------------------------------------------------------------------

# Published in vitro coinfection experiments, used as the PRIMARY validation evidence
# (spec 11.1). Ecology is held at 1.0: these experiments forced coinfection.
#
# `contributes_to_prior` records whether a pair's counts are actually inside
# RPI_WEIGHTS["constellation_counts"]. Only three of the five are:
#     close   {M: 12, S: 12, L: 1} = SNV_NMR11_x_CC107 alone
#     distant {M: 11, S:  0, L: 0} = SNV_x_ANDV (10) + SNV_x_BCCV (1)
# The two Tier-3 pairs (PUUV_x_PHV, DOBV_Aa_x_DOBV_Af) were deliberately NOT folded
# into the weights, so holding either out changes nothing and neither may be listed as
# a contributor. This flag is what makes the leave-one-pair-out subtraction in
# retrospective_reproduction correct rather than approximate.
RPI_INVITRO_BENCHMARK: Dict[str, Dict[str, Any]] = {
    "SNV_NMR11_x_CC107": {
        "source": "Rodriguez et al. 1998, Virology 242:99-106, Table 7 "
                  "(TIER 1 PRIMARY -- counts read off the source PDF)",
        "stratum": "close", "moi": "1.0-1.5", "genotyping": "plaque-level, no amplification",
        "observed": {"M": 12, "S": 12, "L": 1},
        "contributes_to_prior": True,
    },
    "SNV_x_BCCV": {
        "source": "Rodriguez et al. 1998, Virology 242:99-106 "
                  "(TIER 1 PRIMARY -- counts read off the source PDF)",
        "stratum": "distant", "moi": "1.0-1.5", "genotyping": "plaque-level, no amplification",
        "observed": {"M": 1, "S": 0, "L": 0},
        "contributes_to_prior": True,
    },
    "SNV_x_ANDV": {
        "source": "Rizvanov et al. 2004, Virology 327:225-232 "
                  "(TIER 1 PRIMARY -- counts read off the source PDF)",
        "stratum": "distant", "moi": "3", "genotyping": "plaque-level, then passaged",
        "observed": {"M": 10, "S": 0, "L": 0},
        "contributes_to_prior": True,
    },
    "PUUV_x_PHV": {
        "source": "Handke et al. 2010, J Gen Virol 91:2351-2359 (TIER 3 -- read the "
                  "library PDF and promote before relying on this row)",
        "stratum": "distant", "moi": "unrecorded", "genotyping": "unrecorded",
        "observed": {"M": 1, "S": 0, "L": 0},
        "contributes_to_prior": False,
    },
    "DOBV_Aa_x_DOBV_Af": {
        "source": "Kirsanovs et al. 2010, Virus Genes 41:319-328, via Klempa 2018 "
                  "(TIER 3 SECONDHAND)",
        "stratum": "close", "moi": "unrecorded", "genotyping": "unrecorded",
        "observed": {"M": 65, "S": 0, "L": 0},
        "contributes_to_prior": False,
    },
}

RPI_BENCHMARK_CAVEAT = (
    "Divergence is confounded with pair identity AND with MOI: every distant observation "
    "involves a different virus pair from the close observation, and no distant-pair MOI-1 "
    "experiment exists in the literature. This stratification is UNTESTABLE with current "
    "data. Trans-complementation bias applies asymmetrically -- Rizvanov's design passages "
    "plaques before genotyping, Rodriguez's does not."
)


def _composition_distance(prior: Dict[str, float], observed: Dict[str, float]) -> float:
    """Total-variation distance between the prior and the observed composition.

    0.5 * sum(|prior[s] - observed_fraction[s]|) over S, M, L. Ranges 0 (identical) to
    1 (disjoint support).

    This exists because `modal_agrees` SATURATES -- see retrospective_reproduction. A
    boolean "is the observed top segment in the predicted top set?" cannot fail when
    every predicted top set contains M and every observed top set contains M. The
    total-variation distance is a real quantitative comparison and does fail: an
    all-M observation against a prior that weights S equally to M scores ~0.53, which
    the boolean reports as a clean pass.

    Returns nan for an empty observation rather than a flattering 0.0 -- an unmeasured
    pair is not a perfectly-matching pair.
    """
    total = sum(max(0.0, observed.get(segment, 0)) for segment in ("S", "M", "L"))
    if total <= 0:
        return float("nan")
    return 0.5 * sum(
        abs(prior[segment] - max(0.0, observed.get(segment, 0)) / total)
        for segment in ("S", "M", "L")
    )


def _prior_counts_excluding(holdout: Optional[str]) -> Dict[str, Dict[str, float]]:
    """RPI_WEIGHTS' constellation counts with `holdout`'s own observations removed.

    THIS IS WHAT MAKES LEAVE-ONE-PAIR-OUT ACTUALLY HOLD SOMETHING OUT. The previous
    version filtered a `derived_from` list and then called constellation_prior, which
    reads the fixed weight table -- so the held-out pair's counts still generated the
    prediction that was supposed to test it independently, and the prior came back
    byte-identical with and without the holdout.

    Pairs flagged contributes_to_prior=False are not inside the weight table at all,
    so holding them out is a no-op here -- correctly, and visibly so.
    """
    counts = copy.deepcopy(RPI_WEIGHTS["constellation_counts"])
    if holdout is None:
        return counts
    record = RPI_INVITRO_BENCHMARK[holdout]  # KeyError on a typo, deliberately loud
    if not record.get("contributes_to_prior"):
        return counts
    stratum = counts[record["stratum"]]
    for segment, observed_count in record["observed"].items():
        stratum[segment] = max(0, stratum.get(segment, 0) - observed_count)
    return counts


def retrospective_reproduction(holdout: Optional[str] = None) -> Dict[str, Any]:
    """Primary validation: does the index reproduce the published constellation profiles?

    Pass `holdout` to exclude a pair from prior derivation, so it becomes an independent
    test rather than a restatement of its own input. The exclusion is a real count
    SUBTRACTION (_prior_counts_excluding), not a list filter -- Jeffreys smoothing is
    applied afterwards, so no category can hit zero even when the subtraction empties it.

    Two agreement measures are reported per pair:
      `composition_distance` -- PRIMARY. Total-variation distance between the prior and
          the observed composition. Quantitative, and able to fail.
      `modal_agrees`         -- SECONDARY, and it SATURATES. The distant top set is {M},
          the close top set is {M, S}, and every benchmark observation's top set
          contains M -- so no possible input in this benchmark returns False. It is kept
          because it is the comparison the spec's prose describes, but it must never be
          quoted as a pass on its own; read composition_distance beside it.
    """
    alpha = RPI_WEIGHTS["jeffreys_alpha"]
    counts = _prior_counts_excluding(holdout)
    contributors = [name for name, record in RPI_INVITRO_BENCHMARK.items()
                    if record.get("contributes_to_prior")]
    derived_from = [name for name in contributors if name != holdout]
    pairs: Dict[str, Any] = {}
    for name, record in RPI_INVITRO_BENCHMARK.items():
        prior = smoothed_prior(counts[record["stratum"]], alpha)
        observed = record["observed"]

        # Ties must be represented explicitly, not resolved by dict insertion order.
        # The close stratum gives M and S EXACTLY equal weight (12 vs 12), and
        # Rodriguez's own observed counts are tied the same way -- max(d, key=d.get)
        # would silently pick whichever key comes first in each dict's own literal
        # order (constellation_prior orders "S" before "M"; RPI_INVITRO_BENCHMARK's
        # `observed` dicts order "M" before "S"), manufacturing a spurious
        # disagreement out of a perfect tie-for-tie match.
        predicted_top = sorted(seg for seg, val in prior.items() if val == max(prior.values()))
        observed_top = sorted(seg for seg, val in observed.items() if val == max(observed.values()))
        modal = predicted_top[0] if len(predicted_top) == 1 else None
        observed_modal = observed_top[0] if len(observed_top) == 1 else None

        pairs[name] = {
            "source": record["source"],
            "stratum": record["stratum"],
            "moi": record["moi"],
            "genotyping": record["genotyping"],
            "observed": observed,
            "predicted_prior": prior,
            "predicted_modal_segment": modal,
            "predicted_top_segments": predicted_top,
            "observed_modal_segment": observed_modal,
            "observed_top_segments": observed_top,
            # PRIMARY quantitative comparison -- see _composition_distance.
            "composition_distance": _composition_distance(prior, observed),
            # SECONDARY and saturating -- see this function's docstring.
            "modal_agrees": bool(set(predicted_top) & set(observed_top)),
        }
    return {
        "pairs": pairs,
        "derived_from": derived_from,
        "prior_counts_used": counts,
        "contributes_to_prior": contributors,
        "holdout": holdout,
        "holdout_changed_the_prior": bool(
            holdout is not None
            and counts != RPI_WEIGHTS["constellation_counts"]
        ),
        "ecological_factor_held_at": 1.0,
        "primary_agreement_measure": "composition_distance",
        "modal_agrees_saturates": (
            "modal_agrees is SECONDARY: the distant top set is {M}, the close top set "
            "is {M, S}, and every benchmark observation's top set contains M, so no "
            "input in this benchmark can return False. Read composition_distance."
        ),
        "caveat": RPI_BENCHMARK_CAVEAT,
    }


# ---------------------------------------------------------------------------
# Synthetic IsolateRecord builders.
#
# These live HERE, in the production area, and not in the Tests section, because
# _rpi_sensitivity_panel (production) calls _eco_isolate. Under `import
# reassortment_detector` the name resolved at call time regardless of definition
# order, but under script-mode invocation (`python reassortment_detector.py`, which
# stops at `sys.exit(main())` before the Tests section is ever executed) the name did
# not exist -- so any production path reaching the sweep raised NameError. Defining
# them above the CLI section removes that dependency on how the module was entered.
#
# The tests below (Gate3CheckTest, EcologicalOpportunityTest) still use them; they are
# module-level names, so nothing about the test side changes.
# ---------------------------------------------------------------------------

def _make_ecology_isolate(species: str, country: str, host_order: str, year: str,
                          lab: bool = False, vaccine: bool = False) -> "ah.IsolateRecord":
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


def _eco_isolate(species: str, country: str, year: str, host_order: str = "Rodentia",
                 host_family: str = "Cricetidae",
                 host_name: str = "Peromyscus maniculatus") -> "ah.IsolateRecord":
    """IsolateRecord has host_name (a binomial), host_order and host_family -- but NO
    host_genus and NO observed_host. Genus is derived from host_name."""
    isolate = _make_ecology_isolate(species, country, host_order, year)
    isolate.host_family = host_family
    isolate.host_name = host_name
    return isolate


def _rpi_sensitivity_panel() -> List[Dict[str, Any]]:
    """Fixed synthetic candidate panel for sensitivity_sweep.

    Exists so the sweep can Spearman-correlate ACTUAL compute_rpi() output under
    perturbed weights against the same under baseline weights -- rather than
    correlating a weight lookup table against a positive scalar multiple of itself,
    which is provably 1.0 regardless of the real answer (see the review that added
    this function: comparing `RPI_WEIGHTS[key].values()` to `v * (1 - factor)` for the
    SAME v's can never show instability).

    Deliberately spans:
      - every odd_segment (S, M, L);
      - host-proximity categories same_species through different_order, PLUS a
        single-implicated-species candidate (Factor E pinned neutral at 1.0 by
        construction -- immune to host_decay) and an unknown-host-metadata candidate;
      - identities on both sides of the stratum boundary in both perturbation
        directions (0.62 crosses at the low boundary; 0.95 crosses at the high one);
      - one matched same-species pair (candidates A and B) whose baseline
        rpi_generation ratio is engineered into (0.15, 0.225) -- exactly the window
        where a +50% retention_by_odd_segment perturbation flips their rpi_retention
        order, because M's retention ratio is already pinned at retention_factor's
        1.0 ceiling and cannot follow the increase while S's (0.15, well under
        the ceiling) does.
    Reuses (does not reimplement) IsolateRecord construction via
    `_eco_isolate`/`_make_ecology_isolate`, defined directly above -- building the
    ~20-field IsolateRecord by hand a second time here would drift from
    EcologicalOpportunityTest's own fixture. Those builders used to live in the Tests
    section, which made this production function depend on a name that does not exist
    under script-mode invocation; they were moved up for that reason.
    """
    def two(label: str, odd_segment: str, identity: float, m_pack: float, m_pol: float,
            host_name_a: str = "Peromyscus maniculatus", host_family_a: str = "Cricetidae",
            host_order_a: str = "Rodentia", host_name_b: str = "Peromyscus maniculatus",
            host_family_b: str = "Cricetidae", host_order_b: str = "Rodentia") -> Dict[str, Any]:
        isolate_a = _eco_isolate(f"Sp_{label}_a", "USA", "2001", host_order=host_order_a,
                                  host_family=host_family_a, host_name=host_name_a)
        isolate_b = _eco_isolate(f"Sp_{label}_b", "USA", "2001", host_order=host_order_b,
                                  host_family=host_family_b, host_name=host_name_b)
        return {
            "label": label, "implicated_species": [isolate_a.species_name, isolate_b.species_name],
            "isolates": [isolate_a, isolate_b], "odd_segment": odd_segment, "identity": identity,
            "m_pack": m_pack, "m_pol": m_pol, "structural_intact": True,
        }

    panel = [
        two("A_same_species_close_M", "M", 0.90, 0.41, 0.41),
        two("B_same_species_close_S", "S", 0.90, 0.95, 0.95),
        two("C_same_genus_close_L", "L", 0.90, 0.5, 0.5, host_name_b="Peromyscus leucopus"),
        two("D_same_family_distant_M", "M", 0.60, 0.3, 0.3, host_name_b="Neotoma albigula"),
        two("E_same_order_distant_S", "S", 0.60, 0.6, 0.6,
            host_family_b="Muridae", host_name_b="Rattus norvegicus"),
        two("F_different_order_distant_L", "L", 0.60, 0.9, 0.9,
            host_order_b="Soricomorpha", host_family_b="Soricidae", host_name_b="Sorex araneus"),
        two("G_boundary_low_crossing_M", "M", 0.62, 0.9, 0.9),
        two("H_boundary_high_crossing_M", "M", 0.95, 0.9, 0.9),
        # K and L: engineered for the jeffreys_alpha direction, the same way A and B
        # are engineered for retention_by_odd_segment. Both sit at same_species host
        # proximity, same country, same year, so Factor E is 1.0 for both and the only
        # thing separating their rpi_generation is k_value x m_pack x m_pol.
        #
        # k(distant, S) and k(close, L) are BOTH near-zero categories, but they move at
        # different rates with alpha -- k(close,L)/k(distant,S) is 2.28 at alpha 0.25,
        # 1.42 at 0.5 and 1.13 at 0.75, because the close stratum's L count is a real
        # observation (1) while the distant stratum's S count is an unobserved zero that
        # alpha alone populates. Scaling L's molecular product to 0.8 (0.8944^2) puts
        # the pair's baseline ratio at 1.13, inside the window where the +50% direction
        # flips their order and the -50% direction does not. Without a pair like this
        # the panel simply has no candidates close enough for alpha to reorder, and the
        # sweep would report stability it never tested.
        two("K_alpha_flip_distant_S", "S", 0.60, 1.0, 1.0),
        two("L_alpha_flip_close_L", "L", 0.90, 0.8944, 0.8944),
    ]
    single = _eco_isolate("Sp_I_single", "USA", "2001")
    panel.append({
        "label": "I_single_species_close_S", "implicated_species": [single.species_name],
        "isolates": [single], "odd_segment": "S", "identity": 0.90,
        "m_pack": 0.6, "m_pol": 0.6, "structural_intact": True,
    })
    unknown_a = _eco_isolate("Sp_J_a", "USA", "2001")
    unknown_b = _eco_isolate("Sp_J_b", "USA", "2001", host_order="", host_family="", host_name="")
    panel.append({
        "label": "J_unknown_host_close_L",
        "implicated_species": [unknown_a.species_name, unknown_b.species_name],
        "isolates": [unknown_a, unknown_b], "odd_segment": "L", "identity": 0.90,
        "m_pack": 0.7, "m_pol": 0.7, "structural_intact": True,
    })
    return panel


def _score_panel(
    panel: Sequence[Dict[str, Any]], weights: Optional[Dict[str, Any]] = None,
) -> Dict[str, List[float]]:
    """rpi_generation and rpi_retention for every panel candidate under `weights`.

    `weights` must be a caller-owned dict (a deep copy of RPI_WEIGHTS with one key
    perturbed) or None to score under the current global unperturbed. It is PASSED DOWN
    to compute_rpi, not installed anywhere.

    This function used to reassign the module-level RPI_WEIGHTS for the duration of the
    scoring loop and restore it in a finally block. That was correct single-threaded and
    unsafe for the use this index is now built for: the reassignment is PROCESS-GLOBAL,
    so anything scoring concurrently -- the WS-API entry point the downstream team calls
    in a loop -- would read perturbed weights for the length of a sweep and return
    numbers attributable to no declared rpi_weights_version. Nothing announced it,
    because every individual row still looked well-formed.

    Passing the table explicitly is what turns "never adjust tau or any weight inside
    the caller's loop" from a README request into a property of the code. See
    ScorePanelDoesNotMutateGlobalWeightsTest.
    """
    table = _weights(weights)
    generation, retention = [], []
    for candidate in panel:
        # This synthetic panel already models "one isolate stands for this species"
        # directly in candidate["isolates"] -- its own known, real parent, not a
        # production-style unknown -- so it maps straight to parent_isolates/
        # "actual_parents" with no fallback resolution needed.
        parent_isolates = {isolate.species_name: isolate for isolate in candidate["isolates"]}
        row = compute_rpi(
            candidate["implicated_species"], candidate["isolates"],
            parent_isolates, "actual_parents", candidate["odd_segment"],
            candidate["identity"], candidate["m_pack"], "scored", candidate["m_pol"], "scored",
            candidate["structural_intact"], weights=table,
        )
        generation.append(row["rpi_generation"])
        retention.append(row["rpi_retention"])
    return {"generation": generation, "retention": retention}


def sensitivity_sweep(factor: float = 0.5) -> Dict[str, Any]:
    """Vary every ASSUMED parameter by +/- `factor` and report rank stability.

    Required before any RPI result is reported (spec 12). The headline question is:
    do the top-ranked candidates stay top-ranked? Answered against a fixed synthetic
    candidate panel (_rpi_sensitivity_panel), scored end-to-end through compute_rpi, so
    the answer reflects real candidate reordering -- not a weight table compared to a
    scalar multiple of itself.

    `rank_correlation` for host_decay, retention_by_odd_segment and
    stratum_identity_boundary is the WORSE (lower) of the two perturbation directions
    (+factor, -factor) -- the sweep's headline question is whether ranking survives,
    so the more conservative direction is reported, not whichever happens to look
    stabler.

    Operates on deep copies -- never mutates RPI_WEIGHTS. `_score_panel` restores the
    original RPI_WEIGHTS object even if scoring raises partway through.
    """
    varied: Dict[str, Any] = {}
    panel = _rpi_sensitivity_panel()

    baseline_generation = _score_panel(panel)["generation"]

    # --- jeffreys_alpha: routed through _score_panel like every other parameter.
    #
    # The ORIGINAL version Spearman-correlated [close_M_prior, distant_M_prior] against
    # itself under perturbation. distant-M exceeds close-M at every alpha, so that rank
    # order is fixed BY CONSTRUCTION and the statistic returned 1.0 unconditionally --
    # a tautology, exactly the defect the host_decay/retention entries were already
    # fixed for. alpha is reachable via RPI_WEIGHTS["jeffreys_alpha"], so a deep-copied
    # perturbed weights dict routes it through compute_rpi unchanged and the answer is
    # now about real candidate reordering. The two-strata `samples` are retained as
    # descriptive context, but they are NOT what rank_correlation is computed from.
    low_alpha = RPI_WEIGHTS["jeffreys_alpha"] * (1 - factor)
    high_alpha = RPI_WEIGHTS["jeffreys_alpha"] * (1 + factor)
    alpha_low_weights = copy.deepcopy(RPI_WEIGHTS)
    alpha_low_weights["jeffreys_alpha"] = low_alpha
    alpha_high_weights = copy.deepcopy(RPI_WEIGHTS)
    alpha_high_weights["jeffreys_alpha"] = high_alpha
    alpha_gen_low = _score_panel(panel, alpha_low_weights)["generation"]
    alpha_gen_high = _score_panel(panel, alpha_high_weights)["generation"]
    varied["jeffreys_alpha"] = {
        "samples": [
            {"alpha": low_alpha,
             "m_prior": [constellation_prior(s, alpha=low_alpha)["M"] for s in ("close", "distant")]},
            {"alpha": high_alpha,
             "m_prior": [constellation_prior(s, alpha=high_alpha)["M"] for s in ("close", "distant")]},
        ],
        "rank_correlation": _worst_rank_correlation(
            _spearman(baseline_generation, alpha_gen_low),
            _spearman(baseline_generation, alpha_gen_high)),
        "note": "Correlated against rpi_generation on the fixed synthetic panel "
                "(_rpi_sensitivity_panel), the same route as every other parameter. "
                "'samples' is descriptive only.",
    }

    # --- host_decay: affects Factor E (rpi_generation), via ecological_opportunity's
    # host term. Correlate against the panel's rpi_generation, not the weight table.
    host_decay_low_weights = copy.deepcopy(RPI_WEIGHTS)
    host_decay_low_weights["host_decay"] = _perturbed_host_decay(1 - factor)
    host_decay_high_weights = copy.deepcopy(RPI_WEIGHTS)
    host_decay_high_weights["host_decay"] = _perturbed_host_decay(1 + factor)
    host_decay_gen_low = _score_panel(panel, host_decay_low_weights)["generation"]
    host_decay_gen_high = _score_panel(panel, host_decay_high_weights)["generation"]
    varied["host_decay"] = {
        "scaled_by": (1 - factor, 1 + factor),
        "values_low": host_decay_low_weights["host_decay"],
        "values_high": host_decay_high_weights["host_decay"],
        "excluded_from_perturbation": sorted(RPI_HOST_DECAY_NOT_PERTURBED),
        "rank_correlation": _worst_rank_correlation(
            _spearman(baseline_generation, host_decay_gen_low),
            _spearman(baseline_generation, host_decay_gen_high)),
        "note": "Correlated against rpi_generation on the fixed synthetic panel "
                "(_rpi_sensitivity_panel), not against the weight table itself. The "
                "neutral unknown_host_metadata sentinel is held fixed and the "
                "perturbed weights are clamped -- see _perturbed_host_decay.",
    }

    # --- retention_by_odd_segment: affects Factor F (rpi_retention), NOT
    # rpi_generation -- compute_rpi never multiplies F into rpi_generation. Correlate
    # against the panel's rpi_retention, the quantity this weight actually scales.
    baseline_retention = _score_panel(panel)["retention"]
    retention_low_weights = copy.deepcopy(RPI_WEIGHTS)
    retention_low_weights["retention_by_odd_segment"] = {
        k: v * (1 - factor) for k, v in RPI_WEIGHTS["retention_by_odd_segment"].items()}
    retention_high_weights = copy.deepcopy(RPI_WEIGHTS)
    retention_high_weights["retention_by_odd_segment"] = {
        k: v * (1 + factor) for k, v in RPI_WEIGHTS["retention_by_odd_segment"].items()}
    retention_low = _score_panel(panel, retention_low_weights)["retention"]
    retention_high = _score_panel(panel, retention_high_weights)["retention"]
    varied["retention_by_odd_segment"] = {
        "scaled_by": (1 - factor, 1 + factor),
        "values_low": retention_low_weights["retention_by_odd_segment"],
        "values_high": retention_high_weights["retention_by_odd_segment"],
        "rank_correlation": _worst_rank_correlation(
            _spearman(baseline_retention, retention_low),
            _spearman(baseline_retention, retention_high)),
        "note": "Correlated against rpi_retention (the factor this weight actually "
                "scales), not rpi_generation, on the same fixed panel.",
    }

    # --- stratum_identity_boundary: reclassification is a STEP change (spec 7), not a
    # smooth one. Report which panel candidates actually move between strata at each
    # swept boundary -- observed on this run, never a hardcoded constant -- alongside
    # the resulting rpi_generation rank correlation.
    boundary = RPI_WEIGHTS["stratum_identity_boundary"]
    low_boundary = boundary * (1 - factor)
    high_boundary = min(1.0, boundary * (1 + factor))
    moved_low = [c["label"] for c in panel if divergence_stratum(c["identity"])
                 != divergence_stratum(c["identity"], boundary=low_boundary)]
    moved_high = [c["label"] for c in panel if divergence_stratum(c["identity"])
                  != divergence_stratum(c["identity"], boundary=high_boundary)]
    boundary_low_weights = copy.deepcopy(RPI_WEIGHTS)
    boundary_low_weights["stratum_identity_boundary"] = low_boundary
    boundary_high_weights = copy.deepcopy(RPI_WEIGHTS)
    boundary_high_weights["stratum_identity_boundary"] = high_boundary
    boundary_gen_low = _score_panel(panel, boundary_low_weights)["generation"]
    boundary_gen_high = _score_panel(panel, boundary_high_weights)["generation"]
    varied["stratum_identity_boundary"] = {
        "low": low_boundary, "high": high_boundary,
        "moved_low": moved_low, "moved_high": moved_high,
        "rank_correlation": _worst_rank_correlation(
            _spearman(baseline_generation, boundary_gen_low),
            _spearman(baseline_generation, boundary_gen_high)),
        "note": "Reclassification between strata is a step change, not a smooth one; "
                "'moved_low'/'moved_high' name which panel candidates ACTUALLY changed "
                "stratum at each swept boundary on this run -- an empty list is itself "
                "an observed result (this boundary choice did not move that candidate "
                "in that direction), never a hardcoded assumption of stability.",
    }
    return {
        "factor": factor,
        # The panel is SYNTHETIC and FIXED (_rpi_sensitivity_panel) -- twelve constructed
        # candidates chosen to span the parameter space, not real pipeline rows. This
        # marker is top-level so no number below can be quoted as a statement about
        # real candidates without the reader seeing what it was measured on.
        "panel": "synthetic_fixed",
        # Reported rather than written into prose so the count cannot drift out of date
        # as candidates are added to _rpi_sensitivity_panel.
        "panel_size": len(panel),
        "panel_labels": [candidate["label"] for candidate in panel],
        "varied": varied,
    }


def _tied_ranks(values: Sequence[float]) -> "np.ndarray":
    """Ranks with TIED values sharing their average rank (midrank).

    The obvious double-argsort idiom, np.argsort(np.argsort(v)), does NOT do this: it
    breaks ties by position, so a constant vector comes back as [0, 1, 2, ...] -- a
    perfectly ordered ranking manufactured out of an input with no order in it at all.
    That defeats the zero-variance guard in _spearman (the guard tested the ranks, and
    the ranks were never constant), and it also lets any tie in a real score panel
    contribute a spurious ordering to the correlation.

    Midranks are the standard Spearman treatment of ties and make both problems go
    away: a constant vector ranks flat, so the guard fires, and tied candidates in the
    panel neither gain nor lose order relative to each other.
    """
    array = np.asarray(values, dtype=float)
    # mergesort for a stable, reproducible grouping of equal values.
    order = np.argsort(array, kind="mergesort")
    sorted_values = array[order]
    ranks = np.empty(array.shape, dtype=float)
    start = 0
    while start < len(array):
        stop = start
        while stop + 1 < len(array) and sorted_values[stop + 1] == sorted_values[start]:
            stop += 1
        ranks[order[start:stop + 1]] = (start + stop) / 2.0
        start = stop + 1
    return ranks


def _spearman(a: Sequence[float], b: Sequence[float]) -> float:
    """Spearman rank correlation, or nan where the statistic is UNDEFINED.

    Undefined means: fewer than two observations, or one side's ranks are constant (a
    correlation with a zero-variance variable has no value, and numpy would emit
    0/0 = nan for it anyway). Returning 1.0 there -- as this used to -- reports PERFECT
    RANK STABILITY from an input that measured nothing, which is the same class of
    defect as a constant factor dressed as a measurement. nan matches
    rank_consistency's own discipline for its undefined AUC.

    Callers must surface the nan rather than fold it away; see
    _worst_rank_correlation, which is nan-contagious for exactly that reason.
    """
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    rank_a, rank_b = _tied_ranks(a), _tied_ranks(b)
    if np.std(rank_a) == 0 or np.std(rank_b) == 0:
        return float("nan")
    return float(np.corrcoef(rank_a, rank_b)[0, 1])


def _worst_rank_correlation(*values: float) -> float:
    """The WORSE (lower) of the perturbation directions -- with nan contagious.

    The sweep's headline question is whether ranking survives, so the more conservative
    direction is reported rather than whichever looks stabler. But builtin min() is
    order-dependent around nan (min(nan, 0.5) is nan; min(0.5, nan) is 0.5), so an
    undefined direction could be silently discarded depending on argument order. If any
    direction is undefined, the answer is undefined.
    """
    if any(math.isnan(value) for value in values):
        return float("nan")
    return min(values)


# host_decay keys held FIXED during the sensitivity sweep. unknown_host_metadata is a
# NEUTRAL SENTINEL (1.0 = "unknown, not separated"), not an estimated decay weight.
# Scaling it up turned missing host metadata into a positive bonus and saturated
# compute_rpi's min(1.0, ...) clamp; scaling it down turned unknown into a penalty,
# which is the exact failure host_rank_proximity's top-of-function guard exists to
# prevent. Neither direction is a meaningful perturbation of an assumption.
RPI_HOST_DECAY_NOT_PERTURBED: FrozenSet[str] = frozenset({"unknown_host_metadata"})


def _perturbed_host_decay(scale: float) -> Dict[str, float]:
    """RPI_WEIGHTS['host_decay'] scaled by `scale`, sentinel held fixed and clamped.

    Clamped to [eco_floor, 1.0]:
      - 1.0 above, because host proximity is a discount on ecological opportunity and
        never a bonus -- an unclamped +50% pushed same_species to 1.5, which
        compute_rpi's own min(1.0, ...) then silently absorbed, making the whole
        upward direction of the sweep a no-op for the candidates it should have moved.
      - eco_floor below rather than a literal 0.0, per the no-literal-zeros rule: at
        factor >= 1.0 the low direction would otherwise zero every weight and drive
        Factor E's host term to an exact 0.
    """
    return {
        key: value if key in RPI_HOST_DECAY_NOT_PERTURBED
        else min(1.0, max(RPI_WEIGHTS["eco_floor"], value * scale))
        for key, value in RPI_WEIGHTS["host_decay"].items()
    }


# The production pipeline carries exactly 8 supported_reassortment_candidate labels
# (spec 0.6). 8 positives cannot support fitting or validating anything -- this
# threshold is set well ABOVE that count, not merely above it, so a dataset the same
# size as today's cannot clear "underpowered" by the coincidence of a boundary drawn
# exactly at the current label count. See spec 11.2.
RPI_RANK_MIN_POSITIVES = 30


def rank_consistency(
    scored_rows: Sequence[Dict[str, Any]], verdict_key: str = "verdict",
) -> Dict[str, Any]:
    """SANITY CHECK ONLY (spec 11.2) -- never headline validation.

    With 8 positives against ~3,000 negatives, the interval on any rank statistic is very
    wide. The interval is therefore returned alongside the point estimate and the caller
    must report both.
    """
    positives = [r["rpi_generation"] for r in scored_rows
                 if r.get(verdict_key) == "supported_reassortment_candidate"]
    negatives = [r["rpi_generation"] for r in scored_rows
                 if r.get(verdict_key) == "no_conflict"]
    n_pos, n_neg = len(positives), len(negatives)
    if not n_pos or not n_neg:
        # Genuinely undefined with no positives or no negatives -- stays nan/nan rather
        # than letting max(0.0, nan)/min(1.0, nan) silently collapse to a full-width
        # (0.0, 1.0) interval that looks like a real (if uninformative) measurement.
        auc = float("nan")
        confidence_interval = (float("nan"), float("nan"))
    else:
        wins = sum(1 for p in positives for n in negatives
                   if p > n) + 0.5 * sum(1 for p in positives for n in negatives if p == n)
        auc = wins / (n_pos * n_neg)
        # Wilson score interval, n = n_pos: the 8 positives are the limiting sample
        # size (n_neg ~= 3,000 contributes negligible extra uncertainty here), and a
        # naive Wald interval (z*sqrt(p(1-p)/n)) collapses to EXACTLY ZERO width at
        # p in {0, 1} -- which 8 positives can reach by chance alone -- falsely
        # reporting perfect certainty from the smallest possible sample. Wilson does
        # not collapse at the boundary.
        z = 1.96
        n = n_pos
        denominator = 1 + z**2 / n
        centre = (auc + z**2 / (2 * n)) / denominator
        half_width = z * math.sqrt(auc * (1 - auc) / n + z**2 / (4 * n**2)) / denominator
        confidence_interval = (max(0.0, centre - half_width), min(1.0, centre + half_width))
    return {
        "auc": auc,
        "n_positive": n_pos,
        "n_negative": n_neg,
        "confidence_interval": confidence_interval,
        "underpowered": n_pos < RPI_RANK_MIN_POSITIVES,
        "interpretation": (
            "Sanity check, not validation. The primary evidence is the retrospective "
            "reproduction of the published in vitro experiments (spec 11.1). Report the "
            "confidence interval, never a bare pass/fail."
        ),
    }


# Threshold below which rpi_e_total counts as "near zero" for the spec 11.3 cross-layer
# check. Defaults to eco_floor (0.01), and that is a reasoned choice, not a round
# number: E = h_host * g_geo * t_time with each term in (0, 1], so ANY one term sitting
# at its floor drags rpi_e_total to <= eco_floor regardless of the other two. "At least
# one ecological term is floored" is precisely the continuous analogue of gate3's
# categorical "no ecological route", which makes it the right comparison point.
RPI_E_NEAR_ZERO_THRESHOLD = RPI_WEIGHTS["eco_floor"]


def veto_vs_continuous_e(
    scored_rows: Sequence[Dict[str, Any]],
    near_zero_threshold: float = RPI_E_NEAR_ZERO_THRESHOLD,
    gate3_key: str = "gate3_status",
    e_key: str = "rpi_e_total",
    id_key: str = "candidate_id",
) -> Dict[str, Any]:
    """Spec 11.3 cross-layer check: hard gate3 veto vs near-zero continuous Factor E.

    Spec 10.2 is why the RPI is computed for vetoed candidates at all -- without this
    comparison that delivered behaviour has no consumer. Divergence between the two
    layers flags a defect in ONE OF THEM and is informative either way; it is not
    automatically a bug in the RPI.

    Disagreement is EXPECTED, per spec 10.2: gate3_check reads implicated_species[:2]
    and has no temporal dimension, while Factor E covers every implicated species and
    includes t_time. The counts below are the measurement, not a pass/fail.

    Rows lacking either key (e.g. junk_sequence rows, which never reach compute_rpi)
    are counted as `skipped` rather than silently scored as agreement.

    Returns agreement counts plus the ids of the disagreeing candidates in both
    directions.
    """
    agree_veto_near_zero = agree_route_nonzero = 0
    veto_but_e_not_near_zero: List[Any] = []
    near_zero_e_but_no_veto: List[Any] = []
    skipped: List[Any] = []

    for row in scored_rows:
        gate3_status = row.get(gate3_key)
        e_total = row.get(e_key)
        if gate3_status is None or e_total is None:
            skipped.append(row.get(id_key))
            continue
        vetoed = gate3_status == "no_ecological_route"
        near_zero = float(e_total) <= near_zero_threshold
        if vetoed and near_zero:
            agree_veto_near_zero += 1
        elif not vetoed and not near_zero:
            agree_route_nonzero += 1
        elif vetoed:
            veto_but_e_not_near_zero.append(row.get(id_key))
        else:
            near_zero_e_but_no_veto.append(row.get(id_key))

    n_agree = agree_veto_near_zero + agree_route_nonzero
    n_disagree = len(veto_but_e_not_near_zero) + len(near_zero_e_but_no_veto)
    n_compared = n_agree + n_disagree
    return {
        "near_zero_threshold": near_zero_threshold,
        "n_compared": n_compared,
        "n_skipped": len(skipped),
        "skipped_candidates": skipped,
        "n_agree": n_agree,
        "n_disagree": n_disagree,
        "agree_veto_and_near_zero_e": agree_veto_near_zero,
        "agree_route_and_nonzero_e": agree_route_nonzero,
        "disagree_veto_but_e_not_near_zero": veto_but_e_not_near_zero,
        "disagree_near_zero_e_but_no_veto": near_zero_e_but_no_veto,
        "agreement_fraction": (n_agree / n_compared) if n_compared else float("nan"),
        "interpretation": (
            "Spec 11.3 / 10.2. gate3_check reads implicated_species[:2] and has no "
            "temporal term; Factor E covers all implicated species and includes "
            "t_time. Disagreement is expected and informative -- it is the point of "
            "the check, not a bug. Investigate the named candidates in both "
            "directions; it flags a defect in one layer or the other, and which one "
            "is a judgement this function does not make."
        ),
    }


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

    result = {
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

    # RPI is computed for EVERY candidate that reaches this point, including ones gate3
    # vetoed -- deliberate: it is what makes the veto-vs-continuous-E cross-check in
    # Task 8 possible. It is assembled from the gate results above but never feeds BACK
    # into verdict/combine_verdict -- see VerdictUnchangedRegressionTest.
    #
    # M_PACK and M_POL (D3, Task 10/11, WS3.3/WS3.4) BOTH now score for real, and both
    # share the SAME odd_segment default and the SAME donor/backbone species assignment
    # (rpi_species_assignment) -- see score_m_pack's and score_m_pol's docstrings for
    # the full method, the status vocabularies, and why the donor==backbone (concordant)
    # case is the correct neutral answer rather than the self-comparison bug both used
    # to have (M_pol via `polymerase_compatibility(terminal_windows, terminal_windows)`;
    # M_pack via the hardcoded 1.0).
    #
    # WS3.1's calibration measured M_pack's ddG as REAL but SMALL signal (Cliff's
    # delta = +0.222; 38% of the calibration panel has ddG <= 0, where M_pack saturates
    # at exactly 1.0). User-acknowledged 2026-08-29: proceed as scoped, let M_pol carry
    # more of the b2 cross-species discrimination weight. See
    # data/rpi_calibration.json -> "discrimination" for the measurement this call site
    # now depends on.
    #
    # m_pol_status is now whatever polymerase_compatibility (or "no_parental_sequence")
    # actually reports, no longer the forced "fixed_constant" lie. It will NOT report
    # "fixed_constant" itself -- that comes only from m_pol_variance_check, a RUN-LEVEL
    # check over a representative sample of lineage pairs (spec 6.2), not a
    # per-candidate one; wiring that in is future work, not Task 11's stated scope. See
    # score_m_pol's docstring.
    m_pack_value, m_pack_status = score_m_pack(
        candidate, reference, gate1, gate1.discordant_segment or "M")
    m_pol_value, m_pol_status = score_m_pol(
        candidate_cds, reference, gate1, gate1.discordant_segment or "M")
    # Factor E isolate-level fix (design doc 2026-08-30-rpi-e-geo-isolate-level-fix):
    # resolve ONE isolate per implicated species -- the candidate's real known parent
    # when the caller supplied one (CandidateInput.known_segment_isolates), otherwise
    # the k-mer-nearest wild reference isolate proxy. `isolate_basis` travels into the
    # row as rpi_e_isolate_basis so no row's E is ever silently on a proxy.
    parent_isolates, isolate_basis = resolve_e_parent_isolates(
        gate1.implicated_species, reference.isolates, gate1, gate0.orf_proteins,
        known_isolates=candidate.known_segment_isolates)
    result.update(compute_rpi(
        implicated_species=gate1.implicated_species,
        isolates=reference.isolates,
        parent_isolates=parent_isolates,
        isolate_basis=isolate_basis,
        # DEFAULT "M" when gate1 named no discordant segment. Not arbitrary: M is the
        # modal odd segment in both strata of the constellation prior, so it is the
        # least surprising assumption -- and it keeps compute_rpi total for candidates
        # that reach here with an unresolved segment (gate1 patterns other than a
        # single-segment discordance) rather than raising mid-run. The consequence is
        # that such a row's rpi_k_value and rpi_f are scored as an M-exchange without
        # gate1 having identified one. The substitution is visible in the output: an
        # empty gate1_discordant_segment on a row with M-valued rpi_f means the default
        # fired. Check that column before quoting rpi_retention for such a row.
        #
        # MUST MATCH the odd_segment passed to score_m_pack above -- both read
        # gate1.discordant_segment with the same "M" default, so rpi_m_pack, rpi_k_value
        # and rpi_f are always scored against the same segment identity.
        odd_segment=gate1.discordant_segment or "M",
        # SCAFFOLDING (K stratum): gate0.best_kmer_identity is PROTEIN k-mer novelty
        # against the nearest reference, maximised over segments -- NOT the parental
        # nucleotide divergence spec section 7 stratifies on, and on a different scale
        # (~0.35-0.46 here vs the boundary's 0.548 / 0.707 / 0.85 calibration points).
        # Every row is therefore marked rpi_k_stratum_provisional -- see
        # RPI_K_STRATUM_PROVISIONAL. Computing real parental identity is a documented
        # follow-up, deliberately out of scope for this wave.
        identity=gate0.best_kmer_identity,
        m_pack=m_pack_value,
        m_pack_status=m_pack_status,
        m_pol=m_pol_value,
        m_pol_status=m_pol_status,
        structural_intact=not (
            gate2a.tmis_conserved is False
            or gate2a.ph_switch_conserved is False
            or gate2a.waasa_intact is False
        ),
    ))
    return result


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
    parser.add_argument(
        "--rpi-report", action="store_true",
        help=("Also write rpi_validation_report.json (spec 11.1/11.2/11.3 and the "
              "spec 12 sensitivity sweep) into the output directory. Off by default: "
              "existing runs are byte-identical without it."))
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


# Filename for the --rpi-report artifact. Named here rather than inline so tests and
# any downstream consumer agree on it.
RPI_REPORT_FILENAME = "rpi_validation_report.json"


def _json_safe(value: Any) -> Any:
    """Replace non-finite floats with None, recursively, so the result is valid JSON.

    nan and +/-inf have no JSON spelling. null is the honest one for "undefined"; 0.0
    would be a claim. Dict keys are left alone -- they are strings throughout this
    report.

    NUMPY SCALARS ARE NORMALISED HERE TOO. Testing `isinstance(value, float)` alone was
    not enough: np.float32 is NOT a float subclass (only np.float64 is), so every other
    numpy scalar fell through, and json's `default=str` then intercepted it before
    `allow_nan=False` could fire. That mis-serialised in two ways, not one --
    np.float32('nan') became the JSON *string* "nan" (a value, where null was meant),
    and FINITE np.float32(1.5) / np.int64(3) became "1.5" / "3", handing a downstream
    reader strings where the schema promises numbers.
    """
    # bool FIRST: it is an int subclass, so any numeric check below would otherwise
    # normalise True to 1 and change the report's schema
    # (rpi_k_stratum_provisional is a bool in this very report).
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        value = float(value)
        return value if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def write_rpi_report(scored_rows: Sequence[Dict[str, Any]], output_dir: Path) -> Path:
    """Write the RPI validation report to `output_dir` and return its path.

    Collects the validation entry points that previously had no caller outside the
    tests: retrospective reproduction (spec 11.1) including one leave-one-pair-out run
    per contributing pair, the rank-consistency sanity check (11.2), the veto-vs-
    continuous-E cross-layer check (11.3), and the sensitivity sweep (12).

    Written as JSON rather than CSV because these results are nested and heterogeneous;
    flattening them into the candidate table would mix per-candidate observations with
    whole-run diagnostics.
    """
    report: Dict[str, Any] = {
        "rpi_weights_version": RPI_WEIGHTS_VERSION,
        "rpi_k_stratum_provisional": RPI_K_STRATUM_PROVISIONAL,
        "retrospective_reproduction": retrospective_reproduction(),
        "leave_one_pair_out": {
            name: retrospective_reproduction(holdout=name)
            for name, record in RPI_INVITRO_BENCHMARK.items()
            if record.get("contributes_to_prior")
        },
        "rank_consistency": rank_consistency(scored_rows),
        "veto_vs_continuous_e": veto_vs_continuous_e(scored_rows),
        "sensitivity_sweep": sensitivity_sweep(),
    }
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / RPI_REPORT_FILENAME
    with path.open("w", encoding="utf-8") as handle:
        # allow_nan=False, with non-finite values mapped to null FIRST. Python's json
        # emits bare `NaN`/`Infinity` tokens by default, which are not JSON: jq and
        # JavaScript's JSON.parse both reject the file outright. Python's own json.loads
        # accepts them, which is exactly why the round-trip test did not notice.
        #
        # This is a routine path, not an edge case: rank_consistency DELIBERATELY
        # returns nan for auc and its confidence interval whenever a run has no
        # positives or no negatives -- true of the documented single-candidate
        # invocation -- and _worst_rank_correlation and _composition_distance are
        # deliberately nan-returning too. null is the right JSON spelling of "this
        # statistic is undefined"; it must not become 0.0, which would be a value.
        json.dump(_json_safe(report), handle, indent=2, sort_keys=True,
                  default=str, allow_nan=False)
    return path


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

    if args.rpi_report:
        report_path = write_rpi_report(rows, output_dir)
        print(f"RPI validation report: {report_path}")

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
import unittest.mock


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

    def test_orf_running_off_the_end_is_kept_not_discarded(self):
        """CHANGED BEHAVIOUR, deliberately. This used to assert "" -- an ORF with no
        terminal stop was thrown away.

        That is wrong for real data: GenBank records are routinely trimmed to (or just
        inside) the CDS, so the stop codon is simply absent from the deposit. Tatenale
        virus Upton Heath M (NC_055637.1) is a RefSeq record annotated coding-complete
        with an 1118 aa glycoprotein, deposited at 3355 nt -- 5 nt short of the 3360 it
        would need to carry its own stop. The real ORF was discarded, a spurious 62 aa
        ORF elsewhere won, and gate1 then reported "too_divergent" for an ordinary
        unmodified isolate."""
        self.assertEqual(longest_orf_protein("ATG" + "GCT" * 5), "A" * 5)

    def test_trailing_partial_codon_is_not_translated(self):
        """A dangling 1-2 nt tail is not a codon and must not extend the protein."""
        self.assertEqual(longest_orf_protein("ATG" + "GCT" * 5 + "GC"), "A" * 5)

    def test_terminus_is_reported_so_callers_can_tell_the_two_cases_apart(self):
        """A stop-terminated ORF must still be reported as terminated (and its stop
        stripped from the protein); an open-ended one must be reported as open. Callers
        that need a real CDS terminus can then check, rather than guessing from length."""
        span, stopped = _longest_orf_span_with_terminus("ATG" + "GCT" * 7 + "TAA")
        self.assertTrue(stopped)
        self.assertTrue(span.endswith("TAA"))
        span, stopped = _longest_orf_span_with_terminus("ATG" + "GCT" * 7)
        self.assertFalse(stopped)
        self.assertEqual(len(span) % 3, 0)

    def test_open_ended_orfs_do_not_soften_junk_detection(self):
        """The guard that matters: an open-ended span is only LONG when a frame runs
        stop-free, which random sequence does not do. A stop appears roughly every 21
        codons by chance."""
        import random as _random
        rng = _random.Random(1)
        junk = "".join(rng.choice("ACGT") for _ in range(3600))
        self.assertLess(len(longest_orf_protein(junk)), 200)

    def test_empty_sequence_returns_empty(self):
        self.assertEqual(longest_orf_protein(""), "")

    def test_no_atg_at_all_returns_empty(self):
        self.assertEqual(longest_orf_protein("CCCCCCCCCCCC"), "")


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

    def _genomic_of_length(self, total_nt, aa_len):
        """A segment carrying a real `aa_len` ORF, PADDED to exactly `total_nt`.

        Padded, never clipped: clipping would cut the stop codon off and the segment
        would fail the ORF check instead of the length check under test.
        """
        body = self._valid_orf_genomic(aa_len, "TAGTAGTAGAC", "GTCTACTACTA")
        assert len(body) <= total_nt, (
            f"a {aa_len} aa ORF needs {len(body)} nt, more than the requested {total_nt}")
        return body + "A" * (total_nt - len(body))

    def test_marginally_unusual_genomic_length_warns_instead_of_calling_it_junk(self):
        """Regression: two REAL NCBI RefSeq complete genomes were being discarded as
        junk_sequence for sitting just outside ICTV's segment-identification band --
        Luxi virus LX309 S (NC_038530.1) at 2107 nt against 1200-2100, over by SEVEN
        nucleotides, and Tatenale virus Upton Heath M (NC_055637.1) at 3355 nt against
        3400-3800. junk_sequence is terminal, so a divergent novel hantavirus was
        dropped unscored for being marginally off a band derived from its better-sampled
        relatives."""
        candidate = CandidateInput(
            candidate_id="MARGINAL",
            s_genomic=self._genomic_of_length(2107, 430),
            m_genomic=self._genomic_of_length(3355, 1091),
            l_genomic=self._genomic_of_length(6500, 2150),
        )
        result = gate0_check(candidate, self.protein_len_ranges, self.panhandle,
                             self.reference_profiles)
        self.assertTrue(result.passed_structure, result.failures)
        self.assertTrue(any("within tolerance" in w for w in result.warnings),
                        result.warnings)

    def test_genuinely_truncated_segment_is_still_junk(self):
        """The tolerance band must not swallow real truncation. This is the partial
        PUUV record shape that motivated keeping the check at all: a 1014 nt M is 70%
        below the band. Its longest ORF (43 aa) sits INSIDE the empirical protein range,
        so the protein check alone would pass it -- genomic length is what catches it."""
        candidate = CandidateInput(
            candidate_id="TRUNCATED",
            s_genomic=self._genomic_of_length(1302, 420),
            m_genomic=self._genomic_of_length(1014, 40),
            l_genomic=self._genomic_of_length(665, 40),
        )
        result = gate0_check(candidate, self.protein_len_ranges, self.panhandle,
                             self.reference_profiles)
        self.assertFalse(result.passed_structure)
        self.assertTrue(any(f.startswith("M: genomic length 1014") for f in result.failures),
                        result.failures)

    def test_absurdly_long_segment_is_rejected_as_a_concatenation(self):
        """all_hanta.py's SEGMENT_LENGTH_ABSURD_MAX stays a hard ceiling -- beyond it a
        record 'is almost certainly a multi-segment concatenation or a mis-label'."""
        candidate = CandidateInput(
            candidate_id="CONCATENATED",
            s_genomic=self._genomic_of_length(4000, 430),
            m_genomic=self._genomic_of_length(3600, 1130),
            l_genomic=self._genomic_of_length(6500, 2150),
        )
        result = gate0_check(candidate, self.protein_len_ranges, self.panhandle,
                             self.reference_profiles)
        self.assertFalse(result.passed_structure)
        self.assertTrue(any("plausible maximum" in f for f in result.failures),
                        result.failures)


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

    def test_implicated_species_come_back_in_a_stable_order(self):
        """implicated_species was list({...}) -- a SET. Python randomises string hashing
        per process, so the order varied run to run, and gate3_check reads
        implicated_species[:2]: which two species gate 3 examined was chosen by the hash
        seed. Measured on three real species, 2 of the 6 orderings returned
        'ecological_overlap_plausible' where the others returned 'no_ecological_route'.

        Sorting fixes WHICH pair is examined. It does NOT address the deeper issue that
        gate 3 only ever examines 2 of N -- see the KNOWN LIMITATION note on
        gate3_check, deferred to the Factor E rebuild."""
        reference = {
            "S": [("SpA", "M" * 50), ("SpB", "W" * 50), ("SpC", "Y" * 50)],
            "M": [("SpA", "M" * 50), ("SpB", "W" * 50), ("SpC", "Y" * 50)],
            "L": [("SpA", "M" * 50), ("SpB", "W" * 50), ("SpC", "Y" * 50)],
        }
        result = gate1_check({"S": "M" * 50, "M": "W" * 50, "L": "Y" * 50}, reference)
        # Guard the guard: a fixture that stopped implicating 3 species would make the
        # ordering assertion below trivially true.
        self.assertEqual(len(result.implicated_species), 3, result.implicated_species)
        self.assertEqual(result.implicated_species,
                         sorted(result.implicated_species))

    def test_implicated_species_order_survives_a_changed_hash_seed(self):
        """The assertion above cannot observe the actual defect: PYTHONHASHSEED is fixed
        for the life of a process, so an in-process test sees one arbitrary order and
        cannot tell a sorted list from a luckily-ordered one. This runs the same call in
        subprocesses under different seeds, which is the only way to see it."""
        source = (
            "import json, reassortment_detector as rd\n"
            "ref = {seg: [('SpA','M'*50), ('SpB','W'*50), ('SpC','Y'*50)]\n"
            "       for seg in ('S','M','L')}\n"
            "r = rd.gate1_check({'S':'M'*50,'M':'W'*50,'L':'Y'*50}, ref)\n"
            "print(json.dumps(r.implicated_species))\n"
        )
        seen = set()
        for seed in ("0", "1", "12345"):
            env = dict(os.environ, PYTHONHASHSEED=seed)
            completed = subprocess.run(
                [sys.executable, "-c", source], capture_output=True, text=True,
                cwd=str(Path(__file__).parent), env=env, check=True)
            seen.add(completed.stdout.strip())
        self.assertEqual(len(seen), 1,
                         f"gate1 implicated_species varies with PYTHONHASHSEED: {seen}")


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


def _write_fake_binary(directory: Path, name: str, source: str) -> Path:
    """Write a fake mafft/iqtree stand-in and return a path the OS can actually exec.

    The scripts are Python, and phylogeny.py invokes them through subprocess.run with
    the script path as argv[0]. On POSIX the shebang plus the exec bit is enough. On
    WINDOWS IT IS NOT: CreateProcess does not consult shebangs, so handing it a .py
    raises `OSError: [WinError 193] %1 is not a valid Win32 application`, phylogeny.py
    catches it, and every tree comes back phylo_status="failed".

    That is a defect in this harness, not in the pipeline, and it was invisible for as
    long as these tests died earlier on `import dendropy`. Fixing the import exposed it.
    On Windows we therefore emit a .cmd shim that re-launches the script under the
    CURRENT interpreter (sys.executable, so the shim cannot pick up a different Python
    than the one running the suite) and forwards its arguments.
    """
    script = directory / name
    script.write_text(source, encoding="utf-8")
    script.chmod(0o755)
    if os.name != "nt":
        return script
    shim = script.with_suffix(".cmd")
    shim.write_text(f'@echo off\r\n"{sys.executable}" "{script}" %*\r\n', encoding="utf-8")
    return shim


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
        self.mafft_bin = _write_fake_binary(tmp_path, "fake_mafft.py", self._FAKE_MAFFT)
        self.iqtree_bin = _write_fake_binary(tmp_path, "fake_iqtree.py", self._FAKE_IQTREE)
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
        self.mafft_bin = _write_fake_binary(
            tmp_path, "fake_mafft.py", Gate4CheckEndToEndTest._FAKE_MAFFT)
        self.iqtree_bin = _write_fake_binary(
            tmp_path, "fake_iqtree.py", Gate4CheckEndToEndTest._FAKE_IQTREE)
        # Reuses the exact fake-binary scripts from Gate4CheckEndToEndTest above -- see
        # that class for the source. Both classes live in this same module now, so no
        # cross-module import is needed to reach it.
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


class ConstellationPriorTest(unittest.TestCase):
    def test_no_category_is_zero_after_smoothing(self):
        """Rizvanov's S-odd=0 and L-odd=0 are unobserved events in n=10, not measured
        impossibilities. A literal zero inside a multiplicative index is a hidden hard
        gate -- the same cliff logic rejected for gate3's veto."""
        prior = constellation_prior("distant")
        for segment in ("S", "M", "L"):
            self.assertGreater(prior[segment], 0.0)

    def test_priors_sum_to_one(self):
        for stratum in ("close", "distant"):
            self.assertAlmostEqual(sum(constellation_prior(stratum).values()), 1.0, places=9)

    def test_distant_stratum_favours_m_exchange(self):
        prior = constellation_prior("distant")
        self.assertGreater(prior["M"], prior["S"])
        self.assertGreater(prior["M"], prior["L"])

    def test_close_stratum_m_and_s_are_equal(self):
        """Rodriguez 1998 observed 12 M-odd and 12 S-odd. The close stratum must NOT
        inherit the M-only rule -- Klempa's own Table 2 reproduces all three genotypes."""
        prior = constellation_prior("close")
        self.assertAlmostEqual(prior["M"], prior["S"], places=9)

    def test_close_stratum_l_exchange_is_rare_but_present(self):
        prior = constellation_prior("close")
        self.assertLess(prior["L"], prior["M"])
        self.assertGreater(prior["L"], 0.0)

    def test_smoothing_constant_is_configurable_for_sensitivity_sweep(self):
        strong = constellation_prior("distant", alpha=5.0)
        weak = constellation_prior("distant", alpha=0.5)
        self.assertGreater(strong["S"], weak["S"])

    def test_unknown_stratum_raises(self):
        with self.assertRaises(KeyError):
            constellation_prior("intermediate")


class DivergenceStratumTest(unittest.TestCase):
    def test_rodriguez_within_species_pair_is_close(self):
        """SNV NMR11 x CC107: ~15% N-gene nt divergence (Rodriguez 1998 Table 7)."""
        self.assertEqual(divergence_stratum(0.85), "close")

    def test_rizvanov_cross_species_pair_is_distant(self):
        """SNV x ANDV: 54.8% S-segment nt similarity (Rizvanov 2004)."""
        self.assertEqual(divergence_stratum(0.548), "distant")

    def test_snv_bccv_pair_is_distant(self):
        """70.7% S similarity per Rizvanov; ~34% N divergence per Rodriguez Table 7.
        The two sources disagree on this pair; both land in 'distant'."""
        self.assertEqual(divergence_stratum(0.707), "distant")
        self.assertEqual(divergence_stratum(0.66), "distant")

    def test_boundary_is_configurable_for_sensitivity_sweep(self):
        self.assertEqual(divergence_stratum(0.75, boundary=0.70), "close")
        self.assertEqual(divergence_stratum(0.75, boundary=0.80), "distant")

    def test_identity_outside_unit_interval_raises(self):
        with self.assertRaises(ValueError):
            divergence_stratum(1.5)
        with self.assertRaises(ValueError):
            divergence_stratum(-0.1)


def _parent_map(*isolates: "ah.IsolateRecord") -> Dict[str, "ah.IsolateRecord"]:
    """species -> isolate, for tests that already have exactly one isolate per
    implicated species and want to hand it to ecological_opportunity as its real,
    known parent (basis "actual_parents")."""
    return {isolate.species_name: isolate for isolate in isolates}


class EcologicalOpportunityTest(unittest.TestCase):
    """FIX (2026-08-30): ecological_opportunity now takes ONE isolate per implicated
    species (parent_isolates) instead of fanning out to every isolate of that species in
    the reference -- see design doc 2026-08-30-rpi-e-geo-isolate-level-fix-design.md and
    notes_diagnosis.md for the bug this replaces (a widely-distributed species' full
    country union defeated a "different continents" construction guarantee). Every test
    below supplies parent_isolates explicitly, matching how a real caller
    (resolve_e_parent_isolates) would resolve it."""

    def test_same_host_genus_same_country_same_years_scores_high(self):
        a = _eco_isolate("A", "USA", "2001")
        b = _eco_isolate("B", "USA", "2002")
        result = ecological_opportunity(
            ["A", "B"], [a, b], _parent_map(a, b), "actual_parents")
        self.assertGreater(result["rpi_e_total"], 0.4)
        self.assertEqual(result["rpi_e_geo"], 1.0)
        self.assertEqual(result["rpi_e_time"], 1.0)
        self.assertEqual(result["rpi_e_isolate_basis"], "actual_parents")

    def test_disjoint_countries_drive_geo_term_down(self):
        a = _eco_isolate("A", "USA", "2001")
        b = _eco_isolate("B", "Chile", "2002")
        result = ecological_opportunity(
            ["A", "B"], [a, b], _parent_map(a, b), "actual_parents")
        self.assertLess(result["rpi_e_geo"], 1.0)

    def test_same_species_country_pool_but_different_parent_isolates_still_separates(self):
        """THE FIX ITSELF: a species with reference isolates in MANY countries (Puumala/
        Seoul-style) must not register geo opportunity against another species just
        because the two species' country sets overlap SOMEWHERE in the corpus -- only
        the two ACTUAL parent isolates' own countries matter. Species A has isolates in
        both USA and China; species B has isolates in both Chile and China. The specific
        parent pair scored here (USA x Chile) never co-occurs, even though 'China' is a
        country both species are seen in elsewhere in the reference."""
        a_near = _eco_isolate("A", "China", "2001")
        a_far = _eco_isolate("A", "USA", "2001")
        b_near = _eco_isolate("B", "China", "2001")
        b_far = _eco_isolate("B", "Chile", "2001")
        isolates = [a_near, a_far, b_near, b_far]
        # The OLD species-wide-union behaviour would find China in both species' sets
        # and score geo = 1.0 no matter which isolate pair is actually implicated.
        result = ecological_opportunity(
            ["A", "B"], isolates, _parent_map(a_far, b_far), "actual_parents")
        self.assertLess(result["rpi_e_geo"], 1.0)
        # The near pair (both China) DOES have a real route -- same fixture pool,
        # different parent pair, opposite answer. Confirms the function reads
        # `parent_isolates`, not incidental structure of the wider isolates list.
        near_result = ecological_opportunity(
            ["A", "B"], isolates, _parent_map(a_near, b_near), "actual_parents")
        self.assertEqual(near_result["rpi_e_geo"], 1.0)

    def test_disjoint_decades_drive_time_term_down(self):
        a = _eco_isolate("A", "USA", "1985")
        b = _eco_isolate("B", "USA", "2020")
        result = ecological_opportunity(
            ["A", "B"], [a, b], _parent_map(a, b), "actual_parents")
        self.assertLess(result["rpi_e_time"], 1.0)

    def test_no_factor_is_ever_exactly_zero(self):
        """A zero in a multiplicative index is a hidden hard gate. Ecological
        impossibility must approach zero, never reach it."""
        a = _eco_isolate("A", "USA", "1985", host_order="Rodentia")
        b = _eco_isolate("B", "Chile", "2020", host_order="Soricomorpha",
                         host_family="Soricidae", host_name="Sorex araneus")
        result = ecological_opportunity(
            ["A", "B"], [a, b], _parent_map(a, b), "actual_parents")
        for key in ("rpi_e_host", "rpi_e_geo", "rpi_e_time", "rpi_e_total"):
            self.assertGreater(result[key], 0.0, key)
        self.assertLess(result["rpi_e_total"], 0.01)

    def test_covers_all_implicated_species_not_just_first_two(self):
        """gate3_check reads implicated_species[:2]. Factor E must not."""
        a = _eco_isolate("A", "USA", "2001")
        b = _eco_isolate("B", "USA", "2001")
        c = _eco_isolate("C", "Chile", "1980", host_order="Soricomorpha",
                         host_family="Soricidae", host_name="Sorex araneus")
        isolates = [a, b, c]
        two = ecological_opportunity(
            ["A", "B"], isolates, _parent_map(a, b), "actual_parents")
        three = ecological_opportunity(
            ["A", "B", "C"], isolates, _parent_map(a, b, c), "actual_parents")
        self.assertLess(three["rpi_e_total"], two["rpi_e_total"])

    def test_lab_and_vaccine_strains_are_excluded(self):
        """B has no WILD isolate, so a real resolver never puts it in parent_isolates --
        the floor-host branch must still fire from the unfiltered `isolates` list."""
        wild = _eco_isolate("A", "USA", "2001")
        lab = _eco_isolate("B", "USA", "2001")
        lab.is_lab_derived = True
        result = ecological_opportunity(
            ["A", "B"], [wild, lab], _parent_map(wild), "insufficient_metadata")
        self.assertLess(result["rpi_e_total"], 0.5)

    def test_species_absent_from_reference_entirely_stays_neutral(self):
        """A species with no isolate of any kind (wild, lab, or vaccine) is unknown to
        the reference -- novel or simply unsampled -- not ecologically separated.
        gate3_check already treats this as 'ecology_unknown'; Factor E must agree."""
        a = _eco_isolate("A", "USA", "2001")
        result = ecological_opportunity(
            ["A", "B"], [a], _parent_map(a), "insufficient_metadata")
        self.assertEqual(result["rpi_e_total"], 1.0)

    def test_unknown_host_metadata_stays_neutral_not_maximal_separation(self):
        """When one side's isolates carry no host taxonomy at all -- host_name,
        host_family, and host_order all empty, as happens when resolve_host_taxonomy
        fails and no biosample host is available (all_hanta.py:1560-1565, flagged
        'missing_host') -- that is unknown, not separated. It must not fall through
        to 'different_order' (0.01), the harshest value in the ladder."""
        a = _eco_isolate("A", "USA", "2001")
        b = _eco_isolate("B", "USA", "2001", host_order="", host_family="", host_name="")
        result = ecological_opportunity(
            ["A", "B"], [a, b], _parent_map(a, b), "actual_parents")
        self.assertEqual(result["rpi_e_host"], 1.0)

    def test_single_implicated_species_returns_neutral_one(self):
        a = _eco_isolate("A", "USA", "2001")
        result = ecological_opportunity(["A"], [a], _parent_map(a), "actual_parents")
        self.assertEqual(result["rpi_e_total"], 1.0)
        self.assertEqual(result["rpi_e_isolate_basis"], "actual_parents")


class NearestReferenceIsolateTest(unittest.TestCase):
    """nearest_reference_isolate is the fallback resolve_e_parent_isolates uses when a
    candidate's true parent isolate is unknown (design doc section 1a)."""

    def test_picks_the_isolate_whose_protein_is_closest(self):
        close = _make_isolate("close", "SpA", "MEQRSTVCLCLFLTGLLLWNPSESRN", "M" * 10, "M" * 10)
        far = _make_isolate("far", "SpA", "AAAAAAAAAAAAAAAAAAAAAAAAAA", "M" * 10, "M" * 10)
        found = nearest_reference_isolate(
            "MEQRSTVCLCLFLTGLLLWNPSESRN", "S", [far, close])
        self.assertEqual(found.isolate_key, "close")

    def test_returns_none_for_empty_candidate_protein(self):
        isolate = _make_isolate("i1", "SpA", "MEQ", "MEQ", "MEQ")
        self.assertIsNone(nearest_reference_isolate("", "S", [isolate]))

    def test_returns_none_when_pool_has_no_usable_protein_for_the_segment(self):
        isolate = _make_isolate("i1", "SpA", "", "MEQ", "MEQ")  # empty s_protein_aa
        self.assertIsNone(nearest_reference_isolate("MEQ", "S", [isolate]))

    def test_returns_none_for_an_empty_pool(self):
        self.assertIsNone(nearest_reference_isolate("MEQ", "S", []))


class ResolveEParentIsolatesTest(unittest.TestCase):
    """resolve_e_parent_isolates (design doc section 1a): the caller-facing function
    score_candidate uses to decide, per candidate, whether Factor E gets the real known
    parent isolates or the k-mer-nearest-reference proxy -- and to declare which, via
    rpi_e_isolate_basis, so a proxy row is always flagged, never silent."""

    @staticmethod
    def _gate1(species_by_segment: Dict[str, str]) -> "Gate1Result":
        top_hits = {seg: [SpeciesHit(species=sp, distance=0.1)]
                    for seg, sp in species_by_segment.items()}
        implicated = sorted(set(species_by_segment.values()))
        return Gate1Result(top_hits=top_hits, implicated_species=implicated)

    def test_known_isolates_covering_every_species_gives_actual_parents(self):
        donor = _make_isolate("donor", "Donor", "MEQRSTVCLCLFLTGLLLW", "M" * 10, "M" * 10)
        backbone = _make_isolate("backbone", "Backbone", "M" * 10, "MEQRSTVCLCLFLTGLLLW", "M" * 10)
        gate1 = self._gate1({"S": "Donor", "M": "Backbone", "L": "Backbone"})
        known = {"S": donor, "M": backbone, "L": backbone}
        resolved, basis = resolve_e_parent_isolates(
            gate1.implicated_species, [donor, backbone], gate1,
            {"S": "MEQRSTVCLCLFLTGLLLW", "M": "MEQRSTVCLCLFLTGLLLW", "L": "M" * 10},
            known_isolates=known)
        self.assertEqual(basis, "actual_parents")
        self.assertEqual(resolved, {"Donor": donor, "Backbone": backbone})

    def test_known_isolates_missing_a_species_falls_back_for_the_whole_row(self):
        """WS4 group c's own real isolate is known, but the discordant segment's
        inferred donor species is NOT -- resolve_e_parent_isolates must not silently
        mix "real" and "proxy" isolates within one row; the whole row uses the
        fallback."""
        own = _make_isolate("own", "OwnSpecies", "MEQRSTVCLCLFLTGLLLW", "M" * 10, "M" * 10)
        donor_ref = _make_isolate(
            "donor_ref", "DonorSpecies", "MEQRSTVCLCLFLTGLLLW", "M" * 10, "M" * 10)
        gate1 = self._gate1({"S": "DonorSpecies", "M": "OwnSpecies", "L": "OwnSpecies"})
        known = {"S": own, "M": own, "L": own}  # the real isolate, on every segment
        resolved, basis = resolve_e_parent_isolates(
            gate1.implicated_species, [own, donor_ref], gate1,
            {"S": "MEQRSTVCLCLFLTGLLLW", "M": "M" * 10, "L": "M" * 10},
            known_isolates=known)
        self.assertEqual(basis, "nearest_reference")
        self.assertEqual(set(resolved), {"DonorSpecies", "OwnSpecies"})

    def test_no_known_isolates_falls_back_to_nearest_reference_or_insufficient(self):
        ref = _make_isolate("ref", "SpA", "MEQRSTVCLCLFLTGLLLW", "M" * 10, "M" * 10)
        # SpB has no isolate anywhere in the reference pool passed in.
        gate1 = self._gate1({"S": "SpA", "M": "SpB", "L": "SpB"})
        resolved, basis = resolve_e_parent_isolates(
            gate1.implicated_species, [ref], gate1,
            {"S": "MEQRSTVCLCLFLTGLLLW", "M": "M" * 10, "L": "M" * 10},
            known_isolates=None)
        self.assertEqual(basis, "insufficient_metadata")
        self.assertEqual(set(resolved), {"SpA"})

    def test_proxy_basis_is_flagged_on_the_shipped_row_not_silent(self):
        """The whole point of rpi_e_isolate_basis: a row scored on the proxy must SAY
        so on the row that ships, never present itself as if it used real known
        parents. Asserted on compute_rpi's returned dict, not on an intermediate."""
        ref_a = _make_isolate("ref_a", "SpA", "MEQRSTVCLCLFLTGLLLW", "M" * 10, "M" * 10)
        ref_b = _make_isolate("ref_b", "SpB", "M" * 10, "MEQRSTVCLCLFLTGLLLW", "M" * 10)
        gate1 = self._gate1({"S": "SpA", "M": "SpB", "L": "SpB"})
        resolved, basis = resolve_e_parent_isolates(
            gate1.implicated_species, [ref_a, ref_b], gate1,
            {"S": "MEQRSTVCLCLFLTGLLLW", "M": "MEQRSTVCLCLFLTGLLLW", "L": "M" * 10},
            known_isolates=None)
        self.assertEqual(basis, "nearest_reference")
        row = compute_rpi(
            implicated_species=gate1.implicated_species, isolates=[ref_a, ref_b],
            parent_isolates=resolved, isolate_basis=basis,
            odd_segment="M", identity=0.5, m_pack=1.0, m_pack_status="scored",
            m_pol=1.0, m_pol_status="scored", structural_intact=True)
        self.assertEqual(row["rpi_e_isolate_basis"], "nearest_reference")
        self.assertNotEqual(row["rpi_e_isolate_basis"], "actual_parents")


class BoltzmannCompatibilityTest(unittest.TestCase):
    def test_identical_energies_score_one(self):
        self.assertAlmostEqual(boltzmann_compatibility(-10.0, -10.0), 1.0, places=9)

    def test_weaker_non_cognate_scores_below_one(self):
        # A weaker duplex has a LESS NEGATIVE dG.
        self.assertLess(boltzmann_compatibility(-5.0, -10.0), 1.0)

    def test_stronger_non_cognate_is_clamped_to_one(self):
        """Non-cognate cannot beat cognate in this index -- McDonald's competition
        argument is that cognate wins; clamping keeps the factor a penalty, never a
        bonus."""
        self.assertEqual(boltzmann_compatibility(-20.0, -10.0), 1.0)

    def test_positive_dg_is_well_defined_and_small(self):
        """The old quotient form was meaningless for an unfavourable duplex."""
        value = boltzmann_compatibility(3.0, -10.0)
        self.assertGreater(value, 0.0)
        self.assertLess(value, 1e-5)

    def test_sign_convention_is_not_flipped(self):
        """Guards the single most damaging possible bug in this index: reversing ddG
        would score INCOMPATIBLE constellations as maximally plausible. Spec section 6.1."""
        weaker_non_cognate = boltzmann_compatibility(-5.0, -10.0)
        stronger_non_cognate = boltzmann_compatibility(-10.0, -5.0)
        self.assertLess(weaker_non_cognate, stronger_non_cognate)

    def test_rt_constant_matches_rnacofold_units(self):
        """kcal/mol at RNAcofold's default 37 C = 310.15 K."""
        self.assertAlmostEqual(RPI_RT_KCAL_PER_MOL, 0.6163, places=3)

    def test_never_returns_zero(self):
        self.assertGreater(boltzmann_compatibility(100.0, -10.0), 0.0)

    def test_extreme_ddg_clamps_instead_of_overflowing(self):
        """A far-stronger non-cognate duplex (ddG very negative) must clamp to 1.0, not
        raise OverflowError from math.exp -- the exponent has to be clamped BEFORE
        exponentiating, not just the result afterward."""
        self.assertEqual(boltzmann_compatibility(-1000.0, 0.0), 1.0)

    def test_extreme_positive_ddg_underflows_to_a_floor_not_to_zero(self):
        """The other direction, and it is the one that bit. math.exp does NOT raise on
        underflow the way it raises on overflow -- it silently returns 0.0, so a large
        positive ddG produced exactly the hidden hard gate the docstring promises never
        to produce: boltzmann_compatibility(500.0, 0.0) was 0.0.

        Kept as an adversarial-input guard rather than an expected-input one. The
        original rationale -- whole-segment duplexes with MFEs in the hundreds -- was
        retired by D2: the scored window is the 50-nt panhandle, where ddG runs ~0-22.
        test_never_returns_zero only probed ddG = 110."""
        for ddg in (500.0, 1000.0, 1e6):
            self.assertGreater(boltzmann_compatibility(ddg, 0.0), 0.0, ddg)
            self.assertLess(boltzmann_compatibility(ddg, 0.0), 1.0, ddg)

    def test_non_finite_ddg_raises_rather_than_silently_scoring(self):
        """max(-700.0, nan) returns -700.0, so a nan ddG silently became MAXIMAL penalty
        (9.86e-305). Before the underflow clamp, min(0.0, nan) returned 0.0 and the same
        input silently became NO penalty (1.0). The two answers are opposites and both
        are wrong: a non-finite ddG means RNAcofold returned an unusable energy, and this
        module's no-silent-sentinel rule says that must surface, not propagate.

        The caller catches this and records rpi_m_pack_status = "undefined_energy"
        rather than aborting the run."""
        for non_cognate, cognate in ((float("nan"), 0.0), (0.0, float("nan")),
                                     (float("inf"), 0.0), (float("inf"), float("inf"))):
            with self.assertRaises(ValueError, msg=f"{non_cognate}/{cognate}"):
                boltzmann_compatibility(non_cognate, cognate)

    def test_the_raise_did_not_narrow_the_working_range(self):
        """Guard the guard: a too-eager finiteness check would reject the extreme-but-
        legitimate inputs the two clamp tests above rely on."""
        self.assertAlmostEqual(boltzmann_compatibility(0.0, 0.0), 1.0, places=9)
        self.assertEqual(boltzmann_compatibility(-5000.0, 0.0), 1.0)
        self.assertGreater(boltzmann_compatibility(5000.0, 0.0), 0.0)


class MPackTauRescalingTest(unittest.TestCase):
    """WS3.2 (Task 9): tau REPLACES RT at boltzmann_compatibility's `rt` argument, so
    these tests exercise that exact call shape rather than the function in isolation --
    the four properties below are what M_pack must have once Task 10 wires it in, and
    tau is the number that must preserve them, not RPI_RT_KCAL_PER_MOL.

    tau itself is NOT re-derived here. It is a measurement from data/rpi_calibration.json
    (rpi_calibration.py's ArtifactShapeTest and RescalingBranchTest cover that derivation
    with a synthetic ground truth). This class asks only: given the shipped tau, does
    M_pack behave the way spec section 6.1 / D1 Option A requires?
    """

    TAU = RPI_WEIGHTS["m_pack_tau"]

    def test_tau_is_the_measured_ws3_1_value_not_the_dead_provisional_one(self):
        """D13: tau ~= 4.86 was provisional pending WS3.1 and must never be quoted as
        settled again. If this constant drifts back toward 4.86 -- e.g. someone reverts
        to the old placeholder -- this is the test that catches it."""
        self.assertAlmostEqual(self.TAU, 3.914078624686338, places=6)
        self.assertNotAlmostEqual(self.TAU, 4.864, places=1)

    def test_bounded_in_zero_one_inclusive_of_one(self):
        """Property 1/4 (D1, Task 9): M_pack in (0, 1]. Swept across the panel's
        observed ddG range (data/rpi_calibration.json: min -15.60, max +21.90) plus
        extremes."""
        for ddg in (-50.0, -15.6, -1.0, 0.0, 1.0, 21.9, 50.0, 700.0):
            value = boltzmann_compatibility(ddg, 0.0, rt=self.TAU)
            self.assertGreater(value, 0.0, ddg)
            self.assertLessEqual(value, 1.0, ddg)

    def test_monotonic_decreasing_in_ddg(self):
        """Property 2/4: a weaker non-cognate duplex (larger ddG) must never score
        HIGHER than a stronger one. Checked over the calibration panel's actual range,
        not just two arbitrary points."""
        ddgs = [-15.6, -6.8, -3.0, 0.0, 1.6, 5.0, 11.22, 21.9]
        values = [boltzmann_compatibility(d, 0.0, rt=self.TAU) for d in ddgs]
        for earlier, later in zip(values, values[1:]):
            self.assertGreaterEqual(earlier, later)

    def test_no_exact_zeros_at_extremes(self):
        """Property 3/4: the no-literal-zeros guarantee, at tau's scale specifically.
        A far worse-than-observed ddG must still floor at the -700-exponent clamp's
        ~1e-304 rather than underflow to a literal 0.0 -- see boltzmann_compatibility's
        own docstring for why exp() underflowing silently was the bug this guards."""
        self.assertGreater(boltzmann_compatibility(1e6, 0.0, rt=self.TAU), 0.0)
        self.assertGreater(boltzmann_compatibility(-1e6, 0.0, rt=self.TAU), 0.0)

    def test_sweepable_at_plus_minus_50_percent(self):
        """Property 4/4: D1 requires tau to be trivially sweepable at +/-50% --
        tau in {1.957, 3.914, 5.871} -- for WS4 Task 15's sensitivity analysis. Checked
        here as a property of THIS tau (all three variants still bounded, monotonic,
        non-zero), not merely that the arithmetic of halving/doubling succeeds."""
        for scaled_tau in (self.TAU * 0.5, self.TAU, self.TAU * 1.5):
            values = [boltzmann_compatibility(d, 0.0, rt=scaled_tau)
                     for d in (-15.6, 0.0, 11.22, 21.9)]
            for value in values:
                self.assertGreater(value, 0.0)
                self.assertLessEqual(value, 1.0)
            for earlier, later in zip(values, values[1:]):
                self.assertGreaterEqual(earlier, later)

    def test_window_matches_the_window_tau_was_measured_under(self):
        """tau is only valid at the N it was calibrated under (D2). If
        RPI_PANHANDLE_WINDOW_NT is ever changed without re-running rpi_calibration.py,
        this is the test that should be updated alongside it -- and its failure without
        that update is the signal that someone forgot."""
        self.assertEqual(RPI_WEIGHTS["m_pack_window_nt"], RPI_PANHANDLE_WINDOW_NT)


class CofoldDgTest(unittest.TestCase):
    """Replaces RnacofoldDgParsingTest, which exercised a stdout regex against
    hand-written sample strings. Per D2b the production path uses the ViennaRNA Python
    bindings, so there is no stdout to parse and no CLI to invoke -- the regex and its
    long-standing "ASSUMED FORMAT, NOT VALIDATED" caveat are both gone rather than left
    in place as a second, divergent mechanism."""

    def test_returns_a_finite_negative_energy_for_a_complementary_duplex(self):
        dg = cofold_dg("GCGCGCGCGCAUAUAUAUAU", "AUAUAUAUAUGCGCGCGCGC")
        self.assertTrue(math.isfinite(dg))
        self.assertLess(dg, 0.0)

    def test_a_complementary_pair_binds_more_tightly_than_a_mismatched_one(self):
        """The property M_pack rests on: a cognate-like duplex must score a LOWER (more
        negative) dG than a non-cognate-like one. If this direction were ever inverted,
        the whole index would invert with it."""
        stem = "GGGGGGGGGGCCCCCCCCCC"
        self.assertLess(cofold_dg(stem, "GGGGGGGGGGCCCCCCCCCC"),
                        cofold_dg(stem, "AAAAAAAAAAAAAAAAAAAA"))

    def test_empty_input_raises_rather_than_returning_a_sentinel(self):
        with self.assertRaises(ValueError):
            cofold_dg("", "GCGC")
        with self.assertRaises(ValueError):
            cofold_dg("GCGC", "")

    def test_is_deterministic_across_repeated_calls(self):
        """WS-API requires a deterministic scorer, and the fold is the one step that
        calls into a C library -- so its stability is asserted, not assumed."""
        a, b = "GCAUGCAUGCAUGCAUGCAU", "AUGCAUGCAUGCAUGCAUGC"
        self.assertEqual(len({cofold_dg(a, b) for _ in range(5)}), 1)

    def test_provenance_records_mechanism_and_version(self):
        """D2b: the calibration artifact must record WHICH implementation produced the
        energies, because the whole tau derivation depends on it. A bare version string
        would not distinguish these from CLI-produced energies."""
        provenance = viennarna_provenance()
        self.assertEqual(provenance["mechanism"], "ViennaRNA Python bindings")
        self.assertEqual(provenance["version"], "2.7.2")


class PanhandleWindowTest(unittest.TestCase):
    def test_takes_the_terminal_n_nt_from_the_requested_end(self):
        seq = "AAAAACCCCCGGGGGTTTTT"
        self.assertEqual(panhandle_window(seq, "5", 5), "AAAAA")
        self.assertEqual(panhandle_window(seq, "3", 5), "TTTTT")

    def test_default_window_is_the_swept_constant(self):
        self.assertEqual(RPI_PANHANDLE_WINDOW_NT, 50)
        self.assertEqual(len(panhandle_window("A" * 400, "5")), 50)

    def test_window_longer_than_the_sequence_returns_the_whole_sequence(self):
        """Short segments must not raise or pad -- a padded window would fabricate
        base-pairing that is not in the data."""
        self.assertEqual(panhandle_window("ACGT", "5", 50), "ACGT")
        self.assertEqual(panhandle_window("ACGT", "3", 50), "ACGT")

    def test_rejects_an_end_that_is_not_5_or_3(self):
        with self.assertRaises(ValueError):
            panhandle_window("ACGT", "middle")


class PolymeraseCompatibilityTest(unittest.TestCase):
    TAG = "AUCAUCAUCUGAGG"  # conserved across the genus AND across S/M/L
    SNV_FOOTPRINT = 23  # Mir & Panganiban: SNV N requires the FIRST 23 nt

    def test_conserved_tag_is_excluded_from_scoring(self):
        """The terminal 14 nt are identical across S, M and L genus-wide. If M_pol
        scored them it would saturate at 1.0 and measure nothing (spec 6.2)."""
        same_tag_different_footprint = polymerase_compatibility(
            {"S": self.TAG + "GGGGGGGGG", "L": self.TAG + "GGGGGGGGG"},
            {"S": self.TAG + "CCCCCCCCC", "L": self.TAG + "CCCCCCCCC"},
            footprint=self.SNV_FOOTPRINT,
        )[0]
        self.assertLess(same_tag_different_footprint, 1.0)

    def test_identical_footprints_score_one(self):
        seqs = {"S": self.TAG + "AAAAAAAAA", "L": self.TAG + "AAAAAAAAA"}
        value, status = polymerase_compatibility(seqs, seqs, footprint=self.SNV_FOOTPRINT)
        self.assertAlmostEqual(value, 1.0, places=9)
        self.assertEqual(status, "scored")

    def test_footprint_length_follows_mir_panganiban(self):
        """SNV N requires the first 23 nt; most hantavirus N needs only the terminal 9."""
        self.assertEqual(n_recognition_footprint("Orthohantavirus sinnombreense"), 23)
        self.assertEqual(n_recognition_footprint("Orthohantavirus puumalaense"), 9)

    def test_never_returns_zero(self):
        value, _ = polymerase_compatibility(
            {"S": self.TAG + "AAAAAAAAA", "L": self.TAG + "AAAAAAAAA"},
            {"S": self.TAG + "UUUUUUUUU", "L": self.TAG + "UUUUUUUUU"},
            footprint=self.SNV_FOOTPRINT,
        )
        self.assertGreater(value, 0.0)

    def test_scoring_window_is_measured_from_the_terminus(self):
        """Mir & Panganiban's footprints are measured FROM THE TERMINUS -- "the first
        23 nt", "the terminal 9". The scoring window is therefore the part of
        sequence[0:footprint] that is not the invariant tag, i.e. sequence[14:footprint].

        The previous slice was sequence[14:14+footprint], which for the 9-nt default
        scored nt 15-23 -- entirely OUTSIDE a footprint that ends at nt 9, with no
        source behind it."""
        sequence = self.TAG + "ABCDEFGHI" + "ZZZZZZZZZZ"
        self.assertEqual(_beyond_tag(sequence, 23), "ABCDEFGHI")
        self.assertEqual(_beyond_tag(sequence, 18), "ABCD")

    def test_default_footprint_leaves_no_window_beyond_the_tag(self):
        """The 9-nt default footprint lies ENTIRELY inside the invariant 14-nt tag, so
        there is nothing left to score. That must be reported as its own status, not as
        a scored 1.0: a neutral value from an EMPTY comparison is not the same claim as
        a neutral value from a perfect match. Only the SNV 23-vs-9 asymmetry -- the part
        that actually varies between lineages -- can open a window at all."""
        self.assertEqual(_beyond_tag(self.TAG + "GGGGGGGGG", RPI_N_FOOTPRINT_DEFAULT_NT), "")
        value, status = polymerase_compatibility(
            {"S": self.TAG + "GGGGGGGGG", "L": self.TAG + "GGGGGGGGG"},
            {"S": self.TAG + "CCCCCCCCC", "L": self.TAG + "CCCCCCCCC"},
        )
        self.assertEqual(status, "no_signal_beyond_tag")
        self.assertEqual(value, 1.0)

    def test_empty_input_also_reports_no_signal(self):
        value, status = polymerase_compatibility({}, {}, footprint=self.SNV_FOOTPRINT)
        self.assertEqual(status, "no_signal_beyond_tag")
        self.assertEqual(value, 1.0)

    def test_variance_check_fixes_constant_when_no_signal(self):
        """Spec 6.2: if M_pol shows negligible variance across real pairs, fix it at 1.0
        and SAY SO, rather than pretend it measures something."""
        identical = [(0.5, 0.5)] * 20
        fixed, status = m_pol_variance_check(identical)
        self.assertEqual(status, "fixed_constant")
        self.assertEqual(fixed, 1.0)

    def test_variance_check_keeps_scoring_when_signal_present(self):
        varied = [(v / 20.0, 0.5) for v in range(20)]
        _, status = m_pol_variance_check(varied)
        self.assertEqual(status, "scored")


class HasGenomicPanhandleTerminiTest(unittest.TestCase):
    """D3/WS3.1: the filter species_consensus_panhandle_window relies on before
    building a consensus."""

    def test_accepts_a_sequence_carrying_both_conserved_tags(self):
        self.assertTrue(has_genomic_panhandle_termini(
            HANTAVIRUS_GENOMIC_5P_TAG + "ACGT" * 10 + HANTAVIRUS_GENOMIC_3P_TAG))

    def test_rejects_a_sequence_truncated_at_either_end(self):
        body = "ACGT" * 10
        self.assertFalse(has_genomic_panhandle_termini(body + HANTAVIRUS_GENOMIC_3P_TAG))
        self.assertFalse(has_genomic_panhandle_termini(HANTAVIRUS_GENOMIC_5P_TAG + body))
        self.assertFalse(has_genomic_panhandle_termini(body))

    def test_the_two_tags_are_reverse_complements(self):
        """The filter is only meaningful if the tags describe the two ends of ONE
        panhandle. If these ever drift apart the filter would silently select for
        something that is not a panhandle at all."""
        self.assertEqual(_reverse_complement(HANTAVIRUS_GENOMIC_5P_TAG),
                         HANTAVIRUS_GENOMIC_3P_TAG)


def _genomic_isolate(species: str, key_suffix: str, s_seq: str = "", m_seq: str = "",
                     l_seq: str = "", lab: bool = False, vaccine: bool = False,
                     ) -> "ah.IsolateRecord":
    """Minimal IsolateRecord carrying real accessions, for tests that need
    reference.genomic_index lookups to actually resolve -- _eco_isolate's isolates
    above all have empty accessions, which is fine for ecological_opportunity but
    cannot exercise species_consensus_panhandle_window or score_m_pack."""
    key = f"{species}:{key_suffix}"
    return ah.IsolateRecord(
        isolate_key=key, species_name=species, species_taxid="1",
        species_status="classified", biosample="", isolate=key_suffix, strain="",
        host_name="Test host", host_taxid="1", host_order="Rodentia",
        host_family="Cricetidae", host_role="reservoir_rodent",
        collection_date="2001", year="2001", location="", country="USA", region="",
        grouping_tier="tier_2_species_isolate", grouping_confidence=1.0,
        s_accession=f"{key}_S" if s_seq else "", m_accession=f"{key}_M" if m_seq else "",
        l_accession=f"{key}_L" if l_seq else "",
        s_completeness="complete", m_completeness="complete", l_completeness="complete",
        s_nt_length=len(s_seq), m_nt_length=len(m_seq), l_nt_length=len(l_seq),
        s_protein_aa="", m_protein_aa="", l_protein_aa="", concatenated_protein_aa="",
        segments_present=[], coding_complete_segments=[], triplet_status="",
        has_complete_triplet=False, is_lab_derived=lab, is_vaccine_strain=vaccine,
        has_nss_orf=False, metadata_conflict=False, taxonomic_incongruence=False,
        cross_species_partners=[], kmer_cluster_id="", s_protein_accessions=[],
        m_protein_accessions=[], l_protein_accessions=[], all_accessions=[], warnings=[],
    )


def _genomic_reference(isolates: List["ah.IsolateRecord"],
                       sequences: Dict[str, str]) -> ReferenceData:
    """ReferenceData whose genomic_index resolves every accession in `sequences` (keyed
    by accession, e.g. "SpA:iso1_S") to a FastaRecord carrying that sequence."""
    genomic_index = {
        accession: ah.FastaRecord(accession=accession, accession_base=accession,
                                  header="", sequence=seq, length=len(seq))
        for accession, seq in sequences.items()
    }
    return ReferenceData(isolates=isolates, cds_index={}, genomic_index=genomic_index)


# Test-only, deliberately SHORTER than RPI_PANHANDLE_WINDOW_NT (50): these tests pass a
# custom weights table with a matching m_pack_window_nt so fixture sequences can stay
# small while still exercising has_genomic_panhandle_termini's real 8-nt tags.
_TEST_PANHANDLE_WINDOW_NT = 12


def _tagged_seq(five_body: str, three_body: str) -> str:
    """A genomic sequence carrying real conserved termini around the given bodies,
    short enough for _TEST_PANHANDLE_WINDOW_NT (12 nt: the 8-nt tag plus a few bases of
    the body it is windowed with)."""
    return HANTAVIRUS_GENOMIC_5P_TAG + five_body + three_body + HANTAVIRUS_GENOMIC_3P_TAG


class SpeciesConsensusPanhandleWindowTest(unittest.TestCase):
    def test_majority_vote_across_isolates_of_one_species(self):
        # Three isolates of SpA's M segment; two agree on the base right after the tag,
        # one differs -- majority vote should pick the two-out-of-three base.
        seqs = {
            "SpA:i1_M": _tagged_seq("AA", "TT"),
            "SpA:i2_M": _tagged_seq("AA", "TT"),
            "SpA:i3_M": _tagged_seq("GG", "TT"),
        }
        isolates = [_genomic_isolate("SpA", f"i{n}", m_seq=seqs[f"SpA:i{n}_M"])
                   for n in (1, 2, 3)]
        reference = _genomic_reference(isolates, seqs)
        five = species_consensus_panhandle_window(
            reference, "SpA", "M", "5", n=_TEST_PANHANDLE_WINDOW_NT)
        self.assertEqual(five, HANTAVIRUS_GENOMIC_5P_TAG + "AATT")

    def test_excludes_deposits_missing_true_termini(self):
        """A truncated deposit -- present, but not carrying the conserved tags --
        must not enter the consensus. Its "terminal window" is arbitrary internal
        sequence (WS3.1 finding), and a majority-vote consensus is not exempt from
        that: including it would corrupt the very average this function computes."""
        seqs = {
            "SpA:i1_M": _tagged_seq("AA", "TT"),
            "SpA:i2_M": "GGGGGGGGGGGGGGGGGGGGGGGG",  # no conserved tags at all
        }
        isolates = [_genomic_isolate("SpA", "i1", m_seq=seqs["SpA:i1_M"]),
                   _genomic_isolate("SpA", "i2", m_seq=seqs["SpA:i2_M"])]
        reference = _genomic_reference(isolates, seqs)
        five = species_consensus_panhandle_window(
            reference, "SpA", "M", "5", n=_TEST_PANHANDLE_WINDOW_NT)
        self.assertEqual(five, HANTAVIRUS_GENOMIC_5P_TAG + "AATT")

    def test_excludes_lab_derived_and_vaccine_isolates(self):
        seqs = {"SpA:i1_M": _tagged_seq("AA", "TT")}
        isolates = [_genomic_isolate("SpA", "i1", m_seq=seqs["SpA:i1_M"], lab=True)]
        reference = _genomic_reference(isolates, seqs)
        self.assertEqual(
            species_consensus_panhandle_window(reference, "SpA", "M", "5"), "")

    def test_returns_empty_string_when_species_has_no_usable_sequence(self):
        reference = _genomic_reference([], {})
        self.assertEqual(
            species_consensus_panhandle_window(reference, "Nobody", "M", "5"), "")


class RpiSpeciesAssignmentTest(unittest.TestCase):
    def test_two_implicated_species_gives_donor_and_backbone(self):
        gate1 = Gate1Result(
            top_hits={"M": [SpeciesHit(species="Donor", distance=0.1)]},
            implicated_species=["Backbone", "Donor"], discordant_segment="M",
            pattern="m_only")
        donor, backbone = rpi_species_assignment(gate1, "M")
        self.assertEqual((donor, backbone), ("Donor", "Backbone"))

    def test_concordant_pattern_returns_the_same_species_for_both(self):
        """Gate1 found one species matching every segment -- donor == backbone is the
        CORRECT answer here, not the self-comparison bug D3 forbids. See
        score_m_pack's docstring for why."""
        gate1 = Gate1Result(
            top_hits={"M": [SpeciesHit(species="OnlySpecies", distance=0.1)]},
            implicated_species=["OnlySpecies"], discordant_segment="", pattern="concordant")
        donor, backbone = rpi_species_assignment(gate1, "M")
        self.assertEqual((donor, backbone), ("OnlySpecies", "OnlySpecies"))

    def test_three_or_more_implicated_species_is_unresolved(self):
        """Every segment disagrees -- gate1's majority-vote heuristic does not resolve
        a single backbone, and forcing a choice would be arbitrary."""
        gate1 = Gate1Result(
            top_hits={"M": [SpeciesHit(species="Donor", distance=0.1)]},
            implicated_species=["Backbone1", "Backbone2", "Donor"],
            discordant_segment="M", pattern="m_only")
        self.assertEqual(rpi_species_assignment(gate1, "M"), ("Donor", ""))

    def test_no_hits_for_the_odd_segment_is_unresolved(self):
        gate1 = Gate1Result(top_hits={"M": []}, implicated_species=["Backbone"],
                            discordant_segment="M", pattern="m_only")
        self.assertEqual(rpi_species_assignment(gate1, "M"), ("", ""))


class ScoreMPackTest(unittest.TestCase):
    """D3/Task 10 (WS3.3). Uses a custom weights table with a short
    _TEST_PANHANDLE_WINDOW_NT so fixture sequences stay small while still exercising
    has_genomic_panhandle_termini's real tags -- production always reads
    RPI_WEIGHTS["m_pack_window_nt"] (50), asserted separately by
    MPackTauRescalingTest.test_window_matches_the_window_tau_was_measured_under."""

    def _weights(self) -> Dict[str, Any]:
        table = copy.deepcopy(RPI_WEIGHTS)
        table["m_pack_window_nt"] = _TEST_PANHANDLE_WINDOW_NT
        return table

    def _reference_and_gate1(self):
        """One backbone species (Backbone) with a consistent S/M/L consensus, and one
        donor species (Donor) whose M segment consensus differs -- the ordinary
        single-discordant-segment (M) case."""
        seqs = {
            "Backbone:i1_S": _tagged_seq("AA", "TT"),
            "Backbone:i1_M": _tagged_seq("AA", "TT"),
            "Backbone:i1_L": _tagged_seq("AA", "TT"),
            "Donor:i1_M": _tagged_seq("CC", "GG"),
        }
        isolates = [
            _genomic_isolate("Backbone", "i1", s_seq=seqs["Backbone:i1_S"],
                            m_seq=seqs["Backbone:i1_M"], l_seq=seqs["Backbone:i1_L"]),
            _genomic_isolate("Donor", "i1", m_seq=seqs["Donor:i1_M"]),
        ]
        reference = _genomic_reference(isolates, seqs)
        gate1 = Gate1Result(
            top_hits={"M": [SpeciesHit(species="Donor", distance=0.1)]},
            implicated_species=["Backbone", "Donor"], discordant_segment="M",
            pattern="m_only")
        return reference, gate1

    def test_candidate_matching_the_backbone_scores_near_one(self):
        """The candidate's ACTUAL M segment is identical to what the backbone species'
        own consensus would produce -- ddG should be ~0 and M_pack near 1.0."""
        reference, gate1 = self._reference_and_gate1()
        candidate = CandidateInput(
            candidate_id="C1",
            s_genomic=_tagged_seq("AA", "TT"), m_genomic=_tagged_seq("AA", "TT"),
            l_genomic=_tagged_seq("AA", "TT"))
        value, status = score_m_pack(candidate, reference, gate1, "M",
                                     weights=self._weights())
        self.assertEqual(status, "scored")
        self.assertGreater(value, 0.9)

    def test_never_self_comparison(self):
        """D3: the candidate side and the reference side must come from genuinely
        different sources. Mutating the candidate's own sequence after scoring must
        not retroactively change a result that already used the reference consensus --
        the two are structurally independent, not aliased."""
        reference, gate1 = self._reference_and_gate1()
        candidate = CandidateInput(
            candidate_id="C1",
            s_genomic=_tagged_seq("AA", "TT"), m_genomic=_tagged_seq("AA", "TT"),
            l_genomic=_tagged_seq("AA", "TT"))
        value_before, _ = score_m_pack(candidate, reference, gate1, "M",
                                       weights=self._weights())
        candidate.m_genomic = _tagged_seq("CC", "GG")
        value_after, _ = score_m_pack(candidate, reference, gate1, "M",
                                      weights=self._weights())
        self.assertNotEqual(value_before, value_after)

    def test_no_parental_sequence_when_backbone_is_unresolved(self):
        reference, gate1 = self._reference_and_gate1()
        gate1.implicated_species = ["Backbone1", "Backbone2", "Donor"]
        candidate = CandidateInput(
            candidate_id="C1", s_genomic=_tagged_seq("AA", "TT"),
            m_genomic=_tagged_seq("AA", "TT"), l_genomic=_tagged_seq("AA", "TT"))
        value, status = score_m_pack(candidate, reference, gate1, "M",
                                     weights=self._weights())
        self.assertEqual(status, "no_parental_sequence")
        self.assertEqual(value, 1.0)

    def test_no_parental_sequence_when_species_has_no_reference_sequence(self):
        reference = _genomic_reference([], {})
        gate1 = Gate1Result(
            top_hits={"M": [SpeciesHit(species="Donor", distance=0.1)]},
            implicated_species=["Backbone", "Donor"], discordant_segment="M",
            pattern="m_only")
        candidate = CandidateInput(
            candidate_id="C1", s_genomic=_tagged_seq("AA", "TT"),
            m_genomic=_tagged_seq("AA", "TT"), l_genomic=_tagged_seq("AA", "TT"))
        value, status = score_m_pack(candidate, reference, gate1, "M",
                                     weights=self._weights())
        self.assertEqual(status, "no_parental_sequence")
        self.assertEqual(value, 1.0)

    def test_undefined_energy_status_on_a_non_finite_ddg(self):
        """N4: boltzmann_compatibility raises ValueError on a non-finite ddG rather
        than silently scoring it. score_m_pack must catch that and report
        undefined_energy, log-and-continue, per the global constraint -- never let the
        exception propagate out of score_candidate."""
        reference, gate1 = self._reference_and_gate1()
        candidate = CandidateInput(
            candidate_id="C1", s_genomic=_tagged_seq("AA", "TT"),
            m_genomic=_tagged_seq("AA", "TT"), l_genomic=_tagged_seq("AA", "TT"))
        with unittest.mock.patch(
                __name__ + ".boltzmann_compatibility",
                side_effect=ValueError("non-finite ddG")):
            value, status = score_m_pack(candidate, reference, gate1, "M",
                                         weights=self._weights())
        self.assertEqual(status, "undefined_energy")
        self.assertEqual(value, 1.0)

    def test_never_produces_a_literal_zero_or_a_value_above_one(self):
        """Property this function must preserve regardless of status: compute_rpi's
        no-literal-zeros guarantee starts here. A maximally mismatched candidate M
        segment (every base opposite the backbone consensus) is the most adverse
        real-sequence input this function will see."""
        reference, gate1 = self._reference_and_gate1()
        candidate = CandidateInput(
            candidate_id="C1", s_genomic=_tagged_seq("AA", "TT"),
            m_genomic=_tagged_seq("GG", "AA"), l_genomic=_tagged_seq("AA", "TT"))
        value, _status = score_m_pack(candidate, reference, gate1, "M",
                                      weights=self._weights())
        self.assertGreater(value, 0.0)
        self.assertLessEqual(value, 1.0)


def _orf_seq(distinguishing_run: str) -> str:
    """A minimal genomic sequence whose longest ORF is ATG + 11 A's + `distinguishing_run`
    + 10 G's + TAA. Indices 0-13 of the CDS are the invariant-tag region (14 nt);
    indices 14-22 are exactly `_beyond_tag`'s scoring window at footprint=23 (SNV's
    footprint, per n_recognition_footprint) -- so `distinguishing_run` (9 nt) is the
    ONLY thing polymerase_compatibility can see differ between two sequences built from
    this helper."""
    return "ATG" + "A" * 11 + distinguishing_run + "G" * 10 + "TAA"


class SpeciesConsensusCdsWindowTest(unittest.TestCase):
    def test_majority_vote_across_isolates_of_one_species(self):
        seqs = {
            "SpA:i1_M": _orf_seq("CCCCCCCCC"),
            "SpA:i2_M": _orf_seq("CCCCCCCCC"),
            "SpA:i3_M": _orf_seq("TTTTTTTTT"),
        }
        isolates = [_genomic_isolate("SpA", f"i{n}", m_seq=seqs[f"SpA:i{n}_M"])
                   for n in (1, 2, 3)]
        reference = _genomic_reference(isolates, seqs)
        window = species_consensus_cds_window(reference, "SpA", "M", window=30)
        self.assertEqual(_beyond_tag(window, 23), "CCCCCCCCC")

    def test_isolates_with_no_findable_orf_are_excluded(self):
        seqs = {
            "SpA:i1_M": _orf_seq("CCCCCCCCC"),
            "SpA:i2_M": "GGGGGGGGGGGGGGGGGGGGGGGG",  # no ATG at all -- no ORF
        }
        isolates = [_genomic_isolate("SpA", "i1", m_seq=seqs["SpA:i1_M"]),
                   _genomic_isolate("SpA", "i2", m_seq=seqs["SpA:i2_M"])]
        reference = _genomic_reference(isolates, seqs)
        window = species_consensus_cds_window(reference, "SpA", "M", window=30)
        self.assertEqual(_beyond_tag(window, 23), "CCCCCCCCC")

    def test_not_filtered_on_genomic_panhandle_termini(self):
        """A deposit missing its UTR panhandle tags can still have a perfectly good
        ORF -- has_genomic_panhandle_termini must NOT gate this function the way it
        gates species_consensus_panhandle_window."""
        seq = _orf_seq("CCCCCCCCC")
        self.assertFalse(has_genomic_panhandle_termini(seq))
        isolates = [_genomic_isolate("SpA", "i1", m_seq=seq)]
        reference = _genomic_reference(isolates, {"SpA:i1_M": seq})
        window = species_consensus_cds_window(reference, "SpA", "M", window=30)
        self.assertEqual(_beyond_tag(window, 23), "CCCCCCCCC")

    def test_returns_empty_string_when_species_has_no_usable_sequence(self):
        reference = _genomic_reference([], {})
        self.assertEqual(
            species_consensus_cds_window(reference, "Nobody", "M", window=30), "")


class ScoreMPolTest(unittest.TestCase):
    """D3/Task 11 (WS3.4). Backbone species is the REAL "Orthohantavirus sinnombreense"
    (footprint 23, via the real n_recognition_footprint lookup) so these tests exercise
    the actual species-specific footprint wiring rather than a fake standing in for it."""

    BACKBONE = "Orthohantavirus sinnombreense"
    DONOR = "Orthohantavirus andesense"

    def _reference_and_gate1(self, m_run: str = "AAAAAAAAA"):
        seqs = {
            f"{self.BACKBONE}:i1_S": _orf_seq("AAAAAAAAA"),
            f"{self.BACKBONE}:i1_M": _orf_seq("AAAAAAAAA"),
            f"{self.BACKBONE}:i1_L": _orf_seq("AAAAAAAAA"),
            f"{self.DONOR}:i1_M": _orf_seq(m_run),
        }
        isolates = [
            _genomic_isolate(self.BACKBONE, "i1", s_seq=seqs[f"{self.BACKBONE}:i1_S"],
                            m_seq=seqs[f"{self.BACKBONE}:i1_M"],
                            l_seq=seqs[f"{self.BACKBONE}:i1_L"]),
            _genomic_isolate(self.DONOR, "i1", m_seq=seqs[f"{self.DONOR}:i1_M"]),
        ]
        reference = _genomic_reference(isolates, seqs)
        gate1 = Gate1Result(
            top_hits={"M": [SpeciesHit(species=self.DONOR, distance=0.1)]},
            implicated_species=[self.BACKBONE, self.DONOR], discordant_segment="M",
            pattern="m_only")
        return reference, gate1

    def test_candidate_matching_the_backbone_scores_one(self):
        reference, gate1 = self._reference_and_gate1()
        candidate_cds = {
            "S": _orf_seq("AAAAAAAAA"), "M": _orf_seq("AAAAAAAAA"),
            "L": _orf_seq("AAAAAAAAA"),
        }
        value, status = score_m_pol(candidate_cds, reference, gate1, "M", window=30)
        self.assertEqual(status, "scored")
        self.assertEqual(value, 1.0)

    def test_candidate_with_a_mismatched_odd_segment_scores_below_one(self):
        """The candidate's ACTUAL M footprint (beyond the tag) does not match what the
        backbone species' own M would look like -- this is the signal M_pol exists to
        carry."""
        reference, gate1 = self._reference_and_gate1()
        candidate_cds = {
            "S": _orf_seq("AAAAAAAAA"), "M": _orf_seq("TTTTTTTTT"),
            "L": _orf_seq("AAAAAAAAA"),
        }
        value, status = score_m_pol(candidate_cds, reference, gate1, "M", window=30)
        self.assertEqual(status, "scored")
        self.assertLess(value, 1.0)
        self.assertGreater(value, 0.0)

    def test_never_self_comparison(self):
        """D3: mutating the candidate's own CDS after scoring must not retroactively
        change an already-computed result -- the reference side is structurally
        independent, sourced from reference.genomic_index, not aliased to the
        candidate's own dict the way the pre-WS3.4 bug aliased it."""
        reference, gate1 = self._reference_and_gate1()
        candidate_cds = {"S": _orf_seq("AAAAAAAAA"), "M": _orf_seq("AAAAAAAAA"),
                         "L": _orf_seq("AAAAAAAAA")}
        value_before, _ = score_m_pol(candidate_cds, reference, gate1, "M", window=30)
        candidate_cds["M"] = _orf_seq("TTTTTTTTT")
        value_after, _ = score_m_pol(candidate_cds, reference, gate1, "M", window=30)
        self.assertNotEqual(value_before, value_after)

    def test_uses_the_backbone_species_own_footprint_length(self):
        """Confirms the footprint is actually keyed on the BACKBONE species (23 nt,
        sinnombreense) and not silently defaulted to 9 -- at the 9-nt default this
        candidate's mismatch would fall entirely inside the invariant tag and score
        1.0 regardless of the real difference."""
        reference, gate1 = self._reference_and_gate1()
        candidate_cds = {"S": _orf_seq("AAAAAAAAA"), "M": _orf_seq("TTTTTTTTT"),
                         "L": _orf_seq("AAAAAAAAA")}
        value, _status = score_m_pol(candidate_cds, reference, gate1, "M", window=30)
        self.assertLess(value, 1.0)

    def test_no_parental_sequence_when_backbone_is_unresolved(self):
        reference, gate1 = self._reference_and_gate1()
        gate1.implicated_species = [self.BACKBONE, "SomeOtherSpecies", self.DONOR]
        candidate_cds = {"S": _orf_seq("AAAAAAAAA"), "M": _orf_seq("AAAAAAAAA"),
                         "L": _orf_seq("AAAAAAAAA")}
        value, status = score_m_pol(candidate_cds, reference, gate1, "M", window=30)
        self.assertEqual(status, "no_parental_sequence")
        self.assertEqual(value, 1.0)

    def test_no_parental_sequence_when_species_has_no_reference_sequence(self):
        reference = _genomic_reference([], {})
        gate1 = Gate1Result(
            top_hits={"M": [SpeciesHit(species=self.DONOR, distance=0.1)]},
            implicated_species=[self.BACKBONE, self.DONOR], discordant_segment="M",
            pattern="m_only")
        candidate_cds = {"S": _orf_seq("AAAAAAAAA"), "M": _orf_seq("AAAAAAAAA"),
                         "L": _orf_seq("AAAAAAAAA")}
        value, status = score_m_pol(candidate_cds, reference, gate1, "M", window=30)
        self.assertEqual(status, "no_parental_sequence")
        self.assertEqual(value, 1.0)

    def test_never_produces_a_literal_zero(self):
        reference, gate1 = self._reference_and_gate1()
        candidate_cds = {"S": _orf_seq("TTTTTTTTT"), "M": _orf_seq("TTTTTTTTT"),
                         "L": _orf_seq("TTTTTTTTT")}
        value, _status = score_m_pol(candidate_cds, reference, gate1, "M", window=30)
        self.assertGreater(value, 0.0)
        self.assertLessEqual(value, 1.0)


class RetentionFactorTest(unittest.TestCase):
    def test_m_exchange_is_favoured_over_s_exchange(self):
        """Generation (Rodriguez, in vitro): M:S:L = 12:12:1. Retention (Henderson,
        natural, Tier 3 secondhand): 14:2:1. S-exchange is generated as often as M but
        survives markedly less. That gap IS the selection filter."""
        self.assertGreater(retention_factor("M", True), retention_factor("S", True))

    def test_structural_damage_reduces_retention(self):
        self.assertLess(retention_factor("M", False), retention_factor("M", True))

    def test_never_returns_zero(self):
        for segment in ("S", "M", "L"):
            for intact in (True, False):
                self.assertGreater(retention_factor(segment, intact), 0.0)

    def test_bounded_above_by_one(self):
        for segment in ("S", "M", "L"):
            self.assertLessEqual(retention_factor(segment, True), 1.0)


class RpiAssemblyTest(unittest.TestCase):
    def _row(self, **kwargs):
        defaults = dict(
            implicated_species=["A", "B"],
            isolates=[_eco_isolate("A", "USA", "2001"), _eco_isolate("B", "USA", "2002")],
            odd_segment="M", identity=0.85,
            m_pack=1.0, m_pack_status="scored", m_pol=1.0, m_pol_status="scored",
            structural_intact=True,
        )
        defaults.update(kwargs)
        # This fixture's `isolates` is already exactly one isolate per implicated
        # species -- its own known parent, not a production-style unknown -- so unless
        # a test overrides parent_isolates/isolate_basis itself, derive them the same
        # way _score_panel does.
        defaults.setdefault(
            "parent_isolates", _parent_map(*defaults["isolates"]))
        defaults.setdefault("isolate_basis", "actual_parents")
        return compute_rpi(**defaults)

    def test_emits_every_column_from_spec_section_9(self):
        row = self._row()
        for column in ("rpi_generation", "rpi_retention", "rpi_e_host", "rpi_e_geo",
                       "rpi_e_time", "rpi_e_total", "rpi_m_pack", "rpi_m_pack_status",
                       "rpi_m_pol", "rpi_m_pol_status", "rpi_k_stratum", "rpi_k_value",
                       "rpi_f", "rpi_smoothing", "rpi_weights_version"):
            self.assertIn(column, row, column)

    def test_reported_e_total_is_never_a_literal_zero(self):
        """The no-literal-zeros guarantee has to hold in the ROW, not just in
        ecological_opportunity's raw return. rpi_e_total is the product of three terms
        each floored at eco_floor, so it bottoms out at floor**3 = 1e-6; reporting it at
        4 dp alongside its own components collapsed ordinary candidates to 0.0. Two
        implicated species, different host orders, different countries, three years
        apart gives 3.33e-05 -- which the CSV showed as 0.0.

        test_no_factor_is_ever_exactly_zero asserts on the raw dict and so cannot see
        this."""
        far_a = _eco_isolate("SpA", "USA", "2001", host_order="Rodentia",
                             host_family="Cricetidae", host_name="Peromyscus maniculatus")
        far_b = _eco_isolate("SpB", "Chile", "2004", host_order="Soricomorpha",
                             host_family="Soricidae", host_name="Sorex araneus")
        row = compute_rpi(["SpA", "SpB"], [far_a, far_b], _parent_map(far_a, far_b),
                          "actual_parents", "M", 0.9, 1.0, "fixed_constant",
                          1.0, "fixed_constant", True)
        self.assertGreater(row["rpi_e_total"], 0.0)
        self.assertGreater(row["rpi_generation"], 0.0)
        self.assertGreater(row["rpi_retention"], 0.0)

    def test_reported_retention_is_never_a_literal_zero_at_the_deepest_floor(self):
        """rpi_retention carries one more floored factor than rpi_generation (F), so it
        bottoms out at floor**5 = 1e-10. At 8 dp that rounded to 0.0."""
        far_a = _eco_isolate("SpA", "USA", "2001", host_order="Rodentia",
                             host_family="Cricetidae", host_name="Peromyscus maniculatus")
        far_b = _eco_isolate("SpB", "Chile", "2004", host_order="Soricomorpha",
                             host_family="Soricidae", host_name="Sorex araneus")
        row = compute_rpi(["SpA", "SpB"], [far_a, far_b], _parent_map(far_a, far_b),
                          "actual_parents", "S", 0.5,
                          RPI_WEIGHTS["eco_floor"], "fixed_constant",
                          RPI_WEIGHTS["eco_floor"], "fixed_constant", False)
        self.assertEqual(row["rpi_generation"], RPI_WEIGHTS["eco_floor"] ** 4)
        self.assertGreater(row["rpi_retention"], 0.0)

    def test_rounding_never_reports_e_total_above_its_own_components(self):
        """rpi_e_total = h * g * t with every term <= 1, so it can never exceed any one
        of them. round() preserves that only at EQUAL precision: reporting the total at
        6 dp beside its components at 4 dp put the PRODUCT above its own factor.

        The fixture is chosen to produce a component that 4 dp cannot represent -- same
        host and country (h = g = 1.0) with sampling years three apart (t = 1/3). At
        split precision this reported rpi_e_total = 0.333333 against
        rpi_e_time = 0.3333."""
        near_a = _eco_isolate("SpA", "USA", "2001")
        near_b = _eco_isolate("SpB", "USA", "2004")
        seen_fractional = False
        for odd_segment in ("S", "M", "L"):
            for identity in (0.5, 0.9):
                row = compute_rpi(["SpA", "SpB"], [near_a, near_b],
                                  _parent_map(near_a, near_b), "actual_parents",
                                  odd_segment, identity, 1.0, "fixed_constant", 1.0,
                                  "fixed_constant", True)
                components = (row["rpi_e_host"], row["rpi_e_geo"], row["rpi_e_time"])
                if any(0.0 < c < 1.0 for c in components):
                    seen_fractional = True
                self.assertLessEqual(row["rpi_e_total"], min(components),
                                     f"{odd_segment}/{identity}: {row['rpi_e_total']} "
                                     f"> min{components}")
        # Guard the guard: if every component came back 1.0 the assertion above is
        # trivially true and this test would pass while measuring nothing.
        self.assertTrue(seen_fractional,
                        "fixture stopped producing a fractional ecological term")

    def test_rounding_never_reports_retention_above_generation(self):
        """F <= 1 makes rpi_retention <= rpi_generation an invariant of the model, but
        rounding only preserves order at EQUAL precision. Reporting retention finer than
        generation put 1.57233e-05 above 1.572e-05 -- an inversion created purely by
        formatting."""
        far_a = _eco_isolate("SpA", "USA", "2001", host_order="Rodentia",
                             host_family="Cricetidae", host_name="Peromyscus maniculatus")
        far_b = _eco_isolate("SpB", "Chile", "2004", host_order="Soricomorpha",
                             host_family="Soricidae", host_name="Sorex araneus")
        seen_strict = False
        for odd_segment in ("S", "M", "L"):
            for intact in (True, False):
                row = compute_rpi(["SpA", "SpB"], [far_a, far_b],
                                  _parent_map(far_a, far_b), "actual_parents",
                                  odd_segment, 0.9,
                                  1.0, "fixed_constant", 1.0, "fixed_constant", intact)
                if row["rpi_retention"] < row["rpi_generation"]:
                    seen_strict = True
                self.assertLessEqual(row["rpi_retention"], row["rpi_generation"],
                                     f"{odd_segment}/{intact}")
        # Guard the guard: where F = 1.0 the two columns are trivially EQUAL and the
        # assertion above measures nothing. At least one combination has to have F < 1,
        # i.e. has to actually exercise the scaling this invariant is about. Its N1
        # sibling test_rounding_never_reports_e_total_above_its_own_components has had
        # this shape since it was written; this one did not.
        self.assertTrue(seen_strict,
                        "fixture stopped producing a row where F < 1.0, so the ordering "
                        "assertion above is trivially satisfied by equality")

    def test_split_precision_is_what_would_invert_the_pair(self):
        """Companion to the test above, which cannot fail through _rpi_decimals at all:
        compute_rpi passes floor**5 to BOTH columns, so no change to that function can
        split their precisions, and 1 of its 6 combinations would catch a call-site
        split. That makes it a guard whose subject is unreachable from the guard.

        This one reproduces the inversion directly, so the equal-precision choice in
        compute_rpi is protected by a test that can actually observe the failure:
        1.57233e-05 at 10 dp is 1.57233e-05, and the SAME value at 8 dp is 1.572e-05, so
        reporting retention finer than generation puts retention above a generation it
        can never exceed."""
        raw = 1.572330e-05
        self.assertGreater(round(raw, 10), round(raw, 8),
                           "split precision no longer inverts -- this test's premise is "
                           "stale, not satisfied")
        # And the precision compute_rpi actually uses is the deeper of the two, applied
        # to BOTH columns, which is what makes the inversion unreachable in production.
        self.assertEqual(_rpi_decimals(RPI_WEIGHTS["eco_floor"] ** 5), 10)
        self.assertEqual(round(raw, 10), round(raw, 10))

    def test_k_stratum_carries_a_provisional_marker(self):
        """The `identity` the production call site supplies is gate0.best_kmer_identity
        -- PROTEIN k-mer novelty against the nearest reference, maximised over segments
        -- not the parental NUCLEOTIDE divergence spec section 7 stratifies on, and on a
        different scale (~0.35-0.46 vs the boundary's 0.548 / 0.707 / 0.85 calibration
        points). Without an explicit marker a reader takes rpi_k_stratum for a real
        divergence stratification. Computing real parental identity is a documented
        follow-up; saying so on every row is not."""
        row = self._row()
        self.assertIn("rpi_k_stratum_provisional", row)
        self.assertTrue(row["rpi_k_stratum_provisional"])

    def test_both_scores_are_within_unit_interval(self):
        # Exercise both the trivial F=1.0 default AND a row where F<1.0 -- the default
        # alone can't show rpi_retention staying in-bounds once it's actually scaled
        # down rather than merely equal to rpi_generation.
        for row in (self._row(), self._row(structural_intact=False)):
            for key in ("rpi_generation", "rpi_retention"):
                self.assertGreater(row[key], 0.0)
                self.assertLessEqual(row[key], 1.0)

    def test_retention_never_exceeds_generation(self):
        row = self._row()
        self.assertLessEqual(row["rpi_retention"], row["rpi_generation"])

    def test_retention_strictly_below_generation_when_f_bites(self):
        """The non-strict test above passes trivially at the M/intact default, where
        retention_factor returns exactly 1.0 and rpi_retention == rpi_generation -- a
        regression that dropped F from the retention computation entirely would still
        pass it. Force F < 1.0 (structural damage) so the ceiling actually bites, and
        assert the inequality is STRICT."""
        row = self._row(structural_intact=False)
        self.assertLess(row["rpi_retention"], row["rpi_generation"])

    def test_m_exchange_outscores_s_exchange_in_distant_stratum(self):
        m = self._row(odd_segment="M", identity=0.55)
        s = self._row(odd_segment="S", identity=0.55)
        self.assertGreater(m["rpi_generation"], s["rpi_generation"])

    def test_ecological_impossibility_suppresses_the_score(self):
        """Multiplicative, not additive: a perfect molecular profile must NOT rescue a
        candidate whose parents never co-occurred."""
        impossible = self._row(isolates=[
            _eco_isolate("A", "USA", "1985"),
            _eco_isolate("B", "Chile", "2020", host_order="Soricomorpha",
                         host_family="Soricidae", host_name="Sorex araneus"),
        ])
        plausible = self._row()
        self.assertLess(impossible["rpi_generation"], plausible["rpi_generation"] / 100)

    def test_no_factor_is_ever_exactly_zero(self):
        row = self._row(m_pack=1e-12, m_pol=1e-12)
        self.assertGreater(row["rpi_generation"], 0.0)

    def test_fixed_constant_m_pol_status_is_propagated(self):
        row = self._row(m_pol=1.0, m_pol_status="fixed_constant")
        self.assertEqual(row["rpi_m_pol_status"], "fixed_constant")


class VerdictUnchangedRegressionTest(unittest.TestCase):
    def test_combine_verdict_itself_is_unaffected_by_rpi_inputs(self):
        """Spec 10.1: the verdict is an OBSERVATIONAL REPORT; the RPI is a MECHANISTIC
        MODEL. This exercises combine_verdict directly -- it duplicates part of the
        pre-existing CombineVerdictTest coverage and is retained only to document that
        combine_verdict's signature/behaviour takes no RPI input at all. It does NOT
        exercise score_candidate's own wiring; see
        test_compute_rpi_output_cannot_overwrite_the_verdict_key below for the
        assertion that actually covers what Task 7 built."""
        gate0 = Gate0Result(passed_structure=True, passed_identity_floor=True)
        gate1 = Gate1Result(pattern="single_segment_discordant", implicated_species=["A", "B"])
        gate2a, gate2b, gate2c = Gate2TrackAResult(), Gate2TrackBResult(), Gate2TrackCResult()
        gate3 = Gate3Result(status="ecological_overlap_plausible")
        gate4 = Gate4Result(phylo_status="analyzed",
                            reassortment_call="supported_reassortment_candidate")
        self.assertEqual(
            combine_verdict(gate0, gate1, gate2a, gate2b, gate2c, gate3, gate4),
            "supported_reassortment_candidate",
        )

    def test_compute_rpi_output_cannot_overwrite_the_verdict_key(self):
        """score_candidate wires RPI in via `result.update(compute_rpi(...))`, AFTER
        `result["verdict"]` is already set from combine_verdict (spec 10.2). dict.update
        overwrites any key compute_rpi's return value shares with `result` -- so the
        guarantee that RPI can never silently clobber the verdict rests entirely on
        compute_rpi's key set staying disjoint from {"verdict"}. This is what a future
        change that folds RPI into the reported verdict inside score_candidate's own
        post-processing (i.e. NOT inside combine_verdict, so CombineVerdictTest and the
        test above would both stay green) would actually break. Asserted directly,
        without needing full score_candidate fixtures/reference data."""
        isolate_a, isolate_b = _eco_isolate("A", "USA", "2001"), _eco_isolate("B", "USA", "2002")
        row = compute_rpi(
            implicated_species=["A", "B"],
            isolates=[isolate_a, isolate_b],
            parent_isolates=_parent_map(isolate_a, isolate_b), isolate_basis="actual_parents",
            odd_segment="M", identity=0.85,
            m_pack=1.0, m_pack_status="scored", m_pol=1.0, m_pol_status="scored",
            structural_intact=True,
        )
        self.assertNotIn("verdict", row)
        self.assertTrue(set(row.keys()).isdisjoint({"verdict"}))


class RetrospectiveReproductionTest(unittest.TestCase):
    def test_predicts_m_odd_for_every_distant_pair(self):
        report = retrospective_reproduction()
        for pair in ("SNV_x_ANDV", "SNV_x_BCCV"):
            self.assertEqual(report["pairs"][pair]["predicted_modal_segment"], "M", pair)

    def test_close_pair_does_not_predict_m_only(self):
        """Rodriguez saw 12 M-odd and 12 S-odd. Klempa's abstract says 'only M', but his
        own Table 2 lists an L-odd genotype for this same pair. The close stratum must
        NOT inherit the M-only rule."""
        report = retrospective_reproduction()
        prior = report["pairs"]["SNV_NMR11_x_CC107"]["predicted_prior"]
        self.assertAlmostEqual(prior["M"], prior["S"], places=6)

    def test_leave_one_pair_out_actually_changes_the_prior(self):
        """A pair used to DERIVE K cannot serve as an independent test OF K.

        The previous version only checked that the held-out name was absent from a
        `derived_from` LIST, while the prediction was still produced by
        constellation_prior reading the fixed RPI_WEIGHTS counts -- so the held-out
        pair's own 10 M-odd observations still generated the prediction that supposedly
        tested it, and the prior came back byte-identical either way. Holding out must
        move the number, or it holds nothing out."""
        baseline = retrospective_reproduction()
        held_out = retrospective_reproduction(holdout="SNV_x_ANDV")

        self.assertNotIn("SNV_x_ANDV", held_out["derived_from"])
        self.assertIn("SNV_x_ANDV", held_out["pairs"])
        self.assertTrue(held_out["holdout_changed_the_prior"])

        base_prior = baseline["pairs"]["SNV_x_ANDV"]["predicted_prior"]
        held_prior = held_out["pairs"]["SNV_x_ANDV"]["predicted_prior"]
        self.assertNotEqual(base_prior, held_prior)
        # SNV_x_ANDV contributes 10 of the distant stratum's 11 M-odd observations, so
        # removing it must visibly weaken the M prediction it is being tested against.
        self.assertLess(held_prior["M"], base_prior["M"] - 0.1)

    def test_leave_one_pair_out_keeps_every_category_strictly_positive(self):
        """Jeffreys smoothing is applied AFTER the subtraction, so a category the
        subtraction empties still cannot reach zero -- a literal zero inside a
        multiplicative model is the hidden hard gate spec section 7 forbids."""
        held_out = retrospective_reproduction(holdout="SNV_NMR11_x_CC107")
        prior = held_out["pairs"]["SNV_NMR11_x_CC107"]["predicted_prior"]
        self.assertEqual(held_out["prior_counts_used"]["close"], {"M": 0, "S": 0, "L": 0})
        for segment in ("S", "M", "L"):
            self.assertGreater(prior[segment], 0.0, segment)
        self.assertAlmostEqual(sum(prior.values()), 1.0, places=9)

    def test_tier_three_pairs_are_not_listed_as_prior_contributors(self):
        """PUUV_x_PHV and DOBV_Aa_x_DOBV_Af are NOT inside
        RPI_WEIGHTS['constellation_counts'] -- close is Rodriguez's 12/12/1 alone and
        distant is Rizvanov's 10 plus BCCV's 1. Listing all five as contributors
        overstated the evidence base and made three of the five leave-one-out runs
        silent no-ops that still claimed to hold something out."""
        report = retrospective_reproduction()
        self.assertEqual(sorted(report["derived_from"]),
                         ["SNV_NMR11_x_CC107", "SNV_x_ANDV", "SNV_x_BCCV"])
        for tier3 in ("PUUV_x_PHV", "DOBV_Aa_x_DOBV_Af"):
            self.assertNotIn(tier3, report["derived_from"], tier3)
            self.assertFalse(
                retrospective_reproduction(holdout=tier3)["holdout_changed_the_prior"], tier3)

    def test_composition_distance_exposes_what_modal_agrees_hides(self):
        """modal_agrees CANNOT return False on this benchmark: the distant top set is
        {M}, the close top set is {M, S}, and every observation's top set contains M.
        DOBV_Aa_x_DOBV_Af is the case that proves the boolean is empty -- 65 M / 0 S /
        0 L against a close prior that weights S EXACTLY as heavily as M. That is a
        genuine disagreement, and modal_agrees reports it as a pass."""
        report = retrospective_reproduction()
        row = report["pairs"]["DOBV_Aa_x_DOBV_Af"]
        self.assertTrue(row["modal_agrees"])
        self.assertGreater(row["composition_distance"], 0.5)
        # The pair the prior was actually derived from must, by contrast, sit close.
        rodriguez = report["pairs"]["SNV_NMR11_x_CC107"]
        self.assertLess(rodriguez["composition_distance"], 0.1)
        self.assertEqual(report["primary_agreement_measure"], "composition_distance")

    def test_composition_distance_is_total_variation_over_the_composition(self):
        prior = {"S": 0.5, "M": 0.4, "L": 0.1}
        self.assertAlmostEqual(_composition_distance(prior, {"S": 0, "M": 10, "L": 0}), 0.6)
        self.assertAlmostEqual(_composition_distance(prior, {"S": 5, "M": 4, "L": 1}), 0.0)
        self.assertTrue(math.isnan(_composition_distance(prior, {"S": 0, "M": 0, "L": 0})))

    def test_report_records_ecology_was_held_neutral(self):
        report = retrospective_reproduction()
        self.assertEqual(report["ecological_factor_held_at"], 1.0)

    def test_report_carries_the_confound_warning(self):
        report = retrospective_reproduction()
        self.assertIn("untestable", report["caveat"].lower())

    def test_close_pair_tie_is_represented_explicitly_not_by_dict_order(self):
        """The close stratum's prior gives M and S EXACTLY equal weight (12 vs 12 in the
        underlying tissue-culture data), and Rodriguez's observed counts are also
        exactly tied 12 M / 12 S. max(d, key=d.get) silently breaks ties by dict
        insertion order, which would report predicted='S' (constellation_prior's dict
        order) against observed='M' (RPI_INVITRO_BENCHMARK's dict order) as a
        DISAGREEMENT -- a pure key-ordering artefact, not a real model failure. A tie
        must surface as None, not a guessed letter, and the two tied sets must be
        recognised as agreeing."""
        report = retrospective_reproduction()
        row = report["pairs"]["SNV_NMR11_x_CC107"]
        self.assertIsNone(row["predicted_modal_segment"])
        self.assertEqual(row["predicted_top_segments"], ["M", "S"])
        self.assertTrue(row["modal_agrees"])


class VetoVsContinuousETest(unittest.TestCase):
    """Spec 11.3 cross-layer check. This function shipped once as a traceability-table
    claim with no implementation behind it; it does not also get to ship untested."""

    @staticmethod
    def _row(candidate_id, gate3_status, e_total):
        return {"candidate_id": candidate_id, "gate3_status": gate3_status,
                "rpi_e_total": e_total}

    def test_both_agreement_directions_are_counted_separately(self):
        result = veto_vs_continuous_e([
            self._row("veto_and_floored", "no_ecological_route", RPI_WEIGHTS["eco_floor"]),
            self._row("route_and_scored", "ecology_plausible", 0.6),
        ])
        self.assertEqual(result["agree_veto_and_near_zero_e"], 1)
        self.assertEqual(result["agree_route_and_nonzero_e"], 1)
        self.assertEqual(result["n_disagree"], 0)
        self.assertEqual(result["agreement_fraction"], 1.0)

    def test_disagreement_names_the_candidates_in_both_directions(self):
        """The counts alone are not actionable -- spec 10.2 says disagreement flags a
        defect in one layer or the other, which cannot be investigated without knowing
        WHICH candidate diverged."""
        result = veto_vs_continuous_e([
            self._row("vetoed_but_e_high", "no_ecological_route", 0.7),
            self._row("floored_but_not_vetoed", "ecology_plausible", 0.001),
        ])
        self.assertEqual(result["disagree_veto_but_e_not_near_zero"], ["vetoed_but_e_high"])
        self.assertEqual(result["disagree_near_zero_e_but_no_veto"], ["floored_but_not_vetoed"])
        self.assertEqual(result["n_agree"], 0)

    def test_rows_missing_either_key_are_skipped_not_scored_as_agreement(self):
        """junk_sequence rows never reach compute_rpi, so they carry no rpi_e_total.
        Counting them as agreement would inflate agreement_fraction with rows that were
        never compared."""
        result = veto_vs_continuous_e([
            {"candidate_id": "junk", "gate3_status": "not_run"},
            self._row("scored", "ecology_plausible", 0.6),
        ])
        self.assertEqual(result["n_skipped"], 1)
        self.assertEqual(result["skipped_candidates"], ["junk"])
        self.assertEqual(result["n_compared"], 1)

    def test_agreement_fraction_is_nan_when_nothing_was_compared(self):
        """0/0 is undefined. Reporting 0.0 or 1.0 for it would state a result about a
        comparison that never happened."""
        result = veto_vs_continuous_e([{"candidate_id": "junk"}])
        self.assertTrue(math.isnan(result["agreement_fraction"]))


def _parse_version(text: str) -> Tuple[int, ...]:
    """Numeric version tuple, ignoring any non-numeric suffix.

    Deliberately minimal rather than pulling in `packaging`: adding a dependency to a
    test that exists to police dependencies would be self-defeating, and every pin in
    these two files is a plain dotted numeric.
    """
    parts: List[int] = []
    for chunk in text.strip().split("."):
        digits = "".join(itertools.takewhile(str.isdigit, chunk))
        if not digits:
            break
        parts.append(int(digits))
    return tuple(parts)


def _read_requirement_pins(path: Path) -> Dict[str, Tuple[str, str]]:
    """{package: (operator, version)} from a requirements.txt-style file."""
    pins: Dict[str, Tuple[str, str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        match = re.match(r"^([A-Za-z0-9_.\-]+)\s*(==|>=)\s*([0-9][0-9.]*)$", line)
        if match:
            pins[match.group(1).lower()] = (match.group(2), match.group(3))
    return pins


def _read_conda_pins(path: Path) -> Dict[str, Tuple[str, str]]:
    """{package: (operator, version)} from an environment.yml dependency list.

    Hand-parsed rather than via PyYAML, which is not a project dependency and would have
    to be added purely to run this check.
    """
    pins: Dict[str, Tuple[str, str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.split("#", 1)[0].strip()
        if not stripped.startswith("- "):
            continue
        match = re.match(r"^-\s*([A-Za-z0-9_.\-]+)\s*(==|=|>=)\s*([0-9][0-9.]*)$", stripped)
        if match:
            operator = ">=" if match.group(2) == ">=" else "=="
            pins[match.group(1).lower()] = (operator, match.group(3))
    return pins


class ScorePanelDoesNotMutateGlobalWeightsTest(unittest.TestCase):
    """_score_panel used to reassign the module-level RPI_WEIGHTS while scoring.

    Correct single-threaded, unsafe for the WS-API use: the reassignment is
    process-global, so a scorer running concurrently would read perturbed weights for
    the duration of a sweep and return numbers attributable to no declared weights
    version -- silently, because each row still looked well-formed."""

    def test_global_is_untouched_while_a_perturbed_table_is_scored(self):
        panel = _rpi_sensitivity_panel()
        perturbed = copy.deepcopy(RPI_WEIGHTS)
        perturbed["eco_floor"] = 0.5
        perturbed["jeffreys_alpha"] = 9.0
        before_object, before_snapshot = RPI_WEIGHTS, copy.deepcopy(RPI_WEIGHTS)
        _score_panel(panel, perturbed)
        self.assertIs(RPI_WEIGHTS, before_object, "RPI_WEIGHTS was rebound")
        self.assertEqual(RPI_WEIGHTS, before_snapshot, "RPI_WEIGHTS contents changed")

    def test_the_perturbation_actually_reaches_compute_rpi(self):
        """Guard the guard: if `weights` were ignored rather than threaded, the test
        above would pass for the wrong reason -- the global would be untouched because
        nothing used the perturbed table at all."""
        panel = _rpi_sensitivity_panel()
        perturbed = copy.deepcopy(RPI_WEIGHTS)
        perturbed["retention_by_odd_segment"] = {"M": 0.02, "S": 0.02, "L": 0.02}
        self.assertNotEqual(_score_panel(panel, perturbed)["retention"],
                            _score_panel(panel)["retention"])

    def test_a_concurrent_scorer_is_unaffected_by_an_in_flight_sweep(self):
        """The failure this fix exists to prevent, stated directly: a row scored WHILE a
        sweep is mid-loop must equal the same row scored with no sweep in flight.

        Scoring after the sweep returns would prove nothing -- the old code restored the
        global in a finally block, so it was always correct once the call had finished.
        The bug was only observable DURING the loop, so the observation has to happen
        there. This panel entry runs the nested score when _score_panel reads its last
        key, which is inside the loop body, with the perturbed table installed under the
        old implementation."""
        def unperturbed_row():
            isolate_a, isolate_b = _eco_isolate("A", "USA", "2001"), _eco_isolate("B", "USA", "2002")
            return compute_rpi(["A", "B"], [isolate_a, isolate_b],
                               _parent_map(isolate_a, isolate_b), "actual_parents",
                               "M", 0.85, 1.0, "scored", 1.0, "scored", True)

        reference_row = unperturbed_row()
        observed = []

        class ObservingCandidate(dict):
            def __getitem__(self, key):
                if key == "structural_intact":       # read last by _score_panel
                    observed.append(unperturbed_row())
                return super().__getitem__(key)

        probe = ObservingCandidate(_rpi_sensitivity_panel()[0])
        perturbed = copy.deepcopy(RPI_WEIGHTS)
        perturbed["eco_floor"] = 0.9
        perturbed["host_decay"] = {k: 0.9 for k in RPI_WEIGHTS["host_decay"]}

        _score_panel([probe], perturbed)

        self.assertEqual(len(observed), 1, "probe never fired -- test measures nothing")
        self.assertEqual(observed[0], reference_row,
                         "a row scored mid-sweep differs from the same row scored "
                         "outside one: the sweep leaked its weights process-wide")


class EnvironmentReconciliationTest(unittest.TestCase):
    """environment.yml and requirements-pip.txt describe two different interpreters.

    They are allowed to differ in SCOPE. They are not allowed to pin INCOMPATIBLE
    versions of the same package. Nothing enforced that before this test, and the README
    recorded the contradiction as an open question rather than a checked invariant.
    """

    def setUp(self):
        root = Path(__file__).parent
        self.conda = _read_conda_pins(root / "environment.yml")
        self.pip = _read_requirement_pins(root / "requirements-pip.txt")

    def test_both_files_were_actually_parsed(self):
        """Guard the guard: a parser that silently matched nothing would make every
        assertion below vacuously true."""
        self.assertGreaterEqual(len(self.conda), 5, self.conda)
        self.assertGreaterEqual(len(self.pip), 4, self.pip)
        self.assertIn("viennarna", self.conda)
        self.assertIn("viennarna", self.pip)

    def test_every_shared_package_is_compatible(self):
        shared = sorted(set(self.conda) & set(self.pip))
        self.assertTrue(shared, "no packages in common -- parser drift?")
        for name in shared:
            conda_op, conda_version = self.conda[name]
            _pip_op, pip_version = self.pip[name]
            installed, required = _parse_version(pip_version), _parse_version(conda_version)
            if conda_op == ">=":
                self.assertGreaterEqual(
                    installed, required,
                    f"{name}: requirements-pip.txt pins {pip_version}, below "
                    f"environment.yml's floor {conda_version}")
            else:
                self.assertEqual(
                    installed, required,
                    f"{name}: environment.yml pins {name}=={conda_version} exactly, but "
                    f"requirements-pip.txt says {pip_version}")

    def test_viennarna_exact_pin_matches_in_both_files(self):
        """Singled out because it is the one pin where a mismatch changes RESULTS rather
        than merely risking an incompatibility: RNAcofold's free energies depend on the
        version's energy-parameter set, so a drift here silently rewrites every RPI
        score."""
        self.assertEqual(self.conda["viennarna"], ("==", "2.7.2"))
        self.assertEqual(self.pip["viennarna"], ("==", "2.7.2"))

    def test_mafft_and_iqtree_are_not_claimed_as_pip_installable(self):
        """`pymafft` on PyPI is an unrelated package. Listing either here would send the
        next reader to install the wrong thing."""
        for absent in ("mafft", "iqtree", "pymafft"):
            self.assertNotIn(absent, self.pip)

    def test_installed_versions_satisfy_the_pip_pins(self):
        """The files agree with each other above; this checks they also agree with the
        interpreter actually running the suite. Without it both files could be
        internally consistent and jointly wrong."""
        import importlib
        for package, module_name in (("numpy", "numpy"), ("pandas", "pandas"),
                                     ("biopython", "Bio"), ("dendropy", "dendropy"),
                                     ("viennarna", "RNA")):
            with self.subTest(package=package):
                module = importlib.import_module(module_name)
                self.assertEqual(_parse_version(module.__version__),
                                 _parse_version(self.pip[package][1]),
                                 f"{package}: installed {module.__version__}, "
                                 f"requirements-pip.txt pins {self.pip[package][1]}")


class JsonSafeNumpyScalarTest(unittest.TestCase):
    """_json_safe tested only Python's `float`. np.float32 is NOT a float subclass --
    only np.float64 is -- so every other numpy scalar fell straight through it, and
    json's default=str then intercepted them before allow_nan=False could fire."""

    def test_numpy_float32_nan_becomes_null_not_the_string_nan(self):
        """The recorded half of the hole: np.float32('nan') serialised as the JSON
        STRING "nan" -- a value, where null was meant."""
        payload = _json_safe({"auc": np.float32("nan")})
        self.assertIsNone(payload["auc"])
        self.assertEqual(json.dumps(payload, allow_nan=False), '{"auc": null}')

    def test_finite_numpy_scalars_survive_as_numbers_not_strings(self):
        """The unrecorded half, and it is wider: the same type hole stringified FINITE
        numpy scalars too. np.float32(1.5) serialised as "1.5" and np.int64(3) as "3",
        so a downstream consumer reading the report got strings where the schema
        promises numbers."""
        payload = _json_safe({"f32": np.float32(1.5), "f64": np.float64(2.5),
                              "i64": np.int64(3)})
        self.assertEqual(payload["f32"], 1.5)
        self.assertEqual(payload["f64"], 2.5)
        self.assertEqual(payload["i64"], 3)
        self.assertEqual(json.dumps(payload, sort_keys=True, allow_nan=False),
                         '{"f32": 1.5, "f64": 2.5, "i64": 3}')

    def test_numpy_non_finite_in_nested_containers_is_reached(self):
        """The recursion has to see them too -- rank_consistency returns its nan inside
        a nested dict, not at the top level."""
        payload = _json_safe({"ci": [np.float32("nan"), np.float64("inf")],
                              "inner": {"x": np.float32("-inf")}})
        self.assertEqual(payload["ci"], [None, None])
        self.assertIsNone(payload["inner"]["x"])

    def test_python_types_and_bool_are_unchanged(self):
        """Guard the guard: bool is an int subclass, so a careless numbers.Real check
        would silently turn True into 1 and change the report's schema.
        rpi_k_stratum_provisional is a bool in this very report."""
        payload = _json_safe({"b": True, "nb": np.bool_(False), "i": 3, "s": "x",
                              "n": None, "f": 1.5, "nan": float("nan")})
        self.assertIs(payload["b"], True)
        self.assertIs(payload["nb"], False)
        self.assertEqual(payload["i"], 3)
        self.assertNotIsInstance(payload["i"], bool)
        self.assertEqual(payload["s"], "x")
        self.assertIsNone(payload["n"])
        self.assertEqual(payload["f"], 1.5)
        self.assertIsNone(payload["nan"])


def _reject_json_constant(token: str) -> Any:
    """Make json.loads refuse the bare NaN/Infinity tokens it accepts by default.

    Python's json.loads accepts `NaN`, `Infinity` and `-Infinity`, none of which are
    legal JSON -- jq and JavaScript's JSON.parse both reject a file containing them.
    That leniency is precisely why a round-trip test could not see a non-conforming
    report. Used by the artifact-level assertions below.
    """
    raise AssertionError(f"report contains non-JSON constant {token!r}")


class WriteRpiReportTest(unittest.TestCase):
    """Spec 14 asks for these validation entry points to be RUN, not merely to exist.
    Before --rpi-report they were unreachable from the CLI and wrote no artifact."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.output_dir = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        self.rows = [
            {"candidate_id": "c1", "gate3_status": "ecology_plausible",
             "rpi_e_total": 0.5, "rpi_generation": 0.4,
             "verdict": "supported_reassortment_candidate"},
            {"candidate_id": "c2", "gate3_status": "no_ecological_route",
             "rpi_e_total": 0.001, "rpi_generation": 0.001, "verdict": "no_call"},
        ]

    def test_report_is_written_and_carries_every_spec_14_section(self):
        path = write_rpi_report(self.rows, self.output_dir)
        self.assertEqual(path.name, RPI_REPORT_FILENAME)
        report = json.loads(path.read_text(encoding="utf-8"))
        for section in ("retrospective_reproduction", "leave_one_pair_out",
                        "rank_consistency", "veto_vs_continuous_e", "sensitivity_sweep"):
            self.assertIn(section, report, section)

    def test_leave_one_pair_out_runs_once_per_contributing_pair_only(self):
        """A pair outside the weight table cannot be held out of it. Including such a
        pair here would report a leave-one-out run that held nothing out -- the exact
        defect the count subtraction was added to fix."""
        report = json.loads(
            write_rpi_report(self.rows, self.output_dir).read_text(encoding="utf-8"))
        expected = {name for name, record in RPI_INVITRO_BENCHMARK.items()
                    if record.get("contributes_to_prior")}
        self.assertEqual(set(report["leave_one_pair_out"]), expected)
        for name, run in report["leave_one_pair_out"].items():
            self.assertTrue(run["holdout_changed_the_prior"], name)

    def test_report_is_strictly_valid_json_even_when_statistics_are_undefined(self):
        """Python emits bare `NaN` tokens by default, which jq and JavaScript's
        JSON.parse both reject -- while Python's own json.loads accepts them, which is
        why a round-trip check does not catch it. nan is routine here, not exotic: a
        single-candidate run (the documented default invocation) leaves rank_consistency
        with no negatives, and it deliberately returns nan for auc and its interval."""
        one_row = [{"candidate_id": "solo", "gate3_status": "ecology_plausible",
                    "rpi_e_total": 0.5, "rpi_generation": 0.4,
                    "verdict": "supported_reassortment_candidate"}]
        raw = write_rpi_report(one_row, self.output_dir).read_text(encoding="utf-8")
        self.assertNotIn("NaN", raw)
        self.assertNotIn("Infinity", raw)
        report = json.loads(raw)
        # The undefined statistic must survive as null -- not as 0.0, which is a claim.
        self.assertIsNone(report["rank_consistency"]["auc"])

    def test_written_report_carries_no_stringified_numerics(self):
        """Assert on the row that ships, not on the function that feeds it.

        The two checks above look at _json_safe's return and at parsed JSON; neither can
        see a value that reached the FILE as a quoted string. json.loads also accepts
        the bare NaN/Infinity tokens that are not legal JSON, so parsing alone cannot
        certify the artifact -- parse_constant closes that."""
        raw = write_rpi_report(self.rows, self.output_dir).read_text(encoding="utf-8")
        for illegal in ('"nan"', '"inf"', '"-inf"', '"NaN"', "NaN", "Infinity"):
            self.assertNotIn(illegal, raw, illegal)
        json.loads(raw, parse_constant=_reject_json_constant)

    def test_cli_flag_is_off_by_default(self):
        """Existing invocations must stay byte-identical."""
        args = parse_args(["--reference-input", "all_hanta", "--output", "output",
                           "--candidate", "CAND01:s.fa,m.fa,l.fa"])
        self.assertFalse(args.rpi_report)
        args = parse_args(["--reference-input", "all_hanta", "--output", "output",
                           "--candidate", "CAND01:s.fa,m.fa,l.fa", "--rpi-report"])
        self.assertTrue(args.rpi_report)


class SensitivityAndRankTest(unittest.TestCase):
    def test_sweep_covers_every_assumed_parameter(self):
        """Spec 12: decay weights AND decay form, smoothing constant, retention ratio,
        stratum boundary."""
        report = sensitivity_sweep()
        for key in ("host_decay", "jeffreys_alpha", "retention_by_odd_segment",
                    "stratum_identity_boundary"):
            self.assertIn(key, report["varied"], key)

    def test_sweep_reports_rank_correlation_against_baseline(self):
        report = sensitivity_sweep()
        for entry in report["varied"].values():
            self.assertIn("rank_correlation", entry)

    def test_sweep_does_not_mutate_global_weights(self):
        before = copy.deepcopy(RPI_WEIGHTS)
        sensitivity_sweep()
        self.assertEqual(RPI_WEIGHTS, before)

    def test_rank_consistency_reports_a_confidence_interval(self):
        """Spec 11.2: with 8 positives this is a SANITY CHECK, not validation. It must
        never be reported as a bare pass/fail."""
        rows = [{"verdict": "no_conflict", "rpi_generation": 0.1} for _ in range(200)]
        rows += [{"verdict": "supported_reassortment_candidate",
                  "rpi_generation": 0.6} for _ in range(8)]
        report = rank_consistency(rows)
        self.assertIn("confidence_interval", report)
        self.assertEqual(report["n_positive"], 8)
        self.assertIn("sanity check", report["interpretation"].lower())

    def test_rank_consistency_warns_when_positives_are_few(self):
        rows = [{"verdict": "no_conflict", "rpi_generation": 0.1} for _ in range(200)]
        rows += [{"verdict": "supported_reassortment_candidate", "rpi_generation": 0.6}]
        self.assertTrue(rank_consistency(rows)["underpowered"])

    def test_sweep_reports_genuine_reordering(self):
        """host_decay and retention_by_odd_segment must NOT be a tautology: comparing
        RPI_WEIGHTS[key].values() to a positive scalar multiple of itself is provably
        1.0 and can never report instability. Correlating rpi_generation /
        rpi_retention on the fixed synthetic panel instead can, and does, on at least
        one perturbation direction for every parameter that actually varies scores."""
        report = sensitivity_sweep()
        for key in ("host_decay", "retention_by_odd_segment", "stratum_identity_boundary",
                    "jeffreys_alpha"):
            self.assertLess(report["varied"][key]["rank_correlation"], 1.0, key)

    def test_jeffreys_alpha_is_swept_through_the_panel_not_against_itself(self):
        """The original alpha entry Spearman-correlated [close_M_prior, distant_M_prior]
        against itself under perturbation. distant-M exceeds close-M at EVERY alpha, so
        the rank order was fixed by construction and the statistic returned 1.0
        unconditionally -- a tautology, the same defect host_decay and
        retention_by_odd_segment were already fixed for. alpha is reachable via
        RPI_WEIGHTS['jeffreys_alpha'], so it routes through _score_panel like the other
        three and the answer is now about real candidate reordering."""
        report = sensitivity_sweep()
        entry = report["varied"]["jeffreys_alpha"]
        self.assertIn("panel", entry["note"])
        # A two-element self-comparison can only ever return exactly 1.0 or -1.0; a
        # value strictly between them proves the statistic came from the 12-candidate
        # panel instead.
        self.assertLess(entry["rank_correlation"], 1.0)
        self.assertGreater(entry["rank_correlation"], -1.0)

    def test_sweep_marks_its_panel_as_synthetic(self):
        """Without a top-level marker the sweep's numbers can be quoted as a statement
        about real candidates. They are not: the panel is twelve constructed rows."""
        self.assertEqual(sensitivity_sweep()["panel"], "synthetic_fixed")

    def test_host_decay_sweep_never_turns_missing_data_into_a_bonus(self):
        """unknown_host_metadata is a NEUTRAL SENTINEL (1.0 = 'unknown, not separated'),
        not an estimated decay weight. Scaling it above 1.0 made missing host taxonomy a
        positive bonus, and scaling every weight above 1.0 saturated compute_rpi's
        min(1.0, ...) so the whole upward direction became a no-op for the candidates it
        should have moved."""
        high = _perturbed_host_decay(1.5)
        low = _perturbed_host_decay(0.5)
        self.assertEqual(high["unknown_host_metadata"], RPI_WEIGHTS["host_decay"]
                         ["unknown_host_metadata"])
        self.assertEqual(low["unknown_host_metadata"], RPI_WEIGHTS["host_decay"]
                         ["unknown_host_metadata"])
        for perturbed in (high, low):
            for key, value in perturbed.items():
                self.assertLessEqual(value, 1.0, key)
                self.assertGreater(value, 0.0, key)
        # The upward direction must still move the weights it is allowed to move.
        self.assertGreater(high["same_genus"], RPI_WEIGHTS["host_decay"]["same_genus"])
        self.assertEqual(high["same_species"], 1.0)  # clamped, not 1.5

    def test_host_decay_sweep_reports_what_it_held_fixed(self):
        entry = sensitivity_sweep()["varied"]["host_decay"]
        self.assertIn("unknown_host_metadata", entry["excluded_from_perturbation"])

    def test_spearman_is_nan_where_the_statistic_is_undefined(self):
        """Returning 1.0 for a degenerate input reports PERFECT RANK STABILITY from an
        input that measured nothing -- the same class of defect as a constant factor
        dressed as a measurement. nan matches rank_consistency's own discipline for its
        undefined AUC."""
        self.assertTrue(math.isnan(_spearman([1.0], [2.0])))
        self.assertTrue(math.isnan(_spearman([], [])))
        self.assertTrue(math.isnan(_spearman([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])))
        self.assertTrue(math.isnan(_spearman([1.0, 2.0, 3.0], [5.0, 5.0, 5.0])))
        self.assertAlmostEqual(_spearman([1.0, 2.0, 3.0], [10.0, 20.0, 30.0]), 1.0)

    def test_ties_share_a_midrank_rather_than_being_ordered_by_position(self):
        """np.argsort(np.argsort(v)) -- the idiom this used to use -- breaks ties by
        POSITION, so a constant vector ranks [0, 1, 2]: a perfect ordering manufactured
        out of an input with no order in it. That defeated the zero-variance guard in
        _spearman (the guard inspected ranks that were never constant) and let any tie
        in a real score panel contribute a spurious ordering to the correlation."""
        self.assertEqual(list(_tied_ranks([5.0, 5.0, 5.0])), [1.0, 1.0, 1.0])
        self.assertEqual(list(_tied_ranks([3.0, 1.0, 3.0, 1.0])), [2.5, 0.5, 2.5, 0.5])
        self.assertEqual(list(_tied_ranks([10.0, 20.0, 30.0])), [0.0, 1.0, 2.0])
        # Two candidates tied on one axis and separated on the other: the tie must not
        # be silently resolved into agreement with whichever order the other axis has.
        self.assertLess(_spearman([1.0, 1.0, 2.0, 3.0], [1.0, 2.0, 3.0, 4.0]), 1.0)

    def test_undefined_direction_is_surfaced_not_folded_away(self):
        """builtin min() is order-dependent around nan (min(nan, 0.5) is nan;
        min(0.5, nan) is 0.5), so an undefined perturbation direction could be silently
        discarded depending on argument order. If a direction is undefined, the reported
        answer is undefined."""
        self.assertTrue(math.isnan(_worst_rank_correlation(float("nan"), 0.5)))
        self.assertTrue(math.isnan(_worst_rank_correlation(0.5, float("nan"))))
        self.assertEqual(_worst_rank_correlation(0.9, 0.4), 0.4)

    def test_stratum_identity_boundary_reports_which_candidates_move(self):
        """Spec: reclassification is a step change, not a smooth one -- the report must
        name which pairs move, not assert stability by fiat."""
        report = sensitivity_sweep()
        entry = report["varied"]["stratum_identity_boundary"]
        self.assertIn("moved_low", entry)
        self.assertIn("moved_high", entry)
        self.assertTrue(entry["moved_low"] or entry["moved_high"])

    def test_rank_consistency_interval_has_nonzero_width_at_perfect_separation(self):
        """A naive Wald interval (z*sqrt(p(1-p)/n)) is exactly 0 width at auc=1.0 --
        false certainty from 8 positives. Wilson must not collapse there."""
        rows = [{"verdict": "no_conflict", "rpi_generation": 0.1} for _ in range(200)]
        rows += [{"verdict": "supported_reassortment_candidate",
                  "rpi_generation": 0.6} for _ in range(8)]
        report = rank_consistency(rows)
        self.assertEqual(report["auc"], 1.0)
        lower, upper = report["confidence_interval"]
        self.assertGreater(upper, lower)

    def test_rank_consistency_underpowered_at_production_dataset_size(self):
        """The production pipeline has exactly 8 supported_reassortment_candidate
        labels. The threshold must fire unambiguously at that exact count, not sit at
        a boundary equal to it by coincidence."""
        rows = [{"verdict": "no_conflict", "rpi_generation": 0.1} for _ in range(200)]
        rows += [{"verdict": "supported_reassortment_candidate",
                  "rpi_generation": 0.6} for _ in range(8)]
        self.assertTrue(rank_consistency(rows)["underpowered"])

    def test_rank_consistency_confidence_interval_undefined_without_both_classes(self):
        """No positives (or no negatives) means the AUC is genuinely undefined -- the
        interval must stay nan/nan, not silently collapse to a full-width (0.0, 1.0)
        via max(0.0, nan) / min(1.0, nan)."""
        rows = [{"verdict": "no_conflict", "rpi_generation": 0.1} for _ in range(200)]
        report = rank_consistency(rows)
        self.assertTrue(math.isnan(report["auc"]))
        lower, upper = report["confidence_interval"]
        self.assertTrue(math.isnan(lower))
        self.assertTrue(math.isnan(upper))
