# detector_reference.py
"""Reference dataset loading for the reassortment candidate detector.

load_reference_data replicates the ingest portion of all_hanta.py's run_pipeline (see
that function, lines ~3960-4001) without any of the table-building, ML-matrix, or
output-writing steps this detector doesn't need.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import all_hanta as ah
import phylogeny as ph


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
