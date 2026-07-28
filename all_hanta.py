#!/usr/bin/env python3
"""Orthohantavirus NCBI Datasets production pipeline.

The script ingests an NCBI Datasets package for Orthohantavirus and produces:
- accession-level records
- virus-level group records
- isolate-level records
- host-level records
- ML-ready feature tables
- QC and suppression tables
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, DefaultDict, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import pandas as pd

__version__ = "1.1.1"

ACCESSION_RE = re.compile(r"[A-Z]{1,4}_?\d{4,}(?:\.\d+)?")
BIOSAMPLE_RE = re.compile(r"\bSAM[NED]?[A-Z]?\d+\b", re.I)
YEAR_RE = re.compile(r"(19|20)\d{2}")

COMPLETE_LENGTH_MIN = {"L": 6000, "M": 3200, "S": 1200}

OUTPUT_TABLE_SCHEMAS = {
    "ortho_accession_table.parquet": [
        "accession", "species_name", "species_taxid", "genus_name", "genus_taxid", "segment", "segment_source",
        "completeness", "nt_length", "biosample", "bioproject", "isolate", "strain", "collection_date", "year",
        "location", "country", "region", "host_name", "host_taxid", "host_source", "host_role", "isolation_source",
        "tissue_source", "protein_accessions", "annotation_products", "group_key", "group_confidence", "is_refseq",
        "warning_flags",
    ],
    "ortho_virus_table.parquet": [
        "group_key", "species_name", "species_taxid", "host_name", "host_taxid", "country", "year", "biosample",
        "isolate", "strain", "group_confidence", "group_size", "s_accessions", "m_accessions", "l_accessions",
        "complete_segments", "completeness_class", "all_hosts", "all_accessions",
    ],
    "ortho_isolate_table.parquet": [
        "strict_isolate_key", "species_name", "species_taxid", "species_cluster_id", "biosample", "isolate", "strain",
        "host_name", "host_taxid", "host_role", "collection_date", "location", "country", "region", "grouping_tier",
        "grouping_confidence", "s_accession", "m_accession", "l_accession", "s_completeness", "m_completeness",
        "l_completeness", "s_nt_length", "m_nt_length", "l_nt_length", "concatenated_protein_aa",
        "coding_complete_status", "mosaic_flag", "taxonomic_incongruence_warning", "reassortment_candidate_warning",
        "ns_s_overlapping_orf", "species_divergence_flag", "s_protein_evidence", "m_protein_evidence",
        "l_protein_evidence", "s_protein_accessions", "m_protein_accessions", "l_protein_accessions",
        "missing_protein_warnings", "warnings",
    ],
    "ortho_host_table.parquet": [
        "host_name", "host_taxid", "host_role", "host_source", "representative_biosample", "representative_isolate",
        "accession_count", "isolate_count", "example_species",
    ],
    "ortho_ml_feature_matrix.parquet": [
        "strict_isolate_key", "species_name", "species_taxid", "host_name", "host_taxid", "host_role", "kmer_*",
    ],
    "ortho_ml_segments.parquet": [
        "accession", "segment", "species_name", "species_taxid", "host_name", "host_taxid", "country", "year",
        "biosample", "completeness", "group_key",
    ],
    "ortho_suppressed.parquet": [
        "accession", "species_name", "species_taxid", "segment", "completeness", "biosample", "isolate", "strain",
        "collection_date", "location", "host_name", "host_taxid", "suppression_reason",
    ],
    "ortho_qc_report.parquet": [
        "accession", "group_key", "grouping_tier", "group_confidence", "segment", "completeness", "warnings",
    ],
}

SEGMENT_TEXT_PATTERNS = {
    "L": [r"\bsegment\s*l\b", r"\blarge\b", r"\bpolymerase\b", r"\brdrp\b", r"rna[- ]dependent rna polymerase", r"\bl protein\b"],
    "M": [r"\bsegment\s*m\b", r"\bmedium\b", r"glycoprotein precursor", r"\bglycoprotein\b", r"\bgpc\b", r"\bgn\b", r"\bgc\b", r"\bg1\b", r"\bg2\b"],
    "S": [r"\bsegment\s*s\b", r"\bsmall\b", r"\bnucleocapsid\b", r"\bnucleoprotein\b", r"\bn protein\b", r"\bnp\b"],
}

LAB_HOST_TERMS = ["vero", "cell line", "cell-line", "cell culture", "passage", "passaged", "adapted", "in vitro"]
PROTEIN_PRODUCT_PATTERNS = {
    "N": [r"\bnucleocapsid\b", r"\bnucleoprotein\b", r"\bnp\b", r"\bn protein\b"],
    "GPC": [r"\bglycoprotein precursor\b", r"\bglycoprotein\b", r"\bgpc\b", r"\bgn\b", r"\bgc\b", r"\bg1\b", r"\bg2\b"],
    "NSs": [r"\bns s\b", r"\bnon[- ]structural protein s\b", r"\bns s protein\b", r"\bns s orf\b", r"\bns s overlapping\b", r"\bns s overlap\b"],
    "RdRp": [r"\bpolymerase\b", r"\brdrp\b", r"rna[- ]dependent rna polymerase", r"\bl protein\b", r"\blarge protein\b"],
}

DATA_REPORT = "data_report.jsonl"
ANNOTATION_REPORT = "annotation_report.jsonl"
BIOSAMPLE_REPORT = "biosample_report.jsonl"
GENOMIC_FASTA = "genomic.fna"
CDS_FASTA = "cds.fna"
PROTEIN_FASTA = "protein.faa"
KMER_K = 3
PED_SPECIES_THRESHOLD = 0.1
N_AMBIGUOUS_THRESHOLD = 100
TERMINAL_MOTIF_5 = "UAGUAGUAG"
TERMINAL_MOTIF_3 = "AUCAUCAUCUG"


@dataclass
class InputPaths:
    root: Path
    data_dir: Path
    data_report: Path
    annotation_report: Path
    biosample_report: Path
    dataset_catalog: Path
    genomic_fasta: Optional[Path]
    cds_fasta: Optional[Path]
    protein_fasta: Optional[Path]
    output_dir: Path
    log_file: Path


@dataclass
class FastaRecord:
    accession: str
    accession_base: str
    header: str
    sequence: str
    length: int


@dataclass
class AnnotationEvidence:
    accession_links: List[str] = field(default_factory=list)
    product_names: List[str] = field(default_factory=list)
    protein_accessions: List[str] = field(default_factory=list)
    raw_hits: int = 0


@dataclass
class BioSampleInfo:
    accession: str = ""
    host: Optional[str] = None
    host_taxid: Optional[str] = None
    isolation_source: Optional[str] = None
    tissue: Optional[str] = None
    collection_date: Optional[str] = None
    geographic_location: Optional[str] = None
    strain: Optional[str] = None
    isolate: Optional[str] = None
    raw_attributes: Dict[str, str] = field(default_factory=dict)


@dataclass
class AccessionRecord:
    accession: str
    accession_base: str
    species_name: str
    species_taxid: str
    genus_name: str
    genus_taxid: str
    segment: str
    segment_source: str
    completeness: str
    nt_length: int
    biosample: str
    bioproject: str
    isolate: str
    strain: str
    collection_date: str
    year: str
    location: str
    country: str
    region: str
    host_name: str
    host_taxid: str
    host_source: str
    host_role: str
    isolation_source: str
    tissue_source: str
    protein_accessions: List[str]
    annotation_products: List[str]
    protein_aa: str
    protein_product: str
    nss_overlap_flag: bool
    strict_isolate_key: str
    group_key: str
    grouping_tier: str
    group_confidence: float
    is_refseq: bool
    warning_flags: List[str]


@dataclass
class VirusGroupRecord:
    group_key: str
    species_name: str
    species_taxid: str
    host_name: str
    host_taxid: str
    country: str
    year: str
    biosample: str
    isolate: str
    strain: str
    group_confidence: float
    group_size: int
    s_accessions: List[str]
    m_accessions: List[str]
    l_accessions: List[str]
    complete_segments: List[str]
    completeness_class: str
    all_hosts: List[str]
    all_accessions: List[str]


@dataclass
class IsolateRecord:
    strict_isolate_key: str
    species_name: str
    species_taxid: str
    species_cluster_id: str
    biosample: str
    isolate: str
    strain: str
    host_name: str
    host_taxid: str
    host_role: str
    collection_date: str
    location: str
    country: str
    region: str
    grouping_tier: str
    grouping_confidence: float
    s_accession: str
    m_accession: str
    l_accession: str
    s_completeness: str
    m_completeness: str
    l_completeness: str
    s_nt_length: int
    m_nt_length: int
    l_nt_length: int
    concatenated_protein_aa: str
    coding_complete_status: str
    mosaic_flag: bool
    taxonomic_incongruence_warning: bool
    reassortment_candidate_warning: bool
    ns_s_overlapping_orf: bool
    species_divergence_flag: bool
    s_protein_evidence: str
    m_protein_evidence: str
    l_protein_evidence: str
    s_protein_accessions: List[str]
    m_protein_accessions: List[str]
    l_protein_accessions: List[str]
    missing_protein_warnings: List[str]
    host_source: str
    raw_metadata: Dict[str, Any]
    warnings: List[str]


@dataclass
class HostRecord:
    host_name: str
    host_taxid: str
    host_role: str
    host_source: str
    representative_biosample: str
    representative_isolate: str
    accession_count: int
    isolate_count: int
    example_species: str


def configure_logging(log_file: Path, verbose: bool = False) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_file, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
    )


def clean_text(value: Any) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, (dict, list)):
        text = json.dumps(value, ensure_ascii=False)
    else:
        text = str(value).strip()
    if not text:
        return None
    low = text.lower()
    if low in {"none", "null", "nan", "na", "n/a", "unknown", "missing", "not collected"}:
        return None
    return re.sub(r"\s+", " ", text).strip()


def first_nonempty(*values: Any) -> Optional[str]:
    for value in values:
        if isinstance(value, list):
            for item in value:
                x = clean_text(item)
                if x:
                    return x
        else:
            x = clean_text(value)
            if x:
                return x
    return None


def uniq(values: Iterable[Any]) -> List[Any]:
    seen = set()
    out = []
    for value in values:
        key = json.dumps(value, ensure_ascii=False, sort_keys=True) if isinstance(value, (dict, list, tuple, set)) else str(value)
        if key not in seen:
            seen.add(key)
            out.append(value)
    return out


def accession_base(accession: str) -> str:
    return accession.split(".")[0] if accession else ""


def flatten_json(obj: Any, prefix: str = "") -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            p = f"{prefix}.{k}" if prefix else str(k)
            out.update(flatten_json(v, p))
    elif isinstance(obj, list):
        if all(not isinstance(x, (dict, list)) for x in obj):
            out[prefix] = obj
        else:
            out[prefix] = json.dumps(obj, ensure_ascii=False)
    else:
        out[prefix] = obj
    return out


def iter_jsonl(path: Path) -> Iterator[Dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"JSON decode error in {path.name} line {line_no}: {exc}") from exc
            if isinstance(obj, dict):
                yield obj


def md5sum(path: Path) -> str:
    if not path.exists():
        return ""
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_dataset_catalog(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def extract_year(text: Optional[str]) -> str:
    if not text:
        return ""
    m = YEAR_RE.search(text)
    return m.group(0) if m else ""


def split_geo(geo: str) -> Tuple[str, str]:
    if not geo:
        return "", ""
    if ":" in geo:
        country, region = geo.split(":", 1)
        return country.strip(), region.strip()
    return geo.strip(), ""


def safe_write_parquet(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        pd.DataFrame().to_parquet(path, engine="pyarrow", index=False)
        return
    pd.DataFrame(rows).to_parquet(path, engine="pyarrow", index=False)
    logging.info("Wrote %s rows to %s", len(rows), path)


def safe_write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def resolve_inputs(root: Path, output: Optional[Path] = None) -> InputPaths:
    root = root.expanduser().resolve()
    data_dir = root / "all_hanta" / "ncbi_dataset" / "data"
    if not data_dir.exists():
        raise FileNotFoundError(f"Could not find data directory: {data_dir}")

    data_report = data_dir / DATA_REPORT
    annotation_report = data_dir / ANNOTATION_REPORT
    biosample_report = data_dir / BIOSAMPLE_REPORT
    dataset_catalog = data_dir / "dataset_catalog.json"

    for required in [data_report, annotation_report, biosample_report, dataset_catalog]:
        if not required.exists():
            raise FileNotFoundError(f"Missing required file: {required}")

    output_dir = (output.expanduser().resolve() if output else root / "output")
    output_dir.mkdir(parents=True, exist_ok=True)

    return InputPaths(
        root=root,
        data_dir=data_dir,
        data_report=data_report,
        annotation_report=annotation_report,
        biosample_report=biosample_report,
        dataset_catalog=dataset_catalog,
        genomic_fasta=(data_dir / GENOMIC_FASTA) if (data_dir / GENOMIC_FASTA).exists() else None,
        cds_fasta=(data_dir / CDS_FASTA) if (data_dir / CDS_FASTA).exists() else None,
        protein_fasta=(data_dir / PROTEIN_FASTA) if (data_dir / PROTEIN_FASTA).exists() else None,
        output_dir=output_dir,
        log_file=output_dir / "ortho_pipeline.log",
    )


def parse_fasta(path: Optional[Path]) -> Dict[str, FastaRecord]:
    if path is None or not path.exists():
        return {}
    index: Dict[str, FastaRecord] = {}
    header: Optional[str] = None
    parts: List[str] = []

    def commit(raw_header: str, seq_parts: List[str]) -> None:
        desc = raw_header[1:].strip()
        acc = desc.split()[0]
        base = accession_base(acc)
        seq = "".join(seq_parts)
        rec = FastaRecord(accession=acc, accession_base=base, header=desc, sequence=seq, length=len(seq))
        index[acc] = rec
        if base and base not in index:
            index[base] = rec

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    commit(header, parts)
                header = line
                parts = []
            else:
                parts.append(line.strip())
        if header is not None:
            commit(header, parts)
    return index


def extract_accessions_from_blob(obj: Any) -> List[str]:
    return uniq(ACCESSION_RE.findall(json.dumps(obj, ensure_ascii=False)))


def build_annotation_index(records: Sequence[Dict[str, Any]]) -> DefaultDict[str, List[Dict[str, Any]]]:
    idx: DefaultDict[str, List[Dict[str, Any]]] = defaultdict(list)
    for rec in records:
        for acc in extract_accessions_from_blob(rec):
            idx[acc].append(rec)
            base = accession_base(acc)
            if base:
                idx[base].append(rec)
    return idx


def extract_biosample_candidates(obj: Any) -> List[str]:
    found: List[str] = []

    def walk(value: Any) -> None:
        if isinstance(value, dict):
            for k, v in value.items():
                if str(k).lower() in {"biosample", "biosampleaccession", "biosample_acc"} and isinstance(v, str):
                    m = BIOSAMPLE_RE.search(v)
                    if m:
                        found.append(m.group(0).upper())
                walk(v)
        elif isinstance(value, list):
            for item in value:
                walk(item)
        elif isinstance(value, str):
            m = BIOSAMPLE_RE.search(value)
            if m:
                found.append(m.group(0).upper())

    walk(obj)
    return uniq(found)


def normalize_attr_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")


def parse_biosample_record(rec: Dict[str, Any]) -> BioSampleInfo:
    flat = flatten_json(rec)
    accession = first_nonempty(rec.get("accession"), flat.get("accession")) or ""
    info = BioSampleInfo(accession=accession)
    attrs = rec.get("attributes", [])
    attr_map: Dict[str, str] = {}
    if isinstance(attrs, list):
        for item in attrs:
            if not isinstance(item, dict):
                continue
            name = normalize_attr_name(str(item.get("name", "")))
            value = clean_text(item.get("value"))
            if name and value:
                attr_map[name] = value
    info.raw_attributes = attr_map
    info.host = first_nonempty(attr_map.get("host"), flat.get("host"))
    info.host_taxid = first_nonempty(attr_map.get("host_taxid"), flat.get("host.taxId"), flat.get("host.taxid"))
    info.isolation_source = first_nonempty(attr_map.get("isolation_source"), attr_map.get("source"))
    info.tissue = first_nonempty(attr_map.get("tissue"), attr_map.get("host_tissue_sampled"))
    info.collection_date = first_nonempty(attr_map.get("collection_date"), flat.get("collectionDate"))
    info.geographic_location = first_nonempty(attr_map.get("geo_loc_name"), flat.get("geoLocName"), flat.get("geographicLocation"))
    info.strain = first_nonempty(attr_map.get("strain"))
    info.isolate = first_nonempty(attr_map.get("isolate"))
    return info


def build_biosample_index(records: Sequence[Dict[str, Any]]) -> Dict[str, BioSampleInfo]:
    out: Dict[str, BioSampleInfo] = {}
    for rec in records:
        info = parse_biosample_record(rec)
        for candidate in extract_biosample_candidates(rec) + ([info.accession] if info.accession else []):
            if candidate:
                out[candidate.upper()] = info
    return out


def extract_annotation_evidence(annotation_records: Sequence[Dict[str, Any]], protein_index: Dict[str, FastaRecord]) -> AnnotationEvidence:
    evidence = AnnotationEvidence(raw_hits=len(annotation_records))
    products: List[str] = []
    links: List[str] = []
    prot_accs: List[str] = []
    for rec in annotation_records:
        links.extend(extract_accessions_from_blob(rec))
        blob = json.dumps(rec, ensure_ascii=False)
        for product_match in re.findall(r'"name"\s*:\s*"([^"]+)"', blob):
            products.append(product_match)
        for acc in extract_accessions_from_blob(rec):
            if acc in protein_index or accession_base(acc) in protein_index:
                prot_accs.append(acc)
    evidence.accession_links = uniq(links)
    evidence.product_names = uniq(products)
    evidence.protein_accessions = uniq(prot_accs)
    return evidence


def segment_from_explicit_fields(flat_contexts: Sequence[Dict[str, Any]]) -> Optional[str]:
    for flat in flat_contexts:
        value = first_nonempty(flat.get("segment"), flat.get("virus.segment"), flat.get("annotation.segment"))
        if not value:
            continue
        v = value.strip().upper()
        if v in {"S", "M", "L"}:
            return v
        if v == "SMALL":
            return "S"
        if v == "MEDIUM":
            return "M"
        if v == "LARGE":
            return "L"
    return None


def segment_from_text_evidence(texts: Sequence[str]) -> Tuple[Optional[str], List[str]]:
    warnings: List[str] = []
    scores = {"S": 0, "M": 0, "L": 0}
    for text in texts:
        low = text.lower()
        if "nss" in low or "ns s" in low:
            warnings.append("contains_nss_label")
        for seg, patterns in SEGMENT_TEXT_PATTERNS.items():
            for pat in patterns:
                if re.search(pat, low, flags=re.I):
                    scores[seg] += 1
    best = max(scores, key=scores.get)
    if scores[best] == 0:
        return None, uniq(warnings)
    non_zero = [k for k, v in scores.items() if v > 0]
    if len(non_zero) > 1:
        warnings.append("multiple_segment_text_hits")
    return best, uniq(warnings)


def choose_segment(flat_contexts: Sequence[Dict[str, Any]], annotation: AnnotationEvidence, nt_len: int) -> Tuple[str, str, List[str]]:
    explicit = segment_from_explicit_fields(flat_contexts)
    if explicit:
        return explicit, "explicit_metadata", []
    texts: List[str] = []
    for flat in flat_contexts:
        for k, v in flat.items():
            lk = k.lower()
            if any(token in lk for token in ["segment", "description", "definition", "title", "product", "protein"]):
                if isinstance(v, list):
                    texts.extend([str(x) for x in v if x is not None])
                else:
                    texts.append(str(v))
    texts.extend(annotation.product_names)
    seg, warnings = segment_from_text_evidence(texts)
    if seg:
        return seg, "annotation_text", warnings
    if nt_len:
        if nt_len >= COMPLETE_LENGTH_MIN["L"]:
            return "L", "length_heuristic", ["segment_from_length_only"]
        if nt_len >= COMPLETE_LENGTH_MIN["M"]:
            return "M", "length_heuristic", ["segment_from_length_only"]
        if nt_len >= COMPLETE_LENGTH_MIN["S"]:
            return "S", "length_heuristic", ["segment_from_length_only"]
    return "", "unresolved", ["segment_unresolved"]


def count_ambiguous_n(sequence: str) -> int:
    return sequence.upper().count("N") if sequence else 0


def has_terminal_motif(sequence: str) -> bool:
    if not sequence:
        return False
    seq = sequence.upper().replace("T", "U")
    if seq.startswith(TERMINAL_MOTIF_5) and seq.endswith(TERMINAL_MOTIF_3):
        return True
    if seq.startswith(TERMINAL_MOTIF_5) or seq.endswith(TERMINAL_MOTIF_3):
        return True
    return sequence.upper().startswith(TERMINAL_MOTIF_5.replace("U", "T")) or sequence.upper().endswith(TERMINAL_MOTIF_3.replace("U", "T"))


def determine_completeness(flat_contexts: Sequence[Dict[str, Any]], segment: str, nt_len: int, raw_blob: str, nt_sequence: Optional[str] = None) -> str:
    structured = None
    for flat in flat_contexts:
        structured = first_nonempty(flat.get("completeness"), flat.get("nucleotideCompleteness"), flat.get("nucleotide.completeness"))
        if structured:
            break
    if structured:
        low = structured.lower()
        if low == "complete":
            return "complete"
        if "partial" in low:
            return "partial"
    low_blob = raw_blob.lower()
    if "complete genome" in low_blob or "complete sequence" in low_blob:
        return "complete"
    if "complete cds" in low_blob and "partial" not in low_blob:
        return "coding_complete_only"
    if count_ambiguous_n(nt_sequence or "") > N_AMBIGUOUS_THRESHOLD:
        return "partial_or_unknown"
    if segment and nt_len and nt_len >= COMPLETE_LENGTH_MIN.get(segment, 10**9):
        return "coding_complete" if has_terminal_motif(nt_sequence or "") else "likely_complete"
    if nt_len:
        return "partial_or_unknown"
    return "unknown"


def is_complete_status(status: str) -> bool:
    return status in {"complete", "coding_complete_only", "likely_complete", "coding_complete"}


def choose_host(flat_contexts: Sequence[Dict[str, Any]], biosample: Optional[BioSampleInfo]) -> Tuple[str, str, str]:
    if biosample and biosample.host:
        return biosample.host, clean_text(biosample.host_taxid) or "", "biosample"
    for flat in flat_contexts:
        host = first_nonempty(flat.get("host.organismName"), flat.get("host.scientificName"), flat.get("hostName"), flat.get("host"))
        if host:
            host_taxid = first_nonempty(flat.get("host.taxId"), flat.get("host.taxid")) or ""
            return host, host_taxid, "data_report"
    return "", "", "missing"


def build_group_key(species_taxid: str, biosample: str, isolate: str, strain: str, collection_date: str, location: str, host_name: str, accession: str) -> Tuple[str, str, float]:
    isolate_or_strain = isolate or strain
    if species_taxid and biosample and isolate_or_strain:
        return f"species_biosample_isolate|{species_taxid}|{biosample}|{isolate_or_strain}", "tier_1", 1.0
    if species_taxid and biosample:
        return f"species_biosample|{species_taxid}|{biosample}", "tier_2", 0.92
    if species_taxid and isolate_or_strain and collection_date and location:
        return f"species_isolate_date_loc|{species_taxid}|{isolate_or_strain}|{collection_date}|{location}", "tier_3", 0.80
    if species_taxid and isolate_or_strain and location:
        return f"species_isolate_loc|{species_taxid}|{isolate_or_strain}|{location}", "tier_4", 0.72
    if species_taxid and strain and host_name and location:
        return f"species_strain_host_loc|{species_taxid}|{strain}|{host_name}|{location}", "tier_5", 0.62
    return f"accession_only|{accession}", "fallback", 0.10


def is_lab_host_record(raw_blob: str) -> bool:
    low = raw_blob.lower()
    return any(term in low for term in LAB_HOST_TERMS)


def expected_length_warning(segment: str, nt_len: int) -> bool:
    if not segment or not nt_len:
        return False
    if segment == "S" and nt_len > 4000:
        return True
    if segment == "M" and not (1500 <= nt_len <= 5000):
        return True
    if segment == "L" and nt_len < 3000:
        return True
    return False


def determine_host_role(host_name: str, host_taxid: str, host_source: str, isolation_source: str, raw_blob: str) -> str:
    low_name = (host_name or "").lower()
    low_source = (host_source or "").lower()
    low_iso = (isolation_source or "").lower()
    low_blob = (raw_blob or "").lower()
    if host_taxid == "9606" or "human" in low_name or "homo sapiens" in low_name:
        return "spillover_host"
    if any(term in low_name for term in LAB_HOST_TERMS) or any(term in low_source for term in LAB_HOST_TERMS) or any(term in low_iso for term in LAB_HOST_TERMS) or any(term in low_blob for term in LAB_HOST_TERMS):
        return "laboratory_host"
    excreta_terms = ["urine", "feces", "faeces", "saliva", "stool", "fecal", "urinary", "salivary", "excreta"]
    if any(term in low_iso for term in excreta_terms) or any(term in low_source for term in excreta_terms) or any(term in low_blob for term in excreta_terms):
        return "reservoir_candidate"
    reservoir_terms = ["salivary gland", "kidney", "liver", "lung"]
    if any(term in low_iso for term in reservoir_terms) or any(term in low_source for term in reservoir_terms):
        return "potential_host"
    serology_terms = ["serum", "serology", "antibody", "igg", "igm", "sero"]
    if any(term in low_iso for term in serology_terms) or any(term in low_source for term in serology_terms) or any(term in low_blob for term in serology_terms):
        return "spillover_host"
    if host_name:
        return "potential_host"
    return "unknown"


def compile_product_patterns() -> Dict[str, List[re.Pattern]]:
    return {key: [re.compile(pat, flags=re.I) for pat in patterns] for key, patterns in PROTEIN_PRODUCT_PATTERNS.items()}


def matches_product(text: str, patterns: List[re.Pattern]) -> bool:
    return any(pattern.search(text.lower()) for pattern in patterns)


def extract_protein_sequence_from_annotation(record: Dict[str, Any], protein_index: Dict[str, FastaRecord]) -> Optional[str]:
    for acc in extract_accessions_from_blob(record):
        rec = protein_index.get(acc) or protein_index.get(accession_base(acc))
        if rec and rec.sequence:
            return rec.sequence
    return None


def resolve_segment_protein_annotation(ann_records: Sequence[Dict[str, Any]], protein_index: Dict[str, FastaRecord], segment: str, virus_name: str) -> Tuple[str, str, bool]:
    patterns = compile_product_patterns()
    best_seq = ""
    best_label = ""
    nss_flag = False
    for rec in ann_records:
        text = json.dumps(rec, ensure_ascii=False).lower()
        if segment == "S":
            if matches_product(text, patterns["NSs"]):
                nss_flag = True
            if matches_product(text, patterns["N"]):
                seq = extract_protein_sequence_from_annotation(rec, protein_index)
                if seq and len(seq) > len(best_seq):
                    best_seq = seq
                    best_label = "N"
        elif segment == "M":
            if matches_product(text, patterns["GPC"]):
                seq = extract_protein_sequence_from_annotation(rec, protein_index)
                if seq and len(seq) > len(best_seq):
                    best_seq = seq
                    best_label = "GPC"
        elif segment == "L":
            if matches_product(text, patterns["RdRp"]):
                seq = extract_protein_sequence_from_annotation(rec, protein_index)
                if seq and len(seq) > len(best_seq):
                    best_seq = seq
                    best_label = "RdRp"
    if segment == "S" and any(name in virus_name.lower() for name in ["puumala", "tula"]) and nss_flag:
        return best_seq, best_label or "N", True
    return best_seq, best_label, nss_flag


def compute_kmer_counts(sequence: str, k: int = KMER_K) -> Dict[str, float]:
    counts: Dict[str, float] = {}
    if not sequence or len(sequence) < k:
        return counts
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        counts[kmer] = counts.get(kmer, 0.0) + 1.0
    total = sum(counts.values())
    if total:
        for kmer in counts:
            counts[kmer] /= total
    return counts


def pearson_correlation(xs: Sequence[float], ys: Sequence[float]) -> float:
    if not xs or not ys or len(xs) != len(ys):
        return 0.0
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    denom_x = sum((x - mean_x) ** 2 for x in xs)
    denom_y = sum((y - mean_y) ** 2 for y in ys)
    if denom_x == 0 or denom_y == 0:
        return 0.0
    return num / ((denom_x * denom_y) ** 0.5)


def select_correlation_features(matrix: Sequence[Dict[str, float]], labels: Sequence[str], max_features: int = 200) -> List[str]:
    if not matrix or not labels:
        return []
    feature_names = sorted({k for row in matrix for k in row.keys()})
    label_values = {label: idx for idx, label in enumerate(sorted(set(labels)))}
    label_vector = [label_values[label] for label in labels]
    scores: Dict[str, float] = {}
    for feature in feature_names:
        values = [row.get(feature, 0.0) for row in matrix]
        mean_val = sum(values) / len(values) if values else 0.0
        between = sum((label_vector[i] - sum(label_vector) / len(label_vector)) * (values[i] - mean_val) for i in range(len(values)))
        variance = sum((v - mean_val) ** 2 for v in values)
        scores[feature] = abs(between) / (variance + 1e-9)
    selected: List[str] = []
    for feature in sorted(scores, key=lambda k: scores[k], reverse=True):
        if len(selected) >= max_features:
            break
        feature_values = [row.get(feature, 0.0) for row in matrix]
        if any(abs(pearson_correlation(feature_values, [row.get(existing, 0.0) for row in matrix])) > 0.92 for existing in selected):
            continue
        selected.append(feature)
    return selected


def normalize_matrix_columns(rows: Sequence[Dict[str, Any]], feature_names: Sequence[str]) -> None:
    for feature in feature_names:
        values = [float(row.get(feature, 0.0)) for row in rows]
        mean = sum(values) / len(values) if values else 0.0
        std = (sum((v - mean) ** 2 for v in values) / len(values)) ** 0.5 if values else 0.0
        if std <= 0:
            continue
        for row in rows:
            row[feature] = (float(row.get(feature, 0.0)) - mean) / std


def compute_protein_distance(seq1: str, seq2: str) -> float:
    if not seq1 or not seq2:
        return 1.0
    len1, len2 = len(seq1), len(seq2)
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    mismatches += abs(len1 - len2)
    return mismatches / max(len1, len2)


def cluster_species_by_ped(isolates: Sequence[IsolateRecord]) -> Dict[str, int]:
    clusters: Dict[str, int] = {}
    species_groups: DefaultDict[str, List[IsolateRecord]] = defaultdict(list)
    for isolate in isolates:
        species_groups[isolate.species_taxid].append(isolate)
    for species_taxid, members in species_groups.items():
        if len(members) == 1:
            clusters[members[0].strict_isolate_key] = 0
            continue
        component: Dict[str, int] = {}
        next_id = 0
        for isolate in members:
            if isolate.strict_isolate_key in component:
                continue
            stack = [isolate]
            component[isolate.strict_isolate_key] = next_id
            while stack:
                current = stack.pop()
                for other in members:
                    if other.strict_isolate_key in component:
                        continue
                    if compute_protein_distance(current.concatenated_protein_aa, other.concatenated_protein_aa) <= PED_SPECIES_THRESHOLD:
                        component[other.strict_isolate_key] = next_id
                        stack.append(other)
            next_id += 1
        clusters.update(component)
    return clusters


def transform_one_record(
    record: Dict[str, Any],
    annotation_index: Dict[str, List[Dict[str, Any]]],
    biosample_index: Dict[str, BioSampleInfo],
    genomic_index: Dict[str, FastaRecord],
    cds_index: Dict[str, FastaRecord],
    protein_index: Dict[str, FastaRecord],
) -> Tuple[Optional[AccessionRecord], Optional[Dict[str, Any]]]:
    try:
        flat = flatten_json(record)
        accession = first_nonempty(record.get("accession"), flat.get("accession"), flat.get("nucleotideAccession"))
        if not accession:
            candidates = extract_accessions_from_blob(record)
            accession = candidates[0] if candidates else "UNKNOWN_ACCESSION"
        acc_base = accession_base(accession)

        ann_matches = uniq(annotation_index.get(accession, []) + annotation_index.get(acc_base, []))
        annotation = extract_annotation_evidence(ann_matches, protein_index)

        biosample_candidates = extract_biosample_candidates(record)
        biosample_acc = first_nonempty(biosample_candidates, flat.get("bioSampleAccession"), flat.get("biosample.accession"), flat.get("biosample")) or ""
        biosample_info = biosample_index.get(biosample_acc.upper()) if biosample_acc else None

        flat_contexts = [flat] + [flatten_json(x) for x in ann_matches]
        if biosample_info:
            flat_contexts.append(flatten_json(asdict(biosample_info)))

        virus_name = first_nonempty(flat.get("virus.organismName"), flat.get("organismName"), flat.get("organism.organismName")) or ""
        virus_taxid = first_nonempty(flat.get("virus.taxId"), flat.get("taxId"), flat.get("organism.taxId")) or ""
        genus_name = ""
        genus_taxid = ""
        lineage = record.get("virus", {}).get("lineage", []) if isinstance(record.get("virus"), dict) else []
        if isinstance(lineage, list):
            for item in lineage:
                if isinstance(item, dict) and item.get("name") == "Orthohantavirus":
                    genus_name = clean_text(item.get("name")) or ""
                    genus_taxid = str(item.get("taxId") or "")

        isolate = first_nonempty(flat.get("isolate.name"), flat.get("virus.isolate"), biosample_info.isolate if biosample_info else None) or ""
        strain = first_nonempty(flat.get("strain"), flat.get("virus.strain"), biosample_info.strain if biosample_info else None) or ""
        collection_date = first_nonempty(flat.get("isolate.collectionDate"), flat.get("collectionDate"), biosample_info.collection_date if biosample_info else None) or ""
        geographic = first_nonempty(flat.get("location.geographicLocation"), biosample_info.geographic_location if biosample_info else None) or ""
        country, region = split_geo(geographic)
        host_name, host_taxid, host_source = choose_host(flat_contexts, biosample_info)
        isolation_source = first_nonempty(flat.get("isolate.source"), biosample_info.isolation_source if biosample_info else None) or ""
        tissue = first_nonempty(flat.get("tissue"), biosample_info.tissue if biosample_info else None) or ""
        host_role = determine_host_role(host_name, host_taxid, host_source, isolation_source, json.dumps([record, biosample_info.__dict__] if biosample_info else [record], ensure_ascii=False))
        year = extract_year(collection_date)

        genomic_rec = genomic_index.get(accession) or genomic_index.get(acc_base)
        cds_rec = cds_index.get(accession) or cds_index.get(acc_base)
        nt_len = (genomic_rec.length if genomic_rec else 0) or (cds_rec.length if cds_rec else 0)
        segment, segment_source, segment_warnings = choose_segment(flat_contexts, annotation, nt_len)
        nt_sequence = genomic_rec.sequence if genomic_rec and genomic_rec.sequence else cds_rec.sequence if cds_rec and cds_rec.sequence else ""
        completeness = determine_completeness(flat_contexts, segment, nt_len, json.dumps([record] + ann_matches, ensure_ascii=False), nt_sequence=nt_sequence)

        protein_accessions: List[str] = []
        for p_acc in annotation.protein_accessions:
            rec = protein_index.get(p_acc) or protein_index.get(accession_base(p_acc))
            protein_accessions.append(rec.accession if rec else p_acc)
        protein_accessions = uniq(protein_accessions)

        protein_aa = ""
        protein_product = ""
        nss_overlap_flag = False
        if segment in {"S", "M"}:
            protein_aa, protein_product, nss_overlap_flag = resolve_segment_protein_annotation(ann_matches, protein_index, segment, virus_name)
            if protein_aa:
                qc_id = f"protein_{protein_product}"
                if qc_id not in protein_accessions:
                    protein_accessions.append(qc_id)

        group_key, grouping_tier, group_confidence = build_group_key(
            species_taxid=virus_taxid,
            biosample=biosample_acc,
            isolate=isolate,
            strain=strain,
            collection_date=collection_date,
            location=geographic,
            host_name=host_name,
            accession=accession,
        )

        warning_flags = uniq([
            *segment_warnings,
            *( ["missing_segment"] if segment == "" else []),
            *( ["segment_from_heuristic"] if segment_source == "length_heuristic" else []),
            *( ["missing_host"] if not host_name else []),
            *( ["missing_biosample"] if not biosample_acc else []),
            *( ["missing_isolate_and_strain"] if not isolate and not strain else []),
            *( ["lab_host_flag"] if is_lab_host_record(json.dumps(record, ensure_ascii=False)) else []),
            *( ["missing_genomic_sequence"] if genomic_rec is None else []),
            *( ["missing_protein_link"] if len(protein_accessions) == 0 else []),
            *( ["unexpected_length_for_segment"] if expected_length_warning(segment, nt_len) else []),
        ])

        rec = AccessionRecord(
            accession=accession,
            accession_base=acc_base,
            species_name=virus_name,
            species_taxid=virus_taxid,
            genus_name=genus_name,
            genus_taxid=genus_taxid,
            segment=segment,
            segment_source=segment_source,
            completeness=completeness,
            nt_length=nt_len,
            biosample=biosample_acc,
            bioproject=first_nonempty(flat.get("bioprojects"), flat.get("bioProjectAccession"), flat.get("bioproject.accession")) or "",
            isolate=isolate,
            strain=strain,
            collection_date=collection_date,
            year=year,
            location=geographic,
            country=country,
            region=region,
            host_name=host_name,
            host_taxid=host_taxid,
            host_source=host_source,
            host_role=host_role,
            isolation_source=isolation_source,
            tissue_source=tissue,
            protein_accessions=protein_accessions,
            annotation_products=annotation.product_names,
            protein_aa=protein_aa,
            protein_product=protein_product,
            nss_overlap_flag=nss_overlap_flag,
            strict_isolate_key="",
            group_key=group_key,
            grouping_tier=grouping_tier,
            group_confidence=group_confidence,
            is_refseq=accession.startswith(("NC_", "NM_", "NP_", "YP_", "XM_", "XP_", "XR_", "NR_")),
            warning_flags=[item for item in warning_flags if item],
        )

        qc_row = {
            "accession": accession,
            "group_key": group_key,
            "grouping_tier": grouping_tier,
            "group_confidence": group_confidence,
            "segment": segment,
            "completeness": completeness,
            "warnings": ";".join(rec.warning_flags),
        }
        return rec, qc_row
    except Exception as exc:
        logging.exception("Failed to transform accession record")
        return None, {"accession": "", "group_key": "", "segment": "", "completeness": "", "warnings": f"transform_error:{exc}"}


def group_accessions(accessions: Sequence[AccessionRecord]) -> List[VirusGroupRecord]:
    grouped: DefaultDict[str, List[AccessionRecord]] = defaultdict(list)
    for rec in accessions:
        grouped[rec.group_key].append(rec)
    out: List[VirusGroupRecord] = []
    for key, members in grouped.items():
        species_name = Counter([m.species_name for m in members if m.species_name]).most_common(1)
        species_taxid = Counter([m.species_taxid for m in members if m.species_taxid]).most_common(1)
        host_name = Counter([m.host_name for m in members if m.host_name]).most_common(1)
        host_taxid = Counter([m.host_taxid for m in members if m.host_taxid]).most_common(1)
        country = Counter([m.country for m in members if m.country]).most_common(1)
        year = Counter([m.year for m in members if m.year]).most_common(1)
        biosample = Counter([m.biosample for m in members if m.biosample]).most_common(1)
        isolate = Counter([m.isolate for m in members if m.isolate]).most_common(1)
        strain = Counter([m.strain for m in members if m.strain]).most_common(1)
        s_accessions = sorted({m.accession for m in members if m.segment == "S"})
        m_accessions = sorted({m.accession for m in members if m.segment == "M"})
        l_accessions = sorted({m.accession for m in members if m.segment == "L"})
        complete_segments = sorted({m.segment for m in members if m.segment in {"S", "M", "L"} and is_complete_status(m.completeness)})
        if complete_segments == ["L", "M", "S"]:
            completeness_class = "complete_triplet"
        elif len(complete_segments) >= 1:
            completeness_class = "partial_group_with_complete_segment"
        else:
            completeness_class = "no_complete_segment"
        out.append(VirusGroupRecord(
            group_key=key,
            species_name=species_name[0][0] if species_name else "",
            species_taxid=species_taxid[0][0] if species_taxid else "",
            host_name=host_name[0][0] if host_name else "",
            host_taxid=host_taxid[0][0] if host_taxid else "",
            country=country[0][0] if country else "",
            year=year[0][0] if year else "",
            biosample=biosample[0][0] if biosample else "",
            isolate=isolate[0][0] if isolate else "",
            strain=strain[0][0] if strain else "",
            group_confidence=max(m.group_confidence for m in members),
            group_size=len(members),
            s_accessions=s_accessions,
            m_accessions=m_accessions,
            l_accessions=l_accessions,
            complete_segments=complete_segments,
            completeness_class=completeness_class,
            all_hosts=sorted({m.host_name for m in members if m.host_name}),
            all_accessions=sorted({m.accession for m in members}),
        ))
    return out


def choose_best_segment_record(records: Sequence[AccessionRecord]) -> Optional[AccessionRecord]:
    if not records:
        return None
    rank = {"complete": 4, "coding_complete_only": 4, "likely_complete": 3, "partial_or_unknown": 2, "unknown": 1}
    return sorted(records, key=lambda rec: (rank.get(rec.completeness, 0), rec.nt_length), reverse=True)[0]


def build_strict_isolate_records(accessions: Sequence[AccessionRecord]) -> List[IsolateRecord]:
    groups: DefaultDict[str, List[AccessionRecord]] = defaultdict(list)
    for rec in accessions:
        groups[rec.group_key].append(rec)

    isolates: List[IsolateRecord] = []
    for key, members in groups.items():
        s_rec = choose_best_segment_record([r for r in members if r.segment == "S"])
        m_rec = choose_best_segment_record([r for r in members if r.segment == "M"])
        l_rec = choose_best_segment_record([r for r in members if r.segment == "L"])
        if not (s_rec or m_rec or l_rec):
            continue

        protein_warnings: List[str] = []
        concatenated_protein = ""
        if s_rec and m_rec and s_rec.protein_aa and m_rec.protein_aa:
            concatenated_protein = s_rec.protein_aa + m_rec.protein_aa
        else:
            protein_warnings.append("missing_n_or_gpc_sequence")

        taxid_counts = Counter([r.species_taxid for r in [s_rec, m_rec, l_rec] if r and r.species_taxid])
        most_common_taxid = taxid_counts.most_common(1)[0][0] if taxid_counts else ""
        taxonomic_incongruence = any(r and r.species_taxid and r.species_taxid != most_common_taxid for r in [s_rec, m_rec, l_rec])
        segment_complete = {"S": is_complete_status(s_rec.completeness) if s_rec else False, "M": is_complete_status(m_rec.completeness) if m_rec else False, "L": is_complete_status(l_rec.completeness) if l_rec else False}
        complete_count = sum(segment_complete.values())
        if complete_count == 3:
            coding_complete_status = "coding_complete_triplet"
        elif complete_count == 2:
            coding_complete_status = "near_complete_triplet"
        elif complete_count == 1:
            coding_complete_status = "partial_triplet"
        else:
            coding_complete_status = "incomplete_triplet"

        mosaic_flag = bool(
            (s_rec and m_rec and s_rec.collection_date != m_rec.collection_date)
            or (s_rec and l_rec and s_rec.collection_date != l_rec.collection_date)
            or (m_rec and l_rec and m_rec.collection_date != l_rec.collection_date)
            or (s_rec and m_rec and s_rec.location != m_rec.location)
            or (s_rec and l_rec and s_rec.location != l_rec.location)
            or (m_rec and l_rec and m_rec.location != l_rec.location)
        )

        strict_key = ""
        if s_rec and s_rec.biosample and (s_rec.isolate or s_rec.strain):
            strict_key = f"{s_rec.species_taxid}|{s_rec.biosample}|{s_rec.isolate or s_rec.strain}"
        elif m_rec and m_rec.biosample and (m_rec.isolate or m_rec.strain):
            strict_key = f"{m_rec.species_taxid}|{m_rec.biosample}|{m_rec.isolate or m_rec.strain}"
        elif l_rec and l_rec.biosample and (l_rec.isolate or l_rec.strain):
            strict_key = f"{l_rec.species_taxid}|{l_rec.biosample}|{l_rec.isolate or l_rec.strain}"
        else:
            strict_key = key

        host_name = (s_rec.host_name if s_rec else "") or (m_rec.host_name if m_rec else "") or (l_rec.host_name if l_rec else "")
        host_taxid = (s_rec.host_taxid if s_rec else "") or (m_rec.host_taxid if m_rec else "") or (l_rec.host_taxid if l_rec else "")
        host_source = (s_rec.host_source if s_rec else "") or (m_rec.host_source if m_rec else "") or (l_rec.host_source if l_rec else "")
        isolation_source = (s_rec.isolation_source if s_rec else "") or (m_rec.isolation_source if m_rec else "") or (l_rec.isolation_source if l_rec else "")
        collection_date = (s_rec.collection_date if s_rec else "") or (m_rec.collection_date if m_rec else "") or (l_rec.collection_date if l_rec else "")
        location = (s_rec.location if s_rec else "") or (m_rec.location if m_rec else "") or (l_rec.location if l_rec else "")
        country = (s_rec.country if s_rec else "") or (m_rec.country if m_rec else "") or (l_rec.country if l_rec else "")
        region = (s_rec.region if s_rec else "") or (m_rec.region if m_rec else "") or (l_rec.region if l_rec else "")
        biosample = (s_rec.biosample if s_rec else "") or (m_rec.biosample if m_rec else "") or (l_rec.biosample if l_rec else "")
        isolate = (s_rec.isolate if s_rec else "") or (m_rec.isolate if m_rec else "") or (l_rec.isolate if l_rec else "")
        strain = (s_rec.strain if s_rec else "") or (m_rec.strain if m_rec else "") or (l_rec.strain if l_rec else "")
        warning_flags = uniq([*(s_rec.warning_flags if s_rec else []), *(m_rec.warning_flags if m_rec else []), *(l_rec.warning_flags if l_rec else [])])

        isolates.append(IsolateRecord(
            strict_isolate_key=strict_key,
            species_name=(s_rec.species_name if s_rec else "") or (m_rec.species_name if m_rec else "") or (l_rec.species_name if l_rec else "") or "",
            species_taxid=most_common_taxid,
            species_cluster_id="",
            biosample=biosample,
            isolate=isolate,
            strain=strain,
            host_name=host_name,
            host_taxid=host_taxid,
            host_role=determine_host_role(host_name, host_taxid, host_source, isolation_source, json.dumps([asdict(r) for r in [s_rec, m_rec, l_rec] if r], ensure_ascii=False)),
            collection_date=collection_date,
            location=location,
            country=country,
            region=region,
            grouping_tier=members[0].grouping_tier,
            grouping_confidence=max(m.group_confidence for m in members),
            s_accession=s_rec.accession if s_rec else "",
            m_accession=m_rec.accession if m_rec else "",
            l_accession=l_rec.accession if l_rec else "",
            s_completeness=s_rec.completeness if s_rec else "unknown",
            m_completeness=m_rec.completeness if m_rec else "unknown",
            l_completeness=l_rec.completeness if l_rec else "unknown",
            s_nt_length=s_rec.nt_length if s_rec else 0,
            m_nt_length=m_rec.nt_length if m_rec else 0,
            l_nt_length=l_rec.nt_length if l_rec else 0,
            concatenated_protein_aa=concatenated_protein,
            coding_complete_status=coding_complete_status,
            mosaic_flag=mosaic_flag,
            taxonomic_incongruence_warning=taxonomic_incongruence,
            reassortment_candidate_warning=taxonomic_incongruence,
            ns_s_overlapping_orf=s_rec.nss_overlap_flag if s_rec else False,
            species_divergence_flag=False,
            s_protein_evidence=s_rec.protein_product if s_rec else "",
            m_protein_evidence=m_rec.protein_product if m_rec else "",
            l_protein_evidence=l_rec.protein_product if l_rec else "",
            s_protein_accessions=s_rec.protein_accessions if s_rec else [],
            m_protein_accessions=m_rec.protein_accessions if m_rec else [],
            l_protein_accessions=l_rec.protein_accessions if l_rec else [],
            missing_protein_warnings=uniq(protein_warnings),
            host_source=host_source,
            raw_metadata={"group_key": key, "grouping_tier": members[0].grouping_tier, "group_confidence": max(m.group_confidence for m in members), "biosample": biosample},
            warnings=uniq([*protein_warnings, *warning_flags]),
        ))

    if isolates:
        cluster_map = cluster_species_by_ped(isolates)
        for isolate in isolates:
            isolate.species_cluster_id = f"{isolate.species_taxid}.cluster{cluster_map.get(isolate.strict_isolate_key, 0)}"
            isolate.species_divergence_flag = cluster_map.get(isolate.strict_isolate_key, 0) != 0
    return isolates


def build_host_table(accessions: Sequence[AccessionRecord], isolates: Sequence[IsolateRecord]) -> List[Dict[str, Any]]:
    hosts: Dict[Tuple[str, str], HostRecord] = {}
    for rec in accessions:
        key = (rec.host_name or "", rec.host_taxid or "")
        if not any(key):
            continue
        entry = hosts.get(key)
        host_role = determine_host_role(rec.host_name, rec.host_taxid, rec.host_source, rec.isolation_source, json.dumps(rec.annotation_products, ensure_ascii=False))
        if entry is None:
            hosts[key] = HostRecord(host_name=rec.host_name, host_taxid=rec.host_taxid, host_role=host_role, host_source=rec.host_source, representative_biosample=rec.biosample, representative_isolate=rec.isolate, accession_count=1, isolate_count=0, example_species=rec.species_name)
        else:
            entry.accession_count += 1
            entry.host_role = entry.host_role or host_role
            entry.representative_biosample = entry.representative_biosample or rec.biosample
            entry.representative_isolate = entry.representative_isolate or rec.isolate
            entry.example_species = entry.example_species or rec.species_name
    for isolate in isolates:
        key = (isolate.host_name or "", isolate.host_taxid or "")
        if not any(key):
            continue
        entry = hosts.get(key)
        if entry is None:
            hosts[key] = HostRecord(host_name=isolate.host_name, host_taxid=isolate.host_taxid, host_role=isolate.host_role, host_source="isolate", representative_biosample=isolate.biosample, representative_isolate=isolate.isolate, accession_count=0, isolate_count=1, example_species=isolate.species_name)
        else:
            entry.isolate_count += 1
            entry.host_role = entry.host_role or isolate.host_role
            entry.example_species = entry.example_species or isolate.species_name
    return [asdict(rec) for rec in hosts.values()]


def build_ml_feature_matrix(isolates: Sequence[IsolateRecord]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    matrices: List[Dict[str, float]] = []
    labels: List[str] = []
    for isolate in isolates:
        if not isolate.concatenated_protein_aa:
            continue
        counts = compute_kmer_counts(isolate.concatenated_protein_aa, KMER_K)
        if not counts:
            continue
        matrices.append(counts)
        labels.append(isolate.host_taxid or isolate.species_taxid)
        rows.append({"strict_isolate_key": isolate.strict_isolate_key, "species_name": isolate.species_name, "species_taxid": isolate.species_taxid, "host_name": isolate.host_name, "host_taxid": isolate.host_taxid, "host_role": isolate.host_role})
    selected_features = select_correlation_features(matrices, labels)
    for row, counts in zip(rows, matrices):
        for feature in selected_features:
            row[f"kmer_{feature}"] = counts.get(feature, 0.0)
    normalize_matrix_columns(rows, [f"kmer_{feature}" for feature in selected_features])
    return rows


def build_suppressed_rows(accessions: Sequence[AccessionRecord], isolates: Sequence[IsolateRecord]) -> List[Dict[str, Any]]:
    gold_accessions = {iso.s_accession for iso in isolates} | {iso.m_accession for iso in isolates} | {iso.l_accession for iso in isolates}
    rows: List[Dict[str, Any]] = []
    for rec in accessions:
        if rec.accession in gold_accessions:
            continue
        rows.append({
            "accession": rec.accession,
            "species_name": rec.species_name,
            "species_taxid": rec.species_taxid,
            "segment": rec.segment,
            "completeness": rec.completeness,
            "biosample": rec.biosample,
            "isolate": rec.isolate,
            "strain": rec.strain,
            "collection_date": rec.collection_date,
            "location": rec.location,
            "host_name": rec.host_name,
            "host_taxid": rec.host_taxid,
            "suppression_reason": ";".join(rec.warning_flags) if rec.warning_flags else "not_best_representative",
        })
    return rows


def build_isolate_rows(isolates: Sequence[IsolateRecord]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for isolate in isolates:
        rows.append({
            "strict_isolate_key": isolate.strict_isolate_key,
            "species_name": isolate.species_name,
            "species_taxid": isolate.species_taxid,
            "species_cluster_id": isolate.species_cluster_id,
            "biosample": isolate.biosample,
            "isolate": isolate.isolate,
            "strain": isolate.strain,
            "host_name": isolate.host_name,
            "host_taxid": isolate.host_taxid,
            "host_role": isolate.host_role,
            "collection_date": isolate.collection_date,
            "location": isolate.location,
            "country": isolate.country,
            "region": isolate.region,
            "s_accession": isolate.s_accession,
            "m_accession": isolate.m_accession,
            "l_accession": isolate.l_accession,
            "s_completeness": isolate.s_completeness,
            "m_completeness": isolate.m_completeness,
            "l_completeness": isolate.l_completeness,
            "s_nt_length": isolate.s_nt_length,
            "m_nt_length": isolate.m_nt_length,
            "l_nt_length": isolate.l_nt_length,
            "concatenated_protein_aa": isolate.concatenated_protein_aa,
            "coding_complete_status": isolate.coding_complete_status,
            "grouping_tier": isolate.grouping_tier,
            "grouping_confidence": isolate.grouping_confidence,
            "s_protein_evidence": isolate.s_protein_evidence,
            "m_protein_evidence": isolate.m_protein_evidence,
            "l_protein_evidence": isolate.l_protein_evidence,
            "s_protein_accessions": ";".join(isolate.s_protein_accessions),
            "m_protein_accessions": ";".join(isolate.m_protein_accessions),
            "l_protein_accessions": ";".join(isolate.l_protein_accessions),
            "missing_protein_warnings": ";".join(isolate.missing_protein_warnings),
            "mosaic_flag": int(isolate.mosaic_flag),
            "reassortment_candidate_warning": int(isolate.reassortment_candidate_warning),
            "taxonomic_incongruence_warning": int(isolate.taxonomic_incongruence_warning),
            "ns_s_overlapping_orf": int(isolate.ns_s_overlapping_orf),
            "species_divergence_flag": int(isolate.species_divergence_flag),
            "warnings": ";".join(isolate.warnings),
        })
    return rows


def build_ml_segments(accessions: Sequence[AccessionRecord]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for rec in accessions:
        if rec.segment in {"S", "M", "L"} and rec.host_name:
            rows.append({
                "accession": rec.accession,
                "segment": rec.segment,
                "species_name": rec.species_name,
                "species_taxid": rec.species_taxid,
                "host_name": rec.host_name,
                "host_taxid": rec.host_taxid,
                "country": rec.country,
                "year": rec.year,
                "biosample": rec.biosample,
                "completeness": rec.completeness,
                "group_key": rec.group_key,
            })
    return rows


def summarize(accessions: Sequence[AccessionRecord], groups: Sequence[VirusGroupRecord], qc_rows: Sequence[Dict[str, Any]], paths: InputPaths) -> Dict[str, Any]:
    catalog = load_dataset_catalog(paths.dataset_catalog)
    return {
        "input_root": str(paths.root),
        "data_dir": str(paths.data_dir),
        "records_total": len(accessions),
        "groups_total": len(groups),
        "records_with_biosample": sum(1 for r in accessions if r.biosample),
        "records_with_host": sum(1 for r in accessions if r.host_name),
        "segment_counts": dict(Counter(r.segment or "unresolved" for r in accessions)),
        "completeness_counts": dict(Counter(r.completeness for r in accessions)),
        "group_completeness_counts": dict(Counter(g.completeness_class for g in groups)),
        "qc_warning_counts": dict(Counter(w for row in qc_rows for w in row.get("warnings", "").split(";") if w)),
        "provenance": {
            "generated_at": datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
            "dataset_catalog_api_version": catalog.get("apiVersion") or catalog.get("api_version") or "unknown",
            "dataset_catalog_md5": md5sum(paths.dataset_catalog),
            "input_md5": {
                "data_report.jsonl": md5sum(paths.data_report),
                "annotation_report.jsonl": md5sum(paths.annotation_report),
                "biosample_report.jsonl": md5sum(paths.biosample_report),
                "dataset_catalog.json": md5sum(paths.dataset_catalog),
            },
            "input_files": [str(paths.data_report), str(paths.annotation_report), str(paths.biosample_report), str(paths.dataset_catalog)],
        },
        "output_files": {
            "accession_table": str(paths.output_dir / "ortho_accession_table.parquet"),
            "virus_table": str(paths.output_dir / "ortho_virus_table.parquet"),
            "isolate_table": str(paths.output_dir / "ortho_isolate_table.parquet"),
            "host_table": str(paths.output_dir / "ortho_host_table.parquet"),
            "ml_feature_matrix": str(paths.output_dir / "ortho_ml_feature_matrix.parquet"),
            "suppressed_table": str(paths.output_dir / "ortho_suppressed.parquet"),
            "ml_segments": str(paths.output_dir / "ortho_ml_segments.parquet"),
            "qc_report": str(paths.output_dir / "ortho_qc_report.parquet"),
            "summary_json": str(paths.output_dir / "ortho_summary.json"),
            "log_file": str(paths.log_file),
        },
        "output_schema": {name: {"columns": columns} for name, columns in OUTPUT_TABLE_SCHEMAS.items()},
    }


def run_pipeline(paths: InputPaths) -> None:
    logging.info("Loading JSONL reports")
    data_records = list(iter_jsonl(paths.data_report))
    annotation_records = list(iter_jsonl(paths.annotation_report))
    biosample_records = list(iter_jsonl(paths.biosample_report))

    logging.info("Building indexes")
    annotation_index = build_annotation_index(annotation_records)
    biosample_index = build_biosample_index(biosample_records)

    logging.info("Parsing FASTA files")
    genomic_index = parse_fasta(paths.genomic_fasta)
    cds_index = parse_fasta(paths.cds_fasta)
    protein_index = parse_fasta(paths.protein_fasta)

    logging.info("Transforming accession records")
    accession_records: List[AccessionRecord] = []
    qc_rows: List[Dict[str, Any]] = []
    for rec in data_records:
        acc_rec, qc_row = transform_one_record(rec, annotation_index, biosample_index, genomic_index, cds_index, protein_index)
        if acc_rec is not None:
            accession_records.append(acc_rec)
        if qc_row is not None:
            qc_rows.append(qc_row)

    logging.info("Grouping accessions")
    virus_groups = group_accessions(accession_records)

    logging.info("Building isolate records")
    isolate_records = build_strict_isolate_records(accession_records)

    logging.info("Building host table")
    host_rows = build_host_table(accession_records, isolate_records)

    logging.info("Building ML feature matrix")
    ml_feature_rows = build_ml_feature_matrix(isolate_records)

    logging.info("Building suppressed accession table")
    suppressed_rows = build_suppressed_rows(accession_records, isolate_records)

    accession_rows = []
    for rec in accession_records:
        accession_rows.append({
            "accession": rec.accession,
            "species_name": rec.species_name,
            "species_taxid": rec.species_taxid,
            "genus_name": rec.genus_name,
            "genus_taxid": rec.genus_taxid,
            "segment": rec.segment,
            "segment_source": rec.segment_source,
            "completeness": rec.completeness,
            "nt_length": rec.nt_length,
            "biosample": rec.biosample,
            "bioproject": rec.bioproject,
            "isolate": rec.isolate,
            "strain": rec.strain,
            "collection_date": rec.collection_date,
            "year": rec.year,
            "location": rec.location,
            "country": rec.country,
            "region": rec.region,
            "host_name": rec.host_name,
            "host_taxid": rec.host_taxid,
            "host_source": rec.host_source,
            "host_role": rec.host_role,
            "isolation_source": rec.isolation_source,
            "tissue_source": rec.tissue_source,
            "protein_accessions": ";".join(rec.protein_accessions),
            "annotation_products": ";".join(rec.annotation_products),
            "group_key": rec.group_key,
            "group_confidence": rec.group_confidence,
            "is_refseq": int(rec.is_refseq),
            "warning_flags": ";".join(rec.warning_flags),
        })

    virus_rows = []
    for group in virus_groups:
        virus_rows.append({
            "group_key": group.group_key,
            "species_name": group.species_name,
            "species_taxid": group.species_taxid,
            "host_name": group.host_name,
            "host_taxid": group.host_taxid,
            "country": group.country,
            "year": group.year,
            "biosample": group.biosample,
            "isolate": group.isolate,
            "strain": group.strain,
            "group_confidence": group.group_confidence,
            "group_size": group.group_size,
            "s_accessions": ";".join(group.s_accessions),
            "m_accessions": ";".join(group.m_accessions),
            "l_accessions": ";".join(group.l_accessions),
            "complete_segments": ";".join(group.complete_segments),
            "completeness_class": group.completeness_class,
            "all_hosts": ";".join(group.all_hosts),
            "all_accessions": ";".join(group.all_accessions),
        })

    logging.info("Writing outputs")
    safe_write_parquet(paths.output_dir / "ortho_accession_table.parquet", accession_rows)
    safe_write_parquet(paths.output_dir / "ortho_virus_table.parquet", virus_rows)
    safe_write_parquet(paths.output_dir / "ortho_isolate_table.parquet", build_isolate_rows(isolate_records))
    safe_write_parquet(paths.output_dir / "ortho_host_table.parquet", host_rows)
    safe_write_parquet(paths.output_dir / "ortho_ml_feature_matrix.parquet", ml_feature_rows)
    safe_write_parquet(paths.output_dir / "ortho_ml_segments.parquet", build_ml_segments(accession_records))
    safe_write_parquet(paths.output_dir / "ortho_suppressed.parquet", suppressed_rows)
    safe_write_parquet(paths.output_dir / "ortho_qc_report.parquet", qc_rows)

    summary = summarize(accession_records, virus_groups, qc_rows, paths)
    safe_write_json(paths.output_dir / "ortho_summary.json", summary)

    logging.info("Done")
    logging.info("Accession records: %s", len(accession_records))
    logging.info("Virus groups: %s", len(virus_groups))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Orthohantavirus NCBI Datasets production pipeline")
    parser.add_argument("--input", required=False, default=".", help="Path to the package root containing all_hanta/ncbi_dataset/data")
    parser.add_argument("--output", required=False, help="Path to write outputs; defaults to <input>/output")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    paths = resolve_inputs(Path(args.input), output=Path(args.output) if args.output else None)
    configure_logging(paths.log_file, verbose=args.verbose)
    logging.info("Starting Orthohantavirus production pipeline")
    logging.info("Input root: %s", paths.root)
    logging.info("Output dir: %s", paths.output_dir)
    run_pipeline(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())