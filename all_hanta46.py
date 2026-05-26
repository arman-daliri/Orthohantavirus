"""
Orthohantavirus Master Dataset Pipeline (production-ready)
NCBI unzipped package -> grouped ML-ready JSON + host-prediction CSV + summary JSON

Key fixes in this revision
--------------------------
1. Correctly populates protein_accessions and protein_lengths_aa from protein.faa and metadata cross-links.
2. Fixes the common inflation artifact where groups with only one complete segment are counted by records
   instead of by groups. Summary now reports both group-level and record-level pattern counts separately.
3. Cross-checks metadata across data_report.jsonl, annotation_report.jsonl, biosample_report.jsonl,
   genomic.fna, cds.fna, and protein.faa.
4. Uses BioSample-first grouping, then exact fallback on isolate_or_strain + collection_date + location + host.
5. Keeps accession-level provenance while also exporting one row per group for host prediction.
6. sequencehash is preserved as an empty string exactly as requested.
7. Fixes BioSample extraction for hierarchical NCBI virus JSON where BioSample may appear as nested arrays/objects.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, DefaultDict, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

INPUT_ROOT = Path("/home/victorfrankestein/Documents/hanta_revised5/all_hanta")
DATA_DIR = INPUT_ROOT / "ncbi_dataset" / "data"
OUTPUT_DIR = Path("/home/victorfrankestein/Documents/hanta_revised5/output")

DATA_REPORT = DATA_DIR / "data_report.jsonl"
ANNOTATION_REPORT = DATA_DIR / "annotation_report.jsonl"
BIOSAMPLE_REPORT = DATA_DIR / "biosample_report.jsonl"
GENOMIC_FASTA = DATA_DIR / "genomic.fna"
CDS_FASTA = DATA_DIR / "cds.fna"
PROTEIN_FASTA = DATA_DIR / "protein.faa"

OUT_JSON = OUTPUT_DIR / "orthohantavirus_ml_ready.json"
OUT_CSV = OUTPUT_DIR / "orthohantavirus_host_prediction_master.csv"
OUT_SUMMARY = OUTPUT_DIR / "orthohantavirus_summary.json"

MIN_COMPLETE_NT = {"L": 5000, "M": 2500, "S": 1000}
ACCESSION_RE = re.compile(r"[A-Z]{1,4}_?\d{4,}(?:\.\d+)?")
BIOSAMPLE_RE = re.compile(r"\bSAM[NED]?[A-Z]?\d+\b", flags=re.I)
REFSEQ_PREFIXES = ("NC_", "NM_", "NP_", "YP_", "XM_", "XP_", "XR_", "NR_")

LAB_PASSAGE_TERMS = [
    "vero",
    "passage",
    "passaged",
    "cell culture",
    "cell line",
    "cell-line",
    "adapted",
    "in vitro",
]

SEGMENT_RULES = {
    "S": [
        r"\bnucleocapsid\b",
        r"\bnucleoprotein\b",
        r"\bnp\b",
        r"\bn protein\b",
        r"\bsmall\b",
        r"\bsegment\s*s\b",
    ],
    "M": [
        r"\bglycoprotein precursor\b",
        r"\bglycoprotein\b",
        r"\bgpc\b",
        r"\bgn\b",
        r"\bgc\b",
        r"\bg1\b",
        r"\bg2\b",
        r"\bgp\b",
        r"\bmedium\b",
        r"\bsegment\s*m\b",
    ],
    "L": [
        r"\bpolymerase\b",
        r"\brdrp\b",
        r"\brna[- ]dependent rna polymerase\b",
        r"\bl protein\b",
        r"\blarge\b",
        r"\bsegment\s*l\b",
    ],
}


def ensure_inputs() -> None:
    missing = [
        p for p in [
            DATA_REPORT,
            ANNOTATION_REPORT,
            BIOSAMPLE_REPORT,
            GENOMIC_FASTA,
            CDS_FASTA,
            PROTEIN_FASTA,
        ] if not p.exists()
    ]
    if missing:
        msg = "\n".join(str(p) for p in missing)
        raise FileNotFoundError(f"Missing required inputs:\n{msg}")


def safe_jsonable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): safe_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [safe_jsonable(v) for v in obj]
    return obj


def clean_text(value: Any) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, (dict, list)):
        text = json.dumps(value, ensure_ascii=False)
    else:
        text = str(value).strip()
    if not text:
        return None
    if text.lower() in {"none", "null", "nan", "na", "n/a", "unknown", "missing", "not collected"}:
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
    out: List[Any] = []
    for value in values:
        key = json.dumps(safe_jsonable(value), ensure_ascii=False, sort_keys=True) if isinstance(value, (dict, list, tuple, set)) else str(value)
        if key not in seen:
            seen.add(key)
            out.append(value)
    return out


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


def accession_base(acc: str) -> str:
    return acc.split(".")[0] if acc else ""


def extract_accessions_from_blob(obj: Any) -> List[str]:
    return uniq(ACCESSION_RE.findall(json.dumps(obj, ensure_ascii=False)))


def extract_biosample_candidates(obj: Any) -> List[str]:
    found: List[str] = []

    def walk(x: Any) -> None:
        if isinstance(x, dict):
            for k, v in x.items():
                lk = str(k).lower()
                if lk in {"biosampleaccession", "biosample_acc", "biosample"} and isinstance(v, str):
                    m = BIOSAMPLE_RE.search(v.strip())
                    if m:
                        found.append(m.group(0).upper())
                elif lk == "accession" and isinstance(v, str):
                    m = BIOSAMPLE_RE.search(v.strip())
                    if m:
                        found.append(m.group(0).upper())
                walk(v)
        elif isinstance(x, list):
            for item in x:
                walk(item)
        elif isinstance(x, str):
            m = BIOSAMPLE_RE.search(x)
            if m:
                found.append(m.group(0).upper())

    walk(obj)
    return uniq(found)


def parse_fasta(path: Path) -> Dict[str, Dict[str, Any]]:
    index: Dict[str, Dict[str, Any]] = {}
    header: Optional[str] = None
    seq_parts: List[str] = []

    def commit(h: str, parts: List[str]) -> None:
        desc = h[1:].strip()
        acc = desc.split()[0]
        base = accession_base(acc)
        seq = "".join(parts)
        rec = {
            "accession": acc,
            "accession_base": base,
            "header": desc,
            "sequence": seq,
            "length": len(seq),
        }
        index[acc] = rec
        if base and base not in index:
            index[base] = rec

    if not path.exists():
        return index

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    commit(header, seq_parts)
                header = line
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            commit(header, seq_parts)
    return index


def dedupe_records(records: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    seen = set()
    out: List[Dict[str, Any]] = []
    for rec in records:
        key = json.dumps(rec, ensure_ascii=False, sort_keys=True)
        if key not in seen:
            seen.add(key)
            out.append(rec)
    return out


def build_annotation_index(records: Sequence[Dict[str, Any]]) -> DefaultDict[str, List[Dict[str, Any]]]:
    idx: DefaultDict[str, List[Dict[str, Any]]] = defaultdict(list)
    for rec in records:
        for acc in extract_accessions_from_blob(rec):
            idx[acc].append(rec)
            base = accession_base(acc)
            if base:
                idx[base].append(rec)
    return idx


def build_biosample_index(records: Sequence[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    idx: Dict[str, Dict[str, Any]] = {}
    for rec in records:
        for biosample in extract_biosample_candidates(rec):
            idx[biosample] = rec
    return idx


def first_from_contexts(contexts: Sequence[Dict[str, Any]], keys: Sequence[str]) -> Optional[str]:
    for key in keys:
        for ctx in contexts:
            if key in ctx:
                value = clean_text(ctx.get(key))
                if value:
                    return value
    return None


def all_from_contexts(contexts: Sequence[Dict[str, Any]], keys: Sequence[str]) -> List[str]:
    vals: List[str] = []
    for key in keys:
        for ctx in contexts:
            if key in ctx:
                value = ctx.get(key)
                if isinstance(value, list):
                    vals.extend(str(v) for v in value if clean_text(v))
                else:
                    x = clean_text(value)
                    if x:
                        vals.append(x)
    return uniq(vals)


def normalize_host_label(text: Optional[str]) -> Optional[str]:
    if not text:
        return None
    text = re.split(r"[\(\[,:;]", text)[0].strip()
    text = re.sub(r"\s+", " ", text).strip()
    if not text:
        return None
    return text.replace(" ", "_")


def choose_host_candidates(contexts: Sequence[Dict[str, Any]]) -> Tuple[List[str], List[str], List[str]]:
    host_fields = [
        "virus.host.scientificName",
        "host.organismName",
        "host.scientificName",
        "host.organismname",
        "host.scientificname",
        "hostName",
        "host.name",
        "host",
        "biosample.attributes.host",
    ]
    host_taxid_fields = [
        "virus.host.taxId",
        "host.taxId",
        "host.tax_id",
        "biosample.attributes.host_taxid",
    ]
    isolation_source_fields = [
        "biosample.attributes.isolation_source",
        "isolationSource",
        "source",
    ]

    raw_hosts = all_from_contexts(contexts, host_fields)
    host_candidates = uniq([x for x in (normalize_host_label(h) for h in raw_hosts) if x])

    host_tax_ids: List[str] = []
    for value in all_from_contexts(contexts, host_taxid_fields):
        m = re.search(r"\d+", value)
        if m:
            host_tax_ids.append(m.group(0))
    host_tax_ids = uniq(host_tax_ids)

    isolation_sources = all_from_contexts(contexts, isolation_source_fields)
    return host_candidates, host_tax_ids, isolation_sources


def infer_segment_from_text(text: str) -> Tuple[Optional[str], List[str]]:
    warnings: List[str] = []
    low = text.lower()
    if "nss" in low or "ns s" in low:
        warnings.append("contains_nss_label")
        return None, warnings
    for seg, patterns in SEGMENT_RULES.items():
        for pat in patterns:
            if re.search(pat, low, flags=re.I):
                return seg, warnings
    return None, warnings


def choose_segment(contexts: Sequence[Dict[str, Any]], protein_headers: Sequence[str], nt_len: int) -> Tuple[str, str, List[str]]:
    texts: List[str] = []
    warnings: List[str] = []
    for ctx in contexts:
        for k, v in ctx.items():
            lk = k.lower()
            if any(t in lk for t in ["segment", "title", "definition", "description", "product", "protein"]):
                if isinstance(v, list):
                    texts.extend(str(x) for x in v if x is not None)
                else:
                    texts.append(str(v))
    texts.extend(protein_headers)

    for text in texts:
        seg, local_warnings = infer_segment_from_text(text)
        warnings.extend(local_warnings)
        if seg:
            return seg, "text_metadata", uniq(warnings)

    if nt_len:
        if nt_len >= MIN_COMPLETE_NT["L"]:
            return "L", "length_heuristic", uniq(warnings)
        if nt_len >= MIN_COMPLETE_NT["M"]:
            return "M", "length_heuristic", uniq(warnings)
        if nt_len >= MIN_COMPLETE_NT["S"]:
            return "S", "length_heuristic", uniq(warnings)

    warnings.append("segment_unresolved")
    return "", "unresolved", uniq(warnings)


def determine_completeness(segment: str, contexts: Sequence[Dict[str, Any]], raw_blob: str, nt_len: int) -> str:
    structured = first_from_contexts(
        contexts,
        ["completeness", "nucleotideCompleteness", "nucleotide.completeness", "annotation.completeness"],
    )
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
    if segment and nt_len and nt_len >= MIN_COMPLETE_NT.get(segment, 10**9):
        return "likely_complete"
    if nt_len:
        return "partial_or_unknown"
    return "unknown"


def is_complete_status(status: str) -> bool:
    return status in {"complete", "coding_complete_only", "likely_complete"}


def is_refseq(acc: str) -> bool:
    return acc.startswith(REFSEQ_PREFIXES)


def extract_year(text: Optional[str]) -> str:
    if not text:
        return ""
    m = re.search(r"(19|20)\d{2}", text)
    return m.group(0) if m else ""


def clean_for_group(text: Optional[str]) -> str:
    if not text:
        return "unknown"
    text = clean_text(text) or "unknown"
    text = text.casefold()
    text = re.sub(r"\s+", " ", text).strip()
    return text


def build_group_key(biosample: str, isolate: str, strain: str, collection_date: str, location: str, host: str, accession: str) -> Tuple[str, str, float]:
    if biosample:
        return f"biosample|{biosample}", "biosample", 1.0

    isolate_or_strain = isolate or strain
    parts = [
        clean_for_group(isolate_or_strain),
        clean_for_group(collection_date),
        clean_for_group(location),
        clean_for_group(host),
    ]
    non_unknown = sum(1 for x in parts if x != "unknown")
    if non_unknown:
        return "exact|" + "|".join(parts), "isolate_or_strain_collection_date_location_host", round(non_unknown / 4.0, 3)

    return f"accession|{accession}", "accession_fallback", 0.0


def choose_best_group_member_for_segment(members: Sequence[Dict[str, Any]], segment: str) -> Optional[Dict[str, Any]]:
    candidates = [m for m in members if m["segment"] == segment]
    if not candidates:
        return None
    candidates.sort(
        key=lambda r: (
            1 if is_complete_status(r["completeness_status"]) else 0,
            1 if is_refseq(r["accession"]) else 0,
            int(r.get("nucleotide_length", 0)),
            r["accession"],
        ),
        reverse=True,
    )
    return candidates[0]


def segment_completion_pattern(complete_segments: Sequence[str]) -> str:
    if not complete_segments:
        return "none_complete"
    ordered = [seg.lower() for seg in ["S", "M", "L"] if seg in complete_segments]
    return "".join(ordered) + "_complete"


def choose_primary_host(host_candidates: Sequence[str], isolation_sources: Sequence[str]) -> Tuple[str, List[str]]:
    observed = uniq([h for h in host_candidates if h])
    if not observed:
        observed = uniq([x for x in (normalize_host_label(v) for v in isolation_sources) if x])
    if not observed:
        return "unknown", []
    return Counter(observed).most_common(1)[0][0], observed


def predict_segment_host(members: Sequence[Dict[str, Any]], segment: str) -> Optional[str]:
    seg_members = [m for m in members if m["segment"] == segment and is_complete_status(m["completeness_status"])]
    if not seg_members:
        return None
    hosts: List[str] = []
    for member in seg_members:
        if member.get("host_candidates"):
            hosts.extend(member["host_candidates"])
        else:
            hosts.extend(member.get("_isolation_sources", []))
    normalized = uniq([x for x in (normalize_host_label(v) for v in hosts) if x])
    if not normalized:
        return None
    return Counter(normalized).most_common(1)[0][0]


def classify_reassortment(hs: Optional[str], hm: Optional[str], hl: Optional[str]) -> Any:
    available = [x for x in [hs, hm, hl] if x]
    if len(available) == 3:
        return 0 if len(set(available)) == 1 else 1
    if len(available) == 2:
        return 1 if len(set(available)) > 1 else "unknown"
    return "unknown"


def extract_biosample_accession_from_record(record: Dict[str, Any], ann_matches: Sequence[Dict[str, Any]]) -> str:
    candidates: List[str] = []

    candidates.extend(extract_biosample_candidates(record))
    for ann in ann_matches:
        candidates.extend(extract_biosample_candidates(ann))

    data_flat = flatten_json(record)
    candidates.extend([
        first_nonempty(
            data_flat.get("bioSampleAccession"),
            data_flat.get("biosample.accession"),
            data_flat.get("virus.bioSampleAccession"),
            data_flat.get("virus.biosample.accession"),
            data_flat.get("assemblyInfo.biosample.accession"),
        ) or ""
    ])

    for ann in ann_matches:
        ann_flat = flatten_json(ann)
        candidates.extend([
            first_nonempty(
                ann_flat.get("bioSampleAccession"),
                ann_flat.get("biosample.accession"),
                ann_flat.get("virus.bioSampleAccession"),
                ann_flat.get("virus.biosample.accession"),
            ) or ""
        ])

    cleaned = []
    for c in candidates:
        c = clean_text(c)
        if not c:
            continue
        m = BIOSAMPLE_RE.search(c)
        if m:
            cleaned.append(m.group(0).upper())
    cleaned = uniq(cleaned)
    return cleaned[0] if cleaned else ""


def main() -> None:
    ensure_inputs()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/7] Loading JSONL reports...")
    data_records = list(iter_jsonl(DATA_REPORT))
    annotation_records = list(iter_jsonl(ANNOTATION_REPORT))
    biosample_records = list(iter_jsonl(BIOSAMPLE_REPORT))

    print("[2/7] Building indexes...")
    annotation_index = build_annotation_index(annotation_records)
    biosample_index = build_biosample_index(biosample_records)
    genomic_index = parse_fasta(GENOMIC_FASTA)
    cds_index = parse_fasta(CDS_FASTA)
    protein_index = parse_fasta(PROTEIN_FASTA)

    accession_rows: List[Dict[str, Any]] = []
    parse_errors = 0
    lab_dropped = 0

    print("[3/7] Parsing and cross-checking accession records...")
    for i, record in enumerate(data_records, start=1):
        try:
            data_flat = flatten_json(record)
            accession = first_nonempty(
                record.get("accession"),
                data_flat.get("accession"),
                data_flat.get("nucleotideAccession"),
                data_flat.get("annotation.accession"),
            )
            if not accession:
                tokens = extract_accessions_from_blob(record)
                accession = tokens[0] if tokens else "UNKNOWN_ACCESSION"
            accession_key = accession_base(accession)

            ann_matches = dedupe_records(annotation_index.get(accession, []) + annotation_index.get(accession_key, []))
            ann_flats = [flatten_json(x) for x in ann_matches]

            biosample_acc = extract_biosample_accession_from_record(record, ann_matches)
            biosample_rec = biosample_index.get(biosample_acc) if biosample_acc else None
            biosample_flat = flatten_json(biosample_rec) if biosample_rec else {}

            contexts = [data_flat] + ann_flats + ([biosample_flat] if biosample_flat else [])
            raw_blob = json.dumps([record] + ann_matches + ([biosample_rec] if biosample_rec else []), ensure_ascii=False)

            genomic_rec = genomic_index.get(accession) or genomic_index.get(accession_key) or {}
            cds_rec = cds_index.get(accession) or cds_index.get(accession_key) or {}

            protein_accessions: List[str] = []
            protein_lengths_aa: List[int] = []
            protein_headers: List[str] = []

            linked_protein_tokens = []
            for token in extract_accessions_from_blob([record] + ann_matches + ([biosample_rec] if biosample_rec else [])):
                if token in protein_index or accession_base(token) in protein_index:
                    linked_protein_tokens.append(token)
            linked_protein_tokens = uniq(linked_protein_tokens)

            if not linked_protein_tokens:
                iso_hint = first_from_contexts(contexts, ["virus.isolate", "isolate.name", "isolate"])
                strain_hint = first_from_contexts(contexts, ["virus.strain", "strain.name", "strain"])
                iso_l = (iso_hint or "").lower()
                strain_l = (strain_hint or "").lower()
                if iso_l or strain_l:
                    for p_acc, p_rec in protein_index.items():
                        if p_acc != p_rec["accession"]:
                            continue
                        header_l = p_rec["header"].lower()
                        if (iso_l and iso_l in header_l) or (strain_l and strain_l in header_l):
                            linked_protein_tokens.append(p_rec["accession"])
                linked_protein_tokens = uniq(linked_protein_tokens)

            for token in linked_protein_tokens:
                p = protein_index.get(token) or protein_index.get(accession_base(token))
                if p:
                    protein_accessions.append(p["accession"])
                    protein_lengths_aa.append(int(p["length"]))
                    protein_headers.append(p["header"])

            protein_accessions = uniq(protein_accessions)
            protein_lengths_aa = [int((protein_index.get(acc) or protein_index.get(accession_base(acc)) or {}).get("length", 0)) for acc in protein_accessions]
            protein_headers = [((protein_index.get(acc) or protein_index.get(accession_base(acc)) or {}).get("header", "")) for acc in protein_accessions]

            organism_name = first_from_contexts(contexts, ["organism.organismName", "organismName", "taxonomy.organismName", "virus.taxonomy.organismName"]) or ""
            tax_id = first_from_contexts(contexts, ["taxonomy.taxId", "taxId", "organism.taxId", "virus.taxonomy.taxId"]) or ""
            isolate = first_from_contexts(contexts, ["virus.isolate", "isolate.name", "virus.isolateName", "biosample.attributes.isolate", "isolate"]) or ""
            strain = first_from_contexts(contexts, ["virus.strain", "strain.name", "virus.strainName", "biosample.attributes.strain", "strain"]) or ""
            collection_date = first_from_contexts(contexts, ["collectionDate", "virus.collectionDate", "biosample.attributes.collection_date", "collection_date", "isolate.collectionDate"]) or ""
            country_full = first_from_contexts(contexts, ["location.geographicLocation", "virus.geographicLocation", "biosample.attributes.geo_loc_name", "geoLocName", "country"]) or ""
            country = country_full.split(":")[0].strip() if country_full else ""
            region = country_full.split(":", 1)[1].strip() if ":" in country_full else ""
            year = extract_year(collection_date)
            isolation_source = first_from_contexts(contexts, ["biosample.attributes.isolation_source", "isolationSource", "source"]) or ""
            tissue_specimen_source = first_from_contexts(contexts, ["biosample.attributes.tissue", "biosample.attributes.host_tissue_sampled", "tissue"]) or ""
            biosample = biosample_acc or ""
            bioproject = first_from_contexts(contexts, ["bioProjectAccession", "bioproject.accession", "bioprojects"]) or ""

            host_candidates, host_tax_ids, isolation_sources = choose_host_candidates(contexts)
            host_predict, all_observed_hosts = choose_primary_host(host_candidates, isolation_sources)

            nt_len = int(genomic_rec.get("length", 0) or cds_rec.get("length", 0) or 0)
            segment, segment_method, segment_warnings = choose_segment(contexts, protein_headers, nt_len)
            completeness_status = determine_completeness(segment, contexts, raw_blob, nt_len)

            low_blob = raw_blob.lower()
            lab_terms = [term for term in LAB_PASSAGE_TERMS if term in low_blob]
            if lab_terms:
                lab_dropped += 1
                continue

            group_string_id, group_key_type, group_key_confidence = build_group_key(
                biosample=biosample,
                isolate=isolate,
                strain=strain,
                collection_date=collection_date,
                location=country_full,
                host=host_candidates[0] if host_candidates else (host_predict or ""),
                accession=accession,
            )

            accession_rows.append(
                {
                    "accession": accession,
                    "sequencehash": "",
                    "group_string_id": group_string_id,
                    "group_key_type": group_key_type,
                    "group_key_confidence": group_key_confidence,
                    "segment": segment,
                    "segment_assignment_method": segment_method,
                    "completeness_status": completeness_status,
                    "organism_name": organism_name,
                    "tax_id": tax_id,
                    "isolate": isolate,
                    "strain": strain,
                    "collection_date": collection_date,
                    "year": year,
                    "country": country,
                    "region": region,
                    "isolation_source": isolation_source,
                    "tissue_specimen_source": tissue_specimen_source,
                    "biosample": biosample,
                    "bioproject": bioproject,
                    "host_candidates": host_candidates,
                    "host_tax_ids": host_tax_ids,
                    "host_predict": host_predict,
                    "all_observed_hosts": all_observed_hosts,
                    "protein_accessions": protein_accessions,
                    "protein_lengths_aa": protein_lengths_aa,
                    "warning_flags": uniq(segment_warnings),
                    "nucleotide_length": nt_len,
                    "group_size": 0,
                    "solo": 0,
                    "segment_completion_pattern": "",
                    "is_reassortant": "unknown",
                    "group_members": [],
                    "predict_host_S": None,
                    "predict_host_M": None,
                    "predict_host_L": None,
                    "_isolation_sources": isolation_sources,
                }
            )

            if i % 3000 == 0:
                print(f"   parsed {i:,} records...")
        except Exception as exc:
            parse_errors += 1
            sys.stderr.write(f"[WARN] parse failure: {exc}\n")

    print("[4/7] Grouping accession records...")
    grouped: DefaultDict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in accession_rows:
        grouped[row["group_string_id"]].append(row)

    print("[5/7] Applying quality control and building group annotations...")
    final_json_rows: List[Dict[str, Any]] = []
    final_csv_rows: List[Dict[str, Any]] = []

    group_pattern_counts = Counter()
    record_pattern_counts = Counter()
    reassortant_group_count = 0

    kept_groups = 0
    for gid, members in grouped.items():
        complete_segments = sorted(
            {m["segment"] for m in members if m["segment"] in {"S", "M", "L"} and is_complete_status(m["completeness_status"])},
            key=lambda x: {"S": 0, "M": 1, "L": 2}[x],
        )

        if set(complete_segments) != {"S", "M", "L"}:
            continue

        kept_groups += 1
        pattern = segment_completion_pattern(complete_segments)
        group_pattern_counts[pattern] += 1

        group_size = len(members)
        group_members = sorted(m["accession"] for m in members)
        solo = 1 if group_size == 1 else 0

        predict_host_s = predict_segment_host(members, "S")
        predict_host_m = predict_segment_host(members, "M")
        predict_host_l = predict_segment_host(members, "L")
        is_reassortant = classify_reassortment(predict_host_s, predict_host_m, predict_host_l)
        if is_reassortant == 1:
            reassortant_group_count += 1

        preferred_s = choose_best_group_member_for_segment(members, "S")
        preferred_m = choose_best_group_member_for_segment(members, "M")
        preferred_l = choose_best_group_member_for_segment(members, "L")

        all_hosts_group = uniq([h for m in members for h in m.get("all_observed_hosts", []) if h])
        primary_host_label = preferred_s.get("host_predict") if preferred_s and preferred_s.get("host_predict") else (all_hosts_group[0] if all_hosts_group else "unknown")

        final_csv_rows.append(
            {
                "group_string_id": gid,
                "primary_host_label": primary_host_label,
                "s_prot_accession": (preferred_s.get("protein_accessions", [""])[0] if preferred_s and preferred_s.get("protein_accessions") else ""),
                "m_prot_accession": (preferred_m.get("protein_accessions", [""])[0] if preferred_m and preferred_m.get("protein_accessions") else ""),
                "l_prot_accession": (preferred_l.get("protein_accessions", [""])[0] if preferred_l and preferred_l.get("protein_accessions") else ""),
                "complete_members": 3,
                "is_reassortant": is_reassortant,
                "segment_completion_pattern": pattern,
                "group_size": group_size,
                "all_observed_hosts": "|".join(all_hosts_group),
                "country": preferred_s.get("country", "") if preferred_s else (members[0].get("country", "") if members else ""),
                "year": preferred_s.get("year", "") if preferred_s else (members[0].get("year", "") if members else ""),
            }
        )

        for member in members:
            rec = dict(member)
            rec["group_size"] = group_size
            rec["solo"] = solo
            rec["segment_completion_pattern"] = pattern
            rec["is_reassortant"] = is_reassortant
            rec["group_members"] = group_members
            rec["predict_host_S"] = predict_host_s
            rec["predict_host_M"] = predict_host_m
            rec["predict_host_L"] = predict_host_l
            final_json_rows.append(rec)
            record_pattern_counts[pattern] += 1

    print("[6/7] Writing outputs...")
    json_columns = [
        "accession",
        "sequencehash",
        "group_string_id",
        "group_key_type",
        "group_key_confidence",
        "segment",
        "segment_assignment_method",
        "completeness_status",
        "organism_name",
        "tax_id",
        "isolate",
        "strain",
        "collection_date",
        "year",
        "country",
        "region",
        "isolation_source",
        "tissue_specimen_source",
        "biosample",
        "bioproject",
        "host_candidates",
        "host_tax_ids",
        "host_predict",
        "all_observed_hosts",
        "protein_accessions",
        "protein_lengths_aa",
        "warning_flags",
        "nucleotide_length",
        "group_size",
        "solo",
        "segment_completion_pattern",
        "is_reassortant",
        "group_members",
        "predict_host_S",
        "predict_host_M",
        "predict_host_L",
    ]

    with OUT_JSON.open("w", encoding="utf-8") as handle:
        json.dump([{col: safe_jsonable(row.get(col)) for col in json_columns} for row in final_json_rows], handle, ensure_ascii=False, indent=2)

    csv_columns = [
        "group_string_id",
        "primary_host_label",
        "s_prot_accession",
        "m_prot_accession",
        "l_prot_accession",
        "complete_members",
        "is_reassortant",
        "segment_completion_pattern",
        "group_size",
        "all_observed_hosts",
        "country",
        "year",
    ]

    with OUT_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=csv_columns)
        writer.writeheader()
        for row in final_csv_rows:
            writer.writerow({col: row.get(col, "") for col in csv_columns})

    summary = {
        "total_records": len(final_json_rows),
        "groups": kept_groups,
        "solo_records": sum(1 for r in final_json_rows if r["solo"] == 1),
        "segment_counts": dict(Counter(r["segment"] or "unknown" for r in final_json_rows)),
        "completeness_counts": dict(Counter(r["completeness_status"] for r in final_json_rows)),
        "segment_completion_patterns_group_level": dict(group_pattern_counts),
        "segment_completion_patterns_record_level": dict(record_pattern_counts),
        "group_key_type_distribution": dict(Counter(r["group_key_type"] for r in final_json_rows)),
        "reassortant_group_count": reassortant_group_count,
        "reassortant_record_count": sum(1 for r in final_json_rows if r["is_reassortant"] == 1),
        "top_host_predictions": dict(Counter(r["host_predict"] for r in final_json_rows if r.get("host_predict")).most_common(25)),
        "output_files": [str(OUT_JSON), str(OUT_CSV), str(OUT_SUMMARY)],
        "final_json_columns": json_columns,
        "final_csv_columns": csv_columns,
        "records_before_grouping": len(accession_rows),
        "groups_before_filtering": len(grouped),
        "parse_errors": parse_errors,
        "lab_dropped": lab_dropped,
    }

    with OUT_SUMMARY.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)

    print("[7/7] Verification report")
    print("=" * 72)
    print(f"Records before grouping: {len(accession_rows):,}")
    print(f"Groups before filtering: {len(grouped):,}")
    print(f"Groups after filtering:  {kept_groups:,}")
    print(f"Final JSON records:      {len(final_json_rows):,}")
    print(f"Lab-dropped records:     {lab_dropped:,}")
    print(f"Parse errors:            {parse_errors:,}")
    print(f"Reassortant groups:      {reassortant_group_count:,}")
    print(f"JSON output:             {OUT_JSON}")
    print(f"CSV output:              {OUT_CSV}")
    print(f"Summary output:          {OUT_SUMMARY}")
    print("=" * 72)


if __name__ == "__main__":
    main()