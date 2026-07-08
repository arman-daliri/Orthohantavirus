# Orthohantavirus NCBI Dataset Preparation Pipeline


`all_hanta.py` — version 1.0.0

A reproducible pipeline for downloading **NCBI Orthohantavirus** genomes, integrating metadata, reconstructing complete viral isolates (S/M/L segments), filtering laboratory-derived sequences, and generating machine learning–ready datasets for host prediction and reassortment analysis.

Parses an NCBI Datasets Orthohantavirus package and produces accession-level,
isolate-level, virus-level, host-level, machine-learning, and QC outputs.
The pipeline classifies every accession without discarding incomplete records,
groups accessions into isolates conservatively, and emits ten Parquet/JSON files
suitable for downstream phylogenetic, epidemiological, and ML analysis.

---

---

# Step 0 — Download Orthohantavirus Data from NCBI

Download all available **Orthohantavirus** genome records together with their associated metadata, coding sequences (CDS), proteins, and BioSample information.

```bash
datasets download virus genome taxon 1980442   
--include genome,cds,protein,annotation,biosample
--filename all_hanta.zip
```

Extract the downloaded package:

```bash
unzip all_hanta.zip -d all_hanta
```

---

## Pipeline Overview

The pipeline integrates genomic sequences, annotations, proteins, and BioSample metadata into a
unified dataset. It performs metadata reconciliation, segment identification, sequence quality
assessment, isolate reconstruction, laboratory sequence flagging, representative sequence selection,
host role classification, and reassortment candidate detection.

**Design principles:**

- `data_report.jsonl` is the primary accession-level metadata source.
- `biosample_report.jsonl` is the preferred host metadata source.
- `annotation_report.jsonl` supports QC but does not override explicit segment metadata.
- Accessions are grouped conservatively using species + BioSample + isolate/strain first, with
  controlled fallbacks.
- Incomplete viruses are retained and classified, not discarded.
- Reassortment candidates are flagged via taxonomic incongruence, not host label comparison.

---

---

## Step 1 — Prepare the NCBI Package

Download the Orthohantavirus dataset from NCBI Datasets and unzip it so the
directory structure matches:

```
<root>/
└── all_hanta/
    └── ncbi_dataset/
        └── data/
            ├── data_report.jsonl          # required
            ├── annotation_report.jsonl    # required
            ├── biosample_report.jsonl     # required
            ├── dataset_catalog.json       # required
            ├── genomic.fna                # optional
            ├── cds.fna                    # optional
            └── protein.faa               # optional
```

The current script implementation expects this nested layout exactly:
`<root>/all_hanta/ncbi_dataset/data`. This is the behavior of `resolve_inputs()`
even though the script docstring may refer more generally to an unzipped package
root.

The four JSONL/JSON files are **required**. The pipeline raises `FileNotFoundError`
if any of them is absent. The three FASTA files are **optional**; when absent they
are silently treated as empty indexes. Sequence-based completeness checks (length
thresholds, terminal motif detection, ambiguous-N counting) are skipped for any
accession with no FASTA entry — those accessions fall back to text-based and
structured-metadata completeness only.

---

## Step 2 — Run Python

```bash
python all_hanta.py --input /path/to/root [--output /path/to/output] [--threads N] [--verbose]
```

### CLI Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--input` | No | `.` (current directory) | Path to the unzipped NCBI package root containing `all_hanta/ncbi_dataset/data` |
| `--output` | No | `<root>/output/` | Directory to write all output files |
| `--threads` | No | auto (1–8, based on CPU count) | Number of parallel worker threads for accession transformation |
| `--verbose` | No | off | Enable DEBUG-level logging to console and log file |

Thread safety note: each accession record is transformed independently. All shared
data structures (annotation index, biosample index, FASTA indexes) are built once
before the thread pool and used read-only. Results are collected in submission order,
so output is deterministic regardless of thread count.

---

## Step 3 — Outputs

All outputs are written to the output directory. The pipeline produces ten files:

```
output/
├── ortho_accession_table.parquet
├── ortho_virus_table.parquet
├── ortho_isolate_table.parquet
├── ortho_host_table.parquet
├── ortho_ml_feature_matrix.parquet
├── ortho_ml_segments.parquet
├── ortho_suppressed.parquet
├── ortho_qc_report.parquet
├── ortho_summary.json
└── ortho_pipeline.log
```

---

## Data Flow

```
data_report.jsonl ──────────────────────────────────────────────────┐
annotation_report.jsonl ──► build indexes ──► transform_one_record ─┤
biosample_report.jsonl ─────────────────────(parallel, 1 per acc.) ─┤
genomic.fna / cds.fna / protein.faa (optional) ─────────────────────┘
                                                       │
                                       ┌───────────────▼──────────────────┐
                                       │      AccessionRecord[]            │
                                       │  (all input accessions retained)  │
                                       └──┬──────────┬────────────┬────────┘
                                          │          │            │
                              group_      │   build_ │   build_   │
                              accessions  │   strict_ │  host_    │
                                          │   isolate │  table    │
                                          ▼   _records▼           ▼
                                    VirusGroup[] IsolateRecord[] HostRecord[]
                                          │          │
                                          │   build_ml_feature_matrix (amino acid 3-mers)
                                          │   build_ml_segments
                                          │   build_suppressed_rows
                                          │
                                          └─► 10 output files
```

An accession is always in `ortho_accession_table.parquet`. It appears in
`ortho_isolate_table.parquet` only if its group contains at least one accession
with a resolved segment (S, M, or L). Groups where no accession resolves to any
segment are excluded entirely from the isolate table. Any accession not chosen as
the best segment representative for its isolate appears in `ortho_suppressed.parquet`.

---

## Output File Reference

### `ortho_accession_table.parquet`

One row per NCBI accession. This is the primary flat table containing all
metadata resolved from `data_report.jsonl`, `biosample_report.jsonl`, and
`annotation_report.jsonl`.

| Column | Description |
|---|---|
| `accession` | Full versioned accession (e.g. `NC_005218.1`) |
| `species_name` | Virus organism name |
| `species_taxid` | NCBI taxonomy ID for the virus species |
| `genus_name` | Genus name (always `Orthohantavirus`) |
| `genus_taxid` | Genus taxonomy ID |
| `segment` | Resolved segment: `S`, `M`, `L`, or empty if unresolved |
| `segment_source` | How the segment was determined: `explicit_metadata`, `annotation_text`, `length_heuristic`, or `unresolved` |
| `completeness` | Completeness classification (see Completeness section) |
| `nt_length` | Nucleotide sequence length in bp (0 if no FASTA entry) |
| `biosample` | BioSample accession (SAM…) |
| `bioproject` | BioProject accession |
| `isolate` | Isolate name from metadata or BioSample |
| `strain` | Strain name |
| `collection_date` | Raw collection date string |
| `year` | Four-digit year extracted from collection date |
| `location` | Full geographic location string |
| `country` | Country extracted from location |
| `region` | Sub-country region extracted from location |
| `host_name` | Host organism name |
| `host_taxid` | Host NCBI taxonomy ID |
| `host_source` | Where host was resolved from: `biosample`, `data_report`, or `missing` |
| `host_role` | Host role classification (see Host Role section) |
| `isolation_source` | Raw isolation source text from BioSample |
| `tissue_source` | Tissue type from BioSample |
| `protein_accessions` | Semicolon-separated protein accessions linked via annotation |
| `annotation_products` | Semicolon-separated product names from annotation |
| `group_key` | Grouping key used for isolate assembly (see Grouping section) |
| `group_confidence` | Confidence score of the grouping tier (0.10–1.00) |
| `is_refseq` | 1 if accession starts with `NC_`, `NM_`, `NP_`, `YP_`, etc. |
| `warning_flags` | Semicolon-separated QC warning flags (see QC Flags section) |

---

### `ortho_virus_table.parquet`

One row per virus group (grouped accession set). Provides a virus-level view
with segment lists and completeness class. Fields include `group_key`,
`species_name`, `species_taxid`, `host_name`, `host_taxid`, `country`, `year`,
`biosample`, `isolate`, `strain`, `group_confidence`, `group_size`,
`s_accessions`, `m_accessions`, `l_accessions` (all semicolon-separated),
`complete_segments`, `completeness_class`, `all_hosts`, `all_accessions`.

`completeness_class` takes three values:

| Value | Meaning |
|---|---|
| `complete_triplet` | S, M, and L all have complete-status segments |
| `partial_group_with_complete_segment` | At least one complete-status segment |
| `no_complete_segment` | No segment meets completeness criteria |

---

### `ortho_isolate_table.parquet`

One row per strict isolate — the best S, M, and L representative chosen per
group. **A group is included in this table only if at least one of its
accessions has a resolved segment (S, M, or L).** Groups whose accessions all
have `segment = ""` are skipped entirely. For each retained group, the best
representative for each segment is chosen by sorting on
`(completeness_rank, nt_length)` descending, where the rank order is:
`complete` = `coding_complete_only` = 4, `likely_complete` = 3,
`partial_or_unknown` = 2, `unknown` = 1. A segment slot is left empty (`""`)
for any group that has no accession for that segment.

| Column | Description |
|---|---|
| `strict_isolate_key` | `species_taxid\|biosample\|isolate_or_strain`, or group key if those fields are absent |
| `species_name` | Consensus species name across segments |
| `species_taxid` | Consensus species taxonomy ID |
| `species_cluster_id` | Protein-distance cluster ID, e.g. `2697049.cluster0` |
| `biosample` | BioSample accession |
| `isolate` | Isolate name |
| `strain` | Strain name |
| `host_name` | Host organism name |
| `host_taxid` | Host taxonomy ID |
| `host_role` | Host role classification (see Host Role section) |
| `collection_date` | Collection date |
| `location` | Full geographic location |
| `country` | Country |
| `region` | Sub-country region |
| `grouping_tier` | Grouping tier used (e.g. `tier_1`) |
| `grouping_confidence` | Confidence score of the grouping tier |
| `s_accession` | Best S-segment accession for this isolate |
| `m_accession` | Best M-segment accession |
| `l_accession` | Best L-segment accession |
| `s_completeness` | Completeness of the S representative (`unknown` if no S) |
| `m_completeness` | Completeness of the M representative |
| `l_completeness` | Completeness of the L representative |
| `s_nt_length` | Nucleotide length of S representative (0 if no S) |
| `m_nt_length` | Nucleotide length of M representative |
| `l_nt_length` | Nucleotide length of L representative |
| `concatenated_protein_aa` | N protein + GPC protein amino acid sequences concatenated; empty if either is absent |
| `coding_complete_status` | Triplet completeness summary (see below) |
| `grouping_tier` | Tier used for grouping |
| `grouping_confidence` | Numeric confidence for the grouping tier |
| `s_protein_evidence` | Product label matched for S (`N`, or empty) |
| `m_protein_evidence` | Product label matched for M (`GPC`, or empty) |
| `l_protein_evidence` | Product label matched for L (`RdRp`, or empty) |
| `s_protein_accessions` | Semicolon-separated protein accessions for S |
| `m_protein_accessions` | Semicolon-separated protein accessions for M |
| `l_protein_accessions` | Semicolon-separated protein accessions for L |
| `missing_protein_warnings` | `missing_n_or_gpc_sequence` if N or GPC is absent |
| `mosaic_flag` | 1 if segments disagree on collection date or location |
| `taxonomic_incongruence_warning` | 1 if segments disagree on species taxid |
| `reassortment_candidate_warning` | 1 when taxonomic incongruence is detected |
| `ns_s_overlapping_orf` | 1 if NSs overlapping ORF was found in S annotation |
| `species_divergence_flag` | 1 if isolate belongs to a non-zero PED cluster |
| `warnings` | All combined warnings for this isolate |

`coding_complete_status` values:

| Value | Meaning |
|---|---|
| `coding_complete_triplet` | All three segments pass `is_strict_complete_segment()` |
| `near_complete_triplet` | Two of three segments pass |
| `partial_triplet` | Exactly one segment passes |
| `incomplete_triplet` | No segment passes |

---

### `ortho_host_table.parquet`

One row per unique `(host_name, host_taxid)` pair, aggregated across all
accessions and isolates. Columns: `host_name`, `host_taxid`, `host_role`,
`host_source`, `representative_biosample`, `representative_isolate`,
`accession_count`, `isolate_count`, `example_species`.

---

### `ortho_ml_feature_matrix.parquet`

One row per isolate that has a non-empty `concatenated_protein_aa`. Features
are amino acid 3-mers (tripeptides, k=3) computed over the concatenated N + GPC
protein sequence. Up to 200 features are selected by a correlation-based ranking:
each 3-mer is scored by covariance-over-variance relative to encoded host/species
labels, and features with Pearson |r| > 0.92 against already-selected features are
pruned (`select_correlation_features`, `CORRELATION_MAX_FEATURES = 200`). Selected
features are z-score normalised column-wise. Columns include the six metadata fields
`strict_isolate_key`, `species_name`, `species_taxid`, `host_name`, `host_taxid`,
`host_role`, followed by up to 200 `kmer_XYZ` columns.

---

### `ortho_ml_segments.parquet`

One row per accession that has a resolved segment (`S`, `M`, or `L`) **and** a
non-empty `host_name`. Accessions without a host name are excluded. Columns:
`accession`, `segment`, `species_name`, `species_taxid`, `host_name`, `host_taxid`,
`country`, `year`, `biosample`, `completeness`, `group_key`. Designed for
per-segment ML and classification tasks.

---

### `ortho_suppressed.parquet`

All accessions that were **not** chosen as the best representative for any
isolate segment. An accession ends up here either because a higher-quality
representative was chosen for the same segment in its group, or because its
group was excluded from the isolate table entirely (no resolved segments).
Columns: `accession`, `species_name`, `species_taxid`, `segment`, `completeness`,
`biosample`, `isolate`, `strain`, `collection_date`, `location`, `host_name`,
`host_taxid`, `suppression_reason`. The `suppression_reason` is
`not_best_representative` when no QC flags are present, or a semicolon-separated
list of QC flag names otherwise.

---

### `ortho_qc_report.parquet`

One row per accession with QC metadata. Columns: `accession`, `group_key`,
`grouping_tier`, `group_confidence`, `segment`, `completeness`, `warnings`.
Malformed accession records that raise an exception during transformation produce
a minimal error row with `warnings = "transform_error:<message>"`.

---

### `ortho_summary.json`

A JSON file containing pipeline-level statistics and provenance:

| Field | Description |
|---|---|
| `records_total` | Total number of accession records processed |
| `groups_total` | Number of virus groups formed |
| `records_with_biosample` | Accessions with a resolved BioSample |
| `records_with_host` | Accessions with a resolved host name |
| `segment_counts` | Count per segment label (`S`, `M`, `L`, `unresolved`) |
| `completeness_counts` | Count per completeness status |
| `group_completeness_counts` | Count per `completeness_class` |
| `qc_warning_counts` | Count per individual QC flag |
| `provenance.generated_at` | ISO 8601 UTC timestamp |
| `provenance.dataset_catalog_api_version` | API version from `dataset_catalog.json` |
| `provenance.dataset_catalog_md5` | MD5 of `dataset_catalog.json` |
| `provenance.input_md5` | MD5 checksums of all four required input files |
| `provenance.input_files` | Absolute paths of all four required input files |
| `output_files` | Absolute paths to each of the ten output files |
| `output_schema` | Column lists for every output table |

---

### `ortho_pipeline.log`

Timestamped log written to the output directory. Contains INFO-level progress
messages by default. Use `--verbose` to set the logger to DEBUG level. The
current implementation does not emit detailed per-record DEBUG messages unless
additional `logging.debug(...)` statements are added. Malformed JSON in any JSONL
file raises a `ValueError` that includes the file name and line number; it is not
separately logged before the pipeline aborts.

---

## Segment Identification

Segment assignment follows a three-step priority cascade (`choose_segment()`):

1. **Explicit metadata** — `segment`, `virus.segment`, or `annotation.segment`
   fields in `data_report.jsonl` or `annotation_report.jsonl`. Values `S`,
   `M`, `L`, `SMALL`, `MEDIUM`, `LARGE` are all accepted.
2. **Annotation text** — product names and field values matching the keys
   `segment`, `description`, `definition`, `title`, `product`, `protein` are
   scored against regular-expression patterns for each segment. The segment with
   the highest score wins. If multiple segments score non-zero, a
   `multiple_segment_text_hits` warning is added.
3. **Length heuristic** — if no text evidence resolves a segment, the nucleotide
   length is compared to minimum thresholds: L ≥ 6000 bp, M ≥ 3200 bp, S ≥ 1200 bp.
   Triggers `segment_from_length_only` and `segment_from_heuristic` warnings.

Records that pass none of these steps have `segment = ""` and receive a
`segment_unresolved` warning. These records appear in `ortho_accession_table`
and `ortho_suppressed` but **not** in `ortho_isolate_table`.

---

## Completeness Classification

Seven completeness statuses are assigned in priority order (`determine_completeness()`):

| Status | Criteria |
|---|---|
| `complete` | Structured metadata field says "complete", or text contains "complete genome" / "complete sequence" |
| `partial` | Structured metadata field says "partial" |
| `coding_complete_only` | Text contains "complete CDS" without "partial" |
| `coding_complete` | Length ≥ segment minimum **and** terminal motif detected (5′ `UAGUAGUAG` or 3′ `AUCAUCAUCUG`, RNA or DNA equivalent) |
| `likely_complete` | Length ≥ segment minimum but no terminal motif |
| `partial_or_unknown` | Length data are present but either ambiguous N count > 100, or the accession does not satisfy structured/text completeness, terminal-motif, or segment-length criteria |
| `unknown` | No length data available |

The four statuses accepted as "complete" for isolate assembly and
`completeness_class` scoring are: `complete`, `coding_complete_only`,
`likely_complete`, and `coding_complete`.

---

## Host Resolution

Host name and taxid are resolved by `choose_host()` in strict priority order:

1. **BioSample** — `biosample_report.jsonl` `host` and `host_taxid` attributes.
   Used first if a BioSample accession links to a record with a non-empty host.
2. **Data report** — `host.organismName`, `host.scientificName`, `hostName`, or
   `host` fields from `data_report.jsonl`. Used as fallback.

The resolved `host_source` field records which source was used (`biosample`,
`data_report`, or `missing`).

---

## Isolate Grouping Tiers

Accessions are grouped into isolates using a six-tier key cascade
(`build_group_key()`). Higher tiers have higher confidence:

| Tier | Key | Confidence |
|---|---|---|
| `tier_1` | species + biosample + isolate/strain | 1.00 |
| `tier_2` | species + biosample | 0.92 |
| `tier_3` | species + isolate/strain + date + location | 0.80 |
| `tier_4` | species + isolate/strain + location | 0.72 |
| `tier_5` | species + strain + host + location | 0.62 |
| `fallback` | accession only | 0.10 |

Within each group, the best representative for each segment is chosen by
`choose_best_segment_record()`, sorting on `(completeness_rank, nt_length)`
descending. Completeness ranks: `complete` = `coding_complete_only` = 4,
`likely_complete` = 3, `partial_or_unknown` = 2, `unknown` = 1.

---

## Host Role Classification

Each accession and isolate is assigned one of six host roles by
`determine_host_role()`, evaluated in strict priority order:

| Role | Assignment criteria |
|---|---|
| `spillover_host` | Host taxid = 9606, or name contains "human" / "Homo sapiens" |
| `laboratory_host` | Name, source, isolation text, or raw blob contains any lab keyword: `vero`, `cell line`, `cell-line`, `cell culture`, `passage`, `passaged`, `adapted`, `in vitro` |
| `reservoir_candidate` | Isolation source or raw blob contains excreta terms: `urine`, `feces`, `faeces`, `saliva`, `stool`, `fecal`, `urinary`, `salivary`, `excreta` |
| `potential_host` | Isolation source contains organ terms: `salivary gland`, `kidney`, `liver`, `lung`; or any named host not matching the criteria above |
| `spillover_host` (serology) | Isolation source or raw blob contains serology terms: `serum`, `serology`, `antibody`, `igg`, `igm`, `sero` (evaluated after organ check) |
| `unknown` | No host name present |

---

## QC Warning Flags

The following flags appear in `warning_flags` (accession table) and `warnings`
(isolate table):

| Flag | Meaning |
|---|---|
| `missing_segment` | Segment could not be resolved |
| `segment_from_heuristic` | Segment was assigned by length alone |
| `segment_annotation_conflict` | Annotation and metadata disagree on segment |
| `contains_nss_label` | NSs annotation term found in text evidence |
| `multiple_segment_text_hits` | Text evidence matches more than one segment |
| `segment_from_length_only` | Set alongside `segment_from_heuristic` |
| `segment_unresolved` | All three resolution methods failed |
| `missing_host` | No host name found |
| `missing_biosample` | No BioSample accession found |
| `missing_isolate_and_strain` | Neither isolate nor strain present |
| `lab_host_flag` | Record contains lab-host keyword in raw JSON |
| `missing_genomic_sequence` | `genomic.fna` has no entry for this accession |
| `missing_protein_link` | No protein accessions linked via annotation |
| `unexpected_length_for_segment` | Length is outside the expected range for the assigned segment (S > 4000 bp; M < 1500 or > 5000 bp; L < 3000 bp) |

---

## Isolate-Level Flags

| Flag | Meaning |
|---|---|
| `mosaic_flag` | Segments disagree on collection date or location |
| `taxonomic_incongruence_warning` | Segments disagree on species taxid |
| `reassortment_candidate_warning` | Set when taxonomic incongruence is detected |
| `ns_s_overlapping_orf` | NSs overlapping ORF detected in S-segment annotation |
| `species_divergence_flag` | Isolate is in a non-zero protein-distance cluster within its species |

---

## Protein Distance Clustering

Isolates within the same species are clustered by `cluster_species_by_ped()`,
which applies a union-find graph walk over pairwise `compute_protein_distance()`
scores. The distance is a character-level Hamming approximation on unaligned
concatenated amino acid sequences (N protein + GPC). This is a fast heuristic
for coarse grouping — it is **not** a substitution-model-corrected phylogenetic
distance and should be reported as an approximation in any publication. The
threshold is `PED_SPECIES_THRESHOLD = 0.1`. Isolates in cluster 0 are canonical;
non-zero cluster membership sets `species_divergence_flag = 1`.
`species_cluster_id` is formatted as `<species_taxid>.cluster<N>`, e.g.
`2697049.cluster0`.

---

## Known Limitations

- **Protein distance approximation** — `compute_protein_distance()` uses unaligned
  character-level Hamming distance on amino acid sequences. Results are suitable for
  coarse clustering only and must not be reported as phylogenetic distances.
- **No reassortment inference** — `reassortment_candidate_warning` is set only when
  S, M, and L segments of the same group disagree on species taxid. It is not inferred
  from host labels or geography.
- **FASTA dependency for sequence-based completeness** — `coding_complete` and
  `likely_complete` statuses require `genomic.fna` or `cds.fna`. Without these files,
  accessions that would qualify are downgraded to `partial_or_unknown` or `unknown`.
- **Single-source biosample** — if two BioSample accessions in the same record refer
  to different entries in `biosample_report.jsonl`, only the first resolved one is used.

---

## Design Principles

1. `data_report.jsonl` is the primary accession-level metadata source.
2. `biosample_report.jsonl` is the preferred host metadata source; it takes
   precedence over `data_report.jsonl` for host fields.
3. `annotation_report.jsonl` is supporting evidence for QC and product
   identification, not the first source of segment identity.
4. Accessions are grouped conservatively — biosample + isolate/strain first,
   then controlled fallbacks.
5. Incomplete records are kept and classified; none are discarded.
6. Reassortment is not inferred from host labels alone.
7. Sequence files (`genomic.fna`, `cds.fna`, `protein.faa`) are used as
   attachment and validation layers, not as primary metadata sources.

---

## Reproducibility

Every run writes MD5 checksums of all four required input files to
`ortho_summary.json` under `provenance.input_md5`. Results are deterministic
given the same input files and Python version. Thread count (`--threads`) does
not affect output content because accession records are transformed independently
using read-only shared indexes and results are collected in submission order.
