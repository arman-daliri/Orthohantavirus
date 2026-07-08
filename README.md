# Orthohantavirus Dataset Preparation Pipeline

A reproducible pipeline for downloading **NCBI Orthohantavirus** genomes, integrating metadata, reconstructing complete viral isolates (S/M/L segments), filtering laboratory-derived sequences, and generating machine learning–ready datasets for host prediction and reassortment analysis.

---

# Step 1 — Download Orthohantavirus Data from NCBI

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

## Step 2 — Run Pthyon to Build the ML-Ready Dataset

Run the pipeline by passing the path to the unzipped package root:

```bash
python3 all_hanta.py
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

## Required Input Files

```
<input>/all_hanta/ncbi_dataset/data/

data_report.jsonl
annotation_report.jsonl
biosample_report.jsonl
dataset_catalog.json
genomic.fna
cds.fna
protein.faa
```
---

## Data Integration

The pipeline performs the following operations:

- Loads all JSONL reports with strict JSON validation.
- Parses genomic, CDS, and protein FASTA files.
- Matches nucleotide records with annotation and BioSample entries.
- Builds unified metadata for every accession.
- Extracts:
  - Organism name and taxonomy ID
  - Genus name and genus taxonomy ID
  - Isolate and strain
  - Collection date and year
  - Geographic origin (country and sub-region)
  - Host name and host taxonomy ID
  - Host role classification
  - Tissue and isolation source
  - BioSample and BioProject accessions
  - Protein accessions and annotation products
  - RefSeq status

Host information is sourced from BioSample first, then from `data_report.jsonl`. Host labels are
classified into roles using `determine_host_role()` (see Host Role Classification below).

Protein sequences are linked via annotation accession cross-references. Segment-specific protein
products (nucleocapsid N, glycoprotein precursor GPC, RNA-dependent RNA polymerase RdRp) are
resolved and attached to each accession record.

---

## Segment Identification

`choose_segment()` determines whether each genome belongs to the **S**, **M**, or **L** segment,
using a three-step priority cascade.

### 1. Explicit metadata (highest priority)

The segment field is read directly from `data_report.jsonl` or annotation metadata. Values S, M, L,
SMALL, MEDIUM, and LARGE are mapped unambiguously.

### 2. Keyword-based assignment from text evidence

If no explicit field is found, product names and description fields are scored against patterns:

| Segment | Evidence keywords |
| ------- | ---------------------------------------------------------------- |
| S       | segment s, small, nucleocapsid, nucleoprotein, n protein, np     |
| M       | segment m, medium, glycoprotein precursor, glycoprotein, gpc, gn, gc, g1, g2 |
| L       | segment l, large, polymerase, RdRp, RNA-dependent RNA polymerase, l protein |

### 3. Length-based fallback (lowest priority)

If text evidence is unavailable, segment is inferred from genome length. A QC flag
`segment_from_length_only` is attached to any record assigned this way.

| Segment | Minimum Length |
| ------- | -------------: |
| L       | ≥ 6000 bp      |
| M       | ≥ 3200 bp      |
| S       | ≥ 1200 bp      |

Records where segment cannot be resolved by any method receive segment `""` and flag
`segment_unresolved`.

---

## Sequence Completeness Assessment

`determine_completeness()` assigns one of **seven** completeness categories:

| Status               | Description                                                             |
| -------------------- | ----------------------------------------------------------------------- |
| `complete`           | Explicitly reported as complete genome or complete sequence             |
| `partial`            | Explicitly reported as partial                                          |
| `coding_complete_only` | Complete CDS confirmed but genome not complete                        |
| `coding_complete`    | Length exceeds segment threshold **and** conserved terminal motif detected in sequence |
| `likely_complete`    | Length exceeds segment threshold (no terminal motif confirmation)       |
| `partial_or_unknown` | Length available but evidence is inconclusive, or ambiguous-N count exceeds 100 |
| `unknown`            | No completeness evidence available                                      |

For downstream analyses and strict isolate assembly, the following categories are considered
acceptable (complete):

- `complete`
- `coding_complete_only`
- `coding_complete`
- `likely_complete`

---

## Host Role Classification

`determine_host_role()` assigns each accession a biological role based on host name, isolation
source, and metadata context:

| Role                  | Assignment criteria                                                        |
| --------------------- | -------------------------------------------------------------------------- |
| `spillover_host`      | Human host (taxid 9606, "human", "Homo sapiens") or serological sample    |
| `laboratory_host`     | Host or source contains lab-passage terms (see Laboratory Filtering below) |
| `reservoir_candidate` | Isolation source contains excreta terms (urine, feces, saliva, etc.)       |
| `potential_host`      | Host name present but no stronger evidence available                       |
| `unknown`             | No host information                                                        |

---

## Laboratory Sequence Filtering

Records with laboratory-associated metadata are flagged with `lab_host_flag` in the QC report and
excluded from strict isolate assembly.

Terms scanned across all combined metadata fields:

```
vero
cell line
cell-line
cell culture
passage
passaged
adapted
in vitro
```

---

## Isolate Reconstruction and Grouping

Accessions are assigned to a virus group using `build_group_key()`, which applies a six-tier
priority cascade. All tiers include species taxonomy ID as a prefix to prevent cross-species
collisions.

| Tier     | Key fields                                      | Confidence |
| -------- | ----------------------------------------------- | ---------: |
| tier_1   | species + BioSample + isolate or strain         | 1.00       |
| tier_2   | species + BioSample                             | 0.92       |
| tier_3   | species + isolate/strain + collection date + location | 0.80  |
| tier_4   | species + isolate/strain + location             | 0.72       |
| tier_5   | species + strain + host + location              | 0.62       |
| fallback | accession only                                  | 0.10       |

Each accession record carries a `group_key`, `grouping_tier`, and `group_confidence` (0.10–1.00).

---

## Strict Isolate Assembly

`build_strict_isolate_records()` assembles one `IsolateRecord` per group, selecting the best
available accession for each of the S, M, and L segments.

Groups with **no segments at all** are dropped. Groups with one or two segments are retained with
the appropriate `coding_complete_status`:

| Status                   | Meaning                              |
| ------------------------ | ------------------------------------ |
| `coding_complete_triplet`  | All three segments (S, M, L) are complete |
| `near_complete_triplet`    | Two complete segments                |
| `partial_triplet`          | One complete segment                 |
| `incomplete_triplet`       | No segment meets completeness threshold |

### Representative segment selection

For each segment, the best accession is selected by:

1. Highest completeness rank (`complete` / `coding_complete_only` > `likely_complete` > `partial_or_unknown` > `unknown`)
2. Longest nucleotide sequence (tiebreaker)

### Additional isolate-level flags

| Flag                          | Meaning                                                              |
| ----------------------------- | -------------------------------------------------------------------- |
| `mosaic_flag`                 | Segment records within the group differ in collection date or location |
| `reassortment_candidate_warning` | Segment records within the group carry different species taxonomy IDs |
| `taxonomic_incongruence_warning` | Same as above (stored separately for downstream filtering)         |
| `ns_s_overlapping_orf`        | NSs overlapping ORF detected in the S segment (e.g. Puumala, Tula)  |
| `species_divergence_flag`     | Protein distance between isolates exceeds PED species threshold      |

---

## Machine Learning Feature Matrix

`build_ml_feature_matrix()` computes amino-acid tri-mer (k=3) frequency profiles from the
concatenated N + GPC protein sequences of each isolate. Features are then:

1. Ranked by HSIC-Lasso correlation with host taxonomy labels.
2. Filtered for redundancy (Pearson r > 0.92 between selected features is penalised).
3. Truncated to a maximum of 200 features.
4. Z-score normalised across isolates.

Only isolates with both N (S segment) and GPC (M segment) protein sequences contribute to this
matrix.

---

## Output Files

All outputs are written to the directory specified by `--output` (default: `<input>/output`).

### 1. `ortho_accession_table.parquet`

One row per accession. Contains all integrated metadata fields including segment, completeness,
host role, protein accessions, grouping tier, confidence score, and QC warning flags.

### 2. `ortho_virus_table.parquet`

One row per virus group. Contains group-level summary fields: species, host, country, year,
BioSample, isolate, group size, per-segment accession lists, complete segment list, and
`completeness_class`.

### 3. `ortho_isolate_table.parquet`

One row per strict isolate. Contains the assembled S/M/L accessions, per-segment completeness and
length, concatenated protein sequence, coding completeness status, and all isolate-level flags
(`mosaic_flag`, `reassortment_candidate_warning`, `taxonomic_incongruence_warning`,
`ns_s_overlapping_orf`, `species_divergence_flag`).

### 4. `ortho_host_table.parquet`

One row per unique host (host name + taxid). Contains host role, representative BioSample and
isolate, accession count, isolate count, and an example virus species.

### 5. `ortho_ml_feature_matrix.parquet`

One row per isolate that has concatenated protein data. Contains species, host, host role, and
HSIC-Lasso–selected normalised k-mer frequency features (`kmer_<triplet>`).

### 6. `ortho_ml_segments.parquet`

One row per accession that has a known segment and a host name. Lightweight table for per-segment
ML tasks: accession, segment, species, host, country, year, BioSample, completeness, group key.

### 7. `ortho_suppressed.parquet`

Accessions that were not selected as segment representatives for any isolate. Includes the
suppression reason derived from QC warning flags.

### 8. `ortho_qc_report.parquet`

One row per accession. Contains group key, grouping tier, group confidence, segment, completeness,
and semicolon-delimited warning flags.

### 9. `ortho_summary.json`

Pipeline execution summary including:

- Total accession records processed
- Total virus groups assembled
- Records with BioSample and host information
- Segment distribution counts
- Completeness distribution counts
- Group completeness class distribution
- QC warning frequency counts
- Provenance block: generation timestamp, dataset catalog API version, MD5 checksums of all
  input files
- Output file paths and output schema

### 10. `ortho_pipeline.log`

Full execution log written to the output directory. Pass `--verbose` to include DEBUG-level
messages.

---
