# Orthohantavirus Dataset Preparation Pipeline

A reproducible pipeline for downloading **NCBI Orthohantavirus** genomes, integrating metadata, reconstructing complete viral isolates (S/M/L segments), filtering laboratory-derived sequences, and generating machine learning–ready datasets for host prediction and reassortment analysis.

---

# Step 1 — Download Orthohantavirus Data from NCBI

Download all available **Orthohantavirus** genome records together with their associated metadata, coding sequences (CDS), proteins, and BioSample information.

```bash
datasets download virus genome taxon 1980442 \
    --include genome,biosample,cds,protein \
    --filename ortho.zip
```

Extract the downloaded package:

```bash
unzip ortho.zip -d ortho
```

Generate a metadata table:

```bash
dataformat tsv virus-genome \
    --package ortho.zip \
    --fields accession,segment,completeness,length,biosample-acc,'virus-infraspecific-isolate','virus-infraspecific-strain','isolate-collection-date',geo-location,host-name \
    > metadata.tsv
```

---

# Step 2 — Build the ML-Ready Dataset

Run:

```bash
python all_hanta46.py
```

This script transforms the raw NCBI package into curated datasets suitable for downstream machine learning analyses.

---

# Pipeline Overview

The pipeline integrates genomic sequences, annotations, proteins, and BioSample metadata into a unified dataset. It performs metadata reconciliation, segment identification, sequence quality assessment, isolate reconstruction, laboratory sequence filtering, representative sequence selection, host prediction, and reassortment detection.

---

# Required Input Files

The script expects the following files inside its configured working directory:

```
data_report.jsonl
annotation_report.jsonl
biosample_report.jsonl
genomic.fna
cds.fna
protein.faa
```

Before execution, `ensure_inputs()` verifies that all required files exist. The pipeline terminates immediately if any file is missing.

> **Note:** Input paths are currently hard-coded and should be updated if the directory structure changes.

---

# Data Integration

The pipeline performs the following operations:

* Loads all JSONL reports with strict validation.
* Parses genomic, CDS, and protein FASTA files.
* Matches nucleotide records with annotation and BioSample entries.
* Builds unified metadata for every accession.
* Extracts:

  * Organism
  * Taxonomy ID
  * Isolate
  * Strain
  * Collection date
  * Geographic origin
  * Host
  * Tissue / isolation source
  * BioSample accession
  * BioProject accession

Host information is collected from dedicated host fields whenever available and otherwise inferred from isolation-source metadata. Host labels are standardized using `normalize_host_label()`.

Protein sequences are linked primarily by accession number. When no direct match exists, isolate and strain identifiers are used to search protein FASTA headers.

---

# Segment Identification

`choose_segment()` determines whether each genome belongs to the **S**, **M**, or **L** segment.

### Keyword-based assignment

| Segment | Evidence                         |
| ------- | -------------------------------- |
| S       | nucleocapsid, nucleoprotein      |
| M       | glycoprotein-related annotations |
| L       | polymerase-related annotations   |

### Length-based fallback

If annotation keywords are unavailable, segment assignment is inferred from genome length:

| Segment | Minimum Length |
| ------- | -------------: |
| L       |      ≥ 5000 bp |
| M       |      ≥ 2500 bp |
| S       |      ≥ 1000 bp |

---

# Sequence Completeness Assessment

`determine_completeness()` assigns one of six completeness categories:

| Status                 | Description                                           |
| ---------------------- | ----------------------------------------------------- |
| `complete`             | Explicitly reported as complete genome/sequence       |
| `partial`              | Explicitly reported as partial                        |
| `coding_complete_only` | Complete CDS but not complete genome                  |
| `likely_complete`      | Length exceeds the expected threshold for its segment |
| `partial_or_unknown`   | Length available but evidence is inconclusive         |
| `unknown`              | No completeness evidence available                    |

For downstream analyses, the following categories are considered acceptable:

* `complete`
* `coding_complete_only`
* `likely_complete`

---

# Laboratory Sequence Filtering

To minimize experimental bias, likely laboratory-derived sequences are excluded.

Records containing keywords such as:

```
vero
passage
cell line
adapted
in vitro
```

within their combined metadata are removed prior to grouping.

---

# Isolate Reconstruction

Sequences are grouped using the following priority:

1. BioSample accession
2. Isolate / strain + collection date + location + host
3. Accession number (fallback)

A confidence score is calculated from the availability of the four metadata fields (host, location, date, isolate/strain) and reported as a value between **0** and **1**.

---

# Representative Segment Selection

Only isolate groups containing complete **S**, **M**, and **L** segments are retained.

For each segment, the representative sequence is selected using the following priority:

1. Highest completeness
2. RefSeq accession
3. Longest nucleotide sequence

---

# Host Prediction

A representative host is inferred independently for each segment using the most frequently observed host among complete records belonging to that segment.

---

# Reassortment Detection

Each isolate receives one of three reassortment labels:

| Label     | Meaning                                                 |
| --------- | ------------------------------------------------------- |
| `0`       | All predicted segment hosts agree                       |
| `1`       | Segment host predictions differ (possible reassortment) |
| `unknown` | Insufficient evidence                                   |

---

# Output Files

## 1. `orthohantavirus_ml_ready.json`

Machine learning dataset containing one curated record per retained accession.

Each record includes:

* accession
* group ID
* segment
* completeness
* host metadata
* protein accession(s)
* sequence length
* group members
* predicted hosts for S, M, and L

---

## 2. `orthohantavirus_host_prediction_master.csv`

Group-level summary containing:

* group ID
* predicted primary host
* representative protein accessions
* reassortment status
* group size
* observed hosts
* country
* collection year

---

## 3. `orthohantavirus_summary.json`

Pipeline execution summary including:

* total records processed
* retained isolate groups
* segment distribution
* completeness statistics
* reassortment statistics
* parsing errors
* laboratory-filtered records
* output file locations

---

# Notes

### Strengths

* Fully reproducible workflow
* Multi-source metadata integration
* Robust fallback logic
* Automatic isolate reconstruction
* Segment-aware completeness assessment
* Representative sequence selection
* Machine learning–ready outputs

### Current Limitations

* Input directories are hard-coded.
* Laboratory filtering relies on keyword matching.
* `complete_members` in the CSV output is currently fixed to `3`, whereas the actual number of records is correctly stored in `group_size`.
