# Orthohantavirus NCBI Dataset Preparation Pipeline

`all_hanta.py` — version 2.0.0

A reproducible pipeline that ingests an NCBI Datasets **Orthohantavirus** package and produces
accession-, isolate-, species- and host-level tables, pathogenicity labels, per-segment divergence
diagnostics, and ML-ready feature matrices for host-association and risk-classification work.

---

## Step 0 — Download the data

```bash
datasets download virus genome taxon 1980442 \
  --include genome,cds,protein,annotation,biosample \
  --filename all_hanta.zip

unzip all_hanta.zip -d all_hanta
```

## Step 1 — Stop git from corrupting the package

**Do this before anything else.** The package ships an `md5sum.txt` containing checksums of the
exact bytes NCBI published. On Windows, `git config core.autocrlf=true` rewrites LF to CRLF on
checkout for anything git considers text, which changes those bytes, breaks every checksum, and
makes results differ between Linux and Windows checkouts.

The repo's `.gitattributes` marks the package binary to prevent this. If your working copy was
already mangled, repair it once:

```bash
rm -f all_hanta/ncbi_dataset/data/* all_hanta/md5sum.txt all_hanta/README.md
git checkout -- all_hanta/
```

The pipeline verifies before doing any work, so you will know immediately either way.

## Step 2 — Run

```bash
python all_hanta.py --input . --output output
```

| Argument | Default | Purpose |
|---|---|---|
| `--input` | `.` | Package root containing `all_hanta/ncbi_dataset/data` |
| `--output` | `<input>/output` | Output directory |
| `--allow-checksum-mismatch` | off | Proceed despite failed integrity check (mismatch still recorded) |
| `--cluster-threshold` | `0.02` | Cosine k-mer distance for within-species single-linkage clustering |
| `--max-features` | `200` | Number of top-variance k-mer features in the ML matrix |
| `--min-species-size` | `8` | Minimum assessable isolates before divergence diagnostics are computed |
| `--threads` | cpu count | Accepted for compatibility; the transform is GIL-bound and runs serially |
| `--verbose` | off | Debug logging |

Requires Python 3.11+ (uses `datetime.UTC`) with `pandas`, `pyarrow`, `numpy`. No alignment or
phylogenetics dependencies.

---

## Reproducibility model

Integrity is established against NCBI's own published checksums, never against digests this
script invents.

**File level.** `md5sum.txt` is parsed and every shipped file verified. On mismatch the pipeline
also tests the CRLF-normalized form purely to *diagnose* the cause, then aborts:

```
CHECKSUM MISMATCH data_report.jsonl
    expected (NCBI): 9c228361ec1892882d118876fcdd263e
    observed (disk): 24dcee9c3c910a59c67093ab818135e7
    diagnosis      : on-disk bytes are CRLF but NCBI published the LF form ...
```

**Content level.** `nucleotide.sequenceHash` (CRC32, emitted per record by NCBI) is used for
sequence identity and de-duplication. It is independent of line endings and platform, so it
behaves identically everywhere even on a mangled checkout. NCBI strips leading zeros from these
values — observed widths in this package are 5–8 characters — so they are **zero-padded to 8**
before any comparison. They also expose 3,393 exact resubmissions (12,598 distinct hashes over
15,991 records) that would otherwise inflate group sizes.

`hashlib` is still used, but only to hash the pipeline's own outputs and configuration.

Runs are deterministic: two consecutive runs produce byte-identical table contents, and
`ortho_summary.json` contains no absolute paths.

---

## Taxonomy resolution

The organism-level `virus.taxId` is **not** usable as a species key. In this package it takes
**553 distinct values**, mixing strain-level taxa (`Hantaan virus Lee` 11601, `Four Corners
hantavirus` 31621) with legacy/current taxid pairs for the same species (`Puumala orthohantavirus`
1980486 vs `Orthohantavirus puumalaense` 3052493).

The node immediately below genus *Orthohantavirus* (taxId 1980442) in `virus.lineage` **is** the
ICTV species. Canonicalising on its binomial name merges the legacy/current pairs, collapsing
553 organism values to **38 species** plus `unclassified Orthohantavirus` (2,339 records).

The remaining `virus.lineage` nodes supply every higher rank for free, giving the full grouping
hierarchy that terminates naturally at the top:

```
accession → isolate → species → subfamily Mammantavirinae → family Hantaviridae
          → order Elliovirales → class Bunyaviricetes → phylum Negarnaviricota
```

Records that are unclassified or genus-only are retained with `species_status` set accordingly —
never silently merged into a named species.

## Isolate grouping

BioSample is present on only **87 of 15,991** records (0.54%), so a BioSample-first cascade
almost never fires. `isolate.name` is present on 15,491 records (96.9%) and is the field
submitters actually use to tie an isolate's S, M and L segments together.

| Tier | Key | Confidence |
|---|---|---|
| 1 | species + biosample + isolate | 1.00 |
| 2 | species + `isolate.name` | 0.95 |
| 3 | species + normalized isolate/strain label | 0.90 |
| 4 | species + strain + collection date + location | 0.75 |
| 5 | species + strain + location | 0.65 |
| — | accession only | 0.10 |

Result on this package: **1,777 groups with all three S, M and L segments present**, and 222
where all three are coding-complete. The coding-complete count is bounded by the data — only
1,339 of 15,991 records are `COMPLETE` at all — not by the grouping.

## Segment resolution

The `segment` field is dirtier than a three-value check assumes. Values observed here include
`Partial L`, `S segment`, `M-segment`, `M; medium`, `small (S)`, `G1`, `G2`, `ANONYMOUS`, `s`
and `N`, plus 3,558 records with no `segment` field at all. Evidence is applied in order:

| `segment_source` | Basis | Confidence | Count |
|---|---|---|---|
| `explicit` | normalized `segment` field (`G1`/`G2` → M, `N` → S) | 1.00 | 12,431 |
| `annotation` | CDS product-name majority vote | 0.90 | 3,117 |
| `defline` | genomic FASTA description | 0.75 | 104 |
| `length_heuristic` | falls in exactly one ICTV size band | 0.40 | 129 |
| `unresolved` | — | 0.00 | 210 |

Unresolved segments drop from 405 to 210.

## Completeness

NCBI's structured `completeness` is present on 100% of records and is authoritative
(PARTIAL 14,652 / COMPLETE 1,339). Length is used only to qualify a COMPLETE call as
coding-complete, never to override NCBI. ICTV segment sizes: S 1.2–2.1 kb, M 3.4–3.8 kb,
L 6.5–6.6 kb.

Length anomalies are flagged only where they contradict those ranges: partial sequences are
*expected* to be short. This drops length warnings from 5,856 to 9.

## Host classification

Host role is derived from `host.lineage` rank, not substring matching:

| Observed | Count | Role |
|---|---|---|
| Rodentia | 11,043 | `reservoir_rodent` |
| *Homo sapiens* | 2,091 | `spillover_host_human` |
| Eulipotyphla | 1,035 | `reservoir_shrew_mole` |
| Chiroptera | 83 | `reservoir_bat` |
| `isLabHost` true | 338 | `laboratory_host` |

`isLabHost`, `labHost` and `isVaccineStrain` are authoritative NCBI fields and are used directly;
lab-derived and vaccine records are excluded from the ML matrix.

**Host jumps.** Each record's host order is compared against its species' expected reservoir
order. 33 records are pathogenic viruses recovered outside their reservoir order — SEOV in
Chinese shrews and bats, ANDV-clade sequences in Brazilian bats. These are documented
host-sharing observations, not data errors, and are surfaced as
`host_jump_status = host_order_mismatch` and weighted in the risk score.

## Pathogenicity labels

`SPECIES_RISK` covers 38 species. Every entry carries its source URL. Species **absent** from the
table are labelled `unknown` — never silently "non-pathogenic".

CFR estimates genuinely conflict between sources, so the reported range is stored verbatim
alongside an ordinal band for modelling. Do not collapse them to a point estimate.

| Class | Species | CFR band |
|---|---|---|
| `pathogenic_hcps` | ANDV, SNV, CHOV, LANV, RIOMV, BAYV, BCCV | up to `high_gt_30pct` |
| `pathogenic_hfrs` | HTNV, DOBV, SEOV, PUUV | `low_lt_1pct` – `moderate_1_15pct` |
| `rare_or_opportunistic` | TULV | `low_lt_1pct` |
| `suspected` | THAIV, SANGV | `unknown` |
| `no_known_human_disease` | shrew-, mole- and bat-borne species | `none_known` |
| `unknown` | anything not curated | `unknown` |

HTNV is reported as both ~1% and 5–15% depending on source; both appear in `cfr_reported_range`.
ANDV is the only orthohantavirus with documented person-to-person transmission. TULV is held out
of `pathogenic_hfrs` deliberately — lumping its handful of mild case reports with HTNV/DOBV would
put 2,123 records into the same class as the viruses that actually drive HFRS mortality.

**Curation audit.** Curated reservoirs are checked against the modal observed host order on every
run. Six of the original 38 entries disagreed with the package's own host annotations at the order
level (asikkalaense, fugongense, lankaense, luxiense, sagercreekense, tigrayense) and were
corrected against the data. The audit runs every invocation so this class of error cannot silently
return.

## Reassortment — not called, and why

**This pipeline does not call reassortment.** A classifier that did was built and then withdrawn
after it failed validation against the literature.

The literature reports reassortment as frequent in PUUV (19.1–32% in Finnish studies), SNV, DOBV
and TULV, and absent or sparse in ANDV. Measured on this dataset, the classifier produced the
inverse ranking:

| Species | n | flagged | literature |
|---|---|---|---|
| ANDV | 68 | 72.1% | absent or sparse |
| CHOV | 8 | 100.0% | not reported |
| DOBV | 62 | 38.7% | frequent |
| PUUV | 416 | 7.7% | frequent (19.1–32%) |
| TULV | 313 | 3.5% | frequent |

Sweeping the discordance threshold from 0.25 to 0.85 scaled all rates down without ever changing
the ordering. The cause is intrinsic to the metric: cosine distance on 3-mer amino-acid
composition tracks within-species sequence diversity and sample size, not tree topology. Species
represented by many near-identical outbreak sequences get stable neighbour orderings and few
flags; species represented by a handful of divergent, scattered sequences get noisy orderings and
many flags. CHOV at 100% with n=8 is the clearest symptom.

What is emitted instead, in `ortho_segment_divergence.parquet`, is the honest underlying
measurement — per-segment nearest neighbours, k-mer distances, and normalized rank separations
(`sl_rank_separation`, `m_rank_offset`) — with every assessable isolate marked
`requires_phylogenetic_confirmation`.

**To detect reassortment properly**, build per-segment alignments (MAFFT) and trees (IQ-TREE) and
test topological incongruence between the S, M and L trees. The constraints worth testing against:
M-segment exchange dominates, S and L cosegregate from the same parent, and events occur within a
reservoir host species or between sympatric lineages.

The one trustworthy reassortment-adjacent signal here is `taxonomic_incongruence`: one isolate
label whose segments NCBI assigned to more than one species (70 isolates). That is metadata, not
inference, and it does feed the risk score. It is computed across groups, not within — the
grouping key contains the species, so a within-group check can never fire.

## Species clustering

Within-species single-linkage clustering on cosine distance over 3-mer profiles of the
concatenated N+GPC proteins. Columns are named `kmer_cluster_id` and `kmer_distance_threshold`,
deliberately **not** `ped_*`.

The ICTV species criterion for Hantaviridae is DEmARC with a pairwise evolutionary distance cutoff
of **0.1** under the WAG substitution model, which requires multiple sequence alignment and
maximum-likelihood inference. These clusters are an alignment-free proxy and are not PED values.
`ortho_summary.json` restates this.

---

## Outputs

| File | Rows × Cols | Contents |
|---|---|---|
| `ortho_accession_table.parquet` | 15,991 × 55 | One row per accession: taxonomy ranks, segment, completeness, host, pathogenicity, host-jump status, QC flags |
| `ortho_species_table.parquet` | 39 × 32 | One row per canonical ICTV species: merged taxids, legacy synonyms, risk profile, segment and host counts |
| `ortho_isolate_table.parquet` | 9,904 × 52 | Reconstructed isolates with S/M/L representatives, triplet status, concatenated protein, divergence diagnostics |
| `ortho_host_table.parquet` | 253 × 13 | One row per host taxon with risk prior and virus species carried |
| `ortho_risk_table.parquet` | 9,904 × 32 | Pathogenicity class, CFR band, citation, and an auditable additive `surveillance_priority_score` |
| `ortho_segment_divergence.parquet` | 9,904 × 18 | Per-segment nearest neighbours and rank separations — **diagnostics, not reassortment calls** |
| `ortho_ml_feature_matrix.parquet` | 7,529 × 228 | Raw features + labels, with `group_id` and `year` for leakage-free CV |
| `ortho_ml_segments.parquet` | 15,781 × 19 | Segment-level table for per-segment models |
| `ortho_suppressed.parquet` | 1,021 × 16 | Non-representative accessions with the reason |
| `ortho_qc_report.parquet` | 15,991 × 11 | Per-accession QC including failed transforms with traceback |
| `ortho_feature_scaling.json` | — | Per-feature mean/std, to be fitted **inside** your CV fold |
| `ortho_summary.json` | — | Provenance, integrity results, counts, methods and stated limitations |

### Risk score

`surveillance_priority_score` is a transparent sum, with every term emitted as its own column so
it can be audited or re-weighted:

| Term | Weight |
|---|---|
| `priority_cfr_term` | 4.0 / 2.5 / 1.0 / 0.5 / 0.0 by CFR band |
| `priority_human_spillover_term` | 2.0 if human-derived |
| `priority_taxonomic_incongruence_term` | 1.5 |
| `priority_host_jump_term` | 1.0 |
| `priority_triplet_term` | 1.0 if coding-complete triplet |
| `priority_nss_term` | 0.5 if an NSs ORF is annotated |

This is a triage aid for surveillance, not a validated epidemiological risk model. No term derives
from the k-mer divergence diagnostics — weighting a withdrawn classifier would have propagated an
artefact of sample size into a ranking meant to guide surveillance.

### Using the ML tables

Features are emitted **raw**. Fit scaling inside your cross-validation fold using
`ortho_feature_scaling.json`; the statistics there were computed on the full table and must not be
applied before splitting. Feature selection is unsupervised (top-variance k-mers).

**Split by `group_id` (species) or by `year`.** A random split badly overstates accuracy on this
data, because segments from the same outbreak are near-identical.

---

## Known limitations

- **No reassortment calls.** See above. Requires MAFFT + IQ-TREE.
- **k-mer distances are not PED.** The ICTV 0.1 cutoff needs aligned sequences and ML inference.
- **222 coding-complete triplets** is a ceiling set by the data (1,339 COMPLETE records total),
  not by the pipeline. 1,777 groups have all three segments present at some completeness.
- **2,386 records are `unclassified Orthohantavirus` or genus-only** and carry no species-level
  risk label. They are retained, not discarded.
- **CFR bands are literature ranges, not computed from this data.** Sources conflict; both values
  are recorded where they do.
- **1,719 records have no host** and 3,703 have no collection date.
- Grouping tiers 4 and 5 rarely fire, because `strain` is only populated from BioSample.
- Only the first BioSample accession in a record is resolved.

## Citations

- ICTV *Hantaviridae* / *Orthohantavirus* — species demarcation, segment sizes:
  <https://ictv.global/report/chapter/hantaviridae/hantaviridae/mammantavirinae/orthohantavirus>
- Zoonotic Hantaviridae with global public health significance — pathogenicity, CFR, reservoirs:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC10459939/>
- Reassortment events in the evolution of hantaviruses — frequencies, segment constraints:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6153690/>
- ECDC factsheet on orthohantavirus infections:
  <https://www.ecdc.europa.eu/en/infectious-disease-topics/hantavirus-infection/factsheet-orthohantavirus-infections>
- NCBI Datasets virus data report schema:
  <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/virus/>
