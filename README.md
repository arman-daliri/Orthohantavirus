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
| `--verbose` | off | Debug logging |

The transform is GIL-bound and runs serially; a full run takes ~16s on the 15,991-record package.

Requires Python 3.11+ (uses `datetime.UTC`) with `pandas`, `pyarrow`, `numpy`. No alignment or
phylogenetics dependencies for a default run -- those only load if Stage 2 is run (see below).

## Pipeline stages

The pipeline is three independent scripts, each runnable any time after the raw NCBI package
exists (none requires another to have already run):

| Stage | Script | Needs | Produces |
|---|---|---|---|
| 1 | `all_hanta.py` | `pandas`/`pyarrow`/`numpy` only | The core tables above, including `ortho_ml_sequence_table.csv` |
| 2 | `phylogeny.py` | `mafft`/`iqtree` on PATH, Biopython, `dendropy` | `ortho_reassortment_table.csv`, `phylogenetics/<species>/<segment>/` |
| 3 | `reassortment_detector.py` | `mafft`/`iqtree` on PATH, Biopython, `dendropy` | Score a new, not-yet-published candidate genome against the 5-gate framework |

`ortho_risk_table.csv` (Stage 1) is always phylo-blind. Combine it with Stage 2's evidence by
joining `ortho_reassortment_table.csv` on `isolate_key` — there's no automatic recombination
between the two scripts by design (see `docs/superpowers/specs/2026-08-13-three-stage-pipeline-
split-design.md` §4 for why).

### Stage 2: `phylogeny.py` (MAFFT + IQ-TREE reassortment analysis)

```bash
python3 phylogeny.py --input . --output output
```

| Argument | Default | Purpose |
|---|---|---|
| `--mafft-bin` | PATH lookup | Path to `mafft` |
| `--iqtree-bin` | PATH lookup | Path to `iqtree` |
| `--phylo-min-taxa` | `4` | Minimum candidate isolates in a species before building segment trees |
| `--bootstrap` | `1000` | IQ-TREE ultrafast bootstrap replicates |
| `--min-bootstrap-support` | `95.0` | Minimum ufboot support to call a topology conflict *supported* |
| `--phylo-seed` | `42` | Fixed IQ-TREE random seed, for reproducible trees |
| `--phylo-threads` | cores − 1 | Threads for `mafft` and `iqtree` |

Requires `mafft` and `iqtree` (IQ-TREE 3.x; `iqtree2` builds work too, pass `--iqtree-bin`) plus
Biopython and `dendropy`. Set up once with the pinned environment file rather than an ad hoc
`conda create`:

```bash
conda config --set channel_priority strict
conda env create -f environment.yml
```

Channel order matters here, not just as style: `-c bioconda -c conda-forge` (bioconda first)
resolves Biopython from bioconda's own channel, whose newest build is `1.70` and drags in a
`Python 3.6` / `numpy 1.12` toolchain to satisfy it -- silently, with no error. `environment.yml`
lists `conda-forge` first and `bioconda` second (the order bioconda's own docs recommend) so
Biopython, pandas, numpy and pyarrow resolve from conda-forge's current builds, while `mafft` and
`iqtree` -- which only exist in bioconda -- still resolve there. `mafft`/`iqtree` are pinned to
exact versions (`7.525` / `3.1.3`); the Python-side packages are floors, not patches, since exact
patch availability on conda-forge isn't something this repo can assert.

Then either `conda activate hanta_phylo` before running this script, or leave your normal
environment active and pass `--mafft-bin`/`--iqtree-bin` pointing at that env's `bin/`. **Don't
just prepend that env's `bin/` to `PATH`** without activating fully -- an unpinned env resolved
this way has been observed to install a Python 3.6 interpreter as `python3`, which would shadow
the one actually running this script.

Runtime is dominated by IQ-TREE's ModelFinder + ultrafast bootstrap and scales with both the
number of candidate isolates per species and CDS length (L segments are the slowest, ~6.4 kb).
A handful of taxa finishes in seconds; species with hundreds of qualifying isolates (PUUV, TULV,
HTNV, SEOV in this package) can take much longer. There is no built-in subsampling -- every
qualifying isolate is analyzed.

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

Result on this package: **1,795 groups with all three S, M and L segments present**, and 223
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
L 6.5–6.6 kb. The pipeline applies S and M as stated and widens L to 6.3–6.8 kb, because
submitted L segments sit just outside the ICTV figure at both ends; see the note on
`SEGMENT_LENGTH_RANGE` in `all_hanta.py`.

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

`isLabHost`, `labHost` and `isVaccineStrain` are authoritative NCBI fields and are used directly.

**Lab-derived isolates are flagged, not filtered, in the ML matrix.** `isLabHost` marks the
immediate recorded host, not sequence contamination: of 338 lab-host records (164 isolates), 21
isolates have `host_taxid` `9606` — lab-acquired/occupational human infections concentrated in the
highest-CFR species (HTNV, ANDV, PUUV, DOBV) — and most of the rest are samples from
colony-maintained reservoir animals, not culture-adapted sequences. Blanket-excluding them from a
pathogenicity-prediction task discarded real positive evidence of human/reservoir infection; a
study comparing human-derived and reservoir Orthohantavirus sequences for zoonotic-risk assessment
found the human sequences are themselves the spillover-risk baseline, not noise to remove
(<https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-025-02059-0>).
`ortho_ml_feature_matrix.csv` carries an `is_lab_derived` column so a model can be trained
with or without them and compared, rather than losing the isolates outright.

**Vaccine strains remain excluded from the ML matrix.** A vaccine strain is deliberately
attenuated, so labelling its sequence with its parent species' wild-type pathogenicity class would
inject a real signal error, not a hypothetical one — unlike `isLabHost`, this is a case where
exclusion is the correct default.

**Stage 2 (`phylogeny.py`) is unaffected by this and keeps excluding both lab-derived and vaccine-strain
isolates** — a passage/attenuation history is a tree-topology risk (branch-length distortion)
regardless of whether the sequence is informative as an ML feature.

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

## Reassortment

**By default, this pipeline does not call reassortment.** A classifier that did was built and
then withdrawn after it failed validation against the literature. What follows is that
withdrawal, and the rigorous replacement that now exists as an opt-in stage.

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

The k-mer diagnostic is retained as an informational signal (`reassortment_status`,
`sl_rank_separation`, `m_rank_offset`, still emitted in `ortho_isolate_table.csv` and
`ortho_risk_table.csv`); only the standalone `ortho_segment_divergence.csv` file was dropped,
since it duplicated data already in those two tables. Stage 2 replaces the *call* (a real
reassortment verdict), not the underlying measurement. **To detect reassortment properly**, build
per-segment alignments (MAFFT) and trees (IQ-TREE) and test topological incongruence between the
S, M and L trees. The constraints worth testing against: M-segment exchange dominates, S and L
cosegregate from the same parent, and events occur within a reservoir host species or between
sympatric lineages.

The one trustworthy reassortment-adjacent signal here is `taxonomic_incongruence`: one isolate
label whose segments NCBI assigned to more than one species (70 isolates). That is metadata, not
inference, and it does feed the risk score. It is computed across groups, not within — the
grouping key contains the species, so a within-group check can never fire.

### Stage 2 in detail: the rigorous reassortment detection

Stage 2 (`phylogeny.py`) runs the real thing this section always pointed to: per-segment CDS
nucleotide alignments (MAFFT) and maximum-likelihood trees with ultrafast bootstrap support
(IQ-TREE 3.x), compared for topological incongruence between S, M and L. See the argument table
above for setup and the runtime caveat.

**Candidate selection.** Three groups, unioned (each is independently a reason an isolate matters
for the reassortment/pathogenicity goal), then intersected with a hard requirement — CDS
available for all of S, M and L, since a topology comparison needs every isolate placed in all
three trees:

- **human-affecting-or-likely** — `pathogenicity_class` in `{pathogenic_hcps, pathogenic_hfrs,
  rare_or_opportunistic, suspected}`, OR a not-yet-curated species observed directly in a human
  host in this package (`host_role == spillover_host_human`) — the "likely" half a static curated
  table alone cannot express.
- **divergence candidates** — isolates the k-mer diagnostic already flagged
  `requires_phylogenetic_confirmation`.
- **complete triplets** — `has_complete_triplet`.

Lab-derived, vaccine-strain, and `unclassified`/genus-only isolates are excluded — the latter
because pooling them under `species_name` would build a tree across genuinely unrelated lineages
that only share a placeholder label.

**Per species** (≥ `--phylo-min-taxa`, sorted for determinism): CDS nucleotide sequences per
segment are extracted via `select_primary_cds` (role-matched from `cds.fna`, longest-wins among
matches — correctly prefers a combined G1+G2 precursor CDS over either half when both are
annotated), aligned with `mafft --auto`, and a tree is built with `iqtree -m MFP -B <bootstrap>
--seqtype DNA --seed <fixed>` (model selection + ultrafast bootstrap; DNA type is asserted, not
auto-detected, because heavily-N'd partial CDS can otherwise confuse the heuristic). Tip labels
are synthetic (`I0000`, `I0001`, ...) rather than `isolate_key` directly, to avoid any risk of
Newick-incompatible characters corrupting a tree.

**Per isolate**, from the three real trees: patristic (branch-length) nearest-neighbour and
rank-separation diagnostics — the *same* statistical shape as `sl_rank_separation`/
`m_rank_offset` above, now computed from real phylogenetic distance instead of k-mer composition
distance, which is the actual fix for the failure mode documented earlier in this section. A
candidate is only called `supported_reassortment_candidate` if the rank pattern matches the
literature's M-segment-exchange signature (S and L agree; M's neighbour differs) **and** the
relevant M-tree clade clears `--min-bootstrap-support` (default 95) — anything below that is
`weak_signal_low_bootstrap_support`, not a finding.

**Per species**, Robinson-Foulds distance (via `dendropy`) between each pair of segment trees on
the shared taxon set — an overall measure of how much the three trees disagree, independent of
any single isolate's placement.

**Nothing is silently dropped.** Every candidate isolate appears exactly once in
`ortho_reassortment_table.csv`, `phylo_status` one of `analyzed` / `failed` (with the specific
MAFFT/IQ-TREE error) / `species_below_min_taxa` — a species/segment failure never aborts the run,
it just marks that species' isolates `failed` and moves on, the same per-unit error isolation
`transform_one_record` uses for accessions. Alignments and trees are written under
`<output>/phylogenetics/<species>/<segment>/` for inspection or reuse.

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
| `ortho_accession_table.csv` | 15,991 × 55 | One row per accession: taxonomy ranks, segment, completeness, host, pathogenicity, host-jump status, QC flags |
| `ortho_species_table.csv` | 39 × 32 | One row per canonical ICTV species: merged taxids, legacy synonyms, risk profile, segment and host counts |
| `ortho_isolate_table.csv` | 9,852 × 52 | Reconstructed isolates with S/M/L representatives, triplet status, concatenated protein, divergence diagnostics |
| `ortho_host_table.csv` | 253 × 13 | One row per host taxon with risk prior and virus species carried |
| `ortho_risk_table.csv` | 9,852 × 33 | Pathogenicity class, CFR band, citation, and an auditable additive `surveillance_priority_score` |
| `ortho_ml_feature_matrix.csv` | 7,651 × 236 | Raw features + labels, with `group_id` and `year` for leakage-free CV |
| `ortho_ml_sequence_table.csv` | 206 × 6 | One row per isolate with a coding-complete S+M+L triplet (lab-derived and vaccine-strain isolates excluded): species name, isolate name, raw S/M/L genomic sequences, CFR range |
| `ortho_suppressed.csv` | 1,032 × 16 | Non-representative accessions with the reason |
| `ortho_qc_report.csv` | 15,991 × 11 | Per-accession QC including failed transforms with traceback |
| `ortho_feature_scaling.json` | — | Per-feature mean/std, to be fitted **inside** your CV fold |
| `ortho_summary.json` | — | Provenance, integrity results, counts, methods and stated limitations |
| `ortho_reassortment_table.csv` *(Stage 2, `phylogeny.py`)* | one row per candidate isolate | Real per-segment tree diagnostics, RF distances, bootstrap-gated `reassortment_call` — join on `isolate_key` with `ortho_risk_table.csv`, which is always phylo-blind |
| `phylogenetics/<species>/<segment>/` *(Stage 2, `phylogeny.py`)* | — | MAFFT alignment + IQ-TREE outputs per species/segment, for inspection or reuse |

All tables are plain CSV (UTF-8, LF line endings). `species_taxid`, `observed_taxid`, `host_taxid`
and `year` are all-digit strings by design (they are identifiers, not quantities, and some rows
leave them blank) — a bare `pd.read_csv(...)` will infer them as numeric and either error on the
blank rows or silently turn them into floats. Read these columns with `dtype=str` if you need them
back as they were written.

### Column reference

Every column in every table, grouped by file. `int()`-coded booleans (`is_lab_derived`,
`has_complete_triplet`, etc.) are `0`/`1`; `;`-joined columns (`sra_accessions`,
`protein_accessions`, `warning_flags`, etc.) hold zero or more values in one cell.

<details>
<summary><code>ortho_accession_table.csv</code> — one row per GenBank accession</summary>

**Identity & taxonomy** — `accession`, `sequence_hash` (CRC32 from NCBI, zero-padded); `species_name`/`species_taxid`/`species_status` (canonical ICTV species, see [Taxonomy resolution](#taxonomy-resolution)); `observed_organism`/`observed_taxid` (NCBI's raw, pre-canonicalization label); `genus`/`subfamily`/`family`/`order`/`class`/`phylum` (from `virus.lineage`); `virus_abbrev`/`pathogenicity_class`/`expected_reservoir_order` (curated, from `SPECIES_RISK`); `host_jump_status` (`host_order_mismatch` if this record's host order differs from its species' expected reservoir order).

**Segment & completeness** — `segment` (S/M/L or empty if unresolved); `segment_source`/`segment_confidence` (which evidence tier resolved it, see [Segment resolution](#segment-resolution)); `completeness`/`completeness_source` (NCBI's `PARTIAL`/`COMPLETE` field, authoritative); `is_coding_complete` (length-qualified completeness, the basis for `has_complete_triplet` downstream); `nt_length`.

**Sample / isolate metadata** — `biosample`, `bioproject`, `sra_accessions`; `isolate`/`strain` (submitter-provided labels, used to group S/M/L into one isolate); `collection_date`/`year`; `location`/`country`/`region`.

**Host** — `host_name`/`host_taxid`/`host_order`/`host_family`; `host_source` (`data_report` or `biosample`, whichever supplied it); `host_role` (`reservoir_*`, `spillover_host_human`, or `laboratory_host`, from lineage rank not string matching); `isolation_source`/`tissue_source`; `is_lab_host`/`lab_host` (NCBI's `isLabHost`/`labHost` fields — the immediate recorded host, not sequence provenance, see [Host classification](#host-classification)); `is_vaccine_strain` (NCBI's `isVaccineStrain`).

**Protein / annotation** — `protein_accessions`, `annotation_products`, `protein_product`, `has_nss_orf` (encodes the NSs virulence-factor ORF).

**Grouping / provenance / QC** — `group_key` (the isolate this accession was grouped into, see [Isolate grouping](#isolate-grouping)); `grouping_tier`/`group_confidence` (which of the 5 grouping tiers fired, 1.00 down to 0.10); `is_refseq`; `source_database`; `warning_flags` (e.g. `missing_host`, `lab_host`, `vaccine_strain`, `missing_genomic_sequence`, `no_annotation`, length anomalies).
</details>

<details>
<summary><code>ortho_species_table.csv</code> — one row per canonical ICTV species (rollup)</summary>

`species_name`/`species_status`; `observed_taxids`/`canonical_taxids` (every raw taxid that canonicalized to this species); `legacy_synonyms` (superseded organism labels merged in); `subfamily`/`family`/`order`; every key from `species_risk_profile` (`virus_abbrev`, `disease_syndrome`, `cfr_band`, `cfr_band_label`, `cfr_reported_range`, `reservoir_host`, `endemic_region`, `pathogenicity_class`, `risk_citation`, `reservoir_order` — see [Pathogenicity labels](#pathogenicity-labels)); `reported_reassortment` (literature summary, informational only, not computed); `accession_count`/`isolate_count`/`complete_triplet_count`; `s_count`/`m_count`/`l_count`; `distinct_hosts`/`top_host`; `human_derived_count`; `distinct_countries`/`top_country`; `year_min`/`year_max`.
</details>

<details>
<summary><code>ortho_isolate_table.csv</code> — one row per reconstructed isolate (the core biological unit)</summary>

**Identity** — `isolate_key` (the grouping key, see [Isolate grouping](#isolate-grouping)); `species_name`/`species_taxid`/`species_status`; `virus_abbrev`/`pathogenicity_class`; `kmer_cluster_id` (within-species cluster, see [Species clustering](#species-clustering)); `biosample`/`isolate`/`strain`.

**Host & sample metadata** — `host_name`/`host_taxid`/`host_order`/`host_family`/`host_role`; `collection_date`/`year`; `location`/`country`/`region`; `grouping_tier`/`grouping_confidence`.

**Segments** — `s_accession`/`m_accession`/`l_accession` (the chosen representative accession per segment); `s_completeness`/`m_completeness`/`l_completeness`; `s_nt_length`/`m_nt_length`/`l_nt_length`; `segments_present`/`coding_complete_segments`; `triplet_status` (`coding_complete_triplet` / `all_segments_present_not_all_complete` / `partial_triplet_with_complete_segment` / `incomplete_triplet`); `has_complete_triplet` (`triplet_status == coding_complete_triplet` — the flag `ortho_ml_sequence_table.csv` filters on); `concatenated_protein_aa`/`concatenated_protein_length` (S+M protein, the ML feature matrix's raw input); `s_protein_accessions`/`m_protein_accessions`/`l_protein_accessions`; `all_accessions` (every accession folded into this isolate, including suppressed duplicates).

**Flags & diagnostics** — `is_lab_derived`/`is_vaccine_strain`/`has_nss_orf`; `metadata_conflict` (segments disagreed on collection date or location); `taxonomic_incongruence` (segments assigned to more than one species); `cross_species_partners`; `divergence_assessable`/`reassortment_status`/`sl_rank_separation`/`m_rank_offset` (k-mer divergence diagnostics — **not** reassortment calls, see [Reassortment](#reassortment)); `warnings`.
</details>

<details>
<summary><code>ortho_host_table.csv</code> — one row per host taxon (rollup)</summary>

`host_name`/`host_taxid`/`host_order`/`host_family`/`host_role` (modal role across its records); `host_risk_prior` (curated, from `HOST_ORDER_RISK_PRIOR`); `accession_count`/`isolate_count`; `distinct_virus_species`/`top_virus_species`/`all_virus_species`; `representative_biosample`/`representative_isolate` (first non-empty example, for spot-checking).
</details>

<details>
<summary><code>ortho_risk_table.csv</code> — one row per isolate, surveillance triage ranking (sorted by score, descending)</summary>

**Risk profile** — `isolate_key`/`species_name`/`virus_abbrev`/`pathogenicity_class`/`disease_syndrome`; `cfr_band`/`cfr_band_label`/`cfr_reported_range` (see [Pathogenicity labels](#pathogenicity-labels)); `reservoir_host`/`endemic_region`/`risk_citation`.

**This isolate** — `host_role`/`host_order`/`expected_reservoir_order`/`host_jump_status`/`observed_host`; `country`/`year`; `triplet_status`/`has_complete_triplet`/`has_nss_orf`; `reassortment_status` (k-mer diagnostic, informational); `taxonomic_incongruence`/`cross_species_partners`; `phylo_reassortment_call` (from Stage 2, empty string if that stage wasn't run or didn't reach this isolate).

**Score** — `priority_cfr_term`, `priority_human_spillover_term`, `priority_taxonomic_incongruence_term`, `priority_host_jump_term`, `priority_triplet_term`, `priority_nss_term`, `priority_reassortment_term` (one column per additive term — always present and `0.0` if Stage 2 wasn't run), and `surveillance_priority_score` (their sum). See [Risk score](#risk-score) below for the weight of each term.
</details>

<details>
<summary><code>ortho_ml_feature_matrix.csv</code> — one row per ML-usable isolate (k-mer feature matrix + labels)</summary>

Usable = has a concatenated S+M protein and is not `is_vaccine_strain` (`is_lab_derived` isolates are *kept*, flagged — see [Host classification](#host-classification)).

**Identity / grouping** — `isolate_key`, `group_id` (= `species_name`, split your CV on this), `kmer_cluster_id`, `year`, `species_name`/`species_taxid`, `host_name`/`host_taxid`/`host_order`/`host_role`, `country`.

**Segment / quality features** — `s_nt_length`/`m_nt_length`/`l_nt_length`, `n_segments_present`, `n_coding_complete`, `has_complete_triplet`, `has_nss_orf`, `protein_length`, `divergence_assessable`, `sl_rank_separation`, `m_rank_offset`, `taxonomic_incongruence`, `host_jump_status`, `is_lab_derived` (retained as a feature to ablate/stratify on, not filtered).

**Effort / confound controls** (prefix `effort_`, one block per species — sampling intensity, not biology; see the `build_ml_tables` docstring) — `effort_n_accessions`, `effort_n_isolates`, `effort_n_countries`, `effort_n_years`, `effort_first_year`, `effort_years_since_first`, `effort_n_complete_triplets`.

**Labels** — `label_pathogenicity_class`, `label_disease_syndrome`, `label_cfr_band`, `label_is_human_derived`.

**Features** — `kmer_<3-mer>` × up to `--max-features` (default 200) columns: raw counts of the top-variance amino-acid 3-mers over the concatenated S+M protein, unscaled (see `ortho_feature_scaling.json` and [Using the ML tables](#using-the-ml-tables)).
</details>

<details>
<summary><code>ortho_ml_sequence_table.csv</code> — one row per qualifying isolate, raw S/M/L sequences for sequence-based ML</summary>

Built by `build_ml_sequence_rows`. Row filter (all three, isolate-level): `has_complete_triplet` (coding-complete S **and** M **and** L), **and not** `is_lab_derived`, **and not** `is_vaccine_strain` — stricter than `ortho_ml_feature_matrix.csv`, which keeps `is_lab_derived` isolates flagged rather than filtered; a raw sequence fed to a model as "wild-type" input is exactly the case passage history would mislead, unlike a k-mer composition feature. An isolate is also silently dropped if its genomic FASTA sequence is missing for a segment despite being coding-complete (logged as a count, not per-row).

**How 15,991 accessions become 206 rows**, on this package:

| # | Step | Function | Count | What's removed and why |
|---|---|---|---|---|
| 1 | Raw accessions | `iter_jsonl(data_report)` | 15,991 | — (starting point) |
| 2 | Reconstructed into isolates | `build_isolate_records` | 9,852 | Accessions sharing an isolate/strain/biosample label are grouped into one S/M/L isolate (see [Isolate grouping](#isolate-grouping)) — this is a many→fewer collapse, not a filter |
| 3 | Isolates with a full-complete S+M+L triplet | `has_complete_triplet` | 223 | −9,629: missing one or more segments outright, or a present segment isn't NCBI-`complete` + length-qualified (see [Completeness](#completeness)) |
| 4 | Minus laboratory-host isolates | `is_lab_derived` | 206 | −17: host recorded as a laboratory host on at least one segment (`isLabHost`) |
| 5 | Minus vaccine-strain isolates | `is_vaccine_strain` | 206 | −0 in this package (none of the 223 complete triplets are vaccine strains) |
| 6 | Minus missing genomic FASTA sequence | `genomic_index` lookup | **206** | −0 in this package (every remaining isolate's S/M/L genomic sequence was present in `genomic.fna`) |

Step 3 is where almost all of the reduction happens — **not** the lab/vaccine filter (steps 4–5), which only removes 17 isolates on top. Most reconstructed isolates simply never had all three segments sequenced at all (see [Isolate grouping](#isolate-grouping): only 1,795 of 9,852 isolates have all three segments present at *any* completeness, before the coding-complete requirement even applies). Host *species* (human vs. reservoir animal) is never a filter at any step — see the row counts by `host_role` in the table itself for the actual mix.

- `species_name` — the isolate's canonical ICTV species (`isolate.species_name`, see [Taxonomy resolution](#taxonomy-resolution)). Never blank — every isolate in this table came through species resolution before grouping. This is also what disambiguates `isolate_name`, see below.
- `isolate_name` — `isolate.isolate`, falling back to `isolate.strain`, falling back to `isolate_key` if both are blank. Note this is the **submitter-provided label verbatim** — it can be a bare number like `"3"` or `"136"` when that's literally what was submitted to NCBI (`isolate.name`/`source` in `data_report.jsonl`); it is not unique across species (two different species can each have an isolate named `"3"`), only unique in combination with `species_name`.
- `s_sequence`/`m_sequence`/`l_sequence` — full genomic segment sequence (not the CDS-only sequence Stage 2 uses), straight from the genomic FASTA, keyed by `s_accession`/`m_accession`/`l_accession`.
- `cfr_reported_range` — `species_risk_profile(species_name)["cfr_reported_range"]`, the literature CFR range string (e.g. `"0.1-0.4%"`), verbatim, not the ordinal `cfr_band`.
</details>

<details>
<summary><code>ortho_suppressed.csv</code> — non-representative accessions, with why</summary>

`accession`/`sequence_hash`; `species_name`/`species_taxid`; `segment`/`completeness`/`nt_length`; `biosample`/`isolate`/`strain`; `collection_date`/`location`; `host_name`/`host_taxid`; `group_key`; `suppression_reason` (`duplicate_sequence_of:<accession>`, `segment_unresolved`, one of `warning_flags`, or `not_best_representative_for_group`).
</details>

<details>
<summary><code>ortho_qc_report.csv</code> — one row per accession, transform pass/fail</summary>

`accession`; `group_key`/`grouping_tier`/`group_confidence`; `species_name`; `segment`/`segment_source`; `completeness`/`nt_length`; `transform_ok` (`False` only if the accession threw during transform — a single bad record never aborts the run); `warnings`; `traceback` (present only on `transform_ok == False` rows).
</details>

<details>
<summary><code>ortho_feature_scaling.json</code> — fit-outside-your-fold reference stats for the ML feature matrix</summary>

`k`/`alphabet`/`selection` (method used to pick the `kmer_*` columns); `n_features`/`features` (the exact column names selected); `mean`/`std` per feature, computed on the **full** table — a caution/`note`/`caveat` field spells out why these must be refit inside your own CV fold, not applied before splitting; `recommended_cv`; `n_label_bearing_groups` (distinct species with a real, non-`unknown` CFR band).
</details>

<details>
<summary><code>ortho_summary.json</code> — provenance, integrity, counts, methods, limitations</summary>

Top-level sections: `pipeline_version`/`generated_at_utc`/`dataset_catalog_api_version`; `provenance` (taxon id, `datasets` CLI version, package release-date range); `integrity` (md5 + sequence-hash verification results, per-file); `counts` (accessions, species, isolates, complete triplets, host/biosample/human/lab/vaccine record counts); a set of `*_counts` breakdowns (`segment_counts`, `segment_source_counts`, `completeness_counts`, `grouping_tier_counts`, `triplet_status_counts`, `host_order_counts`, `host_role_counts`, `pathogenicity_counts`, `host_jump_counts`, `reassortment_status_counts`, `qc_warning_counts`); `taxonomic_incongruence_count`; `inputs` (package dir name and every file found in it); and `methods` — prose mirroring this README, with nested `reassortment` (withdrawal rationale and numbers) and `reservoir_curation_audit` objects, a `LIMITATION` string (alignment-free k-mer distances vs. true ICTV PED), and `citations`.
</details>

<details>
<summary><code>ortho_reassortment_table.csv</code> — Stage 2 only, one row per candidate isolate per requested scope</summary>

Every field is always present (even for a not-yet-analyzed or failed row) so the table has one stable, joinable schema — see `_phylo_row`'s docstring. Full methodology in [Stage 2 in detail](#stage-2-in-detail-the-rigorous-reassortment-detection).

**Identity / selection** — `isolate_key`/`species_name` (the isolate's true species, never the phylogroup label); `candidate_tags` (why it was selected: `human_affecting_or_likely` / `divergence_candidate` / `complete_triplet`, can be more than one); `analysis_scope` (`species` or `phylogroup` — a candidate can appear once per scope under `--phylo-scope=both`, the default); `phylogroup_assignment_basis` (`curated_reservoir` or `observed_host_order`, only meaningful at `analysis_scope == phylogroup`).

**Tree status** — `n_taxa_in_species_tree` (taxa in *this row's* tree, despite the name — imprecise for phylogroup-scope rows on purpose, see the code comment); `phylo_status` (`analyzed` / `failed` / `species_below_min_taxa`); `phylo_error`.

**Per-isolate tree diagnostics** (real phylogenetic distance, the same statistical shape as the k-mer diagnostics but the actual fix for their failure mode) — `tree_s_nearest_isolate`/`tree_m_nearest_isolate`/`tree_l_nearest_isolate`; `tree_s_distance`/`tree_m_distance`/`tree_l_distance` (patristic/branch-length); `tree_s_rank_offset`/`tree_m_rank_offset`/`tree_l_rank_offset`; `tree_sm_rank_separation`/`tree_sl_rank_separation`/`tree_ml_rank_separation`; `discordant_segment`; `matches_canonical_m_exchange` (S and L agree, M's neighbour differs — the literature's M-exchange signature).

**Bootstrap gate & call** — `m_clade_bootstrap_support`/`discordant_clade_bootstrap_support` (`None`, not `0.0`, when no discordant clade exists to measure — a real 0% support must stay distinguishable from "not applicable"); `bootstrap_support_source` (`mrca` / `nearest_supported_ancestor` / `unsupported_root`, empty for not-yet-analyzed/`no_conflict` rows); `reassortment_call` (`not_analyzed` / `no_conflict` / `weak_signal_low_bootstrap_support` / `supported_reassortment_candidate` — **only** the last one contributes to `priority_reassortment_term` in `ortho_risk_table.csv`).

**Per-tree (species or phylogroup) topology summary** — `species_rf_s_vs_m`/`species_rf_m_vs_l`/`species_rf_s_vs_l` (Robinson-Foulds distance between segment tree pairs, via `dendropy`) and their `_normalized` counterparts; `s_alignment_path`/`m_alignment_path`/`l_alignment_path`/`s_tree_path`/`m_tree_path`/`l_tree_path` (on-disk locations under `<output>/phylogenetics/<species>/<segment>/`).
</details>

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
artefact of sample size into a ranking meant to guide surveillance. `reassortment_call` from
Stage 2 is likewise **not yet wired into this score** — it lives in a separate,
opt-in output table so the default risk table's schema and values stay unchanged whether or not
Stage 2 was run. Join on `isolate_key` to combine them for your own weighting.

### Using the ML tables

Features are emitted **raw**. Fit scaling inside your cross-validation fold using
`ortho_feature_scaling.json`; the statistics there were computed on the full table and must not be
applied before splitting. Feature selection is unsupervised (top-variance k-mers).

**Split by `group_id` (species) or by `year`.** A random split badly overstates accuracy on this
data, because segments from the same outbreak are near-identical.

---

## Known limitations

- **No reassortment calls by default.** See above; run Stage 2 for the real
  MAFFT + IQ-TREE analysis. Even then, it only covers the three named candidate groups with full
  S+M+L CDS available (1,464 of 15,991 accessions' worth of isolates in this package) — not every
  isolate, and species below `--phylo-min-taxa` are reported but not analyzed.
- **Stage 2 (`phylogeny.py`) calls are a single ML tree search per segment, bootstrap-gated, not a formal
  incongruence test** (e.g. no SH-test/AU-test on constrained topologies). `supported_reassortment_
  candidate` means the rank pattern and bootstrap support are consistent with reassortment, not a
  statistically tested certainty.
- **k-mer distances are not PED.** The ICTV 0.1 cutoff needs aligned sequences and ML inference —
  Stage 2 builds real per-segment ML trees but does not compute a DEmARC-style
  concatenated-PED value; it is aimed at the reassortment question, not species re-demarcation.
- **223 coding-complete triplets** is a ceiling set by the data (1,339 COMPLETE records total),
  not by the pipeline. 1,795 groups have all three segments present at some completeness.
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
